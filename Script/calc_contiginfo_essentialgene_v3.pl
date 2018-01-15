#!/usr/bin/perl
use warnings;
use strict;
###########################################################

=head1 SYNOPSIS

   Basic usage : perl calc_contiginfo_essentialgene_v2.pl -i fasta -1 reads1 -2 reads2 [-options]

   Example :
            1. $ perl calc_contiginfo_essentialgene_v2.pl -i scaffolds.fasta -1 R1_paired.fq -2 R2_paired.fq
            2. $ perl calc_contiginfo_essentialgene_v2.pl -i scaffolds.fasta -1 R1_paired.fq -2 R2_paired.fq -N 100 -X 1200 -p phred33 -b "--end-to-end --fast"
   
   Option Description :
            -h  --help_message
                Print USAGE, DESCRIPTION and ARGUMENTS;
            -r  --reference_file
                Input contigs/scaffolds;
            -1  --reads_1
            -2  --reads_2
            -m  --min_contig_length <int>
                Default: 500;
            -b  --bowtie2_presets_options
                Default = '--end-to-end --sensitive';
            -p  --phred
                Input quals [phred33/phred64] (Default: 'phred33')
            -I  --minins <int>
                Minimum insert size for paired-end alignment (Default: 0)
            -X  --maxins <int>
                Maximum insert size for paired-end alignment (Default: 500)
            -t  --threads <int>
                Number of alignment threads to launch (Default: 20)


=head1 AUTHOR

  liqi , liqi@ihb.ac.cn

=head1 COPYRIGHT
    
  Copyright - Lab of Algal Genomics

=cut

###########################################################

use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use threads;

my ($help,$contigs,$reads1,$reads2);
my $bowtie2_presets_options = "--end-to-end --sensitive";
my $min_contig_length = 500;
my $phred = "phred33";
my $threads = 20;
my $minins = 0;
my $maxins = 500;

GetOptions(
    'h!'     => \$help,
    'r=s{1}' => \$contigs,
    '1=s{1}' => \$reads1,
    '2=s{1}' => \$reads2,
    'p=s{1}' => \$phred,
    't=i{1}' => \$threads,
    'b=s{1}' => \$bowtie2_presets_options,
    'm=i{1}' => \$min_contig_length,
    'I=i{1}' => \$minins,
    'X=i{1}' => \$maxins,
);

pod2usage if $help;

(warn ( "Warning:\n\tMissing a required option!\n\tYou must appoint three input files : contigs.fasta, reads1.fastq, reads2.fastq\n\n") and  pod2usage)unless ($contigs && $reads1 && $reads2);

my $bowtie2_options = "-I $minins -X $maxins --$phred $bowtie2_presets_options";

my $in  = Bio::SeqIO->new( -file => "$contigs", -format => 'fasta' );
open(NewSeq,">Newseq.fasta");
while ( my $seq = $in->next_seq() ) {
    my $length = $seq->length;
    if ( $length >= $min_contig_length ) {
        my $id = $seq->id;
        my $seqence = $seq->seq;
        print NewSeq ">$id\n$seqence\n";
    }
}
close NewSeq;

#print "$contigs\t$reads1\t$reads2\t$min_contig_length\t$phred\n$bowtie2_options\n";

my $thr1 = threads->new( \&contig_info, "Newseq.fasta","$reads1","$reads2");
my $thr2 = threads->new( \&essential_gene, "$contigs");
$thr1->join;
$thr2->join;

mkdir "image" || die "can't create directory: image\n";
print "\nComplete!\n\n";
exit;

sub contig_info{
    my ($contig,$read1,$read2) = @_;
    system("\$CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct $threads -p '$bowtie2_options' $read1 $read2 pair $contig metagenome Rm_dup_Generate_cov");
    unlink glob "*.bt2";
    unlink "Newseq.fasta.fai";

    my %hash;
    open(COV,"Rm_dup_Generate_cov/coverage.percontig") || die "Error:\n\tCannot open the file: $!\n";
    while(<COV>){
        chomp;
        my($name,$coverage) = $_ =~ /(\S+)\s+(\S+)/;
        $hash{$name} = $coverage;
    }
    close COV;

    my $in  = Bio::SeqIO->new( -file => "$contig", -format => 'fasta' );
    open(OUT,">contig.info");
    print OUT "name\tlength\tgc\tcoverage\n";

    while ( my $seq = $in->next_seq() ) {
        my $length = $seq->length;
        my $id = $seq->id;
        my $seqence = $seq->seq;
        my $num_g = $seqence =~ tr/Gg//;
        my $num_c = $seqence =~ tr/Cc//;
        my $num_a = $seqence =~ tr/Aa//;
        my $num_t = $seqence =~ tr/Tt//;
        my $gc = sprintf("%.2f",($num_g + $num_c) / ($num_g + $num_c + $num_a + $num_t) * 100);
        print OUT "$id\t$length\t$gc\t$hash{$id}\n" if ($hash{$id} > 0);
    }
    close OUT;
}

sub essential_gene{
    my ($contig) = @_;
    system("prodigal -a temp_prodigal_orf.faa -i $contig -m -o temp_prodigal_stdout -p meta -q ");
    system(`cut -f1 -d " " temp_prodigal_orf.faa > temp_orf.faa`);
    system("hmmsearch --tblout hmm_orf_alignment.txt --cut_tc --notextw --cpu $threads /home/liqi/metagenome_analysis_microcystis/essential.hmm temp_orf.faa > temp_hmm_stdout");

    open(File,"hmm_orf_alignment.txt");
    open(Temp,">temp.txt");
    while(<File>){
        next if (/^#/);
        my ($name) = $_ =~ /^(\S+)_\d+\s+/;
        print Temp "$name\n";
    }
    close File;
    close Temp;
    system("sort temp.txt | uniq -c > essential_gene.info");
    unlink glob "temp*";
}
