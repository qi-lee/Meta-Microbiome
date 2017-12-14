#!/usr/bin/perl
use strict;
use Pod::Usage;

###########################################################

=head1 SYNOPSIS

    Basic usage : perl Rm_fastq_which_mapping_fasta_V2.pl -i fasta -1 reads1 -2 reads2 -s singlereads

    Option Description :
              -h  --help_message
                  Print USAGE, DESCRIPTION and ARGUMENTS;
              -i  --input_fasta
                  Input fasta file;
              -1  --forward_reads
                  File with forward reads for paired-end / mate-pair library;
              -2  --reverse_reads
                  File with reverse reads for paired-end / mate-pair library;
              -s  --unpaired_reads
                  File with unpaired reads for single reads library;

 Author: liqi (lab of algal genomics)
 Created Time: Wed 26 Mar 2014 04:28:32 PM CST

=head1 AUTHOR

    liqi , liqi@ihb.ac.cn

=head1 COPYRIGHT

    Copyright - Lab of Algal Genomics

=cut

###########################################################


use threads;
use Getopt::Long;
use File::Basename;
use Bio::SeqIO;

my($help,$fasta,$forward_reads,$reverse_reads,$unpaired_reads);

GetOptions(
    'h!'     => \$help,
    'i=s{1}' => \$fasta,
    '1=s{1}' => \$forward_reads,
    '2=s{1}' => \$reverse_reads,
    's=s{1}' => \$unpaired_reads,
);

pod2usage if $help;

my @suffixlist = qw(.fastq .fq);
my @thr;

if($forward_reads){
    my $forward_reads_name = basename($forward_reads, @suffixlist);
    $thr[0] = threads->new(\&format,$forward_reads,$forward_reads_name,"1");
}
if($reverse_reads){
    my $reverse_reads_name = basename($reverse_reads, @suffixlist);
    $thr[1] = threads->new(\&format,$reverse_reads,$reverse_reads_name,"2");
}
if($unpaired_reads){
    my $unpaired_reads_name = basename($unpaired_reads, @suffixlist);
    $thr[2] = threads->new(\&format,$unpaired_reads,$unpaired_reads_name,"s");
}

for(my $j=0;$j<=$#thr;$j++){
    $thr[$j]->join;
}

sub format{
    my ($file,$filename,$lable) = @_;
    open(File,"$file");
    open(OUT,">$filename.format.fastq");
    while(<File>){
        if($.%4 == 1){
            s/\s+.*//;
            s/\/.?$//;
            s/^(\S+)/$1\/$lable/;
        }
        if($.%4 == 3){
            s/.*/+/;
        }
        chomp;
        print OUT "$_\n";
    }
    close OUT;
    close File;
}

