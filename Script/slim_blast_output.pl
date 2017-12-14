#!/usr/bin/perl
use warnings;
use strict;
###########################################################

=head1 SYNOPSIS

   Basic usage : perl slim_blast_output.pl -i blast.out [-options]

   Example :
            1. $ perl slim_blast_output.pl -i blast.out -p 90
            2. $ perl slim_blast_output.pl -i blast.out -e 1e-5 -p 85.00 -s 70
   
   Option Description :
            -h  --help_message
                Print USAGE, DESCRIPTION and ARGUMENTS;
            -i  --input_file
                Input blast output file;
            -e  --min_evalue
                Default = '1e-5';
            -p  --min_PctIdentity(0~100)
                Minimum percentage of identical matches;
                Default = '75.00';
            -s  --min_FORSimilarity(0~100)
                Minimum full overlapping region similarity;
                Default = '50.00';


=head1 AUTHOR

  liqi , liqi@ihb.ac.cn

=head1 COPYRIGHT
    
  Copyright - Lab of Algal Genomics

=cut

###########################################################

use Getopt::Long;
use Pod::Usage;

my ($help,$blastoutfile);
my $minPctIdentity = 75.00;
my $min_FORSimilarity = 50.00;
my $min_evalue = 1e-5;

GetOptions(
    'h!'     => \$help,
    'i=s{1}' => \$blastoutfile,
    'p=f{1}' => \$minPctIdentity,
    's=f{1}' => \$min_FORSimilarity,
    'e=f{1}' => \$min_evalue,
);

pod2usage if $help;
pod2usage if !$blastoutfile;

open(File,"$blastoutfile");

while(<File>){
    chomp;
    my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$qlen,$sstart,$send,$slen,$evalue,$bitscore) = split /\s+/;
    next if ($pident < $minPctIdentity);
    next if ($evalue > $min_evalue);
    my $start = $sstart;
    my $end = $send;
    if ( $sstart > $send ) {
        $start = $slen - $sstart + 1;
        $end = $slen - $send + 1;
    }
    my $numLeftMismatch = $qstart < $start ? $qstart : $start;
    my $numRightMismatch = ($qlen - $qend) < ($slen - $end) ? ($qlen - $qend):($slen - $end);
    my $FORSimilarity = ($length / ($length + $numLeftMismatch + $numRightMismatch)) * 100;
    next if ($FORSimilarity <= $min_FORSimilarity);
    print "$qseqid\t$sseqid\t$pident\t$length\t$mismatch\t$gapopen\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$bitscore\n";
    #print "$qseqid\t$sseqid\t$pident\t$length\t$mismatch\t$gapopen\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$bitscore\t$FORSimilarity\n";
}
