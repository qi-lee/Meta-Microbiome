#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;

###########################################################

=head1 SYNOPSIS

    Function:

    Basic usage : perl script.pl [-options]

    Example :
             1. $ perl script.pl
             2. $ perl script.pl

    Option Description :
             -h  --help_message
                 Print USAGE, DESCRIPTION and ARGUMENTS
             -i  --Input file name
             -o  --output file name
             -g  --gi_taxid list file
             -l  --gi list file

=head1 AUTHOR

    liqi , liqi at ihb.ac.cn

=head1 COPYRIGHT

    Copyright - Lab of Algal Genomics
    Created Time: Thu 11 Aug 2016 11:45:32 PM CST

=cut

###########################################################

my ($help,$input,$gi_taxid,$gi_list);
my $output = "blast_add_staxids.out";

GetOptions(
    'h!'     => \$help,
    'i=s{1}' => \$input,
    'o=s{1}' => \$output,
    'g=s{1}' => \$gi_taxid,
    'l=s{1}' => \$gi_list,
);

my %hash;
open(List,"$gi_list");
my %gi = map {chomp; $_ => 1} (<List>);
close List;

open(Gi_Taxid,"$gi_taxid");
while(<Gi_Taxid>){
    my($key,$value) = $_ =~ /(\S+)\s(\S+)/;
    if (defined $gi{$key}) {
        $hash{$key} = $value;
    }
}
close Gi_Taxid;

open(OUT,">$output");
open(File,"$input");
while(<File>){
    chomp;
    my ($gi) = $_ =~ /\s+gi\|(\d+)\|/;
    print OUT "$_\t100\t$hash{$gi}\n" if exists($hash{$gi});
}
close File;
close OUT;
