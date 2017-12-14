#!/usr/bin/perl
use strict;

###############################################
##                ReadMe
###############################################
## Function: 
## Usage $: found_essential_gene_by_group.pl hmm_orf_alignment.txt group_1.list 
## Author: liqi (lab of algal genomics)
## Mail: liqi@ihb.ac.cn
## Created Time: Wed 09 Sep 2015 04:10:53 PM CST
###############################################


open(Lib,$ARGV[0]);
open(List,$ARGV[1]);
my @a=<Lib>;
close Lib;

while(<List>){
    (my $name)=$_=~/^>(\S+)/;
    foreach(@a){
         my @b=split;
         next if /^#/;
         (my $name1)=$b[0]=~/(\S+)_\d+$/;
         if($name eq $name1){
             print "$name\t$b[2]\n";
         }
    }
}

close Lib;
close List;
