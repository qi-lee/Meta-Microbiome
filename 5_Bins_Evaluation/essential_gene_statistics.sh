#!/bin/bash
#############################################
# AUTHOR:liqi , liqi@ihb.ac.cn
# COPYRIGHT:Copyright - Lab of Algal Genomics
#############################################


groupnum=`ls | grep "group_.*.fa$" | wc -l`

echo "group_id    size    essgene_num    contig_snum    orf_num"

for (( i = 0 ; $i < $groupnum; i++ ))
do
    grep ">" group_$i.fa > group_$i.list
    contigsnum=`grep ">" group_$i.fa | wc -l`
    ../Script/found_essential_gene_by_group.pl hmm_orf_alignment.txt group_$i.list > group_$i.essential_gene
    essgenenum=`cut -f 2 group_$i.essential_gene | sort | uniq -c | wc -l`
    prodigal -a group_$i.prodigal_orf.faa -i group_$i.fa -m -o group_$i.prodigal_stdout -p meta -q
    cut -f1 -d " " group_$i.prodigal_orf.faa > group_$i.orf.faa
    orfnum=`grep ">" group_$i.orf.faa | wc -l`
    echo "group_$i    `ls -lh group_$i.fa | awk '{print $5}'`    $essgenenum    $contigsnum    $orfnum"
    rm group_$i.prodigal_orf.faa
    rm group_$i.prodigal_stdout
    rm group_$i.list
done
