#!/bin/bash

#############################################
# AUTHOR:liqi , liqi@ihb.ac.cn
# COPYRIGHT:Copyright - Lab of Algal Genomics
#############################################


ARRAY=($(cut -b 1-4 $1 | sort | uniq))

for i in ${ARRAY[*]}
do
#echo ${i}
total=`grep "${i}" $1 | wc -l`
echo -e "$total\t${i}"
grep "${i}" $1 |cut -f 6 | sort | uniq -c | sort -k1 -nr | perl -ane 'print "$F[0]\t$F[1]\n";exit if ($.>3)'

echo "";

done
