#!/bin/bash

for f in `ls data/sequencing/alignment/STAR/*_amgm350_Log.final.out`
do
    sample_name=$(echo $f | sed 's/_amgm350.*//' | sed 's/.*\///')
    awk -v name=$sample_name 'NR==6{print name"\t"$6}' $f
    done > data/sequencing/library_size.txt
