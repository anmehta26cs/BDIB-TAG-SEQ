#!/bin/bash

# Go into spleen file (JA201.../FASTQ...) and go into each directory and unzip .fastq.gz file with gzip
cd /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/JA20122-167546383/FASTQ_Generation_2020-05-07_09_01_50Z-245884639
for d in *
do
    cd $d
    gzip -d *.gz
    cd ..
done

# Go into spine file (JA212.../FASTQ...) and go into each directory and unzip .fastq.gz file with gzip
cd /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/JA21249-274011746/FASTQ_Generation_2021-06-22_10_01_17Z-430904477
for d in *
do
    cd $d
    gzip -d *.gz
    cd ..
done
