# concatenate L002 onto L001 for every fastq file

cd /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files

for d in * # JA20122 and JA21249
do
    cd $d
    for file in *.fastq # the fastq file (both L001 and L002)
    do
        echo $file > files.txt
    done

    cd ..
done 