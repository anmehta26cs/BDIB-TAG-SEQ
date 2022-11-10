# a script to create the script to concatenate lines one and two for each fastq file for each sample
cd /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files
for d in *
do
    cd $d
    # Concatenate all the fastq files from the same sample into one file
    for f in *.fastq
    do
        # if file contains L001
        if [[ $f == *L001* ]]
        then
            # using L as delimiter to get the prefix of the filename
            prefix=$(echo $f | cut -d'L' -f1)
            # removing underscore at the end of the prefix
            prefix=${prefix::-1}
            # echo command to concatenate L001 and L002 files for each sample into separate file
            file1="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files/"$d"/"$prefix"_L001_R1_001.fastq"
            file2="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files/"$d"/"$prefix"_L002_R1_001.fastq"
            echo "cat "$file1" "$file2" > /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files/"$d"/cat/"$prefix"_concat.fastq"
        fi
    done     

    cd ..
done > /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/bash_scripts/concat_files.sh