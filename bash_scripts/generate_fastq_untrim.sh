cd /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files

for d in *
do 
  
  cd $d"/cat_"$d"/untrim_"$d
    # go into the cat directory
    # for each file in the cat directory run fastqc
    for f in *.fastq
        do
            echo "fastqc /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files/"$d"/cat_"$d"/untrim_"$d"/"$f
            echo "rm /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files/"$d"/cat_"$d"/untrim_"$d"/"$f"_fastqc.html"
            echo "mv /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastq_files/"$d"/cat_"$d"/untrim_"$d"/"$f"_fastqc.zip /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastqc_untrim/"$d"_fastqc"
        done 
  cd ..
  cd ..
  cd ..
done > /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/bash_scripts/fastqc_untrim.sh