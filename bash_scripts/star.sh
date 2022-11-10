# STAR --runThreadN 1 --outFileNamePrefix ${f%_*} --outSAMstrandField intronMotif --genomeDir /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/reference_genome/STAR_genome --outSAMtype BAM Unsorted --readFilesIn $f

path="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/trimmed_data/JA21249/"
for d in $path
do 
  cd $d
  for f in *.fastq.trim
    do
      echo "STAR --runThreadN 1 --outFileNamePrefix ${f%_*} --outSAMstrandField intronMotif --genomeDir /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/reference_genome/STAR_genome --outSAMtype BAM Unsorted --readFilesIn $path$f"
    done
  cd ..
done > /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/star/commands-JA21249.star