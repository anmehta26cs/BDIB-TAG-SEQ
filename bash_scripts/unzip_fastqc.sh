cd /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastqc

for d in *
do 
  cd $d
  for f in *.zip
      do
        rm $f
      done 
  cd ..
done