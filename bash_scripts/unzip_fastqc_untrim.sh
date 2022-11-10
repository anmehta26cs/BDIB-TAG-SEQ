cd /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/fastqc_untrim

for d in *
do 
  cd $d
  for di in *
    do
        cd $di
        for f in *
        do
        unzip $f
        done
        cd ..
    done
done