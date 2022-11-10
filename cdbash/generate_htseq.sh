cd /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/star

for d in *
do 

    cd $d 
        for i in *.bam
            do 
                echo "htseq-count -f bam -m intersection-nonempty -i gene_id /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/star/$d/$i /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/map_to_genome/reference_genome/mouse_GRCm38/gencode.vM19.annotation.gtf > /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/htseq/${i%A*}.counts"
            done
    cd ..
done > /stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/cdbash/htseq.sh