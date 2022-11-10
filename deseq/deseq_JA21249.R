library("DESeq2")
library(ggplot2)

#GET HTSEQ COUNTS AND SET UP SAMPLE TABLE
samples <- readLines("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/htseq_JA21249_counts.txt")
metadata <- read.csv("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/ordered_metadata_JA21249.csv")
sampleTable<-data.frame(sampleName=metadata$FileName, fileName=samples, treatment=metadata$Treatment, sex=metadata$Sex, timepoint=metadata$Time.Point)
sampleTableNoOutlier <- sampleTable[-c(31), ]
write.csv(sampleTableNoOutlier, file = "/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/annotationsTableNoOutlier.csv")

#BUILD A DESEQ2 OBJECT FROM THE HTSEQ COUNTS DATA
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/htseq/htseq_JA21249", design=~treatment)
# design will be annotation table, 
colData(ddsHTSeq)$treatment<-factor(colData(ddsHTSeq)$treatment, levels=c("Control", "Alcohol"))

#ddsHTSeq No Outlier
ddsHTSeqNoOutlier<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTableNoOutlier, directory="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/htseq/htseq_JA21249", design=~treatment)
# design will be annotation table, 
colData(ddsHTSeqNoOutlier)$treatment<-factor(colData(ddsHTSeq)$treatment, levels=c("Control", "Alcohol"))
colData(ddsHTSeqNoOutlier)


#LOOK AT THE DESEQ2 OBJECT WE'VE CREATED BY READING IN HTSEQ COUNT FILES
ddsHTSeq
ddsHTSeqNoOutlier
#assay(ddsHTSeqNoOutlier)

#Optionally, if you want to do variance stabilizing transformation or regularized log transformation on this count matrix and t$
vsd<- vst(ddsHTSeq)
write.csv(assay(vsd), file = "/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/deseq/JA21249_genercounts_variancestabilized.csv")

#vsd no outlier
vsdNoOutlier <- vst(ddsHTSeqNoOutlier)
write.csv(assay(vsd), file = "/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/deseq/JA21249_nooutlier_genercounts_variancestabilized.csv")

# PCA analysis
pca <- plotPCA(vsd, intgroup = c("treatment", "sex"))
#prcomp 

#pca no outlier 
png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/outlier_removed/pca_no_treatment_timepoint.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = c("timepoint", "treatment")) + ggtitle('Treatment, Timepoint')
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/outlier_removed/pca_no_treatment.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = "treatment") + ggtitle("Treatment")
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/outlier_removed/pca_no_timepoint.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = "timepoint") + ggtitle("Timepoint")
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/outlier_removed/pca_no_sex.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = "sex") + ggtitle("Sex")
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/outlier_removed/pca_no_treatment_sex.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = c("sex", "treatment")) + ggtitle("Sex, Treatment")
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/outlier_removed/pca_no_sex_timepoint.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = c("timepoint", "sex")) + ggtitle("Sex, Timepoint")
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/outlier_removed/pca_no_sex_timepoint_treatment.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = c("timepoint", "sex", "treatment")) + ggtitle("Sex, Timepoint, Treatment")
dev.off()

#top 1000
png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/top1000_outlier_removed/pca_1000_no_treatment_timepoint.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = c("timepoint", "treatment"), ntop = 1000) + ggtitle('Treatment, Timepoint')
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/top1000_outlier_removed/pca_1000_no_treatment.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = "treatment", ntop = 1000) + ggtitle('Treatment')
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/top1000_outlier_removed/pca_1000_no_timepoint.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = "timepoint", ntop = 1000) + ggtitle('Timepoint')
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/top1000_outlier_removed/pca_1000_no_sex.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = "sex", ntop = 1000) + ggtitle('Sex')
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/top1000_outlier_removed/pca_1000_no_treatment_sex.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = c("sex", "treatment"), ntop = 1000) + ggtitle('Sex, Treatment')
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/top1000_outlier_removed/pca_1000_no_sex_timepoint.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = c("sex", "timepoint"), ntop = 1000) + ggtitle('Sex, Timepoint')
dev.off()

png(file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA21249/top1000_outlier_removed/pca_1000_no_sex_timepoint_treatment.png",
    width=500, height=500)
plotPCA(vsdNoOutlier, intgroup = c("timepoint", "sex", "treatment"), ntop = 1000) + ggtitle('Sex, Timepoint, Treatment')
dev.off()



ploted<-ggplot(assay(vsdNoOutlier), aes(x = PC1, y = PC2, color = factor(timepoint), shape = factor(treatment))) +
  geom_point(size =3, aes(fill=factor(mouseline))) +
  geom_point(size =3) +
  scale_shape_manual(values=c(21,22)) +
  scale_alpha_manual(values=c("F"=0, "M"=1)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of all genes, no covariate adjusted")
ploted

