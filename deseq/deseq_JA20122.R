library("DESeq2")

#GET HTSEQ COUNTS AND SET UP SAMPLE TABLE
samples <- readLines("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/htseq_JA20122_counts.txt")
metadata <- read.csv("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/ordered_metadata_JA20122.csv")
sampleTable<-data.frame(sampleName=metadata$FileName, fileName=samples, treatment=metadata$Treatment, sex=metadata$Sex, timepoint=metadata$Time.Point)
write.csv(sampleTable, file = "/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/visualizations/JA20122/annotationsTable.csv")

#BUILD A DESEQ2 OBJECT FROM THE HTSEQ COUNTS DATA
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/htseq/htseq_JA20122", design=~treatment)
# design will be annotation table, 
colData(ddsHTSeq)$treatment<-factor(colData(ddsHTSeq)$treatment, levels=c("Control", "Alcohol"))


#LOOK AT THE DESEQ2 OBJECT WE'VE CREATED BY READING IN HTSEQ COUNT FILES
ddsHTSeq

#Optionally, if you want to do variance stabilizing transformation or regularized log transformation on this count matrix and t$
vsd<- vst(ddsHTSeq)
write.csv(assay(vsd), file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/deseq/JA20122_genecounts_variancestabilized.csv")

# PCA analysis
timepoint_treatment <- plotPCA(vsd, intgroup = c("timepoint", "treatment"))
treatment <- plotPCA(vsd, intgroup = "treatment")
timepoint <- plotPCA(vsd, intgroup = "timepoint")
sex <- plotPCA(vsd, intgroup = "sex")
sex_treatment <- plotPCA(vsd, intgroup = c("sex", "treatment"))
sex_timepoint <- plotPCA(vsd, intgroup = c("timepoint", "sex"))
sex_timepoint_treatment <- plotPCA(vsd, intgroup = c("timepoint", "sex", "treatment"))

#top 1000
pca_1000_no_treatment_timepoint <- plotPCA(vsd, intgroup = c("timepoint", "treatment"), ntop = 1000)
pca_1000_no_treatment <- plotPCA(vsd, intgroup = "treatment", ntop = 1000)
pca_1000_no_timepoint <- plotPCA(vsd, intgroup = "timepoint", ntop = 1000)
pca_1000_no_sex <- plotPCA(vsd, intgroup = "sex", ntop = 1000)
pca_1000_no_treatment_sex <- plotPCA(vsd, intgroup = c("sex", "treatment"), ntop = 1000)
pca_1000_no_sex_timepoint <- plotPCA(vsd, intgroup = c("timepoint", "sex"), ntop = 1000)
pca_1000_no_sex_timepoint_treatment <- plotPCA(vsd, intgroup = c("timepoint", "sex", "treatment"), ntop = 1000)

# ggplot pca
library(ggplot2)
pcaData <- plotPCA(vsd, intgroup = c("timepoint", "treatment"), returnData=TRUE, ntop = 1000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=timepoint, shape=treatment)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()

