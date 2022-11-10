library("DESeq2")

#GET HTSEQ COUNTS AND SET UP SAMPLE TABLE
samples <- readLines("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/htseq_JA20122_counts.txt")
metadata <- read.csv("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/ordered_metadata_JA20122.csv")
sampleTable<-data.frame(sampleName=metadata$FileName, fileName=samples, treatment=metadata$Treatment, sex=metadata$Sex, timepoint=metadata$Time.Point)

# make time + treatment column in metadata
metadata$Timepoint_Treatment = paste(metadata$Treatment, metadata$Time.Point, metadata$Sample.ID, sep='_')

#BUILD A DESEQ2 OBJECT FROM THE HTSEQ COUNTS DATA
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/htseq/htseq_JA20122", design=~treatment + timepoint + treatment:timepoint)
# design will be annotation table, 
coldata <- data.frame(row.names = metadata$Timepoint_Treatment, Time = metadata$Time.Point, Treatment = metadata$Treatment)
# colData(ddsHTSeq)$treatment<-factor(colData(ddsHTSeq)$treatment, levels=c("Control", "Alcohol"))

#LOOK AT THE DESEQ2 OBJECT WE'VE CREATED BY READING IN HTSEQ COUNT FILES
ddsHTSeq

# changing design
design(dds) <- formula(~ timepoint + treatment + timepoint:treatment)

reduced = formula(~ timepoint + treatment)

#Optionally, if you want to do variance stabilizing transformation or regularized log transformation on this count matrix and t$
vsd<- vst(ddsHTSeq)
#write.csv(assay(vsd), file="/stor/work/FRI-BigDataBio/ms_alcohol/Melamed_tagseq/deseq/JA20122_genecounts_variancestabilized.csv")

# PCA analysis
#pca <- plotPCA(vsd, intgroup = c("treatment", "sex"))

# RUN STATISTICAL TEST
dds <- DESeq(ddsHTSeq, test="LRT", reduced=~ timepoint + treatment)
res <- results(dds)
res<-res[order(res$padj),]
mcols(res,use.names=TRUE)
summary(res)

# write.csv(as.data.frame(res), file='deseq2_htseq_JA20122_all.csv')

# get genes w/ padj < 0.05 and 0.1
resdf <- as.data.frame(res)
cutoff1_JA20122 <- na.omit(resdf[resdf$padj < 0.05,])
cutoff2_JA20122 <- na.omit(resdf[resdf$padj < 0.1,])
write.csv(as.data.frame(cutoff1_JA20122), file='deseq2_htseq_JA20122_cutoff1.csv')
write.csv(as.data.frame(cutoff2_JA20122), file='deseq2_htseq_JA20122_cutoff2.csv')
