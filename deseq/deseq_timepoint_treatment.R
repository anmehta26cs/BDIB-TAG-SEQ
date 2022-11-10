library("DESeq2")

#GET HTSEQ COUNTS AND SET UP SAMPLE TABLE
directory<-(getwd())
samples0 <- readLines("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/htseq_JA20122_counts.txt")
samples1 <- readLines("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/htseq_JA21249_counts.txt")
metadata0 <- read.csv("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/ordered_metadata_JA20122.csv")
metadata1 <- read.csv("/stor/work/FRI-BigDataBio/ms_alcohol/metadata/ordered_metadata_JA21249.csv")
sampleTable<-data.frame(sampleName=metadata$FileName, fileName=samples, treatment=metadata$Treatment, sex=metadata$Sex, timepoint=metadata$Time.Point)

# ADD TIME/TREAT COLUMN
metadata0$Timepoint_Treatment = paste(metadata0$Treatment, metadata0$Time.Point, metadata0$Sample.ID, sep='_')
metadata1$Timepoint_Treatment = paste(metadata1$Treatment, metadata1$Time.Point, metadata1$Sample.ID, sep='_')

#BUILD A DESEQ2 OBJECT FROM THE HTSEQ COUNTS DATA
ddsHTSeq0<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable0, directory=directory, design=~Timepoint_Treatment + Sex)
ddsHTSeq1<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable1, directory=directory, design=~Timepoint_Treatment + Sex)

#RUN THE STATISTICAL TEST IN ONE GO- NORMALIZATION, ESTIMATE DISPERSION/VARIANCE AND DO TEST FOR DIFFERENTIAL EXPRESSION
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
mcols(res,use.names=TRUE)
summary(res)

#GENERATE MA PLOT
png('MAPlot_htseq.png')
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.off()

#WRITE RESULTS INTO FILE
write.csv(as.data.frame(res),file="deseq2_htseq_C1_vs_C2.csv")