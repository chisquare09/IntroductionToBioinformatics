library("DeSeq2")
library(ggplot2)

countData <- read.csv('pasilla-countData.csv',header=TRUE,sep = ",")
head(countData)

metaData <- read.csv('pasilla-colData.csv',header=TRUE,sep = ",")
head(metaData)

dds <- DeSeqDataSetFromMatrix(countData=countData,
                              colData = metaData,
                              design=~dex,
                              tidy = TRUE)
dds <- DESeq(dds)

res <- results(dds)
head(results(dds,tidy=TRUE))
summary(res)

par(mfrow=c(2,3))
