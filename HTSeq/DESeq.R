# Analysis of Schistosoma mansoni bulk RNAseq

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install('DESeq2')
BiocManager::install('PCAtools')
BiocManager::install('pheatmap')
BiocManager::install('RColorBrewer')
BiocManager::install('gplots')
BiocManager::install('ggplot2')

library(DESeq2)
library(PCAtools)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gplots)

directory <- ('/Users/skb/Desktop/Smansoni_rnaseq/HTSeq')
sampleFiles <- grep ("htseq.o", list.files(directory),value=TRUE)
sampleFiles
sampleCondition <- gsub("_R.htseq.o","",sampleFiles)
sampleCondition

#sampleTable <- data.frame(sampleName=sampleFiles,
                          fileName=sampleFiles,
                          condition=sampleCondition)

#list.files('/Users/skb/Desktop/Smansoni_rnaseq/')
