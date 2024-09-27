#loading required packages
library(tximport)
library(tidyverse)
library(tximportData)
library(rhdf5)
library(DESeq2)
library(pheatmap)
#moved files to system file for tximport for ease of use
path <- system.file("quant_tximport", package = "tximportData")
runs <- read_delim(file.path(path, "SraRunTable.txt"))

#Obtain only runs for which we ran data 
runs<-filter(runs, Run %in%  c("SRR191917", "SRR191918", "SRR191867", "SRR191922", "SRR191385", "SRR191388", "SRR191390", "SRR191861"))
#obtaining path for where salmon files are located
files <- file.path(path, "salmon", runs$Run, "quant.sf")
names(files) <- paste0(c("SRR191917", "SRR191918", "SRR191867", "SRR191922", "SRR191385", "SRR191388", "SRR191390", "SRR191861"))

#converting transcript ids to gene symbols
#calling file with annotation inforamtion
tx2gene<-read_tsv(file.choose())
#only need trancript id and symbols
tx2gene<-tx2gene%>%select(product_accession,symbol )

##importing files from salmon
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

##importing files from salmon
files_kallisto <- file.path(path, "kallisto", runs$Run, "abundance.h5")
names(files_kallisto) <- paste0(c("SRR191917", "SRR191918", "SRR191867", "SRR191922", "SRR191385", "SRR191388", "SRR191390", "SRR191861"))
txi.kallisto <- tximport(files_kallisto, type = "kallisto", tx2gene = tx2gene)

##DESeq2
#unable to call the names of the samples with removing the space between column header, renaming it from sample name to "sample_name" 
names(runs)[24]<-"sample_name"
##Creating a DESeqdata set for analysis of gene expression for salmon
dds.salmon <- DESeqDataSetFromTximport(txi.salmon, colData = runs , design = ~ sample_name)
dds.S<-DESeq(dds.salmon)

###salmon heat map
select.s <- order(rowMeans(counts(dds.S, normalized=TRUE)),
                decreasing=TRUE)[1:20]
heatmap.S <- as.data.frame(colData(dds.salmon)[,c("sample_name", "LibraryLayout")])
ntd.s <- normTransform(dds.S)
pheatmap(assay(ntd.s)[select.s,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=heatmap.S)
##Creating a DESeqdata set for analysis of gene expression for kallisto
dds.kallisto<-DESeqDataSetFromTximport(txi.kallisto, colData = runs , design = ~ sample_name)
dds.K<-DESeq(dds.kallisto)
###kallisto heat map 
select.k <- order(rowMeans(counts(dds.K, normalized=TRUE)),
                  decreasing=TRUE)[1:20]
heatmap.K <- as.data.frame(colData(dds.kallisto)[,c("sample_name", "LibraryLayout")])
ntd.k <- normTransform(dds.K)
pheatmap(assay(ntd.k)[select.k,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=heatmap.K)
