Appendix:
#To obtain data from SRA, you need to download these modules
Zuhaa Ali
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit
#obtain fasta files required for this analysis (focused on short term and long term elevated
temperature exposure
fasterq-dump --fasta SRR19197 SRR191918 SRR191867 SRR191922 SRR191385 SRR191388 SRR191390
SRR191861
#make new directory to move fasta files to
mkdir fasta
mv *.fasta fasta
#obtaining transcriptome data
wget
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/753/865/GCF_013753865.1_Amil_v2.1/GCF_013753
865.1_Amil_v2.1_rna.fna.gz
#transcriptome must first be indexed prior to using it
nano index.sh
creating sbatch script for this
#!/bin/sh
#SBATCH --account=def-lukens
#SBATCH --time=0-00:15:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)
module load salmon
salmon index -t GCF_013753865.1_Amil_v2.1_rna.fna.gz -i indexed_genome
#quantify samples
nano quantify.sh
#script for quantifying samples
#!/bin/sh
#SBATCH --account=def-lukens
#SBATCH --time=0-03:00:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)
#!/bin/bash

# Load the Salmon module
module load salmon

# Set the index and options
index="acropora_index"
lib_type="A"
validate_mappings="--validateMappings"
output_dir="./quant_results"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# List of FASTA files and corresponding output names
declare -A samples=(
    [SRR191867]="/scratch/zee/Project4/fasta/SRR191867.fasta"
    [SRR191918]="/scratch/zee/Project4/fasta/SRR191918.fasta"
    [SRR191917]="/scratch/zee/Project4/fasta/SRR191917.fasta"
    [SRR191922]="/scratch/zee/Project4/fasta/SRR191922.fasta"
    [SRR191385]="/scratch/zee/Project4/fasta/SRR191385.fasta"
    [SRR191388]="/scratch/zee/Project4/fasta/SRR191388.fasta"
    [SRR191390]="/scratch/zee/Project4/fasta/SRR191390.fasta"
    [SRR191861]="/scratch/zee/Project4/fasta/SRR191861.fasta"
)

# Run Salmon quant for each sample
for sample in "${!samples[@]}"; do
    fasta_file="${samples[$sample]}"
    output_dir_path="$output_dir/$sample"
    salmon quant -i "$index" -l "$lib_type" -r "$fasta_file" $validate_mappings -o "$output_dir_path"
done

#script for building index
#!/bin/sh
#SBATCH --account=def-lukens
#SBATCH --time=0-00:15:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
Zuhaa Ali
#SBATCH --mem=8000 # requested memory (in MB)
module load kallisto
kallisto index -i transcripts.idx GCF_013753865.1_Amil_v2.1_rna.fna.gz
#quantifying reads using kallisto
#!/bin/sh
#SBATCH --account=def-lukens
#SBATCH --time=0-03:00:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)
#!/bin/bash

# Load the Kallisto module
module load kallisto

# Set parameters
index="kallisto_index.idx"
output_dir="./kallisto_quant"
num_bootstraps=100
lib_type="--single"
frag_length=150
std_dev=20

# Create output directory if it doesn't exist
mkdir -p $output_dir

# List of FASTA files and corresponding output names
declare -A samples=(
    [SRR191867]="/scratch/zee/Project4/fasta/SRR191867.fasta"
    [SRR191917]="/scratch/zee/Project4/fasta/SRR191917.fasta"
    [SRR191918]="/scratch/zee/Project4/fasta/SRR191918.fasta"
    [SRR191922]="/scratch/zee/Project4/fasta/SRR191922.fasta"
    [SRR191385]="/scratch/zee/Project4/fasta/SRR191385.fasta"
    [SRR191388]="/scratch/zee/Project4/fasta/SRR191388.fasta"
    [SRR191390]="/scratch/zee/Project4/fasta/SRR191390.fasta"
    [SRR191861]="/scratch/zee/Project4/fasta/SRR191861.fasta"
)

# Run Kallisto quant for each sample
for sample in "${!samples[@]}"; do
    fasta_file="${samples[$sample]}"
    kallisto quant -i "$index" -o "$output_dir/$sample" -b "$num_bootstraps" "$lib_type" -l "$frag_length" -s "$std_dev" "$fasta_file"
done

# Move files into the output directory to keep files organized
mv *_kallisto "$output_dir"

#In R
#loading required packages
#BiocManager::install("tximport")
library(tximport)
library(tidyverse)
#BiocManager::install("tximportData")
#BiocManager::install("rhdf5")
BiocManager::install("DESeq2")
#BiocManager::install("rhdf5")
library(tximportData)
library(rhdf5)
library(DESeq2)
#install.packages("pheatmap")
library(pheatmap)
#moved files to system file for tximport for ease of use
path <- system.file("quant_tximport", package = "tximportData")
runs <- read_delim(file.path(path, "SraRunTable.txt"))
#Obtain only runs for which we ran data i.e(SRR191917, SRR191918, SRR191867, SRR191922,
SRR191385, SRR191388, SRR191390, SRR191861)
runs<-filter(runs, Run %in% c("SRR191917", "SRR191918", "SRR191867", "SRR191922", "SRR191385",
"SRR191388", "SRR191390", "SRR191861"))
#obtaining path for where salmon files are located
files <- file.path(path, "salmon", runs$Run, "quant.sf")
Zuhaa Ali
names(files) <- paste0(c("SRR191917", "SRR191918", "SRR191867", "SRR191922", "SRR191385",
"SRR191388", "SRR191390", "SRR191861"))
#converting transcript ids to gene symbols
#calling file with annotation information
tx2gene<-read_tsv(file.choose())
#only need trancript id and gene symbols
tx2gene<-tx2gene%>%select(product_accession,symbol )
##importing files from salmon
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
##importin files from Kallisto
files_kallisto <- file.path(path, "kallisto", runs$Run, "abundance.h5")
names(files_kallisto) <- paste0(c("SRR191917", "SRR191918", "SRR191867", "SRR191922", "SRR191385",
"SRR191388", "SRR191390", "SRR191861"))
txi.kallisto <- tximport(files_kallisto, type = "kallisto", tx2gene = tx2gene)
#DESeq2
#unable to call the names of the samples with removing the space between column header, renaming
it from sample name to “sample_name”
names(runs)[24]<-"sample_name"
##Salmon
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
##Kallisto
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
