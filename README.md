# Salmon-vs-Kallisto
## Introduction

Salmon and Kallisto are two next generation methods for quantifying transcript abundance (Patro et al.,
2017). This process is foundational to genetics and biology as it helps classify diseases as well as monitor
their progression (Patro et al., 2017). There are multiple studies that have been conducted using these
pipelines which include classification of diseases and their subtypes such as what is done by Hoadley et
al., 2014, “Multiplatform Analysis of 12 Cancer Types Reveals Molecular Classification within and across
Tissues of Origin” , understanding developmental change in gene expression by by Li et al., 2014 done
through “Comparison of D. melanogaster and C. elegans developmental stages, tissues, and cells by
modENCODE RNA-seq data” and finally tracking the progression of cancer such as what is done by Shaw
et al., 2013 through “The Cancer Genome Atlas Pan-Cancer analysis project” . Both salmon and kallisto
have been advantageous in bioinformatics realm as they are pseudo-aligners which means that do not
require alignment to an actual genome which in the past was both time consuming and costly creating a
bottle neck for users to perform these analyses.
In this tutorial we will use both salmon and kallisto to determine which software is more user friendly
and efficient (Patro et al., 2017) as well as determine whether similar results as the article can be
achieved using these workflows. Data was taken from the article “Profiling gene expression responses of
coral larvae (Acropora millepora) to elevated temperature and settlement inducers using a novel RNASeq procedure can be achieved” (Meyer et al., 2011). We will document the workflow, challenges, and
downstream analysis using DESeq2 and pheatmap.

## Methods
### Dataset details
This dataset was obtained through NCBI with the SRA accession code SRA029780. Although multiple
variables were tested in the original article, for our analysis we focused on larvae exposed to short term
and long-term temperature elevation, we only took data library 1, using 35 base pair reads, with 2 reads
each and their control. We used the sra-toolkit module to obtain this data, along with run information
the following code was used to obtain these samples
```
fasterq-dump --fasta <input run names>
```
We have included the exact run names for each sample used in the appendix section. We also obtained
the transcriptome (and annotation file) from ncbi(link in appendix) which was available in a .fna format
whereas the experimental data was downloaded as a .fasta format. The transcriptome was obtained
using wget, along with the ftp link which can be found in the appendix.

First the transcriptome was indexed using both salmon and kallisto this allows for quick and efficient
retrieval of information instead of going through the whole transcriptome to match samples (note: all
analyses were done through sbatch to allow for sufficient computational power). In salmon we can use
following code to index our genome
```
salmon index -t <transcriptomename.fa.gz> -i <directoryofindexedgenome>
```
Whereas in kallisto you have to use
```
kallisto index -I <nameofindexedtranscrptome.idx> <transcriptomename.fa.gz>
```
Next, we quantified our samples, for both salmon and kallisto we used the base default settings to
prevent us from over complicating things. In salmon we used the following code. The -i option points
salmon to the indexed transcriptome of the data. The -l A option tells salmon to automatically
determine which library is being used. Since this was a single instead of paired read end library, we used
the -r option, however for paired reads you can use -1 and -2 options. Finally, the -o option tells salmon
where to output the file from this command.
```
salmon quant -i <indexed transcriptome> -l A -r <samples> --validateMappings -o <output>
```
For kallisto quantification we used the following code:
```
kallisto quant -i <indexed kallisto transcriptome.idx>-o <output directory> -b 100 --single -l <fragment
size> -s <fragment standard deviation> <sample>
```
The -i option tells kallisto where to find the indexed transcriptome, the -o option tells kallisto what the
name of the output directory should be. The -b option indicates the number of bootstrap samples, we
set this to 100, however this should be dependant on your analysis. The –single -l -s option specifies
single read ends; with whatever fragment size and standard deviation your fragment was conducted
with. For our analysis we chose arbitrary values based on what literature recommends as fragment size
and standard deviation were not available with supplemental information.

Next these files were taken imported into R with the help of tximport so that correct formatting for
DESeq2 was attained. The files obtained from kallisto and salmon were moved to the tximportData file
in R for ease of use. Run information was first obtained for the variables used in our analysis
(SRR191917", "SRR191918", "SRR191867", "SRR191922", "SRR191385", "SRR191388", "SRR191390",
"SRR191861). Then transcript IDs were converted to gene symbols to quantify gene level expression
rather than transcript level. Finally, a DESeqDataSet was generated which was compatible with DESeq2.
The format can be seen below. Once this was complete DESeq2 was used to quantify differential gene
expression and lastly pheatmap was used to generate a heatmap. The code for this portion of the
analysis is quite large and so it can be seen in the appendix with well documented comments.
### Format of files obtained from salmon and Kallisto:
Salmon

![image](https://github.com/user-attachments/assets/115ead8e-085e-4f31-bb48-599c54daeae8)

Kallisto
![image](https://github.com/user-attachments/assets/ff030353-6d73-4662-9aba-39586a09e463)

Format required for DESeq2
![image](https://github.com/user-attachments/assets/408990d5-e422-4623-a355-a8be0ab21c56)

## Results/Discussion:
The purpose of this analysis was to determine whether one RNA quantifier was more user friendly and
efficient in comparison to the other and whether results similar to the article could be achieved using
these workflows. It was predicted that short term temperature elevation would lead to the upregulation
of heat shock proteins whereas long term temperature elevation would cause a downregulation of these
proteins and an increase in genes controlling ion transport of Ca2+ and carbonate ions (Meyer et al.,
2011). However, once the workflow was developed and results generated, we determined that it was
not possible to achieve similar results as the libraries were in their raw format and had not followed any
QC steps. Due to time constraint, it was not possible to generate a tutorial on this. In addition, we did
not normalize any genes like the article as this step in the vignette was listed under optional commands.
So, this likely led to further dissimilarities. The gene symbols along with their names can be observed in
table 1 and 2. For next steps we can incorporate a tutorial on quality control, as well as try this workflow
with normalization according to the article to see whether similar results are achieved.
We were able to generate a heat map for both salmon (figure 1) and kallisto (figure 2). This indicates
that our analysis did run. In term of user friendliness and efficiency we found salmon be better than
kallisto. Salmons’ options and user manual is relatively straight forward, each option clearly
demonstrated its function with no room for ambiguity whereas the options for kallisto were more
ambiguous such as fragment length which could be confused with read length. Kallisto also requires
more (5) options to run a default base quantification whereas salmon only requires three options
(indexed transcriptome, library type and type of read ends (single vs paired)). Indexing of the
transcriptome and quantification for salmon took less than 10 minutes each whereas kallisto took
relatively longer. So overall it is favourable for beginners to use salmon rather than kallisto

![image](https://github.com/user-attachments/assets/9d139db4-7253-40d3-9aa9-285f59a4a7df)


Figure 1.
Fig1: A heat map of differential gene expression in Acropora millepora larvae (n=20-30) under 4
conditions, long term control, long term elevated temperature, short term control and short-term
elevated temperature analysed using the software salmon.
Table 1. Gene names of differentially expressed genes using the software salmon.
![image](https://github.com/user-attachments/assets/ff21cf66-52cd-48cc-9e2e-fe2285f57edf)
![image](https://github.com/user-attachments/assets/69a8fdd0-e50b-4872-9ee1-bb96d1786341)


![image](https://github.com/user-attachments/assets/93e48279-9966-4dd8-bcf1-01d39a34db5a)

Figure 2.
Fig2: A heat map of differential gene expression in Acropora millepora larvae (n=20-30) under 4
conditions, long term control, long term elevated temperature, short term control and short-term
elevated temperature analysed using the software kallisto.
Table 2. Gene names differentially expressed genes using the software kallisto.
![image](https://github.com/user-attachments/assets/bc2495a4-ae0b-4596-9d7a-b2b2a10e100b)

### Reflection:
The overall experience using the workflow was quite pleasant. Both RNA-seq quantifiers have well
documented, accessible, and in-depth manuals. In general, we found salmon to be more user friendly as
all the options were clearly labelled leaving no room for ambiguity. From our experience salmon was
also quicker than kallisto. In this tutorial we mainly covered the basics with default settings which can be
changed depending on the type of analysis you are conducting. One of the more challenging experiences
was getting kallisto to work as it required both standard deviation and the size of fragment. Since this
data was taken from sequence read archive this information could not be obtained as it was not
provided in the supplemental information. From research, we believe that you can get this information
from the bioanalyzer trace. We estimated these values based on the most common settings people use
when this data is not available. Hence, we do not know for certain whether the difference seen between
the two softwares is due to wrongly assigning values or whether a true difference lies between them.
Due to time constraint, we could not follow up with a different data set, however this should certainly
be considered when using kallisto. The most challenging part of this experience was the downstream
analysis using DESeq2 to work with the files obtained from salmon and kallisto. For these files to work
with DESeq2 you must import them correctly using tximport and although this is thought to be a
relatively user-friendly software, we had difficulties trying to understand the manual. More specifically it
was difficult to determine what each line of the code was doing so; we ran the code with dataset they
had provided and then made changes according to the needs of our analysis. We found this to be much
easier rather than customizing the code beforehand.
This analysis was chosenbecause all through undergrad I have heard about RNA-seq quantification as
well as heard about the lab-based aspect of it however no one has taught the computational side of it. It
was interesting to see that you don’t obtain results right after you have the data but that it goes through
quite a bit of steps before generating a figure that is needed to display those results. In addition, I also
chose this analysis because it is required for the final project in BINF *6999 where I will be using single
cell ATAC seq data obtained from maize to determine the genes responsible for flowering time.
Overall, I believe the workflow is simple and well documented enough for any beginner to follow. For
best practices we would recommend running the code as is, so you get a better understanding of what
each line of code is doing and then customize it according top your needs. We believe that kallistos
fragment size and standard deviation is likely to cause some problems as some users will try to use


downloaded data which may or may not have this information. In this case we recommend the settings
used above. In addition, we believe that tximport for importing data may also cause some issues,
however these can be resolved by running our code as is and then customizing it later. Lastly, we believe
this workflow can be improved by including a tutorial on QC steps as they are required for all libraries
generated and running the normalization step to see what impact it has on results.
### Reference:
Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and biasaware quantification of transcript expression. Nature Methods, 14(4), 417–419.
https://doi.org/10.1038/nmeth.4197
Hoadley, K. A., Yau, C., Wolf, D. M., Cherniack, A. D., Tamborero, D., Ng, S., Leiserson, M. D. M., Niu, B.,
McLellan, M. D., Uzunangelov, V., Zhang, J., Kandoth, C., Akbani, R., Shen, H., Omberg, L., Chu, A.,
Margolin, A. A., Van't Veer, L. J., Lopez-Bigas, N., Laird, P. W., … Stuart, J. M. (2014). Multiplatform
analysis of 12 cancer types reveals molecular classification within and across tissues of
origin. Cell, 158(4), 929–944. https://doi.org/10.1016/j.cell.2014.06.049
Li, J. J., Huang, H., Bickel, P. J., & Brenner, S. E. (2014). Comparison of D. melanogaster and C. elegans
developmental stages, tissues, and cells by modENCODE RNA-seq data. Genome research, 24(7),
1086–1101. https://doi.org/10.1101/gr.170100.113
Cancer Genome Atlas Research Network, Weinstein, J. N., Collisson, E. A., Mills, G. B., Shaw, K. R.,
Ozenberger, B. A., Ellrott, K., Shmulevich, I., Sander, C., & Stuart, J. M. (2013). The Cancer Genome
Atlas Pan-Cancer analysis project. Nature genetics, 45(10), 1113–1120.
https://doi.org/10.1038/ng.2764
Meyer, E., Aglyamova, G. V., & Matz, M. V. (2011). Profiling gene expression responses of coral larvae
(Acropora millepora) to elevated temperature and settlement inducers using a novel RNA-Seq
procedure. Molecular ecology, 20(17), 3599–3616. https://doi.org/10.1111/j.1365-
294X.2011.05205.x
Vignettes and Manuals :
