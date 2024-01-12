# Prepare System

## R update and load librarys

``` r
BiocManager::install(update = TRUE, ask = FALSE)

# library("ggVennDiagram")
# library("ggvenn")
# library(AnnotationHub)
# library(dbplyr)
# library(tidyverse)
# library(ChIPQC)
# library(BiocParallel)
# library(readxl)
# library(writexl)
# library(xlsx)
# library(beepr)
# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# library(TxDb.Mmusculus.UCSC.mm39.knownGene)
# # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# # txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
# library(clusterProfiler)
# library(UpSetR)
# library(ggimage)
# library(ReactomePA)
# library(ggcorrplot)
# library(gplots)
# library(pheatmap)
# library(EnsDb.Mmusculus.v79)
# library(Rsamtools)
# library(VennDiagram)
# library(biomaRt)
# library(biobtreeR)
# library(memes)
# library(magrittr)
# library(ggplot2)
# library(GenomicRanges)
# library(Motif2Site)
# library(updater)
# library(viridis)
# library(BiocManager)
# library(ChIPQC)
# library(rtracklayer)
# library(ChIPpeakAnno)
# library(trackViewer)
# library(Gviz)
# library(rtracklayer)
# library(trackViewer)
# library(GenomicFeatures)
# library(UpSetR)
# library(grid)
# library(plyr)
# library(tidyGenomeBrowser)
# library(rentrez)
# library("GenomicFeatures")
# library("Gviz")
# library(GenomicRanges)
# # library()





# run name
# run_name <- "test1"
```

## R folders

``` r
if(Sys.info()[["sysname"]]== "Linux") {
  s <- "/mnt/s"
} else(
  s <- "S:"
)

ddir <- paste(s,"AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF",sep="/")

wdir <- "E:/Simon/P3026_ChIP_Seq_epiSVF"

setwd(ddir)
getwd()
```

    ## [1] "/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF"

``` r
outdir <- paste(s,"AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF/output/",sep="/")
dir.create(outdir)
```

    ## Warning in dir.create(outdir):
    ## '/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF/output' existiert
    ## bereits

## linux

``` bash
##  linux ####
sudo apt update
sudo apt upgrade
conda update conda
conda update --all
sudo mount -t drvfs S: /mnt/s
```

# 0 Prerequisits

## Folder & Names

``` bash
run="Kelly";

dirdata="/mnt/s/AG/AG-Scholz-NGS/Daten/ChIP_epiSVF_P3026";
dir="/mnt/e/Simon/P3026_ChIP_Seq_epiSVF"; 
dirs="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF"; 

index_GRCh38="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/Genomic_data/Human/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as";
index_GRCm39="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/Genomic_data/Mus/GRCm39/GRCm39";
index_mm10="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/Genomic_data/Mus/mm10/mm10";
index_mm39="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/Genomic_data/Mus/GRCm39/USCS/mm39";
index_cm_mm39="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/Genomic_data/Mus/GRCm39/USCS/mm39.index";

fastqdir="$dirdata/fastq";

ardir="$dirdata/Adapterremoval";
  test -d $ardir && echo "# Folder $ardir exists" ||
  (mkdir $ardir && echo "# Folder $ardir created");
  
fastqcdir="$dirdata/fastqc/";
  test -d $fastqcdir && echo "# Folder $fastqcdir exists" ||
  (mkdir $fastqcdir && echo "# Folder $fastqcdir created");
  
multiqcdir="$dirdata/multiqc/";
  test -d $multiqcdir && echo "# Folder $multiqcdir exists" ||
  (mkdir $multiqcdir && echo "# Folder $multiqcdir created");
  
bdir="$dir/Bowtie2";
  test -d $bdir && echo "# Folder $bdir exists" ||
  (mkdir $bdir && echo "# Folder $bdir created");
  
mdir="$dirs/MACS3";
  test -d $mdir && echo "# Folder $mdir exists" ||
  (mkdir $mdir && echo "# Folder $mdir created");

rgtmdir="$dir/RGT-Motif";
  test -d $rgtmdir && echo "# Folder $rgtmdir exists" ||
  (mkdir $rgtmdir && echo "# Folder $rgtmdir created");

dtoolsdir="$dir/deeptools";
  test -d $dtoolsdir && echo "# Folder $dtoolsdir exists" ||
  (mkdir $dtoolsdir && echo "# Folder $dtoolsdir created");

chromapdir="$dir/chromap";
  test -d $chromapdir && echo "# Folder $chromapdir exists" ||
  (mkdir $chromapdir && echo "# Folder $chromapdir created");

blacklhs="/mnt/e/Simon/ChIPSeq_Kelly/Bowtie2/hg38-blacklist.v2.bed"
blacklmm39="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/Genomic_data/Mus/GRCm39/mm39.excluderanges.bed.gz"
blacklmm10="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/Genomic_data/Mus/mm10/mm10.blacklist.bed.gz"
```

    ## # Folder /mnt/s/AG/AG-Scholz-NGS/Daten/ChIP_epiSVF_P3026/Adapterremoval exists
    ## # Folder /mnt/s/AG/AG-Scholz-NGS/Daten/ChIP_epiSVF_P3026/fastqc/ exists
    ## # Folder /mnt/s/AG/AG-Scholz-NGS/Daten/ChIP_epiSVF_P3026/multiqc/ exists
    ## mkdir: cannot create directory ‘/mnt/e/Simon/P3026_ChIP_Seq_epiSVF/Bowtie2’: No such file or directory
    ## # Folder /mnt/s/AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF/MACS3 exists
    ## mkdir: cannot create directory ‘/mnt/e/Simon/P3026_ChIP_Seq_epiSVF/RGT-Motif’: No such file or directory
    ## mkdir: cannot create directory ‘/mnt/e/Simon/P3026_ChIP_Seq_epiSVF/deeptools’: No such file or directory
    ## mkdir: cannot create directory ‘/mnt/e/Simon/P3026_ChIP_Seq_epiSVF/chromap’: No such file or directory

``` bash

source activate chipseq
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
