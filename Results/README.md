---
title: "ChIP_Wt1 - Results"
author: "Kelterborn"
date: "2024-01-31"
output:
  html_document: 
    toc: true
    keep_md: true
    self_contained: false
    df_print: kable
---

<!-- <style> -->
<!-- .vscroll-plot { -->
<!--     width: 1000px; -->
<!--     height: 1000px; -->
<!--     overflow-y: scroll; -->
<!--     overflow-x: scroll; -->
<!-- } -->
<!-- </style> -->




# R Prepare System
## R update and load librarys
BiocManager::install("")


```r
BiocManager::install(update = TRUE, ask = FALSE)
library(dbplyr)
library(tidyverse)
library(ChIPseeker)
library(rtracklayer)
library(trackViewer)
library(GenomicRanges)
library(IRanges)
library(ChIPpeakAnno)
library(AnnotationHub)
library(ggplot2)
library(viridis)
library(kableExtra)
library(DT)
library(patchwork)
library(gridExtra)

library(TxDb.Mmusculus.UCSC.mm39.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
```




## R folders


# Unix Prepare System
## Unix Update System


## Unix Folder & Names


# 3 Results

![Analysis Overview](../Data/sheme.pdf){height=100%, width=100%}

![Analysis Overview](../Data/sheme.png){height=100%, width=100%}



```r
knitr::include_graphics("https://github.com/DNAborn/ChIPseq_Wt1/blob/main/Results/sheme.png")
knitr::include_graphics("https://github.com/DNAborn/ChIPseq_Wt1/blob/main/Results/sheme.pdf")
# knitr::include_graphics("./sheme.png")
# knitr::include_graphics("./sheme.pdf")
# knitr::include_graphics("../sheme.png")
# knitr::include_graphics("../sheme.pdf")
# knitr::include_graphics("/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF/ChIPseq_Wt1_P3026/Results/sheme.png")
# knitr::include_graphics("/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF/ChIPseq_Wt1_P3026/Results/sheme.pdf")
getwd()
```



```r
print("include figure in r chunk")
knitr::include_graphics("../Data/sheme.pdf") # doesn't show
```


### Generate combined peak list


#### Peak Tables


#### Hists & Tables
<img src="README_files/figure-html/hits_tables-1.png" width="100%" height="1000px" /><img src="README_files/figure-html/hits_tables-2.png" width="100%" height="1000px" /><img src="README_files/figure-html/hits_tables-3.png" width="100%" height="1000px" /><img src="README_files/figure-html/hits_tables-4.png" width="100%" height="1000px" /><img src="README_files/figure-html/hits_tables-5.png" width="100%" height="1000px" />




#### Overlap Peaks
![](README_files/figure-html/venn_overlaps-1.png)<!-- -->

## Annotate Peaks
#### 1 Run
![](README_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

#### All peaks


#### Venns


