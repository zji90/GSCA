GSCA: Gene Set Context Analysis
====

## Overview
Gene Set Context Analysis (GSCA) is an open source software package to help with transforming massive amounts of Publicly available gene Expression Data (PED) into a tool for making new discoveries.  GSCA is constructed based on a large collection of human and mouse gene expression data consisting of 25,000+ consistently normalized samples. Users can interactively visualize and examine transcriptional activities of genes and gene sets in these samples
which represent a broad spectrum of biological contexts such as different cell lines, tissues, diseases, and developmental time points. Given one or multiple gene sets, GSCA can query the expression compendium to systematically identify biological contexts associated with specific gene set activity patterns. 

## GSCA and Data Package Installation

Currently GSCA software can be installed via Github (recommended) and Bioconductor. 
Users should have R installed on their computer before installing GSCA. R can be downloaded here: http://www.r-project.org/. To ensure that GSCA GUI runs stably, users are recommended (but not required) to install Rstudio which is a user interface of R. Rstudio can be downloaded here: http://www.rstudio.com/products/rstudio/download/. 

### Install Data Packages
To run GSCA, users should first install at least one of the four data packages in R, which can be done by running the following commands:
```{r }
#Affyhgu133aExpr (Human GPL96 array, about 316 MB in size) 
source("http://bioconductor.org/biocLite.R")
biocLite("Affyhgu133aExpr")
#Affymoe4302Expr (Mouse GPL1261 array, about 418 MB in size)
source("http://bioconductor.org/biocLite.R")
biocLite("Affymoe4302Expr")
#Affyhgu133Plus2Expr (Human GPL570 array, about 224 MB in size)
source("http://bioconductor.org/biocLite.R")
biocLite("Affyhgu133Plus2Expr")
#Affyhgu133A2Expr (Human GPL571 array, about 9 MB in size)
source("http://bioconductor.org/biocLite.R")
biocLite("Affyhgu133A2Expr")
```

### Install GSCA via Github (Recommended)
To install the latest version of GSCA package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("GSCA","zji90")
```
To invoke user interface after installation, run following commands in R:
```{r }
library(GSCA)
GSCAui()
```

### Install GSCA via Bioconductor
GSCA can also be installed via Bioconductor. Note that the GSCA package is not most up-to-dated on Bioconductor. To install GSCA via Bioconductor, run the following commands in R:
```{r }
source("http://bioconductor.org/biocLite.R")
biocLite("GSCA")
```

## GSCA Online User Interface
GSCA user interface can be directly launched online without installing R or any R package: http://spark.rstudio.com/jzc19900805/GSCA. However, running the GSCA online user interface could be significantly slower than running on a local computer.

## GSCA Demonstration Video
For users who are not familiar with GSCA, here is a demonstration video on Youtube for a quick walk-through: https://www.youtube.com/watch?v=1OeZ1PAUMhw

## Contact the Author
Author: Zhicheng Ji, Hongkai Ji

Maintainer: Zhicheng Ji (zji4@jhu.edu)

## FAQ
### 1. Why GSCA GUI does not react in some web browser?
Due to the potential compatibility issues, in some cases GSCA GUI may not react in some web browser (especially Chrome on Mac). One solution is to download and install Rstudio and launch GSCA GUI within Rstudio. Another solution is to change the default web browser to IE (Windows) or Safari (Mac) and launch GSCA GUI again in R, or change to a computer with another operating system (Windows is recommended). 

### 2. Can I directly use the gene expression profiles in the four data packages?
Yes. The four data packages are built primarily for the use of GSCA, but can also be used for other purposes. Note that the four data packages are stored in hdf5 format.
To build a complete gene expression matrix for Affymoe4302Expr data package for example, run the following commands in R:
```{r }
# Install the rhdf5 package on Bioconductor if you haven't done so
library(rhdf5)
library(Affymoe4302Expr)
path <- system.file("extdata",package="Affymoe4302Expr")
load(paste0(path,"/geneid.rda"))
data(Affymoe4302Exprtab)
geneexpr <- t(h5read(paste0(path,"/data.h5"),"expr"))/1000
row.names(geneexpr) <- geneid
###Use SampleID as the column name
colnames(geneexpr) <- Affymoe4302Exprtab[,1]
###Or use annotated Sample type as the column name
colnames(geneexpr) <- Affymoe4302Exprtab[,3]
```

### 3. How can I reproduce the results for the interactive POI system?
GSCA offers a convenient function for users to save the old POI in both numeric and interactive POI systems. The function appears at the bottom-left region of the GSCA user interface. Users can load the POI file into the GSCA system next time when they want to reproduce the GSCA analysis results with the same POI.

To correctly load the saved POI file, you should have the exact setting as the setting you used when you saved the POI file. The gene sets, scaling options and choice of numeric/interactive POI system should be exactly the same in order to load the POI. In most cases the problem happens when you fail to choose the same scaling option.





