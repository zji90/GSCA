GSCA: Gene Set Context Analysis
====

## Overview
Gene Set Context Analysis (GSCA) is an open source software package to transform massive amounts of Publicly available gene Expression Data (PED) into a tool for making new discoveries.  GSCA is constructed based on a large collection of human and mouse gene expression data consisting of 25,000+ consistently normalized samples. Users can interactively visualize and examine transcriptional activities of genes and gene sets in these samples
which represent a broad spectrum of biological contexts such as different cell lines, tissues, diseases, and developmental time points. Given one or multiple gene sets, GSCA can query the expression compendium to systematically identify biological contexts associated with specific gene set activity patterns. 

## GSCA Online User Interface
GSCA user interface can be directly launched online without installing any software package: https://zhiji.shinyapps.io/GSCA. However, currently the online version only allows one concurrent user and running the online user interface could be slower than running GSCA on a local computer. Users are thus recommended to  install GSCA and data packages on their own computers with following procedures.

## GSCA and Data Package Installation

GSCA software can be installed via Github (recommended) and Bioconductor. 
Users should have R installed on their computer before installing GSCA. R can be downloaded here: http://www.r-project.org/.

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
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("GSCA","zji90")
```
To launch user interface after installation, run following commands in R:
```{r }
library(GSCA)
GSCAui()
```
For users with R programming experience, command line tools are also available in GSCA R package. Please check the manual package included in the package for details. For now, command line tools only support numeric POI.

### Install GSCA via Bioconductor
GSCA can also be installed via Bioconductor. Note that the GSCA package is not most up-to-dated on Bioconductor. To install GSCA via Bioconductor, run the following commands in R:
```{r }
source("http://bioconductor.org/biocLite.R")
biocLite("GSCA")
```

## GSCA Demonstration Video
For users who are not familiar with GSCA, here is a demonstration video on Youtube for a quick walk-through: https://www.youtube.com/watch?v=wqv_dmlxdcI

## GSCA User Manual
The user manual for GSCA GUI can be directly opened in GSCA GUI (top-left corner) or available at https://github.com/zji90/GSCA/blob/master/inst/shiny/www/GSCAmanual.pdf. The user manual for command line tools can be viewed in R and is also available at http://www.bioconductor.org/packages/release/bioc/vignettes/GSCA/inst/doc/GSCA.pdf. The documentation of all R functions in GSCA is available at http://www.bioconductor.org/packages/release/bioc/manuals/GSCA/man/GSCA.pdf.

## Contact the Author
Author: Zhicheng Ji, Hongkai Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zji4@jhu.edu)

Or open a new issue in this Github page

## FAQ
### 1. Why GSCA GUI does not react in my web browser?
With certain extensions installed on Chrome (e.g. Ghostery), in some cases GSCA GUI may not react normally in Chrome. This is because such extensions prevent some Javascript functions from running and will affect all web applications using these Javascript modules. The easiest solution is to uninstall or unable these extensions when you are running GSCA GUI.

There are other two solutions if the problem still cannot be solved. One solution is to download and install Rstudio (http://www.rstudio.com/products/rstudio/download/) and launch GSCA GUI within Rstudio. Another solution is to change the default web browser to IE (Windows) or Safari (Mac) and launch GSCA GUI again in R, or change to a computer with another operating system. 

### 2. Can I upload my own database/compendium?
Yes. GSCA provides the function for users to upload their own compendium which includes a gene expression file and an annotation file. The function is under "Select geneset and compendium" -> "Upload your own compendium". Please read carefully the instrutions before preparing these files as these files should follow some rigorous formats.

### 3. Can I directly use the gene expression profiles in the four data packages?
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

### 4. How can I reproduce the results for the interactive POI system?
GSCA offers a convenient function for users to save the old POI in both numeric and interactive POI systems. The function appears at the bottom-left region of the GSCA user interface. Users can load the POI file into the GSCA system next time when they want to reproduce the GSCA analysis results with the same POI.

To correctly load the saved POI file, you should have the exact setting as the setting you used when you saved the POI file. The gene sets, scaling options and choice of numeric/interactive POI system should be exactly the same in order to load the POI. In most cases the problem happens when you fail to choose the same scaling option.




