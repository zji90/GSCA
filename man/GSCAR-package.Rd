\name{GSCAR-package}
\alias{GSCAR-package}
\alias{GSCAR-package}
\docType{package}
\title{
GSCAR: Geneset Context Analysis R Platform
}
\description{

GSCAR analyzes biological contexts enriched within given patterns of geneset expression activity.
GSCAR takes as input several lists of activated and repressed genes. Though the input genesets could contain any gene which interest users, they are usually dervied from ChIP-chip or ChIP-seq (ChIPx) and gene expression data in one or more biological systems, for example TF target genes (genes that are both TF-bound in the ChIPx data and differentially expressed in the gene expression data). Then GSCAR uses the given genesets to scan through a compendium of gene expression profiles constructed from publicly available gene expression data to search for patterns of geneset expression activity specified by the users. The final output of GSCAR is a ranked table of biological contexts that are significantly enriched with the specified pattern of geneset expression activity.
    After the initial GSCAR analysis, users can further study the predicted biological contexts and related contexts in more detail using the tabSearch function to search for contexts of interest in the human or mouse compendium, and the GSCAReda function to visualize and test for differences in geneset expression activities of the recovered contexts. Further functions to help annotate peaks and construct TF target genes are also provided if users are interested in exploring enriched biological contexts in given TF expression and target gene activity.
    Besides traditional R functions, GSCAR also provides a user-friendly interactive user interface based on R shiny. Users can run GSCARui function to run the UI in the web browser on their own computer (need to install shiny and GSCARdata package) or go to http://spark.rstudio.com/jzc19900805/GSCAR/ to run the UI on shiny server (only a web browser is required, do not need to install GSCAR, GSCARdata or R).
}
\details{
\tabular{ll}{
Package: \tab GSCAR\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-1-13\cr
License: \tab GPL-2\cr
}
}
\author{
Author: Zhicheng Ji, Hongkai Ji
Maintainer: Zhicheng Ji <zji4@jhu.edu>
}
\references{
George Wu, et al. ChIP-PED enhances the analysis of ChIP-seq and ChIP-chip data. Bioinformatics 2013 Apr 23;29(9):1182-1189.
}
\keyword{ package, GSCAR }
