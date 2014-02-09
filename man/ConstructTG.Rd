\name{ConstructTG}
\alias{ConstructTG}
\title{
Construct target genes for a TF using TF-bound genes and differentially expressed genes from ChIP-chip or ChIP-seq and TF perturbation gene expression data.
}
\description{
This function requires users to first analyze their own ChIP-chip and ChIP-seq data to detect significant peaks and then annotate the peaks with their corresponding regulated target genes using the annotatePeaks function in the GSCAR package. Users must also use the limma package to detect differentially expressed genes in their gene expression data (preprocessing and noramzliation can be done with any algorithm the user desires), then the resulting output needs to be annotated into Entrez GeneIDs. Finally, with both inputs ConstructTG will identify the activated and repressed TF target genes.
}
\usage{
ConstructTG(annonPeaksOut, limmaOut)
}
\arguments{
  \item{annonPeaksOut}{
    Output from the annotatePeaks function in the GSCAR package. Contains the genes that correspond to the significant peaks detected from TF ChIP-chip or ChIP-seq data.
}
  \item{limmaOut}{
    Differential expression output from the limma package, and requires the first column of the data.frame to contain the EntrezGeneIDs that match the microarray probeset IDs.
}
}
\details{
This function is designed as one method to allow users to construct target genes after obtaining a list of significant peaks from ChIP-chip or ChIP-seq data and differential expression results from using limma to anaylze their microarray data. It is not designed to be flexible to account for all methods to obtain TF-bound and/or differentially expressed genes. Users can choose to manually intersect their own TF-bound and differentially expressed genes by classifying activated genes as genes, whose expression increases when the TF expression increases and repressed genes as genes, who expression decreases when the TF expression increases. Note, that significant cutoffs for peaks and differentially expressed genes need to be already applied prior to input.
}
\value{
    Returns a list with two items:
    \item{PosTG }{Activated TF target genes}
    \item{NegTG }{Repressed TF target genes}
}
\references{
George Wu, et al. ChIP-PED enhances the analysis of ChIP-seq and ChIP-chip data. Bioinformatics 2013 Apr 23;29(9):1182-1189.
}
\author{
Zhicheng Ji, Hongkai Ji
}
\examples{
### Read in example ChIP-seq analyzed data output from GSE11431
### for Oct4 in ESCs directly downloaded from NCBI GEO
path <- system.file("extdata",package="GSCAR")
chipxfile <- read.delim(paste(path,"GSM288346_ES_Oct4.txt",sep="/"),
                        header=FALSE,stringsAsFactors=FALSE)

### annotate each peak with the corresponding gene target
annon.out <- annotatePeaks(chipxfile,"mm8",10000,5000)

### Read in example limma output from gene expression data obtained
### by analyzing Oct4 RNAi knockdown gene with RMA then limma 
### from the raw CEL files in GSE4189
### The first column contains the Entrez GeneID for each probeset ID
### annotated using the mouse4302.db package in Bioconductor.
gp.out <- read.delim(paste(path,"Pou5f1_E14TG2a_GSE4189_Limma.txt",sep="/"),
                           stringsAsFactors=FALSE)

ConstructTG(annon.out,gp.out)

}
\keyword{ target genes }
