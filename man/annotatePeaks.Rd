\name{annotatePeaks}
\alias{annotatePeaks}
\title{
Annotate ChIPx peaks with genes by Entrez GeneIDs
}
\description{
This function finds all genes that overlap with each peak detected from TF ChIP-chip or ChIP-seq data. Assigned genes are assumed to be genes bound by the TF.
}
\usage{
annotatePeaks(inputfile, genome, up = NULL, down = NULL)
}
%- could also be 'usage' for other objects documented here.
\arguments{
  \item{inputfile}{
    A data.frame where each row corresponds to a peak. The first column is the chromosome on which the peak is found (e.g., chr1) and the second and third columns are the peak starting and ending sites.
}
  \item{genome}{
    Should be one of 'hg19', 'hg18', 'mm9', or 'mm8' genome. More genomes may be supported in future versions of GSCAR.
}
  \item{up}{
    Region upstream of the TSS. A gene will be annotated to a peak if the region upstream to downstream of each gene TSS, as defined by the up and down arguments, overlap with the peak.
}
  \item{down}{
    Region downstream of the TSS. A gene will be annotated to a peak if the region upstream to downstream of each gene TSS, as defined by the up and down arguments, overlap with the peak.
}
}
\details{
A gene will be annotated to a peak if the region upstream to downstream of each gene TSS, as defined by the up and down arguments, overlap with the peak.
}
\value{
    Returns a data.frame with the same columns as the input data.frame and an additional column containing the Enterz GeneIDs for all genes that overlap with the peak. Multiple genes will be separated with ';' and '-9' will be reported if no genes are found.
}
\references{
Chen X, Xu H, Yuan P, Fang F et al. Integration of external signaling pathways with the core transcriptional network in embryonic stem cells. Cell 2008 Jun 13;133(6):1106-17.

George Wu, et al. ChIP-PED enhances the analysis of ChIP-seq and ChIP-chip data. Bioinformatics 2013 Apr 23;29(9):1182-1189.
}
\author{
Zhicheng Ji, Hongkai Ji
}
\examples{
### Read in example ChIP-seq analyzed data output from GSE11431
### for Oct4 in ESCs directly downloaded from NCBI GEO
path <- system.file("extdata",package="GSCAR")
inputfile <- read.delim(paste(path,"GSM288346_ES_Oct4.txt",sep="/"),
                        header=FALSE,stringsAsFactors=FALSE)

### Note that 1st column is chr, 2nd and 3rd columns are starting and ending sites of peaks
### Remaining columns are other output from the peak detection algorithm
head(inputfile)

### annotatePeaks only requires the first 3 columns
annon.out <- annotatePeaks(inputfile,"mm8",10000,5000)
head(annon.out)
}
\keyword{ annotate }
\keyword{ peaks }
