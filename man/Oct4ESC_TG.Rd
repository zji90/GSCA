\name{Oct4ESC_TG}
\alias{Oct4ESC_TG}
\docType{data}
\title{
Oct4 activated (+) and repressed (-) target genes in embryonic stem cells
}
\description{
List of Oct4 target genes derived from ChIP-seq and gene expression data from embryonic stem cells (ESCs). Activated target genes are the first item in the list and repressed target genes are the second item in the list.
}
\usage{data(Oct4ESC_TG)}
\format{
  The format is:
List of 2
 $ : chr [1:519] "100678" "106298" "14609" "12468" ...
 $ : chr [1:337] "246703" "15441" "70579" "20333" ...
}
\details{
Oct4 target genes are defined as genes that are both predicted to be TF-bound in E14 ESCs and differentially expressed after Oct4 knockdown via RNAi in E14TG2a ESCs.
}
\source{
Chen X, Xu H, Yuan P, Fang F et al. Integration of external signaling pathways with the core transcriptional network in embryonic stem cells. Cell 2008 Jun 13;133(6):1106-17.

Loh YH, Wu Q, Chew JL, Vega VB et al. The Oct4 and Nanog transcription network regulates pluripotency in mouse embryonic stem cells. Nat Genet 2006 Apr;38(4):431-40.
}
\references{
http://www.ncbi.nlm.nih.gov/geo/
}
\examples{
data(Oct4ESC_TG)
}
\keyword{datasets}
