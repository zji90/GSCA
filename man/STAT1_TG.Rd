\name{STAT1_TG}
\alias{STAT1_TG}
\docType{data}
\title{
STAT1 activated (+) target genes defined from experimental ChIP-seq data and literature survey.
}
\description{
List of STAT1 target genes derived from ChIP-seq data in Hela cells and further refined by making sure each target gene was further supported by experiments in literature as described in GSE15353. No represed target genes were defined.
}
\usage{data(STAT1_TG)}
\format{
  The format is:
 chr [1:23] "9636" "2537" "2633" "1435" "103" "3433" "3434" ...
}
\details{
STAT1 target genes are defined as TF-bound from Hela ChIP-seq data and then further verified as target genes through literature survey. This procedure is described in GSE15353.
}
\source{
Robertson G, Hirst M, Bainbridge M, Bilenky M et al. Genome-wide profiles of STAT1 DNA association using chromatin immunoprecipitation and massively parallel sequencing. Nat Methods 2007 Aug;4(8):651-7.
}
\references{
http://www.ncbi.nlm.nih.gov/geo/
}
\examples{
data(STAT1_TG)
}
\keyword{datasets}
