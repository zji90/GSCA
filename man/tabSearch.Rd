\name{tabSearch}
\alias{tabSearch}
\title{
Searches through the Affymetrix Human hgu133a Array (GPL96) or Affymetrix Mouse 430 2.0 Array (GPL1261) compendium data for biological contexts of interest.
}
\description{
tabSearch requires users to provide keyword(s), the species, and either 'AND' or 'OR'. Then the function uses grep and the keywords to iteratively search for biological contexts or experiment IDs that match the keywords, where 'AND' requires all recovered contexts to satisfiy all keywords and 'OR' requires all recovered contexts to match at least one keyword.
}
\usage{
tabSearch(keyword, chipdata, option = "OR")
}
\arguments{
  \item{keyword}{
    A character vector of biological context words or experiment IDs. e.g. 'liver' or 'GSE7123'.
}
  \item{chipdata}{
    A character value of either 'hgu133a' or 'moe430'. This argument specifies which compendium to use. Requires the corresponding data package.
}
  \item{option}{
    Either 'AND' or 'OR' to specify whether the recovered contexts need to be found by ALL keywords (AND) or found by at least one keyword (OR).
}
}
\details{
If the users want to search for a specific list of contexts, simply input the contexts as a character vector, where each element is a different context. Alternatively, the contexts can also be a series of keywords in short-hand. The 'AND' option is primarly used when users want to search for contexts from a specific experiment. In most cases 'OR' should be used.

After tabSearch finishes running, it will return a list of contexts that match the inputted keywords and parameters. Users can then further study these contexts for activity of given gensets using the function GSCAeda.
}
\value{
    A data frame consisting of three columns. The 1st column is the experiment ID, the 2nd column is the biological context label, and the 3rd column is the number of samples for each biological context. 
}
\references{
George Wu, et al. ChIP-PED enhances the analysis of ChIP-seq and ChIP-chip data. Bioinformatics 2013 Apr 23;29(9):1182-1189.
}
\author{
Zhicheng Ji, Hongkai Ji
}
\examples{
library(GSCA)
## Search for all contexts in GSE7123 in hgu133a
tabSearch("GSE7123","hgu133a")

## Search for all contexts labeled 'fetal' or 'liver' in moe430
tabSearch(c("Fetal","Liver"),"moe430")

## Search for all contexts labeled 'fetal liver' AND in GSE13044 in moe430
tabSearch(c("Fetal Liver","GSE13044"),"moe430","AND")
} 
\keyword{ search }
