\name{GSCAeda}
\alias{GSCAeda}
\title{
GSCA follow-up exploratory data analysis
}
\description{
GSCAeda is used to further study GSCA significant predictions in more detail to obtain additional insight into biological function. GSCAeda requires users to first run the tabSearch function to identify the biological contexts of interest. By default, GSCAeda will run automatically after an initial GSCA analysis by searching for all contexts related to the experimentID for each significant GSCA prediction. Alternatively, users can use GSCAeda by itself to further study any geneset or biological contexts of interest that are found in the compendium.
   The output of GSCAeda are multiple plots displaying the geneset activity values and genes of interest in the input biological contexts. Also included are the usual GSCA analysis results table showing the enrichment of each contexts for the geneset activity pattern of interest, t-test results (t-statistics and p-values) for all pair-wise combinations of inputted contexts in each geneset, and a summary of raw geneset activity values for each context of interest. Users can then use the raw geneset activity values for further statistical analyses if desired.
}
\usage{
GSCAeda(genedata,pattern,chipdata,SearchOutput,Pval.co=0.05,Ordering="Average",Title=NULL,outputdir=NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{genedata}{
      A data.frame with three columns specifying the input genesets. Each row specifies an activated or repressed gene in a geneset. First column: character value of geneset name specified by the user, could be any name easy to remember e.g. GS1,GS2,...; Second column: numeric value of Entrez GeneID of the gene; Third column: numeric value of 1,-1 indicating whether gene is activated or repressed. 1 for activated gene and -1 for repressed gene. Here, activated gene means that increases in expression of the gene also increases the overall activity of the whole geneset, while increases in expression of the repressed genes will decrease the overall activity of the whole geneset.
}
\item{pattern}{
      A data.frame with four columns indicating the activity patterns corresponding to the given genedata. Each row specifies activity pattern for one geneset. First column: character value of the same geneset name used in genedata, each geneset name in genedata should appear exactly once in this column. Second column: character value of whether high or low activity of the whole geneset is interested. "High" stands for high activity and "Low" stands for low activity. Third column: character value of which cutoff type is going to be used. 3 cutoff types can be specified: "Norm", "Quantile", or "Exprs". If cutoff type is "Norm", then the fourth column should be specified as p-value between 0 and 1, where the geneset expression cutoff will correspond to the specified p-value (one-sided) based on a fitted normal distribution; If cutoff type is "Quantile", then the fourth column should be specified as a desired quantile between 0 and 1, where the geneset expression cutoff will correspond to the specified quantile. Finally, if cutoff type is "Exprs", the geneset expression cutoff will be equal to the value given in the fourth column. Fourth column: numeric value of cutoff value based on different cutoff types specified in the third column.
}

\item{chipdata}{
    A character value of either 'hgu133a' or 'moe430'. This argument specifies which compendium to use. Requires the corresponding data package.
}
  \item{SearchOutput}{
    Output of the tabSearch function. More specifically, a data frame where the 1st column is the ExperimentIDs (GSE ids), the 2nd column is the SampleTypes, and the 3rd column is the sample count for each SampleType.
}
  \item{Pval.co}{
    A numeric value specifying the adjusted p-value cutoff. Only the biological contexts with significant enrichment above the adjusted p-value cutoff will be reported in the final ranked table output.
}
  \item{Ordering}{
    A character value of either one geneset name or 'Average'. If Ordering is one geneset name, the plot of geneset activity values and heatmap of the t-statistics/pvalues will be ordered from the highest to lowest according the Ordering geneset activity value. If Ordering is 'Average', the plots and heatmap will be organized by the average rank across all geneset activity values.
}
  \item{Title}{
    Title of the plot, will appear on the top of the plot.
}
  \item{outputdir}{
    Either null or a character value giving the directory in which GSCAeda will save the output files.
}
}
\details{
GSCAeda is designed to be used in combination with tabSearch after an initial GSCA analysis. GSCAeda is used to further study each predicted biological context in more detail by comparing the functional activity across related contexts through the geneset activities. To do so, GSCAeda requires users to specific genedata, pattern, species, pval cutoff, and the search results from tabSearch containing the list of biological contexts of interest.
Then GSCAeda will calculate the mean and standard deviation of each geneset activity value for each inputted context, and will perform t-tests comparing the mean geneset activity values for all pair-wise combinations of inputted contexts, and test for enrichment of the geneset activity pattern of interest. The results will be shown in several plots and tables(number of files varying with number of given genesets), along with the raw geneset activity values for further followup statistical analyses. Check the value part of this help file to see how GSCAeda saves the outputs.
For information on the GSCA parameters, see the GSCA help file which explains in more detail how functional enrichment of a geneset activity pattern of interest is tested.
}
\value{
    If outputdir is specified, GSCAeda will first produce a boxplot depicting the distribution of all geneset activities in different biological contexts of interest. Then, for each geneset GSCAeda will produce tow heatmaps showing respectively the t-statistics and p-values obtained from the t-tests testing the mean of geneset activity for each pair-wise combination of the input biological contexts. Finally, GSCAeda will output two csv files. The first one contains the raw geneset activity values for each input context and the second one contains the mean and standard deviation of the geneset activity values for each context, the GSCA enrichment test results, and the p-values/t-statistics of the t-tests.
If outputdir is NULL, all plots and the result table will be directly displayed in the R console.
}
\references{
George Wu, et al. ChIP-PED enhances the analysis of ChIP-seq and ChIP-chip data. Bioinformatics 2013 Apr 23;29(9):1182-1189.
}
\author{
Zhicheng Ji, Hongkai Ji
}
\seealso{
GSCA
}
\examples{
library(GSCA)

## Load example STAT1 target genes defined ChIP-seq and literature
data(STAT1_TG)

## Construct genedata and pattern using the same way as GSCA
Statgenenum <- length(STAT1_TG)
Statgenedata <- data.frame(gsname=c("GS1",rep("GS2",Statgenenum)),gene=c(6772,STAT1_TG),actrep=1,stringsAsFactors=FALSE)
Statpattern <- data.frame(gsname=c("GS1","GS2"),acttype="High",cotype="Norm",cutoff=0.1,stringsAsFactors=FALSE)

## Find all contexts in human compendium from GSE7123
GSE7123out <- tabSearch("GSE7123","hgu133a")

## Run GSCAeda
GSCAeda(Statgenedata,Statpattern,"hgu133a",GSE7123out,Pval.co=0.05,Ordering="Average",Title=NULL,outputdir=NULL) 

## To save the results, instead of displaying in R console, specifiy an outputdir argument

GSCAeda(Statgenedata,Statpattern,"hgu133a",GSE7123out,Pval.co=0.05,Ordering="Average",Title=NULL,outputdir=tempdir())
}
\keyword{ GSCAeda }
