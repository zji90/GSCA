\name{GSCA}
\alias{GSCA}
\title{
GSCA
}
\description{
The function takes as input several lists of activated and repressed genes. It then searches through a compendium of publicly available gene expression profiles for biological contexts that are enriched with a specified pattern of gene expression.
}
\usage{
GSCA(genedata,pattern,chipdata,Pval.co=0.05,directory=NULL)
}
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
 \item{Pval.co}{
    A numeric value specifying the adjusted p-value cutoff. Only the biological contexts with significant enrichment above the adjusted p-value cutoff will be reported in the final ranked table output.
}
\item{directory}{
      Either null or a character value giving a directory path. If directory is not null, then additional follow-up GSCA analyses will be performed and stored in the folder specified by directory. If directory is null then no additional follow-up GSCA analyses will be performed.
}
}
\details{
GSCA requires as input user-specified genesets together with their corresponding activity patterns. Each geneset contained the Entrez GeneID of activied and repressed genes. Activated gene means that increases in expression of the gene also increases the overall activity of the whole geneset, while repressed gene means that increases in expression of the gene decreases the overall activity of the whole geneset.

GSCA also requires activity patterns of the genesets. Users can choose either high or low level of activity for each geneset. Cutoffs are given by the users to determine what activity level should be considered high or low. There are three types of cutoffs available: normal, quantile and expression value. For normal cutoff type, a specified p-value (one-sided) based on a fitted normal distribution will be used as cutoff, and all samples having p-value larger(smaller) than this p-value will be considered having high(low) expression activity in a certain geneset.Likewise, for quantile cutoff type, a quantile will be used as cutoff. As for cutoff type of expression level, a numeric value will be directly used as cutoff.

GSCA then searches through the compendium for all samples that exhibit the specified activity pattern of interest. For example, if activity patterns of all genesets are set to be high, then GSCA will find all samples in the compendium that have greater geneset expression levels than the respective cutoffs. Since each of the samples correspond to different biological contexts, the Fisher's exact test will then be used to test the association between each biological context and the geneset activity pattern of interest based on the number of samples in each biological context that exhibits the specified geneset activity pattern of interest.

The final output is a ranked table of biological contexts enriched with the geneset activity pattern of interest. The p-values are also adjusted by the Bonferroni correction.

If directory is not null, then GSCA will peform detail analyses for all contexts in each of the experimental IDs in the final GSCA results table. For each of the experiment IDs, tabSearch will be run to locate all contexts in the compendium for that experiment ID, and then GSCAeda will be run using the same genedata and pattern as input specific to the contexts recovered by tabSearch. See GSCAeda for more details. Note, this automated process could be time-consuming and produce a lot of files and directories.
}
\value{
    Returns a list with
    \item{Ranking }{Data.frame of ranked table of biological contexts significantly enriched with the specified geneset activity pattern. It includes information of Ranking, number of samples exhibiting the given activity pattern, total number of samples, fold change values, adjusted p-values, name of biological context and corresponding experiment ID.}
    \item{Score }{Numeric matrix of geneset expression values for each sample in the compendium. Each row stands for a certain geneset and each column stands for a certain sample.}
    \item{Pattern}{Data.frame of geneset activity pattern. The same as the input value.}
    \item{Cutoff}{Numeric vector of cutoff values calculated for each geneset based on the input pattern.}
    \item{SelectedSample}{Numeric vector of all samples that exhibits the given geneset activity pattern.}
    \item{Totalgene }{Numeric vector of total number of genes use to calculate the geneset activity in each geneset}
    \item{Missinggene }{Numeric vector of number of genes that do not have corresponding expression measuremnts on the platform}
    \item{Species }{Character value of the species analyzed.}

    If directory is not null, then pdf and csv files containing the GSCAeda follow-up analysis results and plots in the directory folder will also be returned.
}
\references{
George Wu, et al. ChIP-PED enhances the analysis of ChIP-seq and ChIP-chip data. Bioinformatics 2013 Apr 23;29(9):1182-1189.
}
\author{
Zhicheng Ji, Hongkai Ji
}
\examples{
## First load the TF target genes derived from Oct4 ChIPx data
## in embryonic stem cells. The data is in the form of a list
## where the first item contains the activated (+) target genes in
## Entrez GeneID format and the second item contains the repressed (-)
## target genes in Entrez GeneID format.
data(Oct4ESC_TG)

## We want to analyze Oct4, so we need to specify the EntrezGeneID for Oct4
## and input the activated (+) and repressed (-) target genes of Oct4.
## Constucting the input genedata required by GSCA. There are two genesets
## one is the TF and another is the TF target genes. Note that constructing genedata
## with many genesets could be laborious, so using the interactive UI is recommended to 
## easily start up the analysis.
activenum <- length(Oct4ESC_TG[[1]])
repressnum <- length(Oct4ESC_TG[[2]])
Octgenedata <- data.frame(gsname=c("GS1",rep("GS2",activenum+repressnum)),gene=c(18999,Oct4ESC_TG[[1]],Oct4ESC_TG[[2]]),actrep=c(rep(1,1+activenum),rep(-1,repressnum)),stringsAsFactors=FALSE)

## We are interested in the pattern that TF and its target genes are all highly expressed.
## We also need to define how high the cutoffs should be such
## that each cutoff corresponds to the p-value of 0.1
## based on fitted normal distributions.
## Constructing pattern required by GSCA, all geneset names in genedata should appear
## exactly once in the first column
Octpattern <- data.frame(gsname=c("GS1","GS2"),acttype="High",cotype="Norm",cutoff=0.1,stringsAsFactors=FALSE)

## Lastly, we specify the chipdata to be "moe430" and the significance of enriched
## biological contexts must be at least 0.05 to be reported.
Octoutput <- GSCA(Octgenedata,Octpattern,"moe430",Pval.co=0.05)

## The first item in the list 'Octoutput[[1]]' contains the ranked table, which
## can then be saved. Additionally, we may be interested in plotting the results
## to visualize the enriched biological contexts within given geneset activity.
## Here, N specifies the top 5 significant biological contexts.
## Since plotfile is NULL, the plot directly shows up in R.
## Check GSCAplot for more details.
GSCAplot(Octoutput,N=5,plotfile=NULL,Title="GSCA plot of Oct4 in ESC")

## If you would like detailed follow-up analyses to be automatically performed
## for the Oct4 analyses in ESCs, just specify a file directory.
## Check GSCAeda for more details.

\dontrun{
Octoutput <- GSCA(Octgenedata,Octpattern,"moe430",Pval.co=0.05,directory="DIR")
}

## All output will be stored in the directory, where DIR should be your
## desired output folder. This process may be time-consuming and 
## generate a lot of files. Alternatively, see GSCAeda for
## more info on manual alternatives.
}
\keyword{ GSCA }
