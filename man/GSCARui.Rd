\name{GSCARui}
\alias{GSCARui}
\title{
Initiating GSCAR interactive UI
}
\description{
GSCARui initiates in the web browser an interactive user interface of GSCAR built using R shiny. This user interface enables users to easily perform nearly all standard GSCAR functions in GSCAR package, and provides more powerful and useful options to specify geneset activity patterns, investigate interested biological contexts and customize output plots and tables. For a complete user manual of GSCARui, please refer to the "Manual" section after the user interface has been launched.
}
\usage{
GSCARui()
}
\details{
The purpose of building GSCAR interactive user interface is to provide an easy way for all users to perform analysis tools offered by GSCAR, even though the users do not have any prior knowledge in computer programming or statistics. 
GSCARui provides users more handy ways to input their original dataset into the program. Users who do not have much experience using R may find themselves having difficulties building genedata and pattern datasets required by standard GSCAR functions. In comparison, GSCARui offers more covenient ways to directly type in gene IDs and specify parameters like cutoff types using pull down menus. Users can also check instantly how many genes they inputted are recorded in a given compendium and decide what geneset to be used in further analysis process.
GSCARui also provides users more flexible and direct means to specify geneset activity patterns. For different number of geneset, GSCARui will automatically generate control panels which are most suitable for users to interactively choose the geneset activity patterns. Users can not only specify geneset activity patterns using traditional GSCAR options, but they can also choose geneset activity pattern on histograms, scatterplots and heatmaps by point and click, which makes the process easier and more explicit.
In addition, GSCARui offers more powerful analysis and plotting options. Both p-value and foldchange cutoffs can be given interactively to select the enriched biological contexts. Besides displaying top ranked enriched biological contexts, users can also select specific biological contexts to be displayed on the plot. Finally, users can specify plotting details like x-axis range and titles of the plots if they want to keep the plots for future use.
Thanks to the shiny server, users can type in the URL: http://spark.rstudio.com/jzc19900805/GSCAR/ to directly launch the UI in their web browser. This does not require any dependent R packages or even R itself installed on users computer. All required is a web browser and the URL.
Please check the "Manual" page in the UI for a complete user manual of the UI.

}
\value{
A UI will be shown in users' default web browser.
}
\references{
George Wu, et al. ChIP-PED enhances the analysis of ChIP-seq and ChIP-chip data. Bioinformatics 2013 Apr 23;29(9):1182-1189.
}
\author{
Zhicheng Ji, Hongkai Ji
}
\seealso{
GSCAR
}
\examples{

## Running this will launch the UI in users' default web browser.
\dontrun{
GSCARui()
}
}
\keyword{ GSCARui }
