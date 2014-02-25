\name{GSCAplot}
\alias{GSCAplot}
\title{
Visualize GSCA output
}
\description{
GSCAplot visualizes the output from GSCA. For one geneset, GSCAplot makes histograms of geneset activities for all samples and samples in each of most significantly enriched biological contexts. For two genesets, GSCAplot makes a scatter plot of sample activities of first geneset versus sample activities of second geneset plotting, and the most significantly enriched biological contexts are highlighted in the plot. For more than two genesets, GSCAplot produces two heatmap according to geneset activities. The first heatmap shows the geneset activities of all samples and indicates which samples belong to enriched biological contexts. The second heatmap shows the geneset activities of samples exhibiting given geneset activity pattern, and the most significantly enriched biological contexts are highlighted.
}
\usage{
GSCAplot(GSCAoutput,N=5,plotfile=NULL,Title=NULL)
}
\arguments{
  \item{GSCAoutput}{
    Exact output from GSCA.
}
  \item{N}{
    N is a numeric value ranging from 1 to 5. It specifies the number of top-ranked biological contexts to plot from the GSCA analysis.
}
  \item{plotfile}{
    A character value specifying the path to save the GSCA plot. If plotfile is null, the plot will not be saved and will appear directly in R console.  
}
  \item{Title}{
    A character value specifying the title of the plot.
}
}
\details{
    GSCAplot is a plotting function that acts as an easy-to-use tool to visualize the GSCA output. For one geneset, GSCAplot uses histogram() to first plot a histogram of geneset activities for all samples in the compendium, then plot N histograms of geneset activities for samples in each of top N most significantly enriched biological contexts. For two genesets, GSCAplots uses plot() to make a scatterplot of all samples in the compendium where x-axis is the activity of the first geneset and y-axis is the activity of the second geneset. Then it highlights the top N most significantly enriched biological contexts in different colors and types. Cutoff of the two genesets will also be represented on the scatterplot as one vertical and one horizontal dotted line. For more than two genesets, GSCAplots uses heatmap.2() from gplots package to plot two heatmaps. In the first heatmap, geneset activities of all samples in the compendium will be shown. A color legend will be drawn on the left upper corner of the heatmap so that users will know the corresponding activity value each color represents. Above the heatmap there is a color bar of light and dark blue indicating which samples exhibit the specific geneset activity pattern. In the second heatmap, geneset activity of all samples which exhibit the specific geneset activity pattern will be shown. A color bar above the heatmap uses different colors to indicate top N most significantly enriched biological contexts. A color legend will also appear in the left upper corner of the heatmap.
    If plotfile is not null, then instead of showing the plots directly in the R console, GSCAplot will save the plots to the designated filepath as a pdf file.
    Note that because there are a lot of samples in both human and mouse compendiums, drawing the first heatmap (and sometimes the second heatmap) could take a lot of time especially a large number of genesets are given. GSCAplot only supports a predefined geneset activity pattern and basic plotting options. Users are encouraged to use the interactive UI if they want to interactively determine the geneset activity pattern, gain more powerful plotting options and further customize their plots.
}
\value{
      A plot consisting of several histograms, a scatterplot or two heatmaps will be returned, depending on numbers of genesets users give.
}
\references{
George Wu, et al. ChIP-PED enhances the analysis of ChIP-seq and ChIP-chip data. Bioinformatics 2013 Apr 23;29(9):1182-1189.
}
\author{
Zhicheng Ji, Hongkai Ji
}
\examples{
## Constructing genedata and pattern.
## Example of mouse gene Gli1,Gli2 and Gli3, all members of GLI-Kruppel family. Their corresponding Entrez GeneID are 14632,14633 and 14634.
gligenedata <- data.frame(gsname=c("Gli1","Gli2","Gli3"),gene=c(14632,14633,14634),weight=1,stringsAsFactors=FALSE)
glipattern <- data.frame(gsname=c("Gli1","Gli2","Gli3"),acttype="High",cotype="Norm",cutoff=0.1,stringsAsFactors=FALSE)

## Case of one geneset: a set of histograms
## Note that for N too large sometimes there is figure margins too large error.
## Decrease N or try to enlarge the plotting area in R console.
oneout <- GSCA(gligenedata[1,],glipattern[1,],"moe4302")
GSCAplot(oneout,N=2)

## Case of two genesets: a scatterplot
twoout <- GSCA(gligenedata[-3,],glipattern[-3,],"moe4302")
GSCAplot(twoout)

## Case of three genesets: two heatmaps, press Enter to switch to the second heatmap
## May take some time, be patient
threeout <- GSCA(gligenedata,glipattern,"moe4302")
GSCAplot(threeout)

## Same plots in designated file path, FILE, which is a pdf file.
## If you want to further customize output plots, for example changing
## range of x-axis, changing titles or altering display of enriched
## biological contexts, please check out the interactive user interface.

GSCAplot(oneout,plotfile=tempfile("plot",fileext=".pdf"),N=2,Title="Demo of one geneset plot")
GSCAplot(twoout,plotfile=tempfile("plot",fileext=".pdf"),Title="Demo of two genesets plot")
GSCAplot(threeout,plotfile=tempfile("plot",fileext=".pdf"),Title="Demo of three genesets plot")

}
\keyword{ GSCAplot}
\keyword{ plot }
