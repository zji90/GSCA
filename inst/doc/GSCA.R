### R code from vignette source 'GSCA.Rnw'

###################################################
### code chunk number 1: GSCA.Rnw:53-57
###################################################
library(GSCA)
data(Oct4ESC_TG)
head(Oct4ESC_TG[[1]]) ##Show some positive target genes of Oct4
head(Oct4ESC_TG[[2]]) ##Show some negative target genes of Oct4


###################################################
### code chunk number 2: GSCA.Rnw:62-69
###################################################
activenum <- length(Oct4ESC_TG[[1]])
repressnum <- length(Oct4ESC_TG[[2]])
Octgenedata <- data.frame(
      gsname=c("TF",rep("TG",activenum+repressnum)),
      gene=c(18999,Oct4ESC_TG[[1]],Oct4ESC_TG[[2]]),
      actrep=c(rep(1,1+activenum),rep(-1,repressnum)),
      stringsAsFactors=FALSE)


###################################################
### code chunk number 3: GSCA.Rnw:73-79
###################################################
Octpattern <- data.frame(
      gsname=c("TF","TG"),
      acttype="High",
      cotype="Norm",
      cutoff=0.1,
      stringsAsFactors=FALSE)


###################################################
### code chunk number 4: GSCA.Rnw:84-87
###################################################
displayoct <- Octoutput <- GSCA(Octgenedata,Octpattern,"moe430",Pval.co=0.05,directory=NULL)
displayoct[[1]]$SampleType <- substr(displayoct[[1]]$SampleType,1,25)
head(displayoct[[1]]) ## Partial results of the ranking table


###################################################
### code chunk number 5: GSCA.Rnw:95-96
###################################################
GSCAplot(Octoutput,N=5,plotfile=NULL,Title="GSCA plot of Oct4 in ESC")


###################################################
### code chunk number 6: GSCA.Rnw:107-108
###################################################
data(STAT1_TG) ### Note, only activated (+) STAT1 target genes were found


###################################################
### code chunk number 7: GSCA.Rnw:112-115
###################################################
Statgenenum <- length(STAT1_TG)
Statgenedata <- data.frame(gsname=c("TF",rep("TG",Statgenenum)),gene=c(6772,STAT1_TG),actrep=1,stringsAsFactors=FALSE)
Statpattern <- data.frame(gsname=c("TF","TG"),acttype="High",cotype="Norm",cutoff=0.1,stringsAsFactors=FALSE)


###################################################
### code chunk number 8: GSCA.Rnw:119-121
###################################################
Statoutput <- GSCA(Statgenedata,Statpattern,"hgu133a",Pval.co=0.05,directory=NULL)
head(Statoutput[[1]]) 


###################################################
### code chunk number 9: GSCA.Rnw:127-129
###################################################
GSE7123out <- tabSearch("GSE7123","hgu133a")
GSE7123out


###################################################
### code chunk number 10: GSCA.Rnw:135-136
###################################################
GSE7123followup <- GSCAeda(Statgenedata,Statpattern,"hgu133a",GSE7123out,Pval.co=0.05,Ordering="Average",Title=NULL,outputdir=NULL) 


###################################################
### code chunk number 11: GSCA.Rnw:139-141
###################################################
GSE7123followup$Tstats
GSE7123followup$Pval


