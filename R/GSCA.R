GSCA <- function(genedata,pattern,chipdata,scaledata=F,Pval.co=0.05,directory=NULL) {

      ###presettings
      genedata[,1] <- as.character(genedata[,1])
      pattern[,1] <- as.character(pattern[,1])
      #setting path of data package
      path <- system.file("extdata",package=paste0("Affy",chipdata,"Expr"))
      #read in gene names
      load(paste0(path,"/geneid.rda"))
      #read in reference tables for given compendium
      if(chipdata == "hgu133a") {
            if (!require(Affyhgu133aExpr)) {
                  stop("Affyhgu133aExpr Package is not found")
            } else {
                  data(Affyhgu133aExprtab)
                  tab <- Affyhgu133aExprtab
            }
      } else if(chipdata == "moe4302"){
            if (!require(Affymoe4302Expr)) {
                  stop("Affymoe4302Expr Package is not found")
            } else {
                  data(Affymoe4302Exprtab)
                  tab <- Affymoe4302Exprtab
            }
      } else if(chipdata == "hgu133A2"){
            if (!require(Affyhgu133A2Expr)) {
                  stop("Affyhgu133A2Expr Package is not found")
            } else {
                  data(Affyhgu133A2Exprtab)
                  tab <- Affyhgu133A2Exprtab
            }
      } else if(chipdata == "hgu133Plus2"){
            if (!require(Affyhgu133Plus2Expr)) {
                  stop("Affyhgu133Plus2Expr Package is not found")
            } else {
                  data(Affyhgu133Plus2Exprtab)
                  tab <- Affyhgu133Plus2Exprtab
            }
      } else {
            stop("Please enter valid name for chipdata. Current Supported chipdata: 'hgu133a', 'moe4302', 'hgu133Plus2', 'hgu133A2'")
      }
      
      tabsamplename <- tab$SampleName
      
      #check whether each geneset has at least one gene on the compendium
      genesetname <- NULL
      for (tmpgenesetname in unique(genedata[,1])) {
            if(sum(geneid %in% genedata[genedata[,1]==tmpgenesetname,2]) == 0) {
                  warning(paste("No matching target genes found on the compendium for gene set",tmpgenesetname))
            } else {
                  genesetname <- c(genesetname, tmpgenesetname)
            }
      }
      if (length(genesetname) == 0)
            stop("No matching target genes found on the compendium for all gene sets")
      
      selectsample <- 1:nrow(tab)
      activity <- matrix(0,nrow = length(genesetname),ncol = nrow(tab))
      rownames(activity) <- genesetname
      genesetcutoff <- genesettotalgenenum <- genesetmissinggene <- rep(0,length(genesetname))
      names(genesetcutoff) <- names(genesettotalgenenum) <- names(genesetmissinggene)  <- genesetname
      for (genesetid in 1:length(genesetname)) {
            ###Scoring geneset activity
            singlegeneset <- genesetname[genesetid]           
            currentgeneset <- genedata[genedata[,1] == singlegeneset & genedata[,2] %in% geneid,]
            score <- rep(0, nrow(tab))
            for (i in 1:nrow(currentgeneset)) {
                  e <- as.vector(h5read(paste0(path,"/data.h5"),"expr",index=list(NULL,which(geneid==currentgeneset[i,2]))))/1000
                  if (scaledata)
                        e <- scale(e)
                  score <- score + currentgeneset[i,3]*e
            }
            score <- score/nrow(currentgeneset)
            
            missinggene <- setdiff(genedata[singlegeneset == genedata[,1],2],geneid)
            genesetmissinggene[genesetid] <- length(missinggene)
            genesettotalgenenum[genesetid] <- length(genedata[,1] == singlegeneset && genedata[,2] %in% geneid)
            activity[genesetid,] <- score
            ###Find samples matching the given pattern
            singlepattern <- pattern[pattern[,1]==singlegeneset,]
            if (singlepattern[,3] == "Norm") {
                  if (singlepattern[,2] == "High") {
                        cutoff <- qnorm(1-singlepattern[,4],mean(score),sd(score))
                        selectsample <- intersect(selectsample,which(score >= cutoff))
                  } else if (singlepattern[,2] == "Low") {
                        cutoff <- qnorm(singlepattern[,4],mean(score),sd(score))
                        selectsample <- intersect(selectsample,which(score < cutoff))
                  } else {
                        stop(paste("Second Column of pattern in",singlegeneset,"is not correctly given"))      
                  }
            } else if (singlepattern[,3] == "Quantile") {
                  if (singlepattern[,2] == "High") {
                        cutoff <- quantile(score,1-singlepattern[,4])
                        selectsample <- intersect(selectsample,which(score >= cutoff))
                  } else if (singlepattern[,2] == "Low") {
                        cutoff <- quantile(score,singlepattern[,4])
                        selectsample <- intersect(selectsample,which(score < cutoff))
                  } else {
                        stop(paste("Second Column of pattern in",singlegeneset,"is not correctly given"))      
                  }
            }  else if (singlepattern[,3] == "Exprs") {
                  if (singlepattern[,2] == "High") {
                        cutoff <- singlepattern[,4]
                        selectsample <- intersect(selectsample,which(score >= cutoff))
                  } else if (singlepattern[,2] == "Low") {
                        cutoff <- singlepattern[,4]
                        selectsample <- intersect(selectsample,which(score < cutoff))
                  } else {
                        stop(paste("Second Column of pattern in",singlegeneset,"is not correctly given"))      
                  }                
            } else {
                  stop(paste("Cutoff type pattern of geneset",singlegeneset,"is not correctly given"))
            }
            genesetcutoff[genesetid] <- cutoff
      }

      ExpID <- tab[selectsample,"ExperimentID"]
      tmpTypes <- tab[selectsample,"SampleType"]

  tabTWO <- table(tab$SampleType)
  tabTWO <- names(tabTWO)[tabTWO > 2]
  tab <- tab[tab$SampleType %in% tabTWO,]
  ExpID <- ExpID[tmpTypes %in% tab$SampleType]
  tmpTypes <- tmpTypes[tmpTypes %in% tab$SampleType]

  ## Fisher's exact test
  ## Let X = sample type i
  if(length(tmpTypes)>0){
  ttT <- table(tmpTypes) ## total of each sample types in enrichment region
  sT <- sum(ttT) ## total samples in enrichment region
  bgT <- table(tab$SampleType)
  bgT <- bgT[names(ttT)]
  ContextN <- sum(table(tab$SampleType)>2)

  SCORE <- matrix(0, nrow=length(ttT), ncol=4)
  ExperimentID <- rep("0",length(bgT))
  for(i in 1:length(bgT)){
    r1c1 <- ttT[i] ### num of samples X in enrichment region
    r1c2 <- sT - ttT[i] ### num != sampletypes x in enrichment region
    r2c1 <- bgT[i] - ttT[i] ### num of sampletypes X not in enrichment region
    r2c2 <- length(tab$SampleType) - r1c1 - r1c2 - r2c1
    tmpmat <- matrix(c(r1c1,r2c1,r1c2,r2c2),ncol=2)
    SCORE[i,1] <- fisher.test(tmpmat,alternative="greater")$p.value
    SCORE[i,2] <- round(((as.numeric(ttT)[i]+sT/length(tab$SampleType)) / 
                        (as.numeric(bgT)[i]+1))/(
                        sT/length(tab$SampleType)),3)
    SCORE[i,3] <- min(SCORE[i,1]*ContextN,1)
    ExperimentID[i] <- paste(unique(unlist(strsplit(
                             paste(ExpID[tmpTypes==names(ttT)[i]],
                             collapse=";"),";"))),collapse=";")
  }

  FIN <- data.frame(as.numeric(ttT),as.numeric(bgT),SCORE[,2],SCORE[,3],
                    rownames(ttT),ExperimentID,stringsAsFactors=F)
  colnames(FIN) <- c("Active","Total","FoldChange","Adj.Pvalue",
                     "SampleType","ExperimentID")
  FIN <- FIN[order(as.numeric(FIN[,"Adj.Pvalue"],
                   -1*as.numeric(FIN[,"Active"]),decreasing=FALSE)),]
  FIN <- FIN[as.numeric(FIN[,"Adj.Pvalue"]) <= Pval.co,]
  FIN[,4] <- signif(FIN[,4],3)

 
      if(!is.null(directory)) {
          expiddex <- unique(unlist(strsplit(as.character(FIN$ExperimentID),";")))
          for(k in 1:length(expiddex)){
              filepath <- paste0(directory,"/",expiddex[k])
              Temp <- tabSearch(expiddex[k],chipdata)
              if(nrow(Temp) > 1){
            dir.create(filepath)
              GSCAeda(genedata,pattern,chipdata=chipdata,
                         SearchOutput=Temp,Pval.co=Pval.co,
                         Ordering="Average",Title=expiddex[k],
                       outputdir=filepath)
              }
          }
      }
  
  if(is.null(dim(FIN)) | nrow(FIN)==0) {
      message("No significant biological contexts found.")
  } else {
      FIN <- cbind(1:nrow(FIN),FIN)
      colnames(FIN)[1] <- "Rank"
      rownames(FIN) <- 1:nrow(FIN)
  }
      colnames(activity) <- tabsamplename
      
  return(list(Ranking=FIN,Score=activity,Pattern=pattern,Cutoff=genesetcutoff,SelectedSample=selectsample,Totalgene=genesettotalgenenum,Missinggene=genesetmissinggene,Chipdata=chipdata))

  } else {
      stop("No samples show the pattern of interest.
       Try relaxing cutoffs.")
  }
}
