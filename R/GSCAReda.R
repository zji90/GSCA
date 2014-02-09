GSCAReda <-
      function(genedata,pattern,chipdata,SearchOutput,Pval.co=0.05,Ordering="Average",Title=NULL,outputdir=NULL) {
            if(nrow(SearchOutput)==0) stop("No SearchOutput specified.")
            require(ggplot2)
            require(reshape2)
            require(RColorBrewer)
            ###presettings
            genedata[,1] <- as.character(genedata[,1])
            pattern[,1] <- as.character(pattern[,1])
            if (!is.null(outputdir)) {
                  if (!grepl("/$",outputdir))
                        outputdir <- paste0(outputdir,"/")
                  plotfilepath <- paste0(outputdir,"GSCAReda_plot.pdf")
                  datafilepath <- paste0(outputdir,"GSCAReda_data.csv")
                  resfilepath <- paste0(outputdir,"GSCAReda_result.csv")
            }
            
            #setting path of package GSCARdata
            path <- system.file("extdata",package=paste0("Affy",chipdata,"Expr"))
            #read in gene names
            compengene <- sub(".rda","",list.files(path))
            #read in reference tables for given compendium
            if(chipdata == "hgu133a") {
                  if (!require(Affyhgu133aExpr)) {
                        stop("Affyhgu133aExpr Package is not found")
                  } else {
                        data(Affyhgu133aExprtab)
                        tab <- Affyhgu133aExprtab
                  }
            } else if(chipdata == "moe430"){
                  if (!require(Affyhgu133aExpr)) {
                        stop("Affymoe430Expr Package is not found")
                  } else {
                        data(Affymoe430Exprtab)
                        tab <- Affymoe430Exprtab
                  }
            } else {
                  stop("Please enter valid name for chipdata. Current Supported chipdata: 'hgu133a' and 'moe430'")
            }
            tabsamplename <- tab$SampleName
            tab$SampleType <- substr(tab$SampleType,1,25)
            SearchOutput$SampleType <- substr(SearchOutput$SampleType,1,25)
            #check whether each geneset has at least one gene on the compendium
            genesetname <- NULL
            for (tmpgenesetname in unique(genedata[,1])) {
                  if(sum(compengene %in% genedata[genedata[,1]==tmpgenesetname,2]) == 0) {
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
            genesetcutoff <- genesettotalTG <- genesetmissingTG <- rep(0,length(genesetname))
            names(genesetcutoff) <- names(genesettotalTG) <- names(genesetmissingTG)  <- genesetname
            for (genesetid in 1:length(genesetname)) {
                  ###Scoring geneset activity
                  singlegeneset <- genesetname[genesetid]
                  posgene <- intersect(compengene,genedata[singlegeneset == genedata[,1]&genedata[,3]==1,2])
                  neggene <- intersect(compengene,genedata[singlegeneset == genedata[,1]&genedata[,3]==-1,2])
                  n.GSup <- length(posgene)
                  GSup <- rep(0, nrow(tab))  
                  for (i in posgene) {
                        load(paste0(path,"/",i,".rda"))
                        GSup <- GSup + Expr
                  }
                  GSup <- GSup / n.GSup
                  n.GSdown <- length(neggene)
                  GSdown <- rep(0, nrow(tab))  
                  for (i in neggene) {
                        load(paste0(path,"/",i,".rda"))
                        GSdown <- GSdown + Expr
                  }
                  GSdown <- GSdown / n.GSdown
                  missinggene <- setdiff(genedata[singlegeneset == genedata[,1],2],compengene)
                  genesetmissingTG[genesetid] <- length(missinggene)
                  genesettotalTG[genesetid] <- n.gene <- length(posgene)+length(neggene)
                  if(length(posgene) > 0 & length(neggene) > 0) {
                        score <- GSup*n.GSup/n.gene-GSdown*n.GSdown/n.gene
                  } else if (length(posgene) > 0) {
                        score <- GSup
                  } else {
                        score <- -1*GSdown
                  }
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
                              cutoff <- quantile(1-singlepattern[,4],mean(score),sd(score))
                              selectsample <- intersect(selectsample,which(score >= cutoff))
                        } else if (singlepattern[,2] == "Low") {
                              cutoff <- quantile(singlepattern[,4],mean(score),sd(score))
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
            
            colnames(activity) <- tabsamplename
            
            ExpID <- tab[selectsample,"ExperimentID"]
            tmpTypes <- tab[selectsample,"SampleType"]           
            
            #### Plotting
            if (!is.null(outputdir)) pdf(plotfilepath)
            
            GSE <- unique(SearchOutput[,"ExperimentID"])
            GSE2 <- unique(unlist(strsplit(SearchOutput[,"ExperimentID"],";")))
            SampleNames <- unique(SearchOutput[,"SampleType"])
            
            #### Expression Plotting
            if(sum(tab$ExperimentID %in% GSE |
                         tab$ExperimentID %in% GSE2 |
                         tab$SampleType %in% SampleNames)==0)
                  stop("No matching samples found.")
            
            tabdex <- tab[(tab$ExperimentID %in% GSE |
                                 tab$ExperimentID %in% GSE2) &
                                tab$SampleType %in% SampleNames,]
            
            eact <- activity[,(tab$ExperimentID %in% GSE |
                                     tab$ExperimentID %in% GSE2) &
                                   tab$SampleType %in% SampleNames]
            
            ExpressionTab <- NULL
            ExpressionTab <- tab[tab$SampleName %in% colnames(eact),]
            ExpressionTab <- data.frame(ExpressionTab,t(eact))
            #### Choice of ordering (Default is average rank)
            if (length(unique(ExpressionTab$SampleType)) == 1) {
                  Expressionorder <- unique(ExpressionTab$SampleType)
            } else {
            if(Ordering %in% genesetname) {
                  Expressionorder <- names(sort(tapply(ExpressionTab[,Ordering],ExpressionTab$SampleType,mean),decreasing=T))
            } else if(Ordering=="Average") {
                  Expressionorder <- rep(0,length(unique(ExpressionTab$SampleType)))
                  for (i in genesetname)
                        Expressionorder <- Expressionorder + scale(tapply(ExpressionTab[,i],ExpressionTab$SampleType,mean))
                  Expressionorder <- names(sort(Expressionorder[,1],decreasing=T))
            } else {
                  stop("No ordering specified. Must be one of the geneset name or 'Average'.")
            }
            }
            Expressionmelt <- melt(ExpressionTab[,-(1:3)],id.vars="SampleType")
            Expressionmelt$SampleType <- factor(Expressionmelt$SampleType,levels=Expressionorder)
            gboxplot <- ggplot(data = Expressionmelt, aes(x=variable, y=value)) + 
                  geom_boxplot(aes(fill=SampleType))+facet_grid(.~SampleType) +
                  labs(title=Title) + 
                  theme_bw() + theme(axis.text.x=element_text(angle=90)) +
                  theme(axis.title.x = element_blank(), axis.title.y=element_blank())
            print(gboxplot)
            
            #### Make Heatmap Plot
            
            if(is.null(outputdir)) par(ask = TRUE)
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
            alltmat <- NULL
            allpmat <- NULL
            for (geneset in genesetname) {
                  pmat <- tmat <- matrix(0,length(Expressionorder),length(Expressionorder))
                  colnames(pmat) <- colnames(tmat) <- rownames(pmat) <- rownames(tmat) <- Expressionorder
                  for (i in 1:length(Expressionorder))
                        for (j in 1:length(Expressionorder)) {
                              if(length(ExpressionTab[ExpressionTab$SampleType==Expressionorder[i],geneset])==1|
                                       length(ExpressionTab[ExpressionTab$SampleType==Expressionorder[j],geneset])==1) {
                                    tmat[i,j] <- NA
                                    pmat[i,j] <- NA
                              } else {
                                    ttestres <- t.test(ExpressionTab[ExpressionTab$SampleType==Expressionorder[i],geneset],ExpressionTab[ExpressionTab$SampleType==Expressionorder[j],geneset])
                                    tmat[i,j] <- ttestres$statistic
                                    pmat[i,j] <- ttestres$p.value
                              }
                        }
                  eval(parse(text=paste0("alltmat$",geneset, "<- tmat")))
                  eval(parse(text=paste0("allpmat$",geneset, "<- pmat")))
                  
                  #T test score heatmap
                  if (!all(is.na(tmat))) {
                  tmatdata <- melt(tmat,value.name="t.stat")
                  tmatdata$Var1 <- factor(tmatdata$Var1,levels=Expressionorder)
                  tmatdata$Var2 <- factor(tmatdata$Var2,levels=Expressionorder)
                  tplot <- ggplot(tmatdata,aes(x = Var1, y = Var2, fill = t.stat)) + 
                        geom_tile() + geom_abline(intercept=0,slope=1,size=1.25) + 
                        scale_fill_gradientn(colours = myPalette(10)) + 
                        scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
                        coord_equal() + labs(title=paste(geneset,"t.stat heatmap")) + 
                        theme_bw() + theme(axis.text.x=element_text(angle=90)) +
                        theme(axis.title.x = element_blank(), axis.title.y=element_blank())
                  print(tplot)
                  
                  #P value Heatmap
                  pmatdata <- melt(pmat,value.name="P.value")
                  pmatdata$Var1 <- factor(pmatdata$Var1,levels=Expressionorder)
                  pmatdata$Var2 <- factor(pmatdata$Var2,levels=Expressionorder)
                  pplot <- ggplot(pmatdata,aes(x = Var1, y = Var2, fill = P.value)) + 
                        geom_tile() + geom_abline(intercept=0,slope=1,size=1.25) + 
                        scale_fill_gradientn(colours = myPalette(10),breaks=c(0.001,0.01,0.05,0.1,0.25,0.5,0.75,1)) + 
                        guides(fill = guide_legend(title.position = "top")) +
                        scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
                        coord_equal() + labs(title=paste(geneset,"p.value heatmap")) + 
                        theme_bw() + theme(axis.text.x=element_text(angle=90)) +
                        theme(axis.title.x = element_blank(), axis.title.y=element_blank())
                  print(pplot)
                  }
            }
            
            #### GSCAR Results with given contexts
            ttT <- table(tmpTypes)
            sT <- sum(ttT)
            bgT <- table(tab$SampleType)
            bgT <- bgT[names(ttT)]
            SCORE <- matrix(0, nrow = length(Expressionorder), ncol = 4)
            ExperimentID <- NULL; FIN <- NULL
            pvaldex <- NULL; fold.change <- NULL; active <- NULL; total <- NULL
            for (i in Expressionorder) {
                  if(sum(names(ttT) %in% i) > 0 &
                           sum(tab$SampleType==i) > 2){
                        r1c1 <- ttT[names(ttT) %in% i]
                        r1c2 <- sT - ttT[names(ttT) %in% i]
                        r2c1 <- bgT[names(ttT) %in% i] - ttT[names(ttT) %in% i]
                        r2c2 <- length(tab$SampleType) - r1c1 - r1c2 - r2c1
                        tmpmat <- matrix(c(r1c1, r2c1, r1c2, r2c2), ncol = 2)
                        pval <- fisher.test(tmpmat, alternative = "greater")$p.value
                        adj.pval <- min(pval * sum(table(tab$SampleType) > 2),1)
                        adj.pval <- signif(adj.pval,3)
                        active <- c(active,as.numeric(ttT)[names(ttT) %in% i])
                        total <- c(total,as.numeric(bgT)[names(bgT) %in% i])
                        fold.change <- c(fold.change,round(
                              ((as.numeric(ttT)[names(ttT) %in% i]+
                                      sT/length(tab$SampleType)) /
                                     (as.numeric(bgT)[names(bgT) %in% i]+1)) /
                                    (sT/length(tab$SampleType)),3)
                        )
                        pvaldex <- c(pvaldex,adj.pval)
                        ExperimentID <- c(ExperimentID,
                                          paste(unique(unlist(strsplit(paste(
                                                tab$ExperimentID[tab$SampleType %in%
                                                                       i], collapse = ";"), ";"))),
                                                collapse = ";"))
                  } else {
                        active <- c(active,0)
                        total <- c(total,sum(tab$SampleType==i))
                        if(sum(tab$SampleType==i) > 2){
                              pvaldex <- c(pvaldex,1)
                        } else { pvaldex <- c(pvaldex,NA) }
                        fold.change <- c(fold.change,1)
                        ExperimentID <- c(ExperimentID,paste(unique(
                              unlist(strsplit(
                                    paste(tab$ExperimentID[tab$SampleType %in%
                                                                 i], collapse = ";"), ";"))),
                                                             collapse = ";"))
                  }
            }
            
            if(!is.null(outputdir)) dev.off()
            if(is.null(outputdir)) par(ask = FALSE)
            
            FIN <- data.frame(active,total,fold.change,pvaldex,Expressionorder,ExperimentID)
            colnames(FIN) <- c("Active","Total","FoldChange","Adj.Pvalue","SampleType","ExperimentID")
            
            ExpressionSum <- NULL
            for(i in Expressionorder){
                  tmpdat <- ExpressionTab[ExpressionTab$SampleType==i,]
                  tmpsummary <- sum(ExpressionTab$SampleType==i)
                  for (j in genedata[,1]) {
                        if (j %in% colnames(tmpdat))
                              tmpsummary <- c(tmpsummary,mean(tmpdat[,j]),sd(tmpdat[,j]))
                  }
                  ExpressionSum <- rbind(ExpressionSum,tmpsummary)
            }
            rownames(ExpressionSum) <- Expressionorder
            expplotcolname <- c("Sample Count")
            for (j in genedata[,1]) {
                  if (j %in% genesetname)
                        expplotcolname <- c(expplotcolname,paste0("Mean_",j),paste0("Sd_",j))
            }
            colnames(ExpressionSum) <- expplotcolname
            
            
            if(!is.null(outputdir)) {
                  write.csv(ExpressionTab,file=datafilepath,row.names=F)
                  write.table(FIN,file=resfilepath,row.names=F,sep=",")
                  cat("\n Summary", file = resfilepath, append = TRUE)
                  write.table(rbind(colnames(ExpressionSum),ExpressionSum),file=resfilepath,row.names=T,sep=",",col.names=F,append=T)
                  for (i in genesetname) {
                        cat(paste(i,"T statistics"), file = resfilepath, append = TRUE)
                        tmpout <- eval(parse(text=paste0("alltmat$",i)))
                        write.table(rbind(colnames(tmpout),tmpout),file=resfilepath,row.names=T,sep=",",col.names=F,append=T)
                        cat(paste(i,"P-values"), file = resfilepath, append = TRUE)
                        tmpout <- eval(parse(text=paste0("allpmat$",i)))
                        write.table(rbind(colnames(tmpout),tmpout),file=resfilepath,row.names=T,sep=",",col.names=F,append=T)     
                  }
                        
            } else {
                  return(list(Ranking=FIN,Data=ExpressionTab,
                              Summary=ExpressionSum,Tstats=alltmat,
                              Pval=allpmat))
            }
      }

