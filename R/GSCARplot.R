GSCARplot <- function(GSCARoutput,N=5,plotfile=NULL,Title=NULL){
      
      require(gplots)
      
      Ranking <- GSCARoutput$Ranking
      Score <- GSCARoutput$Score
      Cutoff <- GSCARoutput$Cutoff
      selectsample <- GSCARoutput$SelectedSample
      chipdata <- GSCARoutput$Chipdata
      STYLES <- c(15:18,3)
      COLORS <- brewer.pal(5,"Set1")

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
            stop("Please enter valid name for chipdata. Current Supported name: 'hgu133a' and 'moe430'")
      }
      
      if(!is.null(plotfile)) 
            pdf(plotfile)

      if (nrow(Score) == 1) {
            ###single geneset plot
            if(!(N %in% 1:5))
                  N <- min(5,nrow(Ranking))
            par(mfrow=c(N+1,1),oma=c(0,0,2,0))
            hist(Score,xlab="Sample Score",xlim=range(Score),main="All Biological contexts")
            abline(v=Cutoff,lty=2)
            for(i in 1:N) {
                  INDEX <- Ranking[i,"SampleType"]
                  hist(Score[tab$SampleType %in% INDEX],xlim=range(Score),main=substr(Ranking[i,"SampleType"],1,25),xlab="Sample Score")
                  abline(v=Cutoff,lty=2)
            }
            title(Title,outer=TRUE)
            par(mfrow=c(1,1),oma=c(0,0,0,0))
      } else if (nrow(Score) == 2) {
            plot(Score[1,],Score[2,],xlab=rownames(Score)[1],ylab=rownames(Score)[2],col="#00000022",pch=20,main=Title)
            if(!(N %in% 1:5)) 
                  N <- min(5,nrow(Ranking))
            toprankingsample <- NULL
            for(i in 1:N) {
                  INDEX <- Ranking[i,"SampleType"]
                  points(Score[1,][tab$SampleType %in% INDEX],Score[2,][tab$SampleType %in% INDEX],
                        col=COLORS[i],pch=STYLES[i],cex=0.75,bg=COLORS[i])
                  toprankingsample <- union(toprankingsample,which(tab$SampleType %in% INDEX))
            }
            abline(v=Cutoff[1], h=Cutoff[2], lty=2)
            points(Score[1,setdiff(selectsample,toprankingsample)],Score[2,setdiff(selectsample,toprankingsample)],pch=20)
            leg.txt <- c(substr(Ranking[1:N,"SampleType"],1,25),"Selected Samples","Not Selected Samples")
            legend("topleft",legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=0.8)
      } else {
            colcolorall <- rep("cyan",ncol(Score))
            colcolorall[selectsample] <- "blue"
            heatmap.2(Score,col=bluered,dendrogram="none",trace="none",Rowv=F,labCol=NA,ColSideColors=colcolorall,main="Score Heatmap for all Samples")
            legend("bottomleft",legend=c("Selected Sample","Unselected Sample"),lwd=1,col=c("blue","cyan"))
            title(Title,outer=TRUE)
            if (is.null(plotfile))
                  par(ask=T)
            colcolorselect <- rep("white",ncol(Score))
            if(!(N %in% 1:5)) 
                  N <- min(5,nrow(Ranking))
            for(i in 1:N) {
                  INDEX <- Ranking[i,"SampleType"]
                  colcolorselect[tab$SampleType %in% INDEX] <- COLORS[i]
            }
            leg.txt <- substr(Ranking[1:N,"SampleType"],1,25)
            heatmap.2(Score[,selectsample],col=bluered,labCol=NA,Rowv=F,dendrogram="none",trace="none",ColSideColors=colcolorselect[selectsample],main="Score Heatmap for selected Samples")
            legend("bottomleft",legend=leg.txt,lwd=1,col=COLORS,cex=0.8)
            if (is.null(plotfile))
                  par(ask=F)
      }
      
      if(!is.null(plotfile)) 
            dev.off()
}


