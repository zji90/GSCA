######################################################
##           GSCA:Gene Set Context Analysis         ##
##             Interactive User Interface           ##
##                     Server File                  ##
##           Author:Zhicheng Ji, Hongkai Ji         ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################

require(gplots)
require(shiny)
require(sp)
require(RColorBrewer)
###two genedata case initiate
polycord <- NULL
resetvalue <- 0
finishpolygonvalue <- 0
addpolygonvalue <- 0
undovalue <- 0
actiontaken <- 0
finishpolygonaction <- 0
polynum <- 1

#three genedata case sample selection initiate
threesampleslidervalue <- matrix(0,nrow=5,ncol=3)

###default styles and colors
STYLES <- rep(c(15:18,3:4),5)
COLORS <- c(brewer.pal(5,"Set1"),brewer.pal(5,"Pastel1"),brewer.pal(10,"Paired"),brewer.pal(10,"BrBG"))


shinyServer(function(input, output, session) {
      
      ######  Mainmethod : Input  ######
      
      Currentdata <- reactiveValues()
      Rawdata <- reactiveValues()
      
      output$InputGenesetcolnumui <- renderUI({
            if (!is.null(Currentdata$rawGenesetFile) && !is.null(ncol(Currentdata$rawGenesetFile)) && ncol(Currentdata$rawGenesetFile) > 2)
                  checkboxInput("InputGenesetcolnumthree","Include an optional third column",value=T)
      })
      
      #Check and automatically modify genesetname if current genesetname is already used
      observe({
            currentgenesetname <- input$InputGenesetname
            if (!is.null(Rawdata$patterndata))
                  if (input$InputGenesetname %in% Rawdata$patterndata[,1]) {
                        i <- 2
                        modifygenesetname <- paste(currentgenesetname,i,sep="_")
                        while (modifygenesetname %in% Rawdata$patterndata[,1]) {
                              i <- i + 1
                              modifygenesetname <- paste(currentgenesetname,i,sep="_")
                        }
                        currentgenesetname <- modifygenesetname
                  }
            Currentdata$genesetname <- currentgenesetname
      })
      
      #Input Geneset: Controls the behavior of specifying the current genedata
      observe({          
            if (input$InputGenesetmethod == "InputspecifyGeneset") {  
                  if (!is.null(input$InputGenesetname) && input$InputGenesetname!="" & (input$InputActGeneID!="" | input$InputRepGeneID!="")) 
                        Currentdata$genedata <- data.frame(Genesetname=Currentdata$genesetname,GeneID=c(strsplit(input$InputActGeneID,";")[[1]],strsplit(input$InputRepGeneID,";")[[1]]),activated_repressed=c(rep(1,length(strsplit(input$InputActGeneID,";")[[1]])),rep(-1,length(strsplit(input$InputRepGeneID,";")[[1]]))),stringsAsFactors=F)
            } else {
                  GenesetFileHandle <- input$InputGenesetFile  
                  if (!is.null(GenesetFileHandle)) {
                        tmpfile <- read.table(GenesetFileHandle$datapath,header=input$InputGenesetheader,sep=input$InputGenesetsep,quote=input$InputGenesetquote,stringsAsFactors=F,blank.lines.skip=TRUE)
                        if (ncol(tmpfile)==1) {
                              Currentdata$genedata <- NULL
                              Currentdata$rawGenesetFile <- NULL
                        } else {
                              tmpfile <- tmpfile[,colSums(is.na(tmpfile)) != nrow(tmpfile)]
                              tmpfile <- tmpfile[rowSums(is.na(tmpfile)) != ncol(tmpfile),]
                              Currentdata$rawGenesetFile <- tmpfile[rowSums(tmpfile=="") != ncol(tmpfile),]
                              if (ncol(Currentdata$rawGenesetFile) > 2 && !is.null(input$InputGenesetcolnumthree) && input$InputGenesetcolnumthree) {
                                    Currentdata$GenesetFile <- Currentdata$rawGenesetFile[,1:3]
                                    Currentdata$genedata <- data.frame(Genesetname=Currentdata$GenesetFile[,3],GeneID=Currentdata$GenesetFile[,1],activated_repressed=Currentdata$GenesetFile[,2],stringsAsFactors=F)
                              } else if (ncol(Currentdata$rawGenesetFile) == 2 || (ncol(Currentdata$rawGenesetFile) > 2 && !is.null(input$InputGenesetcolnumthree) && !input$InputGenesetcolnumthree)) {
                                    Currentdata$GenesetFile <- Currentdata$rawGenesetFile[,1:2]
                                    Currentdata$genedata <- data.frame(Genesetname=Currentdata$genesetname,GeneID=Currentdata$GenesetFile[,1],activated_repressed=Currentdata$GenesetFile[,2],stringsAsFactors=F)
                              }
                        }
                  }
            }
      })
      
      #Input Geneset: Controls the behavior of adding the current genedata to the whole genedata
      observe({   
            if(input$Inputgenesetadd != 0) {
                  isolate({
                        Rawdata$genedata <- rbind(Rawdata$genedata,Currentdata$genedata)
                        if (input$InputGenesetmethod == 'InputuploadGeneset' && ncol(Currentdata$rawGenesetFile) > 2 && !is.null(input$InputGenesetcolnumthree) && input$InputGenesetcolnumthree) {
                              Rawdata$patterndata <- rbind(Rawdata$patterndata,data.frame(Genesetname=unique(Currentdata$genedata$Genesetname),Activity=input$Inputgenesetpatternactivity,cotype=input$Inputgenesetpatterncotype,cutoff=input$Inputgenesetpatternco,stringsAsFactors=F))
                        } else {
                              Rawdata$patterndata <- rbind(Rawdata$patterndata,data.frame(Genesetname=Currentdata$genesetname,Activity=input$Inputgenesetpatternactivity,cotype=input$Inputgenesetpatterncotype,cutoff=input$Inputgenesetpatternco,stringsAsFactors=F))
                        }
                  })
            }
      })
      
      output$OutputCurrentGenedatawarnui <- renderUI({
            if (input$InputGenesetmethod == 'InputuploadGeneset' && is.null(Currentdata$genedata))
                  helpText("Try changing the separator if no Genedata table appears")
      })
      output$OutputCurrentGenedata <- renderDataTable(Currentdata$genedata)
      output$OutputAllGenedata <- renderDataTable(Rawdata$genedata)    
      output$OutputAllPattern <- renderDataTable(Rawdata$patterndata)
      
      output$OutputGenedataname <- renderText({
            if (is.null(Rawdata$genedata)){
                  "No input genedata in GSCA"
            } else {
                  do.call("paste",c(as.list(unique(Rawdata$genedata$Genesetname)),sep=","))
            }
      })
      
      #Delete genedata
      output$Inputgenesetdeleteui <- renderUI({
            if (!is.null(Rawdata$genedata))
                  selectInput("Inputdeletegenesetselect","Hold Control(Win) or Command(Mac) to enable multiple selection",unique(Rawdata$genedata$Genesetname),multiple=T)
      })
      
      observe({   
            if(input$Inputgenesetdelete != 0) {
                  isolate({
                        Rawdata$genedata <- Rawdata$genedata[!Rawdata$genedata$Genesetname %in% input$Inputdeletegenesetselect,]
                        Rawdata$patterndata <- Rawdata$patterndata[!Rawdata$patterndata$Genesetname %in% input$Inputdeletegenesetselect,]
                  })
            }
      })      
      
      observe({
            input$Inputgenesetreset
            isolate({
                  Rawdata$genedata <- NULL
                  Rawdata$patterndata <- NULL
            })
            
      })
      
      ######  Mainmethod : Summary  ######
      
      Maindata <- reactiveValues()
      
      #Sidebar checkbox group for selecting dataset
      output$Summarydataselect <- renderUI({
            if (!is.null(Rawdata$patterndata)) {
                  checkboxGroupInput("Selectedgeneset","Select Genesets",Rawdata$patterndata[,1],selected=Rawdata$patterndata[,1])
            }
      })
      
      #UI for selecting compendium
      output$Summarycompselectui <- renderUI({
            complist <- list()
            if (require(Affyhgu133aExpr)) {
                  complist <- c(complist,"Affymetrix Human hgu133a Array (GPL96)"="hgu133a")
            } 
            if (require(Affymoe4302Expr)) {
                  complist <- c(complist,"Affymetrix Mouse 430 2.0 Array (GPL1261)"="moe4302")
            }
            if (require(Affyhgu133Plus2Expr)) {
                  complist <- c(complist,"Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)"="hgu133Plus2")
            }
            if (require(Affyhgu133A2Expr)) {
                  complist <- c(complist,"Affymetrix Human Genome U133A 2.0 Array (GPL571)"="hgu133A2")
            }
            radioButtons("Summarycompselect","Select Compendium", complist)
      })
      
      output$Summarycompinfo <- renderUI({
            if (!is.null(input$Summarycompselect)) {
                  if(input$Summarycompselect=="moe4302"){
                        helpText("This compendium contains 9444 human profiles on 20630 genes")
                  } else if (input$Summarycompselect=="hgu133a"){
                        helpText("This compendium contains 11778 human profiles on 12495 genes")
                  } else if (input$Summarycompselect=="hgu133Plus2"){
                        helpText("This compendium contains 5153 human profiles on 19944 genes")
                  } else if (input$Summarycompselect=="hgu133A2"){
                        helpText("This compendium contains 313 human profiles on 12494 genes")
                  }    
            }
      })
      
      #Update Maindata information
      observe({   
            if (!is.null(input$Summarycompselect)) {
            if(input$Summarycompselect=="moe4302"){
                  data(Affymoe4302Exprtab)
                  Maindata$tab <- Affymoe4302Exprtab
            } else if (input$Summarycompselect=="hgu133a"){
                  data(Affyhgu133aExprtab)
                  Maindata$tab <- Affyhgu133aExprtab
            } else if (input$Summarycompselect=="hgu133Plus2"){
                  data(Affyhgu133Plus2Exprtab)
                  Maindata$tab <- Affyhgu133Plus2Exprtab
            } else if (input$Summarycompselect=="hgu133A2"){
                  data(Affyhgu133A2Exprtab)
                  Maindata$tab <- Affyhgu133A2Exprtab
            }
            if (!is.null(Rawdata$genedata)) {
                  path <- system.file("extdata",package=paste0("Affy",input$Summarycompselect,"Expr"))
                  compengene <- sub(".rda","",list.files(path))
                  Maindata$genedata <- Rawdata$genedata[Rawdata$genedata[,2] %in% compengene & Rawdata$genedata[,1] %in% input$Selectedgeneset,]
                  Maindata$patterndata <- Rawdata$patterndata[Rawdata$patterndata[,1] %in% Maindata$genedata[,1],]
                  Maindata$dim <- nrow(Maindata$patterndata)
            }       
            }
      })
      
      #Render summary of selected dataset
      output$OutputDataSummary <- renderDataTable({
            if (!is.null(Maindata$patterndata))
                  if(nrow(Maindata$patterndata) != 0) {
                        sumtable <- data.frame(Genesetname=Maindata$patterndata[,1],Originalgene=0,Compendiumgene=0)
                        for (i in 1:nrow(sumtable)) {
                              sumtable[i,2] <- sum(Rawdata$genedata[,1] == sumtable[i,1])
                              sumtable[i,3] <- sum(Maindata$genedata[,1] == sumtable[i,1])
                        }
                        sumtable <- merge(sumtable,Maindata$patterndata)
                        names(sumtable) <- c("Genesetname","Original Number of Genes","Number of Genes in Compendium","Cutoff Activity","Cutoff Type","Cutoff Value")
                        sumtable
                  }
      })    
      
      #Missing genedata report
      output$Outputmissinggenesetreport <- renderUI({
            if (!is.null(Maindata$dim))
                  if (Maindata$dim != length(input$Selectedgeneset)) {
                        tagList(
                              br(h4("Following genedata do not contain any gene in the selected compendium, thus omitted:")),
                              textOutput("Outputmissinggenesetreporttext")
                        )
                  }
      })
      
      output$Outputmissinggenesetreporttext <- renderText({
            do.call("paste",c(as.list(setdiff(input$Selectedgeneset,Maindata$patterndata$Genesetname)),sep=","))
      })
      
      ######  Mainmethod : GSCA  ######
      
      GSCA <- reactiveValues()
      
      #When reentering GSCA page, GSCAmethod should be set to GSCAdefault
      observe({
            if (input$Mainmethod == "Input" | input$Mainmethod == "Select")
                  updateRadioButtons(session,"GSCAmethod",choices=list("Default Enrichment Region Selection"="GSCAdefault","Interactive Enrichment Region Selection"="GSCAinteractive"),selected="Default Enrichment Region Selection")
      })
      
      #Function to calculate enriched samples
      enrichfunc <- function(Score,selectedsample,tab) {
            ExpID <- tab[selectedsample,"ExperimentID"]
            tmpTypes <- tab[selectedsample,"SampleType"]
            tabTWO <- table(tab$SampleType)
            tabTWO <- names(tabTWO)[tabTWO > 2]
            tab <- tab[tab$SampleType %in% tabTWO,]
            ExpID <- ExpID[tmpTypes %in% tab$SampleType]
            tmpTypes <- tmpTypes[tmpTypes %in% tab$SampleType]
            if(length(tmpTypes) == 0){
                  return(NULL)
            } else {
                  ttT <- table(tmpTypes)
                  sT <- sum(ttT)
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
                        SCORE[i,2] <- round(((as.numeric(ttT)[i]+sT/length(tab$SampleType)) / (as.numeric(bgT)[i]+1))/(sT/length(tab$SampleType)),3)
                        SCORE[i,3] <- min(SCORE[i,1]*ContextN,1)
                        ExperimentID[i] <- paste(unique(unlist(strsplit(paste(ExpID[tmpTypes==names(ttT)[i]],collapse=";"),";"))),collapse=";")
                  }      
                  FIN <- data.frame(as.numeric(ttT),as.numeric(bgT),SCORE[,2],SCORE[,3],rownames(ttT),ExperimentID,stringsAsFactors=F)
                  colnames(FIN) <- c("Active","Total","FoldChange","Adj.Pvalue","SampleType","ExperimentID")
                  FIN <- FIN[order(as.numeric(FIN[,"Adj.Pvalue"],-1*as.numeric(FIN[,"Active"]),decreasing=FALSE)),]
                  FIN[,4] <- signif(FIN[,4],3)
                  if(is.null(dim(FIN)) | nrow(FIN)==0) {
                        message("No significant biological contexts found.")
                  } else {
                        FIN <- cbind(1:nrow(FIN),FIN)
                        colnames(FIN)[1] <- "Rank"
                        rownames(FIN) <- 1:nrow(FIN)
                  }
            }
            return(FIN)
      }      
      
      #Calculating GSCA Score, default selected samples, default Ranking, also hclust dendrogram if there are geq three genesets
      observe({
            if (input$Mainmethod == 'GSCA' || input$Mainmethod == 'Download' || input$Mainmethod == 'about') {
                  isolate({
                        if (!is.null(Maindata$patterndata))
                              if(nrow(Maindata$patterndata) != 0) {
                                    waiting$defaultstatus <- 1
                                    scoremat <- matrix(0,nrow=Maindata$dim,ncol=nrow(Maindata$tab))
                                    row.names(scoremat) <- rep("tmp",Maindata$dim)
                                    selectsample <- 1:nrow(Maindata$tab)
                                    cutoffval <- rep(0,Maindata$dim)
                                    for (genesetid in 1:Maindata$dim) {
                                          ###Calculate Score
                                          singlegeneset <- Maindata$patterndata[genesetid,1]
                                          posgene <- Maindata$genedata[Maindata$genedata[,1] == singlegeneset & Maindata$genedata[,3] == "1",2]
                                          neggene <- Maindata$genedata[Maindata$genedata[,1] == singlegeneset & Maindata$genedata[,3] == "-1",2]
                                          n.GSup <- length(posgene)
                                          GSup <- rep(0, nrow(Maindata$tab))
                                          for (i in posgene) {
                                                path <- system.file("extdata",package=paste0("Affy",input$Summarycompselect,"Expr"))
                                                load(paste0(path,"/",i,".rda"))
                                                GSup <- GSup + e
                                          }
                                          GSup <- GSup / n.GSup
                                          n.GSdown <- length(neggene)
                                          GSdown <- rep(0, nrow(Maindata$tab))  
                                          for (i in neggene) {
                                                path <- system.file("extdata",package=paste0("Affy",input$Summarycompselect,"Expr"))
                                                load(paste0(path,"/",i,".rda"))
                                                GSdown <- GSdown + e
                                          }
                                          GSdown <- GSdown / n.GSdown
                                          if(length(posgene) > 0 & length(neggene) > 0) {
                                                score <- GSup*n.GSup/(length(posgene)+length(neggene))-GSdown*n.GSdown/(length(posgene)+length(neggene))
                                          } else if (length(posgene) > 0) {
                                                score <- GSup
                                          } else {
                                                score <- -1*GSdown
                                          }
                                          scoremat[genesetid,] <- score
                                          row.names(scoremat)[genesetid] <- singlegeneset
                                          ###Calculate selected samples                   
                                          if (Maindata$patterndata[genesetid,3] == "Norm") {
                                                if (Maindata$patterndata[genesetid,2] == "High") {       
                                                      cutoff <- qnorm(1-as.numeric(Maindata$patterndata[genesetid,4]),mean(score),sd(score))
                                                      selectsample <- intersect(selectsample,which(score >= cutoff))
                                                      cutoffval[genesetid] <- cutoff
                                                } else {
                                                      cutoff <- qnorm(as.numeric(Maindata$patterndata[genesetid,4]),mean(score),sd(score))
                                                      selectsample <- intersect(selectsample,which(score < cutoff))
                                                      cutoffval[genesetid] <- cutoff
                                                }
                                          } else if (Maindata$patterndata[genesetid,3] == "Quantile") {
                                                if (Maindata$patterndata[genesetid,2] == "High") {
                                                      cutoff <- quantile(1-as.numeric(Maindata$patterndata[genesetid,4]),mean(score),sd(score))
                                                      selectsample <- intersect(selectsample,which(score >= cutoff))
                                                      cutoffval[genesetid] <- cutoff
                                                } else {
                                                      cutoff <- quantile(as.numeric(Maindata$patterndata[genesetid,4]),mean(score),sd(score))
                                                      selectsample <- intersect(selectsample,which(score < cutoff))
                                                      cutoffval[genesetid] <- cutoff
                                                }
                                          } else if (Maindata$patterndata[genesetid,3] == "Exprs") {
                                                if (Maindata$patterndata[genesetid,2] == "High") {
                                                      cutoff <- as.numeric(Maindata$patterndata[genesetid,4])
                                                      selectsample <- intersect(selectsample,which(score >= cutoff))
                                                      cutoffval[genesetid] <- cutoff
                                                } else {
                                                      cutoff <- as.numeric(Maindata$patterndata[genesetid,4])
                                                      selectsample <- intersect(selectsample,which(score < cutoff))
                                                      cutoffval[genesetid] <- cutoff
                                                }
                                          }
                                    }
                                    
                                    Maindata$GSCAscore <- scoremat
                                    Maindata$defaultsample <- selectsample
                                    Maindata$cutoffval <- cutoffval
                                    if (Maindata$dim >= 3) {
                                          Maindata$GSCAclust <- hclust(dist(t(scoremat)))
                                          Rowv <- rowMeans(Maindata$GSCAscore)
                                          hcr <- hclust(dist(Maindata$GSCAscore))
                                          ddr <- as.dendrogram(hcr)
                                          Maindata$GSCArowclust <- reorder(ddr, Rowv)
                                    }
                                    if (threesampleslidervalue[1,3] != ncol(Maindata$GSCAscore)) {
                                          threesampleslidervalue[,3] <<- rep(ncol(Maindata$GSCAscore),5)
                                          threesampleslidervalue[,1] <<- rep(0,5)
                                          threesampleslidervalue[,2] <<- rep(ncol(Maindata$GSCAscore),5)
                                    }
                                    if (Maindata$dim == 2) {
                                          Maindata$twocorr <- cor(Maindata$GSCAscore[1,],Maindata$GSCAscore[2,])
                                          Maindata$twocorrp <- cor.test(Maindata$GSCAscore[1,],Maindata$GSCAscore[2,])$p.value
                                          twolm <- summary(lm(Maindata$GSCAscore[2,]~Maindata$GSCAscore[1,]))
                                          Maindata$twoslope <- twolm$coefficients[2,1]
                                          Maindata$twoslopep <- twolm$coefficients[2,4]
                                    }
                              }
                  })
            }
      })
      
      #GSCA Ranking
      observe({
            Maindata$AllRanking <- enrichfunc(Maindata$GSCAscore,Maindata$selectsample,Maindata$tab)
            Maindata$Ranking <- Maindata$AllRanking[as.numeric(Maindata$AllRanking[,"Adj.Pvalue"]) <= as.numeric(input$Inputpvalco) & as.numeric(Maindata$AllRanking[,"FoldChange"]) >= as.numeric(input$Inputfoldchangeco),]
            if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAdefault') {
                  Maindata$defaultranking <- Maindata$Ranking
            } else if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAinteractive'){
                  Maindata$interactranking <- Maindata$Ranking
            }    
      })
      
      #GSCA sidebar
      output$InputNslider <- renderUI({
            maxval <- min(30,nrow(Maindata$Ranking))
            defaultval <- min(5,nrow(Maindata$Ranking))
            if (maxval > 0)
                  sliderInput("InputN","Number of Top Ranking Context Displayed",min=0,max=maxval,value=defaultval,step=1)
      })
      
      output$InputGSCAsidebar <- renderUI({
            if (!is.null(Maindata$patterndata)) 
                  if(nrow(Maindata$patterndata) != 0) {
                        if (Maindata$dim == 1) {
                              sliderInput(inputId = "GSCAoneslider", label = Maindata$patterndata[,1], min = min(Maindata$GSCAscore), max = max(Maindata$GSCAscore), value = c(min(Maindata$GSCAscore), max(Maindata$GSCAscore)))
                        } else if (Maindata$dim == 2) {
                              tagList(
                                    helpText("Build Polygon to Specify Interested Region"),
                                    p(actionButton('addpolygon', 'Add New Polygon')),
                                    p(actionButton('finishpolygon', 'Finish Drawing Polygon')),
                                    p(actionButton('undo', 'Undo last Operation')),
                                    p(actionButton('reset', 'Reset'))
                              )
                        } else {             
                              tagList(
                                    checkboxInput("Threecutoffvalue","Select sample by exact value"),
                                      conditionalPanel(condition="input.Threecutoffvalue == 1",
                                                       lapply(1:Maindata$dim, function(i) {
                                                             sliderInput(inputId = paste0("GSCAthreevalueslider",i), label = Maindata$patterndata[i,1], min = min(Maindata$GSCAscore[i,]), max = max(Maindata$GSCAscore[i,]), value = c(min(Maindata$GSCAscore[i,]), max(Maindata$GSCAscore[i,])))
                                                       })
                                      )
                              )
                        }
                  }
      })
      
      #Checkbox whether to display enriched context in all area
      output$plotenrichedareaui <- renderUI({
            if (!is.null(Maindata$dim))
                  if (Maindata$dim == 2)
                        checkboxInput("Inputenrichedareaonly","Show Enriched Context only in interested region")
      })
      
      #Checkbox whether to cluster on rows for heatmaps
      output$heatmapthreerowvui <- renderUI({
            if (!is.null(Maindata$dim))
                  if (Maindata$dim > 2)
                        checkboxInput("heatmapthreerowv","Cluster on rows")
      })
      
      #Specify Context
      output$InputGSCAspecifycontextui <- renderUI({
            if (!is.null(Maindata$Ranking))
                  selectInput("GSCAspecifycontext","Hold Control(Win) or Command(Mac) to enable multiple selection",Maindata$Ranking[,"SampleType"],multiple=T)
      })
      
      #Establish GSCA Context
      observe({
            if(!is.null(Maindata$Ranking)) {
                  if (input$Inputcontexttype=='Toprank') {
                        if (nrow(Maindata$Ranking) == 0) {
                              Maindata$GSCAcontext <- NULL
                        } else {
                              Maindata$GSCAcontext <- Maindata$Ranking[1:as.numeric(input$InputN),"SampleType"]
                        }
                  } else {
                        Maindata$GSCAcontext <- input$GSCAspecifycontext 
                  }
            }
            
            if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAdefault') {
                  Maindata$defaultcontext <- Maindata$GSCAcontext
            } else if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAinteractive'){
                  Maindata$interactcontext <- Maindata$GSCAcontext
            }
      })
      
      #Establish Maindata$selectsample
      observe({
            if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAinteractive' & !is.null(Maindata$dim)) {
                  if (Maindata$dim == 1) {
                        Maindata$selectsample <- which(Maindata$GSCAscore > input$GSCAoneslider[1] & Maindata$GSCAscore < input$GSCAoneslider[2])
                  } else if (Maindata$dim == 2) {
                        Maindata$selectsample <- inpoly$tf
                  } else {
                        if (!is.null(input$GSCAinteractivethreeupdate)) 
                              if (input$GSCAinteractivethreeupdate>0) {
                                    isolate({
                                          Maindata$selectsample <- GSCAthreeinfo$selectsample
                                    })
                              }
                  }
                  Maindata$interactsample <- Maindata$selectsample
            } else {
                  Maindata$selectsample <- Maindata$defaultsample
            }
      })
      
      #GSCA default plot
      output$GSCAdefaultplot <- renderUI({
            if (!is.null(Maindata$patterndata)) 
                  if(nrow(Maindata$patterndata) != 0) {
                        if (Maindata$dim == 1) {
                              tagList(
                                    helpText(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                                    plotOutput("GSCAdefaultplotone",height=300*(length(Maindata$GSCAcontext)+1))
                              )
                        } else if (Maindata$dim == 2) {
                              div(align="center",
                                  helpText(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                                  plotOutput("GSCAdefaultplottwo",width=800,height=800),
                                  helpText(paste0("Correlation: ",Maindata$twocorr)),
                                  helpText(paste0("Pearson Correlation Test p-value: ",Maindata$twocorrp)),
                                  helpText(paste0("Regression Slope: ",Maindata$twoslope)),
                                  helpText(paste0("Regression Slope t-test p-value: ",Maindata$twoslopep))
                              )
                        } else {
                              tagList(
                                    plotOutput("GSCAdefaultplotthree"),
                                    helpText(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                                    plotOutput("GSCAdefaultplotthreeplus")
                              )
                        }
                  }      
      })
      
      output$GSCAdefaultplotone <- renderPlot({
            if (Maindata$dim == 1) {
                  par(mfrow=c(length(Maindata$GSCAcontext)+1,1),oma=c(0,0,2,0))
                  hist(Maindata$GSCAscore,xlab="Sample Score",xlim=range(Maindata$GSCAscore),main="All Biological contexts")
                  abline(v=Maindata$cutoffval[1], lty=2)
                  if (!is.null(Maindata$GSCAcontext))
                        for(INDEX in Maindata$GSCAcontext) {
                              hist(Maindata$GSCAscore[Maindata$tab$SampleType %in% INDEX],xlim=range(Maindata$GSCAscore),main=substr(INDEX,1,25),xlab="Sample Score")
                              abline(v=Maindata$cutoffval[1], lty=2)
                        }     
            }
      })
      
      output$GSCAdefaultplottwo <- renderPlot({  
            if (Maindata$dim == 2) {
                  plot(Maindata$GSCAscore[1,],Maindata$GSCAscore[2,],col="#00000022",pch=20,xlab=Maindata$patterndata[1,1],cex=0.7,ylab=Maindata$patterndata[2,1])
                  toprankingsample <- NULL
                  if (!is.null(Maindata$GSCAcontext)) {
                        for(INDEX in Maindata$GSCAcontext) {
                              toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                        }
                  }
                  points(Maindata$GSCAscore[1,setdiff(Maindata$selectsample,toprankingsample)],Maindata$GSCAscore[2,setdiff(Maindata$selectsample,toprankingsample)],cex=0.7,pch=20)
                  abline(v=Maindata$cutoffval[1], h=Maindata$cutoffval[2], lty=2,lwd=2)
                  if (!is.null(Maindata$GSCAcontext)) {
                        i <- 1
                        for(INDEX in Maindata$GSCAcontext) {
                              if (input$Inputenrichedareaonly == TRUE) {
                                    points(Maindata$GSCAscore[1,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$selectsample)],Maindata$GSCAscore[2,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$selectsample)],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i])
                              } else {
                                    points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX],Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i])
                              }
                              i <- i+1
                        }
                  }
                  leg.txt <- c(substr(Maindata$GSCAcontext,1,25),"Selected Samples","Not Selected Samples")
                  if (length(leg.txt)==2) {
                        legend("topleft",legend=leg.txt,pch=c(20,20),pt.bg=c("black","#00000022"),col=c("black","#00000022"),cex=0.8)                        
                  } else {
                        legend("topleft",legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=0.8)     
                  }
            }
      })
      
      waiting <- reactiveValues()
      
      output$GSCAdefaultplotthree <- renderPlot({
            if (Maindata$dim >= 3) {
                  par(oma=c(0.5,0,0.5,max(nchar(Maindata$patterndata[,1]))/1.5))
                  colcolorall <- rep("cyan",ncol(Maindata$GSCAscore))
                  colcolorall[Maindata$defaultsample] <- "blue"
                  tmprowv <- F
                  if (input$heatmapthreerowv)
                        tmprowv <- Maindata$GSCArowclust
                  heatmap.2(Maindata$GSCAscore,col=bluered,Colv=as.dendrogram(Maindata$GSCAclust),dendrogram="none",trace="none",Rowv=tmprowv,labCol=NA,ColSideColors=colcolorall,main="All Sample Heatmap")
                  legend("bottomleft",legend=c("Selected Sample","Unselected Sample"),lwd=1,col=c("blue","cyan"))
            }
      })
      
      output$GSCAdefaultplotthreeplus <- renderPlot({
            if (Maindata$dim >= 3 && length(Maindata$defaultsample) > 0) {
                  colcolorselect <- rep("white",ncol(Maindata$GSCAscore))
                  par(oma=c(0.5,0,0.5,max(nchar(Maindata$patterndata[,1]))/1.5))
                  tmprowv <- F
                  if (input$heatmapthreerowv)
                        tmprowv <- Maindata$GSCArowclust
                  if (!is.null(Maindata$GSCAcontext)) {
                        i <- 1
                        for(INDEX in Maindata$GSCAcontext) {
                              colcolorselect[Maindata$tab$SampleType %in% INDEX] <- COLORS[i]
                              i <- i+1
                        }
                        heatmap.2(Maindata$GSCAscore[,Maindata$defaultsample],col=bluered,labCol=NA,Rowv=tmprowv,dendrogram="none",trace="none",ColSideColors=colcolorselect[Maindata$defaultsample],main="Selected Sample Heatmap")
                        leg.txt <- substr(Maindata$GSCAcontext,1,25)
                        legend("bottomleft",legend=leg.txt,lwd=1,col=COLORS)
                  }
                  else {
                        heatmap.2(Maindata$GSCAscore[,Maindata$defaultsample],col=bluered,labCol=NA,Rowv=tmprowv,dendrogram="none",trace="none",main="Selected Sample Heatmap")
                  }
            }
      })
      
      #GSCA interactive plot
      output$GSCAinteractiveplot <- renderUI({
            if (!is.null(Maindata$patterndata)) 
                  if(nrow(Maindata$patterndata) != 0) {
                        if (Maindata$dim == 1) {
                              tagList(
                                    helpText(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                                    plotOutput("GSCAinteractiveplotone",height=300*(as.numeric(input$InputN)+1))
                              )
                        } else if (Maindata$dim == 2) {
                              div(align="center",
                                  plotOutput("GSCAinteractiveplottwo",clickId="coords",width=800,height=800),
                                  plotOutput("GSCAinteractiveplottwoplus",width=800,height=800),
                                  helpText(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found",""))
                              )      
                        } else {
                              tagList(          
                                    conditionalPanel(condition="input.Threecutoffvalue == 0",
                                                     helpText("Select sample range using sliders. The range would be the union set of multiple sliders."),
                                                     actionButton("GSCAthreesampleaddslider","Add Slider"),
                                                     actionButton("GSCAthreesampledeleteslider","Delete Slider"),
                                                     tags$div(class="row-fluid",
                                                     tags$div(class="span11",lapply(1:ifelse(is.null(GSCAthreeinfo$sampleslidernum),1,GSCAthreeinfo$sampleslidernum) , function(i) {
                                                           sliderInput(inputId = paste0("GSCAthreesampleslider",i),"",min=0,max=ncol(Maindata$GSCAscore),value=c(threesampleslidervalue[i,1],threesampleslidervalue[i,2]),step=1)
                                                     })
                                                     )
                                                     )
                                    ),
                                    tags$div(class="row-fluid",
                                             tags$div(class="span11",p(plotOutput("GSCAinteractiveplotthreecolbar",height=20)))
                                    ),
                                    tags$div(class="row-fluid",
                                          tags$div(class="span11",plotOutput("GSCAinteractiveplotthreeheatmap")),
                                          tags$div(class="span1",plotOutput("GSCARinteractiveplotthreeheatmaprowlab"))
                                    ),
                                    tags$div(class="row-fluid",
                                             tags$div(class="span11",p(plotOutput("GSCAinteractiveplotthreecolbarunder",height=20)))
                                    ),
                                    conditionalPanel(condition="input.Threecutoffvalue == 0",
                                                     checkboxInput("GSCAinteractiveplotthreeheatmapzoomincheck","Heatmap Zoom in & precise sample selection")),
                                    conditionalPanel("input.GSCAinteractiveplotthreeheatmapzoomincheck==1 && input.Threecutoffvalue == 0",
                                                     helpText("Select zoom in sample range"),
                                                     tags$div(class="row-fluid",
                                                            tags$div(class="span1",textInput("GSCAinteractiveplotthreeheatmapzoominrowone","From",value="0")),
                                                            tags$style(type='text/css', "#GSCAinteractiveplotthreeheatmapzoominrowone { width: 70px; }"),
                                                            tags$div(class="span1",textInput("GSCAinteractiveplotthreeheatmapzoominrowtwo","To",value="0")),
                                                            tags$style(type='text/css', "#GSCAinteractiveplotthreeheatmapzoominrowtwo { width: 70px; }")
                                                     ),
                                                     p(actionButton("GSCAinteractiveplotthreeheatmapzoominupdate","Update range to current slider")),
                                                     tags$div(class="row-fluid",
                                                              tags$div(class="span11",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot")),
                                                              tags$div(class="span1",plotOutput("GSCAinteractiveplotthreeheatmapzoominplotlab"))
                                                     )
                                                     ),
                                    helpText(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                                    plotOutput("GSCAinteractiveplotthreeplus"))
                        }
                  }      
      })

      observe({
            if (!is.null(input$GSCAinteractiveplotthreeheatmapzoominupdate) && input$GSCAinteractiveplotthreeheatmapzoominupdate>0) {
                  isolate({
                        threesampleslidervalue[GSCAthreeinfo$sampleslidernum,1] <<- as.numeric(input$GSCAinteractiveplotthreeheatmapzoominrowone)
                        threesampleslidervalue[GSCAthreeinfo$sampleslidernum,2] <<- as.numeric(input$GSCAinteractiveplotthreeheatmapzoominrowtwo)
                        eval(parse(text=paste0("updateSliderInput(session,\"GSCAthreesampleslider",GSCAthreeinfo$sampleslidernum,"\",value=as.numeric(c(input$GSCAinteractiveplotthreeheatmapzoominrowone,input$GSCAinteractiveplotthreeheatmapzoominrowtwo)))")))
                  })
            }
      })
      
      output$testui <- renderUI({
            if (!is.null(Maindata$dim) && Maindata$dim >= 3 && input$GSCAmethod=='GSCAinteractive')
                  p(actionButton("GSCAinteractivethreeupdate","Update Sample Selection"))
      })
      
      GSCAthreeinfo <- reactiveValues()
      
      #initiate number of sample slider, record current slider value
      observe({
            if (is.null(GSCAthreeinfo$sampleslidernum))
                  GSCAthreeinfo$sampleslidernum <- 1
            for (i in 1:GSCAthreeinfo$sampleslidernum) 
                  if(!is.null(eval(parse(text=paste0("input$GSCAthreesampleslider",i,"[1]"))))) {
                        threesampleslidervalue[i,1] <<- eval(parse(text=paste0("input$GSCAthreesampleslider",i,"[1]"))) 
                        threesampleslidervalue[i,2] <<- eval(parse(text=paste0("input$GSCAthreesampleslider",i,"[2]")))
                  }
      })
      
      #set up action for add slider button in three genedata case
      observe({
            if (!is.null(input$GSCAthreesampleaddslider) && input$GSCAthreesampleaddslider>0)
                  isolate({
                        if (GSCAthreeinfo$sampleslidernum < 5) {
                              GSCAthreeinfo$sampleslidernum <- GSCAthreeinfo$sampleslidernum + 1
                        }
                  })
      })
      
      #set up action for delete slider button in three genedata case
      observe({
            if (!is.null(input$GSCAthreesampledeleteslider) && input$GSCAthreesampledeleteslider>0)
                  isolate({
                        if (GSCAthreeinfo$sampleslidernum != 1) {
                              GSCAthreeinfo$sampleslidernum <- GSCAthreeinfo$sampleslidernum - 1
                        }
                  })
      })
      
      output$GSCAinteractiveplotone <- renderPlot({
            if (Maindata$dim == 1) {
                  par(mfrow=c(as.numeric(input$InputN)+1,1),oma=c(0,0,2,0))
                  hist(Maindata$GSCAscore,xlab="Sample Score",xlim=range(Maindata$GSCAscore),main="All Biological contexts")
                  abline(v=c(input$GSCAoneslider[1],input$GSCAoneslider[2]), lty=2)
                  if (!is.null(Maindata$GSCAcontext)) {
                        for(INDEX in Maindata$GSCAcontext) {
                              hist(Maindata$GSCAscore[Maindata$tab$SampleType %in% INDEX],xlim=range(Maindata$GSCAscore),main=substr(INDEX,1,25),xlab="Sample Score")
                              abline(v=c(input$GSCAoneslider[1],input$GSCAoneslider[2]), lty=2)
                        }     
                  }
            }
      })
      
      inpoly <- reactiveValues()
      
      get.coords <- reactive({
            data.frame(x=input$coords$x, y=input$coords$y)
      })
      
      output$GSCAinteractiveplottwo <- renderPlot({  
            if (Maindata$dim == 2) { 
                  inputcoords <- get.coords()
                  plot(Maindata$GSCAscore[1,], Maindata$GSCAscore[2,],pch=20, xlim = range(Maindata$GSCAscore[1,]), ylim = range(Maindata$GSCAscore[2,]),xlab=Maindata$patterndata[1,1],ylab=Maindata$patterndata[2,1])
                  if (input$reset != resetvalue) {
                        polycord <<- NULL
                        polynum <<- 1
                        resetvalue <<- input$reset
                        actiontaken <<- 1
                        inpoly$tf <- rep(F,ncol(Maindata$GSCAscore))
                  } 
                  if (input$undo != undovalue) {
                        if (!is.null(polycord)) {
                              if (polynum == max(polycord[,3]))
                                    polycord <<- polycord[-nrow(polycord),,drop=F]
                              if (nrow(polycord) == 0) {
                                    polynum <<- 1
                                    polycord <<- NULL
                              } else {
                                    polynum <<- max(polycord[,3])
                              }
                        }
                        undovalue <<- input$undo
                        actiontaken <<- 1
                  } 
                  if (input$finishpolygon != finishpolygonvalue)  {
                        isolate({
                              inpoly$tf <- rep(F,ncol(Maindata$GSCAscore))
                              for (j in 1:polynum) {
                                    tmpcord <- polycord[polycord[,3]==j,]
                                    if (is.matrix(tmpcord)) {
                                          inpoly$tf[which(as.logical(point.in.polygon(Maindata$GSCAscore[1,],Maindata$GSCAscore[2,],tmpcord[,1],tmpcord[,2])))] <- T
                                    }
                              }
                        })
                        points(Maindata$GSCAscore[1,inpoly$tf], Maindata$GSCAscore[2,inpoly$tf],col = 'green',pch=20)      
                        finishpolygonvalue <<- input$finishpolygon
                        finishpolygonaction <<- 1
                        actiontaken <<- 1
                  }            
                  if (input$addpolygon != addpolygonvalue) {
                        if (!is.null(polycord))
                              polynum <<- max(polycord[,3]) + 1
                        addpolygonvalue <<- input$addpolygon
                        actiontaken <<- 1
                  }
                  if (actiontaken == 0) 
                        if (!is.null(inputcoords$x))
                              polycord <<- rbind(polycord,c(inputcoords$x,inputcoords$y,polynum))      
                  if (is.matrix(polycord) && nrow(polycord) != 0) {
                        points(polycord[,1], polycord[,2], pch = 19, col = 'red', cex = 1.5)
                        for (j in 1:polynum) {
                              tmpcord <- polycord[polycord[,3]==j,]
                              if (is.matrix(tmpcord) && nrow(tmpcord) > 1) {
                                    for (i in 1:(nrow(tmpcord)-1)) 
                                          lines(c(tmpcord[i,1],tmpcord[i+1,1]),c(tmpcord[i,2],tmpcord[i+1,2]),type="l",col="blue")
                                    if (j==polynum && finishpolygonaction == 0) {
                                          lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lty=2)
                                    } else {
                                          lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue")
                                    }
                              }
                        }
                  }
                  finishpolygonaction <<- 0
                  actiontaken <<- 0
            }
      })
      
      
      output$GSCAinteractiveplottwoplus <- renderPlot({  
            if (Maindata$dim == 2) { 
                  if (sum(inpoly$tf) != 0) {
                        plot(Maindata$GSCAscore[1,], Maindata$GSCAscore[2,], xlim = range(Maindata$GSCAscore[1,]), ylim = range(Maindata$GSCAscore[2,]),xlab=Maindata$patterndata[1,1],ylab=Maindata$patterndata[2,1],pch=20,cex=0.7,col="#00000022")            
                        for (j in 1:polynum) {
                              tmpcord <- polycord[polycord[,3]==j,]
                              if (is.matrix(tmpcord)) {
                                    for (i in 1:(nrow(tmpcord)-1)) 
                                          lines(c(tmpcord[i,1],tmpcord[i+1,1]),c(tmpcord[i,2],tmpcord[i+1,2]),type="l",col="blue",lty=2,lwd=2)
                                    lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lty=2,lwd=2)
                              }
                        }
                        toprankingsample <- NULL
                        if (!is.null(Maindata$GSCAcontext)) {
                              for(INDEX in Maindata$GSCAcontext) {
                                    toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                              }
                        }
                        points(Maindata$GSCAscore[1,setdiff(which(inpoly$tf),toprankingsample)], Maindata$GSCAscore[2,setdiff(which(inpoly$tf),toprankingsample)],cex=0.7,pch=20)
                        if (!is.null(Maindata$GSCAcontext)) {
                              i <- 1
                              for (INDEX in Maindata$GSCAcontext) {    
                                    if (input$Inputenrichedareaonly == TRUE) {
                                          points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0], Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0],col = COLORS[i], pch = STYLES[i], bg = COLORS[i])
                                    } else {
                                          points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX], Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX],col = COLORS[i], pch = STYLES[i],bg = COLORS[i])
                                    }
                                    i <- i+1
                              }
                        }
                        leg.txt <- c(substr(Maindata$GSCAcontext,1,25),"Selected Samples","Not Selected Samples")
                        if (length(leg.txt)==2) {
                              legend("topleft",legend=leg.txt,pch=c(20,20),pt.bg=c("black","#00000022"),col=c("black","#00000022"),cex=0.8)                        
                        } else {
                              legend("topleft",legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=0.8)                        
                        }
                  }
            }
      })
      
      observe({
            if (!is.null(input$Threecutoffvalue)) {
                  if (input$Threecutoffvalue == F) {
                        selectsample <- NULL
                        for (i in 1:GSCAthreeinfo$sampleslidernum) {
                              if(!is.null(eval(parse(text=paste0("input$GSCAthreesampleslider",i,"[1]"))))) 
                                    eval(parse(text=paste0("selectsample <- union(selectsample,Maindata$GSCAclust$order[input$GSCAthreesampleslider",i,"[1]:input$GSCAthreesampleslider",i,"[2]])")))
                        }
                        GSCAthreeinfo$selectsample <- selectsample
                  } else {
                        selectsample <- 1:ncol(Maindata$GSCAscore)
                        for (i in 1:Maindata$dim) {
                              eval(parse(text=paste0("selectsample <- intersect(selectsample,which(input$GSCAthreevalueslider",i,"[1]<Maindata$GSCAscore[",i,",] & input$GSCAthreevalueslider",i,"[2]>Maindata$GSCAscore[",i,",]))")))
                        }
                        GSCAthreeinfo$selectsample <- selectsample
                  }
            }
      })
            
      output$GSCARinteractiveplotthreeheatmaprowlab <- renderPlot({
            par(mar = c(0,0,0,0))
            plot(c(1, 1), c(Maindata$dim, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            pos <- seq(1.5-1/2/Maindata$dim,Maindata$dim-0.5+1/2/Maindata$dim,length.out=Maindata$dim)
            if (input$heatmapthreerowv) {
                  Rowv <- rowMeans(Maindata$GSCAscore)
                  hcr <- hclust(dist(Maindata$GSCAscore))
                  ddr <- as.dendrogram(hcr)
                  ddr <- reorder(ddr, Rowv)
                  rowInd <- order.dendrogram(ddr)
                  names <- row.names(Maindata$GSCAscore)[rowInd]
            } else {
                  names <- rev(row.names(Maindata$GSCAscore))
            }
            for (i in 1:Maindata$dim) {
                  text(1,pos[i],names[i],cex=1.5)
            }
      })
      
      output$GSCAinteractiveplotthreecolbar <- renderPlot({
            if (Maindata$dim >= 3) { 
                  colcolorall <- rep("cyan",ncol(Maindata$GSCAscore))
                  colcolorall[GSCAthreeinfo$selectsample] <- "blue"
                  par(mar=c(0,0,0,0))
                  image(cbind(1:ncol(Maindata$GSCAscore)),col=colcolorall[Maindata$GSCAclust$order],axes=F)
            }
      })
      
      output$GSCAinteractiveplotthreecolbarunder <- renderPlot({
            if (Maindata$dim >= 3) { 
                  colcolorall <- rep("cyan",ncol(Maindata$GSCAscore))
                  colcolorall[GSCAthreeinfo$selectsample] <- "blue"
                  par(mar=c(0,0,0,0))
                  image(cbind(1:ncol(Maindata$GSCAscore)),col=colcolorall[Maindata$GSCAclust$order],axes=F)
            }
      })
      
      output$GSCAinteractiveplotthreeheatmap <- renderPlot({
            if (Maindata$dim >= 3) {
                  par(mar=c(0,0,0,0))
                  if (input$heatmapthreerowv) {
                        rowInd <- order.dendrogram(Maindata$GSCArowclust)
                        image(1:ncol(Maindata$GSCAscore),1:nrow(Maindata$GSCAscore),t(Maindata$GSCAscore[rowInd,Maindata$GSCAclust$order]),col=bluered(100),axes=F)
                  } else {
                        image(1:ncol(Maindata$GSCAscore),1:nrow(Maindata$GSCAscore),t(Maindata$GSCAscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order]),col=bluered(100),axes=F)
                  }
            }
      })

      output$GSCAinteractiveplotthreeheatmapzoominplot <- renderPlot({
            par(mar=c(0,0,0,0))
            rowone <- as.numeric(input$GSCAinteractiveplotthreeheatmapzoominrowone)
            rowtwo <- as.numeric(input$GSCAinteractiveplotthreeheatmapzoominrowtwo)
            if (rowone+rowtwo>0 && rowone < rowtwo && max(rowone,rowtwo) <= ncol(Maindata$GSCAscore)) {
                  if (input$heatmapthreerowv) {
                        rowInd <- order.dendrogram(Maindata$GSCArowclust)
                        image(1:(rowtwo-rowone+1),1:nrow(Maindata$GSCAscore),t(Maindata$GSCAscore[rowInd,Maindata$GSCAclust$order[rowone:rowtwo]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$GSCAscore),max(Maindata$GSCAscore),length.out=101))
                  } else {
                        image(1:(rowtwo-rowone+1),1:nrow(Maindata$GSCAscore),t(Maindata$GSCAscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[rowone:rowtwo]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$GSCAscore),max(Maindata$GSCAscore),length.out=101))
                  }
            }
      })
      
      output$GSCAinteractiveplotthreeheatmapzoominplotlab <- renderPlot({
            par(mar = c(0,0,0,0))
            plot(c(1, 1), c(Maindata$dim, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            pos <- seq(1.5-1/2/Maindata$dim,Maindata$dim-0.5+1/2/Maindata$dim,length.out=Maindata$dim)
            if (input$heatmapthreerowv) {
                  Rowv <- rowMeans(Maindata$GSCAscore)
                  hcr <- hclust(dist(Maindata$GSCAscore))
                  ddr <- as.dendrogram(hcr)
                  ddr <- reorder(ddr, Rowv)
                  rowInd <- order.dendrogram(ddr)
                  names <- row.names(Maindata$GSCAscore)[rowInd]
            } else {
                  names <- rev(row.names(Maindata$GSCAscore))
            }
            for (i in 1:Maindata$dim) {
                  text(1,pos[i],names[i],cex=1.5)
            }
      })
      
      output$GSCAinteractiveplotthreeplus <- renderPlot({
            if (Maindata$dim >= 3 && length(Maindata$selectsample) > 0) {
                  if(input$GSCAinteractivethreeupdate>0) {
                        Maindata$Ranking
                        Maindata$GSCAcontext
                        input$heatmapthreerowv
                        isolate({
                              par(oma=c(0.5,0,0.5,max(nchar(Maindata$patterndata[,1]))/1.5))
                              colcolorselect <- rep("white",ncol(Maindata$GSCAscore))
                              tmprowv <- F
                              if (input$heatmapthreerowv)
                                    tmprowv <- Maindata$GSCArowclust
                              if (!is.null(Maindata$GSCAcontext)) {
                                    i <- 1
                                    for(INDEX in Maindata$GSCAcontext) {
                                          colcolorselect[Maindata$tab$SampleType %in% INDEX] <- COLORS[i]
                                          i <- i+1
                                    }
                                    leg.txt <- substr(Maindata$GSCAcontext,1,25)
                                    heatmap.2(Maindata$GSCAscore[,Maindata$selectsample],col=bluered,labCol=NA,Rowv=tmprowv,dendrogram="none",trace="none",ColSideColors=colcolorselect[Maindata$selectsample],main="Selected Sample Heatmap")
                                    legend("bottomleft",legend=leg.txt,lwd=1,col=COLORS)
                              } else {
                                    heatmap.2(Maindata$GSCAscore[,Maindata$selectsample],col=bluered,labCol=NA,Rowv=tmprowv,dendrogram="none",trace="none",main="Selected Sample Heatmap")
                              }
                        })
                  }
            }
      })
      
      output$GSCArankingtable <- renderDataTable(Maindata$Ranking)
      
      #####   Mainmethod : Download   #####
      
      output$Downloadsidebarui <- renderUI({
            if (!is.null(Maindata$genedata)) {
                  if (Maindata$dim == 1) {
                        tagList(
                              wellPanel(
                                    selectInput("Downloadplotonetype","Select File Type for Plot",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotonefilename","Enter File Name","GSCA Plot"),
                                    textInput("Downloadplotonefilewidth","Enter plot width (inches)",10),
                                    textInput("Downloadplotonefileheight","Enter plot height (inches)",3*(as.numeric(input$InputN)+1)),
                                    p(downloadButton("Downloadplotone","Save Plot"))
                              ),
                              wellPanel(
                                    helpText("Change Plotting Details"),
                                    textInput("Downloadmaintitleone","Enter Main Title","GSCA result"),
                                    textInput("Downloadxlabone","Enter Title for X Axis","SampleScore"),
                                    textInput("Downloadylabone","Enter Title for Y Axis","Frequency"),
                                    textInput("Downloadxlimminone","Enter minimum value of X Axis",min(Maindata$GSCAscore)),
                                    textInput("Downloadxlimmaxone","Enter maximum value of X Axis",max(Maindata$GSCAscore)),
                                    selectInput("Downloadcolone","Select filling color",c("NULL","gray","red","blue","green","yellow","purple"))
                              )
                        )
                  } else if (Maindata$dim == 2) {
                        tagList(
                              wellPanel(
                                    selectInput("Downloadplottwotype","Select File Type for Plot",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplottwofilename","Enter File Name","GSCA Plot"),
                                    textInput("Downloadplottwofilewidth","Enter plot width (inches)",15),
                                    textInput("Downloadplottwofileheight","Enter plot height (inches)",15),
                                    p(downloadButton("Downloadplottwo","Save Plot"))
                              ),
                              wellPanel(
                                    helpText("Change Plotting Details"),
                                    textInput("Downloadmaintitletwo","Enter Main Title","GSCA result"),
                                    textInput("Downloadxlabtwo","Enter Title for X Axis",Maindata$patterndata[1,1]),
                                    textInput("Downloadylabtwo","Enter Title for Y Axis",Maindata$patterndata[2,1]),
                                    textInput("Downloadxlimmintwo","Enter minimum value of X Axis",min(Maindata$GSCAscore[1,])),
                                    textInput("Downloadxlimmaxtwo","Enter maximum value of X Axis",max(Maindata$GSCAscore[1,])),
                                    textInput("Downloadylimmintwo","Enter minimum value of Y Axis",min(Maindata$GSCAscore[2,])),
                                    textInput("Downloadylimmaxtwo","Enter maximum value of Y Axis",max(Maindata$GSCAscore[2,])),
                                    textInput("Downloadlegcextwo","Enter legend size",0.8),
                                    checkboxInput("Downloadcorvaluetwo","Show correlation and p-value")
                              )
                        )
                  } else {
                        tagList(
                              wellPanel(
                                    selectInput("Downloadplotthreeplotonetype","Select File Type for Heatmap One",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotthreeplotonefilename","Enter File Name","GSCA Heatmap One"),
                                    textInput("Downloadplotthreeplotonefilewidth","Enter plot width (inches)",20),
                                    textInput("Downloadplotthreeplotonefileheight","Enter plot height (inches)",10),
                                    p(downloadButton("Downloadplotthreeplotone","Save Heatmap One")),
                                    selectInput("Downloadplotthreeplottwotype","Select File Type for Heatmap Two",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotthreeplottwofilename","Enter File Name","GSCA Heatmap Two"),
                                    textInput("Downloadplotthreeplottwofilewidth","Enter plot width (inches)",20),
                                    textInput("Downloadplotthreeplottwofileheight","Enter plot height (inches)",10),
                                    p(downloadButton("Downloadplotthreeplottwo","Save Heatmap Two"))
                              ),
                              wellPanel(
                                    helpText("Change Plotting Details"),
                                    textInput("Downloadmaintitlethreeplotone","Enter Main Title For Heatmap One","All Sample Heatmap"),
                                    textInput("Downloadmaintitlethreeplottwo","Enter Main Title For Heatmap Two","Selected Sample Heatmap"),
                                    selectInput("Downloadthreecolplotone","Select Palette for Heatmap One",choices=c("Bluered","Heat","Rainbow","Terrain","Topo","CM","gray")),
                                    selectInput("Downloadthreecolplottwo","Select Palette for Heatmap Two",choices=c("Bluered","Heat","Rainbow","Terrain","Topo","CM","gray")),
                                    checkboxInput("Downloadthreecoldendplotone","Display Column Histogram in Heatmap One (Could take some time)"),
                                    checkboxInput("Downloadthreecoldendplottwo","Display Column Histogram in Heatmap Two"),
                                    checkboxInput("Downloadthreerowvplotone","Cluster on rows in Heatmap One"),
                                    uiOutput("Downloadthreerowdendplotoneui"),
                                    checkboxInput("Downloadthreerowvplottwo","Cluster on rows in Heatmap Two"),
                                    uiOutput("Downloadthreerowdendplottwoui"),
                                    checkboxInput("Downloadthreerotatelabplotone","Rotate row label in Heatmap One"),
                                    checkboxInput("Downloadthreerotatelabplottwo","Rotate row label in Heatmap Two"),
                                    textInput("Downloadlegcexthreeone","Enter legend size in Heatmap One",1),
                                    textInput("Downloadlegcexthreetwo","Enter legend size in Heatmap Two",1)
                              )
                        )
                  }
            }
      })
      
      output$Downloadthreerowdendplotoneui <- renderUI({
            if (!is.null(input$Downloadthreerowvplotone) && input$Downloadthreerowvplotone==TRUE)
                  checkboxInput("Downloadthreerowdendplotone","Display Row Histogram in Heatmap One")
      })
      
      output$Downloadthreerowdendplottwoui <- renderUI({
            if (!is.null(input$Downloadthreerowvplottwo) && input$Downloadthreerowvplottwo==TRUE)
                  checkboxInput("Downloadthreerowdendplottwo","Display Row Histogram in Heatmap Two")
      })


      observe({
            if (input$Mainmethod=='Download') {
                  if (input$Downloadregionselect == 'GSCAdefault') {
                        Maindata$downloadsample <- Maindata$defaultsample
                        Maindata$downloadcontext <- Maindata$defaultcontext
                        Maindata$downloadranking <- Maindata$defaultranking
                  } else {
                        Maindata$downloadsample <- Maindata$interactsample
                        Maindata$downloadcontext <- Maindata$interactcontext
                        Maindata$downloadranking <- Maindata$interactranking
                  }
            }
      })
      
      #Download Ranking Table
      
      output$Downloadshowrankingtable <- renderDataTable(Maindata$downloadranking)
      
      output$Downloadranktable <- downloadHandler(
            filename = function() { paste0(input$Downloadranktablefilename,'.',input$Downloadranktabletype) },
            content = function(file) {
                  if (input$Downloadranktabletype == 'txt') {
                        write.table(Maindata$downloadranking,file,row.names=F)
                  } else {
                        write.csv(Maindata$downloadranking,file,row.names=F)
                  }
            }
      )
      
      #Download Plot
      
      output$Downloadshowplotui <- renderUI({
            if (!is.null(Maindata$dim))
                  if (Maindata$dim == 1) {
                        plotOutput("downloadonerenderplot",height=300*(as.numeric(input$InputN)+1))
                  } else if (Maindata$dim == 2) {
                        plotOutput("downloadtworenderplot",height=800)
                  } else {
                        tagList(
                              h4("Heatmap One"),
                              plotOutput("downloadthreeonerenderplot"),
                              h4("Heatmap Two"),
                              plotOutput("downloadthreetworenderplot")
                        )
                  }
      })
      
      output$downloadonerenderplot <- renderPlot(downloadonefunc())
      
      output$downloadtworenderplot <- renderPlot(downloadtwofunc())
      
      output$downloadthreeonerenderplot <- renderPlot(downloadthreeonefunc())
      
      output$downloadthreetworenderplot <- renderPlot(downloadthreetwofunc())
      
      downloadonefunc <- function() {
            colone <- input$Downloadcolone
            if (colone == "NULL")
                  colone <- NULL
            par(mfrow=c(as.numeric(input$InputN)+1,1),oma=c(0,0,2,0))
            hist(Maindata$GSCAscore,xlab=input$Downloadxlabone,ylab=input$Downloadylabone,xlim=as.numeric(c(input$Downloadxlimminone,input$Downloadxlimmaxone)),col=colone,main="All Biological contexts")
            if (input$Downloadregionselect == 'GSCAdefault') {
                  abline(v=Maindata$cutoffval[1], lty=2)      
            } else {
                  abline(v=c(input$GSCAoneslider[1],input$GSCAoneslider[2]), lty=2)
            }
            if (!is.null(Maindata$downloadcontext)) {
                  for(INDEX in Maindata$downloadcontext) {
                        hist(Maindata$GSCAscore[Maindata$tab$SampleType %in% INDEX],xlab=input$Downloadxlabone,ylab=input$Downloadylabone,xlim=as.numeric(c(input$Downloadxlimminone,input$Downloadxlimmaxone)),col=colone,main=substr(INDEX,1,25))
                        if (input$Downloadregionselect == 'GSCAdefault') {
                              abline(v=Maindata$cutoffval[1], lty=2)      
                        } else {
                              abline(v=c(input$GSCAoneslider[1],input$GSCAoneslider[2]), lty=2)
                        }
                  }     
            }
            title(input$Downloadmaintitleone,outer=T)
      }
      
      downloadtwofunc <- function() {
            cortext <- paste0("Correlation: ",round(Maindata$twocorr,3),"; ","Correlation p-value: ",round(Maindata$twocorrp,3),"; ","Slope: ",round(Maindata$twoslope,3),"; ","Slope p-value: ",round(Maindata$twoslopep,3))
            if (input$Downloadregionselect == 'GSCAdefault') {
                  plot(Maindata$GSCAscore[1,],Maindata$GSCAscore[2,],col="#00000022",pch=20,cex=0.7,xlab=input$Downloadxlabtwo,ylab=input$Downloadylabtwo,xlim = as.numeric(c(input$Downloadxlimmintwo,input$Downloadxlimmaxtwo)), ylim = as.numeric(c(input$Downloadylimmintwo,input$Downloadylimmaxtwo)),main=input$Downloadmaintitletwo)
                  toprankingsample <- NULL
                  if (!is.null(Maindata$downloadcontext)) {
                        for(INDEX in Maindata$downloadcontext) {
                              toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                        }
                  }
                  points(Maindata$GSCAscore[1,setdiff(Maindata$downloadsample,toprankingsample)],Maindata$GSCAscore[2,setdiff(Maindata$downloadsample,toprankingsample)],cex=0.7,pch=20)
                  abline(v=Maindata$cutoffval[1], h=Maindata$cutoffval[2], lty=2,lwd=2)
                  if (!is.null(Maindata$downloadcontext)) {
                        i <- 1
                        for(INDEX in Maindata$downloadcontext) {
                              if (input$Inputenrichedareaonly == TRUE) {
                                    points(Maindata$GSCAscore[1,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$downloadsample)],Maindata$GSCAscore[2,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$downloadsample)],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i])
                              } else {
                                    points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX],Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i])
                              }
                              i <- i+1
                        }     
                  }
                  leg.txt <- c(substr(Maindata$downloadcontext,1,25),"Selected Samples","Not Selected Samples")
                  if (length(leg.txt)==2) {
                        legend("topleft",legend=leg.txt,pch=c(20,20),pt.bg=c("black","#00000022"),col=c("black","#00000022"),cex=as.numeric(input$Downloadlegcextwo))                        
                  } else {
                        legend("topleft",legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=as.numeric(input$Downloadlegcextwo))   
                  }
                  if (input$Downloadcorvaluetwo)
                        mtext(cortext)
            } else {
                  if (sum(inpoly$tf) != 0) {
                        plot(Maindata$GSCAscore[1,], Maindata$GSCAscore[2,], xlim = as.numeric(c(input$Downloadxlimmintwo,input$Downloadxlimmaxtwo)), ylim = as.numeric(c(input$Downloadylimmintwo,input$Downloadylimmaxtwo)),xlab=input$Downloadxlabtwo,ylab=input$Downloadylabtwo,col="#00000022",pch=20,cex=0.7,main=input$Downloadmaintitletwo)            
                        for (j in 1:polynum) {
                              tmpcord <- polycord[polycord[,3]==j,]
                              if (is.matrix(tmpcord)) {
                                    for (i in 1:(nrow(tmpcord)-1)) 
                                          lines(c(tmpcord[i,1],tmpcord[i+1,1]),c(tmpcord[i,2],tmpcord[i+1,2]),type="l",col="blue",lty=2,lwd=2)
                                    lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lty=2,lwd=2)
                              }
                        }
                        toprankingsample <- NULL
                        if (!is.null(Maindata$downloadcontext)) {
                              for(INDEX in Maindata$downloadcontext) {
                                    toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                              }
                        }
                        points(Maindata$GSCAscore[1,setdiff(which(inpoly$tf),toprankingsample)], Maindata$GSCAscore[2,setdiff(which(inpoly$tf),toprankingsample)],cex=0.7,pch=20)
                        if (!is.null(Maindata$downloadcontext)) {
                              i <- 1
                              for (INDEX in Maindata$downloadcontext) {    
                                    if (input$Inputenrichedareaonly == TRUE) {
                                          points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0], Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0],col = COLORS[i], pch = STYLES[i], bg = COLORS[i])
                                    } else {
                                          points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX], Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX],col = COLORS[i], pch = STYLES[i], bg = COLORS[i])
                                    }
                                    i <- i+1
                              }     
                        }
                        leg.txt <- c(substr(Maindata$downloadcontext,1,25),"Selected Samples","Not Selected Samples")
                        if (length(leg.txt)==2) {
                              legend("topleft",legend=leg.txt,pch=c(20,20),pt.bg=c("black","#00000022"),col=c("black","#00000022"),cex=as.numeric(input$Downloadlegcextwo))                        
                        } else {
                              legend("topleft",legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=as.numeric(input$Downloadlegcextwo))     
                        }
                        if (input$Downloadcorvaluetwo)
                              mtext(cortext)
                  }
            }
      }    
      
      threepalette <- function(inputcol) {
            switch(inputcol,
                   Bluered=bluered,
                   Heat=heat.colors,
                   Rainbow=rainbow,
                   Terrain=terrain.colors,
                   Topo=topo.colors,
                   CM=cm.colors,
                   gray=gray.colors
            )
      }
      
      downloadthreeonefunc <- function() { 
            if(input$Downloadthreecoldendplotone && input$Downloadthreerowvplotone && !is.null(input$Downloadthreerowdendplotone) && input$Downloadthreerowdendplotone) {
                  threeonedendro <- "both"
            } else if (input$Downloadthreecoldendplotone && !(input$Downloadthreerowvplotone && !is.null(input$Downloadthreerowdendplotone) && input$Downloadthreerowdendplotone)) {
                  threeonedendro <- "column"
            } else if (!input$Downloadthreecoldendplotone && input$Downloadthreerowvplotone && !is.null(input$Downloadthreerowdendplotone) && input$Downloadthreerowdendplotone) {
                  threeonedendro <- "row"
            } else if (!input$Downloadthreecoldendplotone && !(input$Downloadthreerowvplotone && !is.null(input$Downloadthreerowdendplotone) && input$Downloadthreerowdendplotone)) {
                  threeonedendro <- "none"
            }
            par(oma=c(0.5,0,0.5,max(nchar(Maindata$patterndata[,1]))/ifelse(input$Downloadthreerotatelabplotone,2,1.5)))
            colcolorall <- rep("cyan",ncol(Maindata$GSCAscore))
            colcolorall[Maindata$downloadsample] <- "blue"
            tmprowv <- F
            if (input$Downloadthreerowvplotone)
                  tmprowv <- Maindata$GSCArowclust
            heatmap.2(Maindata$GSCAscore,col=threepalette(input$Downloadthreecolplotone),Colv=as.dendrogram(Maindata$GSCAclust),srtRow=ifelse(input$Downloadthreerotatelabplotone,-45,0),dendrogram=threeonedendro,trace="none",Rowv=tmprowv,labCol=NA,ColSideColors=colcolorall,main=input$Downloadmaintitlethreeplotone)
            legend("bottomleft",legend=c("Selected Sample","Unselected Sample"),lwd=1,col=c("blue","cyan"),cex=as.numeric(input$Downloadlegcexthreeone))           
      }
      
      downloadthreetwofunc <- function() {
            colcolorselect <- rep("white",ncol(Maindata$GSCAscore))
            if(input$Downloadthreecoldendplottwo && input$Downloadthreerowvplottwo && !is.null(input$Downloadthreerowdendplottwo) && input$Downloadthreerowdendplottwo) {
                  threetwodendro <- "both"
            } else if (input$Downloadthreecoldendplottwo && !(input$Downloadthreerowvplottwo && !is.null(input$Downloadthreerowdendplottwo) && input$Downloadthreerowdendplottwo)) {
                  threetwodendro <- "column"
            } else if (!input$Downloadthreecoldendplottwo && input$Downloadthreerowvplottwo && !is.null(input$Downloadthreerowdendplottwo) && input$Downloadthreerowdendplottwo) {
                  threetwodendro <- "row"
            } else if (!input$Downloadthreecoldendplottwo && !(input$Downloadthreerowvplottwo && !is.null(input$Downloadthreerowdendplottwo) && input$Downloadthreerowdendplottwo)) {
                  threetwodendro <- "none"
            }
            if (!is.null(Maindata$downloadcontext)) {
                  par(oma=c(0.5,0,0.5,max(nchar(Maindata$patterndata[,1]))/ifelse(input$Downloadthreerotatelabplotone,2,1.5)))
                  i <- 1
                  for(INDEX in Maindata$downloadcontext) {
                        colcolorselect[Maindata$tab$SampleType %in% INDEX] <- COLORS[i]
                        i <- i+1
                  }
                  tmprowv <- F
                  if (input$Downloadthreerowvplottwo)
                        tmprowv <- Maindata$GSCArowclust
                  heatmap.2(Maindata$GSCAscore[,Maindata$downloadsample],col=threepalette(input$Downloadthreecolplottwo),labCol=NA,Rowv=tmprowv,srtRow=ifelse(input$Downloadthreerotatelabplottwo,-45,0),dendrogram=threetwodendro,trace="none",ColSideColors=colcolorselect[Maindata$downloadsample],main=input$Downloadmaintitlethreeplottwo)
                  leg.txt <- substr(Maindata$downloadcontext,1,25)
                  legend("bottomleft",legend=leg.txt,lwd=1,col=COLORS,cex=as.numeric(input$Downloadlegcexthreetwo))
            }
      }
      
      output$Downloadplotone <- downloadHandler(
            filename = function() { paste0(input$Downloadplotonefilename,'.',input$Downloadplotonetype) },
            content = function(file) {
                  if (input$Downloadplotonetype == 'pdf') {
                        pdf(file,width=as.numeric(input$Downloadplotonefilewidth),height=as.numeric(input$Downloadplotonefileheight))
                  } else if (input$Downloadplotonetype == 'ps') {
                        postscript(file,width=as.numeric(input$Downloadplotonefilewidth),height=as.numeric(input$Downloadplotonefileheight),paper="special")
                  } else {
                        eval(parse(text=paste0(input$Downloadplotonetype,"(file,res=72,units=\"in\",width=as.numeric(input$Downloadplotonefilewidth),height=as.numeric(input$Downloadplotonefileheight))")))
                  }
                  downloadonefunc()
                  dev.off()
            }
      )
      
      output$Downloadplottwo <- downloadHandler(
            filename = function() { paste0(input$Downloadplottwofilename,'.',input$Downloadplottwotype) },
            content = function(file) {
                  if (input$Downloadplottwotype == 'pdf') {
                        pdf(file,width=as.numeric(input$Downloadplottwofilewidth),height=as.numeric(input$Downloadplottwofileheight))
                  } else if (input$Downloadplottwotype == 'ps') {
                        postscript(file,width=as.numeric(input$Downloadplottwofilewidth),height=as.numeric(input$Downloadplottwofileheight),paper="special")
                  } else {
                        eval(parse(text=paste0(input$Downloadplottwotype,"(file,res=72,units=\"in\",width=as.numeric(input$Downloadplottwofilewidth),height=as.numeric(input$Downloadplottwofileheight))")))
                  }
                  downloadtwofunc()
                  dev.off()
            }
      )
      
      output$Downloadplotthreeplotone <- downloadHandler(
            filename = function() { paste0(input$Downloadplotthreeplotonefilename,'.',input$Downloadplotthreeplotonetype) },
            content = function(file) {
                  if (input$Downloadplotthreeplotonetype == 'pdf') {
                        pdf(file,width=as.numeric(input$Downloadplotthreeplotonefilewidth),height=as.numeric(input$Downloadplotthreeplotonefileheight))
                  } else if (input$Downloadplotthreeplotonetype == 'ps') {
                        postscript(file,width=as.numeric(input$Downloadplotthreeplotonefilewidth),height=as.numeric(input$Downloadplotthreeplotonefileheight),paper="special")
                  } else {
                        eval(parse(text=paste0(input$Downloadplotthreeplotonetype,"(file,res=72,units=\"in\",width=as.numeric(input$Downloadplotthreeplotonefilewidth),height=as.numeric(input$Downloadplotthreeplotonefileheight))")))                        
                  }
                  downloadthreeonefunc()
                  dev.off()
            }
      )
      
      output$Downloadplotthreeplottwo <- downloadHandler(
            filename = function() { paste0(input$Downloadplotthreeplottwofilename,'.',input$Downloadplotthreeplottwotype) },
            content = function(file) {
                  if (input$Downloadplotthreeplottwotype == 'pdf') {
                        pdf(file,width=as.numeric(input$Downloadplotthreeplottwofilewidth),height=as.numeric(input$Downloadplotthreeplottwofileheight))
                  } else if (input$Downloadplotthreeplottwotype == 'ps') {
                        postscript(file,width=as.numeric(input$Downloadplotthreeplottwofilewidth),height=as.numeric(input$Downloadplotthreeplottwofileheight),paper="special")
                  } else {
                        eval(parse(text=paste0(input$Downloadplotthreeplottwotype,"(file,res=72,units=\"in\",width=as.numeric(input$Downloadplotthreeplottwofilewidth),height=as.numeric(input$Downloadplotthreeplottwofileheight))")))                        
                  }
                  downloadthreetwofunc()
                  dev.off()
            }
      )
      
      
})




