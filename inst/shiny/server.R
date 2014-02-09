######################################################
##     GSCAR:Geneset Context Analysis R Platform    ##
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
      
      output$InputGenesetnameui <- renderUI({
            if (input$InputGenesetmethod == 'InputspecifyGeneset' || (input$InputGenesetmethod == 'InputuploadGeneset' && !is.null(Currentdata$rawGenesetFile) && !is.null(ncol(Currentdata$rawGenesetFile)) &&  (ncol(Currentdata$rawGenesetFile) == 2 || (!is.null(input$InputGenesetcolnumthree) && input$InputGenesetcolnumthree == F))))
                  textInput("InputGenesetname","Input Geneset Name","Geneset 1")
      })
      
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
                  "No input genedata in GSCAR"
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
            if (require(Affymoe430Expr)) {
                  complist <- c(complist,"Affymetrix Mouse 430 2.0 Array (GPL1261)"="moe430")
            }        
            radioButtons("Summarycompselect","Select Compendium", complist)
      })
      
      
      #Update Maindata information
      observe({   
            if (!is.null(input$Summarycompselect)) {
            if(input$Summarycompselect=="moe430"){
                  data(Affymoe430Exprtab)
                  Maindata$tab <- Affymoe430Exprtab
            } else if (input$Summarycompselect=="hgu133a"){
                  data(Affyhgu133aExprtab)
                  Maindata$tab <- Affyhgu133aExprtab
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
      
      ######  Mainmethod : GSCAR  ######
      
      GSCAR <- reactiveValues()
      
      #When reentering GSCAR page, GSCARmethod should be set to GSCARdefault
      observe({
            if (input$Mainmethod == "Input" | input$Mainmethod == "Select")
                  updateRadioButtons(session,"GSCARmethod",choices=list("Default Enrichment Region Selection"="GSCARdefault","Interactive Enrichment Region Selection"="GSCARinteractive"),selected="Default Enrichment Region Selection")
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
      
      #Calculating GSCAR Score, default selected samples, default Ranking, also hclust dendrogram if there are geq three genesets
      observe({
            if (input$Mainmethod == 'GSCAR') {
                  isolate({
                        if (!is.null(Maindata$patterndata))
                              if(nrow(Maindata$patterndata) != 0) {
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
                                                GSup <- GSup + Expr
                                          }
                                          GSup <- GSup / n.GSup
                                          n.GSdown <- length(neggene)
                                          GSdown <- rep(0, nrow(Maindata$tab))  
                                          for (i in neggene) {
                                                path <- system.file("extdata",package=paste0("Affy",input$Summarycompselect,"Expr"))
                                                load(paste0(path,"/",i,".rda"))
                                                GSdown <- GSdown + Expr
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
                                    
                                    Maindata$GSCARscore <- scoremat
                                    Maindata$defaultsample <- selectsample
                                    Maindata$cutoffval <- cutoffval
                                    if (Maindata$dim >= 3) {
                                          Maindata$GSCARclust <- hclust(dist(t(scoremat)))
                                    }
                                    if (threesampleslidervalue[1,3] != ncol(Maindata$GSCARscore)) {
                                          threesampleslidervalue[,3] <<- rep(ncol(Maindata$GSCARscore),5)
                                          threesampleslidervalue[,1] <<- rep(0,5)
                                          threesampleslidervalue[,2] <<- rep(ncol(Maindata$GSCARscore),5)
                                    }
                                    if (Maindata$dim == 2) {
                                          Maindata$twocorr <- cor(Maindata$GSCARscore[1,],Maindata$GSCARscore[2,])
                                          Maindata$twocorrp <- cor.test(Maindata$GSCARscore[1,],Maindata$GSCARscore[2,])$p.value
                                          twolm <- summary(lm(Maindata$GSCARscore[2,]~Maindata$GSCARscore[1,]))
                                          Maindata$twoslope <- twolm$coefficients[2,1]
                                          Maindata$twoslopep <- twolm$coefficients[2,4]
                                    }
                              }
                  })
            }
      })
      
      #GSCAR Ranking
      observe({
            Maindata$AllRanking <- enrichfunc(Maindata$GSCARscore,Maindata$selectsample,Maindata$tab)
            Maindata$Ranking <- Maindata$AllRanking[as.numeric(Maindata$AllRanking[,"Adj.Pvalue"]) <= as.numeric(input$Inputpvalco) & as.numeric(Maindata$AllRanking[,"FoldChange"]) >= as.numeric(input$Inputfoldchangeco),]
            if (input$Mainmethod == 'GSCAR' & input$GSCARmethod == 'GSCARdefault') {
                  Maindata$defaultranking <- Maindata$Ranking
            } else if (input$Mainmethod == 'GSCAR' & input$GSCARmethod == 'GSCARinteractive'){
                  Maindata$interactranking <- Maindata$Ranking
            }    
      })
      
      #GSCAR sidebar
      output$InputNslider <- renderUI({
            maxval <- min(30,nrow(Maindata$Ranking))
            defaultval <- min(5,nrow(Maindata$Ranking))
            if (maxval > 0)
                  sliderInput("InputN","Number of Top Ranking Context Displayed",min=0,max=maxval,value=defaultval,step=1)
      })
      
      output$InputGSCARsidebar <- renderUI({
            if (!is.null(Maindata$patterndata)) 
                  if(nrow(Maindata$patterndata) != 0) {
                        if (Maindata$dim == 1) {
                              sliderInput(inputId = "GSCARoneslider", label = Maindata$patterndata[,1], min = min(Maindata$GSCARscore), max = max(Maindata$GSCARscore), value = c(min(Maindata$GSCARscore), max(Maindata$GSCARscore)))
                        } else if (Maindata$dim == 2) {
                              tagList(
                                    helpText("Build Polygon to Specify Interested Region"),
                                    p(actionButton('addpolygon', 'Add New Polygon')),
                                    p(actionButton('finishpolygon', 'Finish Drawing Polygon')),
                                    p(actionButton('undo', 'Undo last Operation')),
                                    p(actionButton('reset', 'Reset'))
                              )
                        } else {             
                              tagList(radioButtons("ThreeCutoffType","Choose Cutoff Type",
                                                   list("Sample"="Sample","Value"="Value")),
                                      conditionalPanel(condition="input.ThreeCutoffType == 'Value'",
                                                       lapply(1:Maindata$dim, function(i) {
                                                             sliderInput(inputId = paste0("GSCARthreevalueslider",i), label = Maindata$patterndata[i,1], min = min(Maindata$GSCARscore[i,]), max = max(Maindata$GSCARscore[i,]), value = c(min(Maindata$GSCARscore[i,]), max(Maindata$GSCARscore[i,])))
                                                       })
                                      ),
                                      p(actionButton("GSCARinteractivethreeupdate","Update Sample Selection"))
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
      output$InputGSCARspecifycontextui <- renderUI({
            if (!is.null(Maindata$Ranking))
                  selectInput("GSCARspecifycontext","Hold Control(Win) or Command(Mac) to enable multiple selection",Maindata$Ranking[,"SampleType"],multiple=T)
      })
      
      #Establish GSCAR Context
      observe({
            if(!is.null(Maindata$Ranking)) {
                  if (input$Inputcontexttype=='Toprank') {
                        if (nrow(Maindata$Ranking) == 0) {
                              Maindata$GSCARcontext <- NULL
                        } else {
                              Maindata$GSCARcontext <- Maindata$Ranking[1:as.numeric(input$InputN),"SampleType"]
                        }
                  } else {
                        Maindata$GSCARcontext <- input$GSCARspecifycontext 
                  }
            }
            
            if (input$Mainmethod == 'GSCAR' & input$GSCARmethod == 'GSCARdefault') {
                  Maindata$defaultcontext <- Maindata$GSCARcontext
            } else if (input$Mainmethod == 'GSCAR' & input$GSCARmethod == 'GSCARinteractive'){
                  Maindata$interactcontext <- Maindata$GSCARcontext
            }
      })
      
      #Establish Maindata$selectsample
      observe({
            if (input$Mainmethod == 'GSCAR' & input$GSCARmethod == 'GSCARinteractive' & !is.null(Maindata$dim)) {
                  if (Maindata$dim == 1) {
                        Maindata$selectsample <- which(Maindata$GSCARscore > input$GSCARoneslider[1] & Maindata$GSCARscore < input$GSCARoneslider[2])
                  } else if (Maindata$dim == 2) {
                        Maindata$selectsample <- inpoly$tf
                  } else {
                        if (!is.null(input$GSCARinteractivethreeupdate) & !is.null(input$ThreeCutoffType)) 
                              if (input$GSCARinteractivethreeupdate>0) {
                                    isolate({
                                          Maindata$selectsample <- GSCARthreeinfo$selectsample
                                    })
                              }
                  }
                  Maindata$interactsample <- Maindata$selectsample
            } else {
                  Maindata$selectsample <- Maindata$defaultsample
            }
      })
      
      #GSCAR default plot
      output$GSCARdefaultplot <- renderUI({
            if (!is.null(Maindata$patterndata)) 
                  if(nrow(Maindata$patterndata) != 0) {
                        if (Maindata$dim == 1) {
                              tagList(
                                    helpText(ifelse(is.null(Maindata$GSCARcontext),"No significantly enriched biological contexts found","")),
                                    plotOutput("GSCARdefaultplotone",height=300*(length(Maindata$GSCARcontext)+1))
                              )
                        } else if (Maindata$dim == 2) {
                              div(align="center",
                                  helpText(ifelse(is.null(Maindata$GSCARcontext),"No significantly enriched biological contexts found","")),
                                  plotOutput("GSCARdefaultplottwo",width=800,height=800),
                                  helpText(paste0("Correlation: ",Maindata$twocorr)),
                                  helpText(paste0("Pearson Correlation Test p-value: ",Maindata$twocorrp)),
                                  helpText(paste0("Regression Slope: ",Maindata$twoslope)),
                                  helpText(paste0("Regression Slope t-test p-value: ",Maindata$twoslopep))
                              )
                        } else {
                              tagList(
                                    plotOutput("GSCARdefaultplotthree"),
                                    helpText(ifelse(is.null(Maindata$GSCARcontext),"No significantly enriched biological contexts found","")),
                                    plotOutput("GSCARdefaultplotthreeplus")
                              )
                        }
                  }      
      })
      
      output$GSCARdefaultplotone <- renderPlot({
            if (Maindata$dim == 1) {
                  par(mfrow=c(length(Maindata$GSCARcontext)+1,1),oma=c(0,0,2,0))
                  hist(Maindata$GSCARscore,xlab="Sample Score",xlim=range(Maindata$GSCARscore),main="All Biological contexts")
                  abline(v=Maindata$cutoffval[1], lty=2)
                  if (!is.null(Maindata$GSCARcontext))
                        for(INDEX in Maindata$GSCARcontext) {
                              hist(Maindata$GSCARscore[Maindata$tab$SampleType %in% INDEX],xlim=range(Maindata$GSCARscore),main=substr(INDEX,1,25),xlab="Sample Score")
                              abline(v=Maindata$cutoffval[1], lty=2)
                        }     
            }
      })
      
      output$GSCARdefaultplottwo <- renderPlot({  
            if (Maindata$dim == 2) {
                  plot(Maindata$GSCARscore[1,],Maindata$GSCARscore[2,],col="#00000022",pch=20,xlab=Maindata$patterndata[1,1],cex=0.8,ylab=Maindata$patterndata[2,1])
                  toprankingsample <- NULL
                  if (!is.null(Maindata$GSCARcontext)) {
                        for(INDEX in Maindata$GSCARcontext) {
                              toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                        }
                  }
                  points(Maindata$GSCARscore[1,setdiff(Maindata$selectsample,toprankingsample)],Maindata$GSCARscore[2,setdiff(Maindata$selectsample,toprankingsample)],cex=0.8,pch=20)
                  abline(v=Maindata$cutoffval[1], h=Maindata$cutoffval[2], lty=2)
                  if (!is.null(Maindata$GSCARcontext)) {
                        i <- 1
                        for(INDEX in Maindata$GSCARcontext) {
                              if (input$Inputenrichedareaonly == TRUE) {
                                    points(Maindata$GSCARscore[1,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$selectsample)],Maindata$GSCARscore[2,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$selectsample)],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i])
                              } else {
                                    points(Maindata$GSCARscore[1,Maindata$tab$SampleType %in% INDEX],Maindata$GSCARscore[2,Maindata$tab$SampleType %in% INDEX],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i])
                              }
                              i <- i+1
                        }
                  }
                  leg.txt <- c(substr(Maindata$GSCARcontext,1,25),"Selected Samples","Not Selected Samples")
                  if (length(leg.txt)==2) {
                        legend("topleft",legend=leg.txt,pch=c(20,20),pt.bg=c("black","#00000022"),col=c("black","#00000022"),cex=0.8)                        
                  } else {
                        legend("topleft",legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=0.8)     
                  }
            }
      })
      
      output$GSCARdefaultplotthree <- renderPlot({
            if (Maindata$dim >= 3) {
                  par(oma=c(0.5,0,0.5,max(nchar(Maindata$patterndata[,1]))/2))
                  colcolorall <- rep("cyan",ncol(Maindata$GSCARscore))
                  colcolorall[Maindata$defaultsample] <- "blue"
                  heatmap.2(Maindata$GSCARscore,col=bluered,Colv=as.dendrogram(Maindata$GSCARclust),dendrogram="none",srtRow=-45,trace="none",Rowv=input$heatmapthreerowv,labCol=NA,ColSideColors=colcolorall,main="All Sample Heatmap")
                  legend("bottomleft",legend=c("Selected Sample","Unselected Sample"),lwd=1,col=c("blue","cyan"))       
            }
      })
      
      output$GSCARdefaultplotthreeplus <- renderPlot({
            if (Maindata$dim >= 3 && length(Maindata$defaultsample) > 0) {
                  colcolorselect <- rep("white",ncol(Maindata$GSCARscore))
                  par(oma=c(0.5,0,0.5,max(nchar(Maindata$patterndata[,1]))/2))
                  if (!is.null(Maindata$GSCARcontext)) {
                        i <- 1
                        for(INDEX in Maindata$GSCARcontext) {
                              colcolorselect[Maindata$tab$SampleType %in% INDEX] <- COLORS[i]
                              i <- i+1
                        }
                        heatmap.2(Maindata$GSCARscore[,Maindata$defaultsample],col=bluered,labCol=NA,Rowv=input$heatmapthreerowv,dendrogram="none",srtRow=-45,trace="none",ColSideColors=colcolorselect[Maindata$defaultsample],main="Selected Sample Heatmap")
                        leg.txt <- substr(Maindata$GSCARcontext,1,25)
                        legend("bottomleft",legend=leg.txt,lwd=1,col=COLORS,cex=0.8)
                  }
                  else {
                        heatmap.2(Maindata$GSCARscore[,Maindata$defaultsample],col=bluered,labCol=NA,Rowv=F,dendrogram="none",srtRow=-45,trace="none",main="Selected Sample Heatmap")
                  }
            }
      })
      
      #GSCAR interactive plot
      output$GSCARinteractiveplot <- renderUI({
            if (!is.null(Maindata$patterndata)) 
                  if(nrow(Maindata$patterndata) != 0) {
                        if (Maindata$dim == 1) {
                              tagList(
                                    helpText(ifelse(is.null(Maindata$GSCARcontext),"No significantly enriched biological contexts found","")),
                                    plotOutput("GSCARinteractiveplotone",height=300*(as.numeric(input$InputN)+1))
                              )
                        } else if (Maindata$dim == 2) {
                              div(align="center",
                                  plotOutput("GSCARinteractiveplottwo",clickId="coords",width=800,height=800),
                                  plotOutput("GSCARinteractiveplottwoplus",width=800,height=800),
                                  helpText(ifelse(is.null(Maindata$GSCARcontext),"No significantly enriched biological contexts found",""))
                              )      
                        } else {
                              tagList(                                    
                                    conditionalPanel(condition="input.ThreeCutoffType == 'Sample'",
                                                     helpText("Select sample range using sliders. The range would be the union set of multiple sliders"),
                                                     actionButton("GSCARthreesampleaddslider","Add Slider"),
                                                     actionButton("GSCARthreesampledeleteslider","Delete Slider"),
                                                     lapply(1:ifelse(is.null(GSCARthreeinfo$sampleslidernum),1,GSCARthreeinfo$sampleslidernum) , function(i) {
                                                           sliderInput(inputId = paste0("GSCARthreesampleslider",i),"",min=0,max=ncol(Maindata$GSCARscore),value=c(threesampleslidervalue[i,1],threesampleslidervalue[i,2]),step=1)
                                                     })
                                    ),
                                    p(plotOutput("GSCARinteractiveplotthreecolbar",height=20)),
                                    plotOutput("GSCARinteractiveplotthreeheatmap"),
                                    helpText(ifelse(is.null(Maindata$GSCARcontext),"No significantly enriched biological contexts found","")),
                                    plotOutput("GSCARinteractiveplotthreeplus"))
                        }
                  }      
      })
      
      GSCARthreeinfo <- reactiveValues()
      
      #initiate number of sample slider, record current slider value
      observe({
            if (is.null(GSCARthreeinfo$sampleslidernum))
                  GSCARthreeinfo$sampleslidernum <- 1
            for (i in 1:GSCARthreeinfo$sampleslidernum) 
                  if(!is.null(eval(parse(text=paste0("input$GSCARthreesampleslider",i,"[1]"))))) {
                        threesampleslidervalue[i,1] <<- eval(parse(text=paste0("input$GSCARthreesampleslider",i,"[1]"))) 
                        threesampleslidervalue[i,2] <<- eval(parse(text=paste0("input$GSCARthreesampleslider",i,"[2]")))
                  }
      })
      
      #set up action for add slider button in three genedata case
      observe({
            if (!is.null(input$GSCARthreesampleaddslider) && input$GSCARthreesampleaddslider>0)
                  isolate({
                        if (GSCARthreeinfo$sampleslidernum < 5) {
                              GSCARthreeinfo$sampleslidernum <- GSCARthreeinfo$sampleslidernum + 1
                        }
                  })
      })
      
      #set up action for delete slider button in three genedata case
      observe({
            if (!is.null(input$GSCARthreesampledeleteslider) && input$GSCARthreesampledeleteslider>0)
                  isolate({
                        if (GSCARthreeinfo$sampleslidernum != 1) {
                              GSCARthreeinfo$sampleslidernum <- GSCARthreeinfo$sampleslidernum - 1
                        }
                  })
      })
      
      output$GSCARinteractiveplotone <- renderPlot({
            if (Maindata$dim == 1) {
                  par(mfrow=c(as.numeric(input$InputN)+1,1),oma=c(0,0,2,0))
                  hist(Maindata$GSCARscore,xlab="Sample Score",xlim=range(Maindata$GSCARscore),main="All Biological contexts")
                  abline(v=c(input$GSCARoneslider[1],input$GSCARoneslider[2]), lty=2)
                  if (!is.null(Maindata$GSCARcontext)) {
                        for(INDEX in Maindata$GSCARcontext) {
                              hist(Maindata$GSCARscore[Maindata$tab$SampleType %in% INDEX],xlim=range(Maindata$GSCARscore),main=substr(INDEX,1,25),xlab="Sample Score")
                              abline(v=c(input$GSCARoneslider[1],input$GSCARoneslider[2]), lty=2)
                        }     
                  }
            }
      })
      
      inpoly <- reactiveValues()
      
      get.coords <- reactive({
            data.frame(x=input$coords$x, y=input$coords$y)
      })
      
      output$GSCARinteractiveplottwo <- renderPlot({  
            if (Maindata$dim == 2) { 
                  inputcoords <- get.coords()
                  plot(Maindata$GSCARscore[1,], Maindata$GSCARscore[2,],pch=20, xlim = range(Maindata$GSCARscore[1,]), ylim = range(Maindata$GSCARscore[2,]),xlab=Maindata$patterndata[1,1],ylab=Maindata$patterndata[2,1])
                  if (input$reset != resetvalue) {
                        polycord <<- NULL
                        polynum <<- 1
                        resetvalue <<- input$reset
                        actiontaken <<- 1
                        inpoly$tf <- rep(F,ncol(Maindata$GSCARscore))
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
                              inpoly$tf <- rep(F,ncol(Maindata$GSCARscore))
                              for (j in 1:polynum) {
                                    tmpcord <- polycord[polycord[,3]==j,]
                                    if (is.matrix(tmpcord)) {
                                          inpoly$tf[which(as.logical(point.in.polygon(Maindata$GSCARscore[1,],Maindata$GSCARscore[2,],tmpcord[,1],tmpcord[,2])))] <- T
                                    }
                              }
                        })
                        points(Maindata$GSCARscore[1,inpoly$tf], Maindata$GSCARscore[2,inpoly$tf],col = 'green',pch=20)      
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
      
      
      output$GSCARinteractiveplottwoplus <- renderPlot({  
            if (Maindata$dim == 2) { 
                  if (sum(inpoly$tf) != 0) {
                        plot(Maindata$GSCARscore[1,], Maindata$GSCARscore[2,], xlim = range(Maindata$GSCARscore[1,]), ylim = range(Maindata$GSCARscore[2,]),xlab=Maindata$patterndata[1,1],ylab=Maindata$patterndata[2,1],pch=20,cex=0.8,col="#00000022")            
                        for (j in 1:polynum) {
                              tmpcord <- polycord[polycord[,3]==j,]
                              if (is.matrix(tmpcord)) {
                                    for (i in 1:(nrow(tmpcord)-1)) 
                                          lines(c(tmpcord[i,1],tmpcord[i+1,1]),c(tmpcord[i,2],tmpcord[i+1,2]),type="l",col="blue",lty=2)
                                    lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lty=2)
                              }
                        }
                        toprankingsample <- NULL
                        if (!is.null(Maindata$GSCARcontext)) {
                              for(INDEX in Maindata$GSCARcontext) {
                                    toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                              }
                        }
                        points(Maindata$GSCARscore[1,setdiff(which(inpoly$tf),toprankingsample)], Maindata$GSCARscore[2,setdiff(which(inpoly$tf),toprankingsample)],cex=0.8,pch=20)
                        if (!is.null(Maindata$GSCARcontext)) {
                              i <- 1
                              for (INDEX in Maindata$GSCARcontext) {    
                                    if (input$Inputenrichedareaonly == TRUE) {
                                          points(Maindata$GSCARscore[1,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0], Maindata$GSCARscore[2,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0],col = COLORS[i], pch = STYLES[i], bg = COLORS[i])
                                    } else {
                                          points(Maindata$GSCARscore[1,Maindata$tab$SampleType %in% INDEX], Maindata$GSCARscore[2,Maindata$tab$SampleType %in% INDEX],col = COLORS[i], pch = STYLES[i],bg = COLORS[i])
                                    }
                                    i <- i+1
                              }
                        }
                        leg.txt <- c(substr(Maindata$GSCARcontext,1,25),"Selected Samples","Not Selected Samples")
                        if (length(leg.txt)==2) {
                              legend("topleft",legend=leg.txt,pch=c(20,20),pt.bg=c("black","#00000022"),col=c("black","#00000022"),cex=0.8)                        
                        } else {
                              legend("topleft",legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=0.8)                        
                        }
                  }
            }
      })
      
      observe({
            if (!is.null(input$ThreeCutoffType)) {
                  if (input$ThreeCutoffType == 'Sample') {
                        selectsample <- NULL
                        for (i in 1:GSCARthreeinfo$sampleslidernum) {
                              if(!is.null(eval(parse(text=paste0("input$GSCARthreesampleslider",i,"[1]"))))) 
                                    eval(parse(text=paste0("selectsample <- union(selectsample,Maindata$GSCARclust$order[input$GSCARthreesampleslider",i,"[1]:input$GSCARthreesampleslider",i,"[2]])")))
                        }
                        GSCARthreeinfo$selectsample <- selectsample
                        
                  } else if (input$ThreeCutoffType == 'Value') {
                        selectsample <- 1:ncol(Maindata$GSCARscore)
                        for (i in 1:Maindata$dim) {
                              eval(parse(text=paste0("selectsample <- intersect(selectsample,which(input$GSCARthreevalueslider",i,"[1]<Maindata$GSCARscore[",i,",] & input$GSCARthreevalueslider",i,"[2]>Maindata$GSCARscore[",i,",]))")))
                        }
                        GSCARthreeinfo$selectsample <- selectsample
                  }
            }
      })
      
      output$GSCARinteractiveplotthreecolbar <- renderPlot({
            if (Maindata$dim >= 3) { 
                  colcolorall <- rep("cyan",ncol(Maindata$GSCARscore))
                  colcolorall[GSCARthreeinfo$selectsample] <- "blue"
                  par(mar=c(0,0,0,0))
                  image(cbind(1:ncol(Maindata$GSCARscore)),col=colcolorall[Maindata$GSCARclust$order],axes=F)
                  par(mar=c(5.1,4.1,4.1,2.1))
            }
      })
      
      output$GSCARinteractiveplotthreeheatmap <- renderPlot({
            if (Maindata$dim >= 3) {
                  par(mar=c(0,0,0,0))
                  if (input$heatmapthreerowv) {
                        Rowv <- rowMeans(Maindata$GSCARscore)
                        hcr <- hclust(dist(Maindata$GSCARscore))
                        ddr <- as.dendrogram(hcr)
                        ddr <- reorder(ddr, Rowv)
                        rowInd <- order.dendrogram(ddr)
                        image(1:ncol(Maindata$GSCARscore),1:nrow(Maindata$GSCARscore),t(Maindata$GSCARscore[rowInd,Maindata$GSCARclust$order]),col=bluered(100),axes=F)
                  } else {
                        image(1:ncol(Maindata$GSCARscore),1:nrow(Maindata$GSCARscore),t(Maindata$GSCARscore[nrow(Maindata$GSCARscore):1,Maindata$GSCARclust$order]),col=bluered(100),axes=F)
                  }
                  par(mar=c(5.1,4.1,4.1,2.1))
            }
      })
      
      output$GSCARinteractiveplotthreeplus <- renderPlot({
            if (Maindata$dim >= 3 && length(Maindata$selectsample) > 0) {
                  if(input$GSCARinteractivethreeupdate>0) {
                        Maindata$Ranking
                        Maindata$GSCARcontext
                        input$heatmapthreerowv
                        isolate({
                              par(oma=c(0.5,0,0.5,max(nchar(Maindata$patterndata[,1]))/2))
                              colcolorselect <- rep("white",ncol(Maindata$GSCARscore))
                              if (!is.null(Maindata$GSCARcontext)) {
                                    i <- 1
                                    for(INDEX in Maindata$GSCARcontext) {
                                          colcolorselect[Maindata$tab$SampleType %in% INDEX] <- COLORS[i]
                                          i <- i+1
                                    }
                                    leg.txt <- substr(Maindata$GSCARcontext,1,25)
                                    heatmap.2(Maindata$GSCARscore[,Maindata$selectsample],col=bluered,labCol=NA,Rowv=input$heatmapthreerowv,dendrogram="none",srtRow=-45,trace="none",ColSideColors=colcolorselect[Maindata$selectsample],main="Selected Sample Heatmap")
                                    legend("bottomleft",legend=leg.txt,lwd=1,col=COLORS,cex=0.8)
                              } else {
                                    heatmap.2(Maindata$GSCARscore[,Maindata$selectsample],col=bluered,labCol=NA,Rowv=input$heatmapthreerowv,dendrogram="none",srtRow=-45,trace="none",main="Selected Sample Heatmap")
                              }
                        })
                  }
            }
      })
      
      output$GSCARrankingtable <- renderDataTable(Maindata$Ranking)
      
      #####   Mainmethod : Download   #####
      
      output$Downloadsidebarui <- renderUI({
            if (!is.null(Maindata$genedata)) {
                  if (Maindata$dim == 1) {
                        tagList(
                              wellPanel(
                                    selectInput("Downloadplotonetype","Select File Type for Plot",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotonefilename","Enter File Name","GSCAR Plot"),
                                    textInput("Downloadplotonefilewidth","Enter plot width (inches)",10),
                                    textInput("Downloadplotonefileheight","Enter plot height (inches)",3*(as.numeric(input$InputN)+1)),
                                    p(downloadButton("Downloadplotone","Save Plot"))
                              ),
                              wellPanel(
                                    helpText("Change Plotting Details"),
                                    textInput("Downloadmaintitleone","Enter Main Title","GSCAR result"),
                                    textInput("Downloadxlabone","Enter Title for X Axis","SampleScore"),
                                    textInput("Downloadylabone","Enter Title for Y Axis","Frequency"),
                                    textInput("Downloadxlimminone","Enter minimum value of X Axis",min(Maindata$GSCARscore)),
                                    textInput("Downloadxlimmaxone","Enter maximum value of X Axis",max(Maindata$GSCARscore)),
                                    selectInput("Downloadcolone","Select filling color",c("NULL","gray","red","blue","green","yellow","purple"))
                              )
                        )
                  } else if (Maindata$dim == 2) {
                        tagList(
                              wellPanel(
                                    selectInput("Downloadplottwotype","Select File Type for Plot",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplottwofilename","Enter File Name","GSCAR Plot"),
                                    textInput("Downloadplottwofilewidth","Enter plot width (inches)",15),
                                    textInput("Downloadplottwofileheight","Enter plot height (inches)",15),
                                    p(downloadButton("Downloadplottwo","Save Plot"))
                              ),
                              wellPanel(
                                    helpText("Change Plotting Details"),
                                    textInput("Downloadmaintitletwo","Enter Main Title","GSCAR result"),
                                    textInput("Downloadxlabtwo","Enter Title for X Axis",Maindata$patterndata[1,1]),
                                    textInput("Downloadylabtwo","Enter Title for Y Axis",Maindata$patterndata[2,1]),
                                    textInput("Downloadxlimmintwo","Enter minimum value of X Axis",min(Maindata$GSCARscore[1,])),
                                    textInput("Downloadxlimmaxtwo","Enter maximum value of X Axis",max(Maindata$GSCARscore[1,])),
                                    textInput("Downloadylimmintwo","Enter minimum value of Y Axis",min(Maindata$GSCARscore[2,])),
                                    textInput("Downloadylimmaxtwo","Enter maximum value of Y Axis",max(Maindata$GSCARscore[2,])),
                                    textInput("Downloadlegcextwo","Enter legend size",0.8),
                                    checkboxInput("Downloadcorvaluetwo","Show correlation and p-value")
                              )
                        )
                  } else {
                        tagList(
                              wellPanel(
                                    selectInput("Downloadplotthreeplotonetype","Select File Type for Heatmap One",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotthreeplotonefilename","Enter File Name","GSCAR Heatmap One"),
                                    textInput("Downloadplotthreeplotonefilewidth","Enter plot width (inches)",20),
                                    textInput("Downloadplotthreeplotonefileheight","Enter plot height (inches)",10),
                                    p(downloadButton("Downloadplotthreeplotone","Save Heatmap One")),
                                    selectInput("Downloadplotthreeplottwotype","Select File Type for Heatmap Two",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotthreeplottwofilename","Enter File Name","GSCAR Heatmap Two"),
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
                                    textInput("Downloadlegcexthreeone","Enter legend size in Heatmap One",0.8),
                                    textInput("Downloadlegcexthreetwo","Enter legend size in Heatmap Two",0.8)
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
                  if (input$Downloadregionselect == 'GSCARdefault') {
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
            hist(Maindata$GSCARscore,xlab=input$Downloadxlabone,ylab=input$Downloadylabone,xlim=as.numeric(c(input$Downloadxlimminone,input$Downloadxlimmaxone)),col=colone,main="All Biological contexts")
            if (input$Downloadregionselect == 'GSCARdefault') {
                  abline(v=Maindata$cutoffval[1], lty=2)      
            } else {
                  abline(v=c(input$GSCARoneslider[1],input$GSCARoneslider[2]), lty=2)
            }
            if (!is.null(Maindata$downloadcontext)) {
                  for(INDEX in Maindata$downloadcontext) {
                        hist(Maindata$GSCARscore[Maindata$tab$SampleType %in% INDEX],xlab=input$Downloadxlabone,ylab=input$Downloadylabone,xlim=as.numeric(c(input$Downloadxlimminone,input$Downloadxlimmaxone)),col=colone,main=substr(INDEX,1,25))
                        if (input$Downloadregionselect == 'GSCARdefault') {
                              abline(v=Maindata$cutoffval[1], lty=2)      
                        } else {
                              abline(v=c(input$GSCARoneslider[1],input$GSCARoneslider[2]), lty=2)
                        }
                  }     
            }
            title(input$Downloadmaintitleone,outer=T)
      }
      
      downloadtwofunc <- function() {
            cortext <- paste0("Correlation: ",round(Maindata$twocorr,3),"; ","Correlation p-value: ",round(Maindata$twocorrp,3),"; ","Slope: ",round(Maindata$twoslope,3),"; ","Slope p-value: ",round(Maindata$twoslopep,3))
            if (input$Downloadregionselect == 'GSCARdefault') {
                  plot(Maindata$GSCARscore[1,],Maindata$GSCARscore[2,],col="#00000022",pch=20,cex=0.8,xlab=input$Downloadxlabtwo,ylab=input$Downloadylabtwo,xlim = as.numeric(c(input$Downloadxlimmintwo,input$Downloadxlimmaxtwo)), ylim = as.numeric(c(input$Downloadylimmintwo,input$Downloadylimmaxtwo)),main=input$Downloadmaintitletwo)
                  toprankingsample <- NULL
                  if (!is.null(Maindata$downloadcontext)) {
                        for(INDEX in Maindata$downloadcontext) {
                              toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                        }
                  }
                  points(Maindata$GSCARscore[1,setdiff(Maindata$downloadsample,toprankingsample)],Maindata$GSCARscore[2,setdiff(Maindata$downloadsample,toprankingsample)],cex=0.8,pch=20)
                  abline(v=Maindata$cutoffval[1], h=Maindata$cutoffval[2], lty=2)
                  if (!is.null(Maindata$downloadcontext)) {
                        i <- 1
                        for(INDEX in Maindata$downloadcontext) {
                              if (input$Inputenrichedareaonly == TRUE) {
                                    points(Maindata$GSCARscore[1,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$downloadsample)],Maindata$GSCARscore[2,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$downloadsample)],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i])
                              } else {
                                    points(Maindata$GSCARscore[1,Maindata$tab$SampleType %in% INDEX],Maindata$GSCARscore[2,Maindata$tab$SampleType %in% INDEX],
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
                        plot(Maindata$GSCARscore[1,], Maindata$GSCARscore[2,], xlim = as.numeric(c(input$Downloadxlimmintwo,input$Downloadxlimmaxtwo)), ylim = as.numeric(c(input$Downloadylimmintwo,input$Downloadylimmaxtwo)),xlab=input$Downloadxlabtwo,ylab=input$Downloadylabtwo,col="#00000022",pch=20,cex=0.8,main=input$Downloadmaintitletwo)            
                        for (j in 1:polynum) {
                              tmpcord <- polycord[polycord[,3]==j,]
                              if (is.matrix(tmpcord)) {
                                    for (i in 1:(nrow(tmpcord)-1)) 
                                          lines(c(tmpcord[i,1],tmpcord[i+1,1]),c(tmpcord[i,2],tmpcord[i+1,2]),type="l",col="blue",lty=2)
                                    lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lty=2)
                              }
                        }
                        toprankingsample <- NULL
                        if (!is.null(Maindata$downloadcontext)) {
                              for(INDEX in Maindata$downloadcontext) {
                                    toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                              }
                        }
                        points(Maindata$GSCARscore[1,setdiff(which(inpoly$tf),toprankingsample)], Maindata$GSCARscore[2,setdiff(which(inpoly$tf),toprankingsample)],cex=0.8,pch=20)
                        if (!is.null(Maindata$downloadcontext)) {
                              i <- 1
                              for (INDEX in Maindata$downloadcontext) {    
                                    if (input$Inputenrichedareaonly == TRUE) {
                                          points(Maindata$GSCARscore[1,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0], Maindata$GSCARscore[2,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0],col = COLORS[i], pch = STYLES[i], bg = COLORS[i])
                                    } else {
                                          points(Maindata$GSCARscore[1,Maindata$tab$SampleType %in% INDEX], Maindata$GSCARscore[2,Maindata$tab$SampleType %in% INDEX],col = COLORS[i], pch = STYLES[i], bg = COLORS[i])
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
            colcolorall <- rep("cyan",ncol(Maindata$GSCARscore))
            colcolorall[Maindata$downloadsample] <- "blue"
            heatmap.2(Maindata$GSCARscore,col=threepalette(input$Downloadthreecolplotone),Colv=as.dendrogram(Maindata$GSCARclust),srtRow=ifelse(input$Downloadthreerotatelabplotone,-45,0),dendrogram=threeonedendro,trace="none",Rowv=input$Downloadthreerowvplotone,labCol=NA,ColSideColors=colcolorall,main=input$Downloadmaintitlethreeplotone)
            legend("bottomleft",legend=c("Selected Sample","Unselected Sample"),lwd=1,col=c("blue","cyan"),cex=as.numeric(input$Downloadlegcexthreeone))           
      }
      
      downloadthreetwofunc <- function() {
            colcolorselect <- rep("white",ncol(Maindata$GSCARscore))
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
                  heatmap.2(Maindata$GSCARscore[,Maindata$downloadsample],col=threepalette(input$Downloadthreecolplottwo),labCol=NA,Rowv=input$Downloadthreerowvplottwo,srtRow=ifelse(input$Downloadthreerotatelabplottwo,-45,0),dendrogram=threetwodendro,trace="none",ColSideColors=colcolorselect[Maindata$downloadsample],main=input$Downloadmaintitlethreeplottwo)
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




