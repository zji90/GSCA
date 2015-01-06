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
#require(shinyRGL)
#require(rgl)
require(rhdf5)
require(GSCA)
###two genedata case initiate
polycord <- NULL
resetvalue <- 0
finishpolygonvalue <- 0
addpolygonvalue <- 0
undovalue <- 0
actiontaken <- 0
finishpolygonaction <- 0
polynum <- 1

#one and three genedata case sample selection initiate
onesampleslidervalue <- matrix(0,nrow=5,ncol=3)
threesampleslidervalue <- matrix(0,nrow=5,ncol=3)

###default styles and colors
STYLES <- rep(c(15:18,3:4),5)
COLORS <- c(brewer.pal(5,"Set1"),brewer.pal(5,"Pastel1"),brewer.pal(10,"Paired"),brewer.pal(10,"BrBG"))


shinyServer(function(input, output, session) {
      
      #Establishing GSCA status panel. Open to change in the future with more advanced shiny options available.
      GSCAstatus <- reactiveValues()
      
      observe({
            if (is.null(Maindata$genedata))
                  GSCAstatus$genedata <- NULL
            if (is.null(Maindata$GSCAcontext))
                  GSCAstatus$GSCAcontext <- NULL
            if ((input$Mainmethod=='GSCA' && !identical(GSCAstatus$genedata,Maindata$genedata)) 
                | input$Mainmethod=='Download' && !identical(GSCAstatus$GSCAcontext,Maindata$GSCAcontext)) {
                  isolate({
                        GSCAstatus$status <- 1
                        GSCAstatus$genedata <- Maindata$genedata
                        GSCAstatus$GSCAcontext <- Maindata$GSCAcontext
                  }) 
            } else {
                  isolate(GSCAstatus$status <- 0)
            }   
      },priority=10)
      
      output$GSCAstatusui <- renderUI({
            if (!is.null(GSCAstatus$status) && GSCAstatus$status) {
                  tagList(
                        tags$head(
                              tags$link(rel="stylesheet", type="text/css",href="style.css"),
                              tags$script(type="text/javascript", src = "busy.js")
                        ),
                        
                        div(class = "busy",  
                            p("Calculation in progress.."), 
                            img(src="ajaxloaderq.gif")
                        )
                  )        
            } else {
                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                   tags$head(
                                         tags$link(rel="stylesheet", type="text/css",href="style.css"),
                                         tags$script(type="text/javascript", src = "busy.js")
                                   ),                                   
                                   div(class = "busy",  
                                       p("Calculation in progress.."), 
                                       img(src="ajaxloaderq.gif")
                                   ))                  
            }
      })
      
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
            if (!is.null(Rawdata$genedata))
                  if (input$InputGenesetname %in% Rawdata$genedata[,1]) {
                        i <- 2
                        modifygenesetname <- paste(currentgenesetname,i,sep="_")
                        while (modifygenesetname %in% Rawdata$genedata[,1]) {
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
                        Currentdata$genedata <- data.frame(Genesetname=Currentdata$genesetname,GeneID=c(strsplit(input$InputActGeneID,";")[[1]],strsplit(input$InputRepGeneID,";")[[1]]),weight=c(rep(1,length(strsplit(input$InputActGeneID,";")[[1]])),rep(-1,length(strsplit(input$InputRepGeneID,";")[[1]]))),stringsAsFactors=F)
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
                                    if (sum(is.na(suppressWarnings(as.numeric(Currentdata$GenesetFile[,2]))))==0) {
                                          Currentdata$genedata <- data.frame(Genesetname=Currentdata$GenesetFile[,3],GeneID=Currentdata$GenesetFile[,1],weight=as.numeric(Currentdata$GenesetFile[,2]),stringsAsFactors=F)
                                    } else {
                                          Currentdata$genedata <- NULL
                                    } 
                              } else if (ncol(Currentdata$rawGenesetFile) == 2 || (ncol(Currentdata$rawGenesetFile) > 2 && !is.null(input$InputGenesetcolnumthree) && !input$InputGenesetcolnumthree)) {
                                    Currentdata$GenesetFile <- Currentdata$rawGenesetFile[,1:2]
                                    if (sum(is.na(suppressWarnings(as.numeric(Currentdata$GenesetFile[,2]))))==0) {
                                          Currentdata$genedata <- data.frame(Genesetname=Currentdata$genesetname,GeneID=Currentdata$GenesetFile[,1],weight=as.numeric(Currentdata$GenesetFile[,2]),stringsAsFactors=F)
                                    } else {
                                          Currentdata$genedata <- NULL
                                    }
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
                  })
            }
      })
      
      output$OutputCurrentGenedatawarnui <- renderUI({
            if (input$InputGenesetmethod == 'InputuploadGeneset' && is.null(Currentdata$genedata))
                  helpText("If nothing appears, make sure entries in weight column are all numeric and try the followings: change the separator; change the quote; check the header option")
      })
      
      output$OutputCurrentGenedata <- renderDataTable(Currentdata$genedata)
      
      output$OutputAllGenedata <- renderDataTable({
            if (!is.null(Rawdata$genedata) && nrow(Rawdata$genedata)>0) {
                  sumtable <- data.frame(Genesetname=unique(Rawdata$genedata[,1]),Totalgene=0,Activatedgene=0,Repressedgene=0)
                  for (i in 1:nrow(sumtable)) {
                        temp <- Rawdata$genedata[Rawdata$genedata[,1]==sumtable[i,1],]
                        sumtable[i,2] <- nrow(temp)
                        sumtable[i,3] <- sum(temp[,3]>0)
                        sumtable[i,4] <- sum(temp[,3]<0)
                  }
                  sumtable      
            }
      })    
      
      output$Indigenesetnameui <- renderUI({
            if (!is.null(Rawdata$genedata)) {
                  namelist <- list()
                  for (i in unique(Rawdata$genedata$Genesetname))
                        eval(parse(text=paste0("namelist <- c(namelist,`",i,"`=i)")))
                  checkboxGroupInput("Indigenesetname","Choose Gene Set",namelist)
            }
      })
      
      output$Indigeneset <- renderDataTable({
            if (!is.null(Rawdata$genedata)) { 
                  tmp <- Rawdata$genedata[Rawdata$genedata[,1] %in% input$Indigenesetname,]
                  if (nrow(tmp) > 0)
                        tmp
            }
      })
      
      output$OutputGenedataname <- renderText({
            if (is.null(Rawdata$genedata)){
                  "No gene set"
            } else {
                  do.call("paste",c(as.list(unique(Rawdata$genedata$Genesetname)),sep=","))
            }
      })
      
      #Delete genedata
      output$Inputgenesetdeleteui <- renderUI({
            if (!is.null(Rawdata$genedata))
                  selectInput("Inputdeletegenesetselect","",unique(Rawdata$genedata$Genesetname),multiple=T)
      })
      
      observe({   
            if(input$Inputgenesetdelete != 0) {
                  isolate({
                        Rawdata$genedata <- Rawdata$genedata[!Rawdata$genedata$Genesetname %in% input$Inputdeletegenesetselect,]
                  })
            }
      })      
      
      observe({
            input$Inputgenesetreset
            isolate({
                  Rawdata$genedata <- NULL
            })
      })
      
      output$Savegenedatafile <- downloadHandler(
            filename = function() { "Genesetfile.csv" },
            content = function(file) {
                  tmp <- data.frame(ENTREZ_ID=Rawdata$genedata[,2],weight=Rawdata$genedata[,3],Genesetname=Rawdata$genedata[,1])
                  write.csv(tmp,file,row.names=F)     
            }
      )
      
      ######  Mainmethod : Summary  ######
      
      Maindata <- reactiveValues()
      
      #Sidebar checkbox group for selecting dataset
      output$Summarydataselect <- renderUI({
            if (!is.null(Rawdata$genedata)) {
                  checkboxGroupInput("Selectedgeneset","",unique(Rawdata$genedata[,1]),selected=unique(Rawdata$genedata[,1]))
            }
      })
      
      #Change geneset order
      output$Summaryswitchorderui <- renderUI({
            tagList(
                  selectInput("Summaryswitchgs1","Gene Set 1",Rawdata$genedata$Genesetname),
                  selectInput("Summaryswitchgs2","Gene Set 2",Rawdata$genedata$Genesetname)
            )
      })
      
      observe({
            if(input$Summaryswitchbut > 0) {
                  isolate({              
                        if (nchar(input$Summaryswitchgs1) > 0 && nchar(input$Summaryswitchgs2) > 0) {
                              id1 <- which(Rawdata$genedata$Genesetname==input$Summaryswitchgs1)
                              id2 <- which(Rawdata$genedata$Genesetname==input$Summaryswitchgs2)
                              if (id1[1] < id2[1]) {
                                    header <- id1
                                    tailer <- id2
                              } else {
                                    header <- id2
                                    tailer <- id1
                              }
                              if(header[1] != 1) {
                                    bfhead <- 1:(header[1]-1)
                              } else {
                                    bfhead <- NULL
                              }
                              if (tailer[length(tailer)] == nrow(Rawdata$genedata)) {
                                    aftail <- NULL
                              } else {
                                    aftail <- (tailer[length(tailer)] + 1):nrow(Rawdata$genedata)
                              }
                              if (header[length(header)] == tailer[1] - 1) {
                                    midpart <- NULL
                              } else {
                                    midpart <- (header[length(header)]+1):(tailer[1] - 1)
                              }
                              tmp <- c(bfhead,tailer,midpart,header,aftail)
                              Rawdata$genedata <- Rawdata$genedata[tmp,]
                        }
                  })
            }
      })
      
      #UI for selecting compendium
      output$Summarycompselectui <- renderUI({
            complist <- list()
            if (require(Affyhgu133aExpr)) {
                  complist <- c(complist,"Affymetrix Human Genome U133A Array, GPL96 (11778 samples)"="hgu133a")
            } 
            if (require(Affymoe4302Expr)) {
                  complist <- c(complist,"Affymetrix Mouse Genome 430 2.0 Array, GPL1261 (9444 samples)"="moe4302")
            }
            if (require(Affyhgu133Plus2Expr)) {
                  complist <- c(complist,"Affymetrix Human Genome U133 Plus 2.0 Array, GPL570 (5153 samples)"="hgu133Plus2")
            }
            if (require(Affyhgu133A2Expr)) {
                  complist <- c(complist,"Affymetrix Human Genome U133A 2.0 Array, GPL571 (313 samples)"="hgu133A2")
            }
            radioButtons("Summarycompselect","Select Compendium", complist)
      })
      
      output$Summarycompinfo <- renderUI({
            if (!is.null(input$Summarycompselect)) {
                  if(input$Summarycompselect=="moe4302"){
                        p(helpText("This compendium contains 20630 genes"),a("NCBI GEO description",href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1261",target="_blank"))
                  } else if (input$Summarycompselect=="hgu133a"){
                        p(helpText("This compendium contains 12495 genes"),a("NCBI GEO description",href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96",target="_blank"))
                  } else if (input$Summarycompselect=="hgu133Plus2"){
                        p(helpText("This compendium contains 19944 genes"),a("NCBI GEO description",href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570",target="_blank"))
                  } else if (input$Summarycompselect=="hgu133A2"){
                        p(helpText("This compendium contains 12494 genes"),a("NCBI GEO description",href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL571",target="_blank"))
                  }    
            }
      })
      
      #Update Maindata information
      observe({   
            if (input$Summarycompmethod=='available' && !is.null(input$Summarycompselect)) {
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
                  }
                  path <- system.file("extdata",package=paste0("Affy",input$Summarycompselect,"Expr"))
                  load(paste0(path,"/geneid.rda"))
            } else {
                  input$Summaryuploadtabfile
                  if (!is.null(input$Summaryuploadtabfile)) {
                        tmptab <- read.table(input$Summaryuploadtabfile$datapath,stringsAsFactors=F,blank.lines.skip=TRUE)
                        Maindata$tab <- data.frame(SampleID=tmptab[,1],ExperimentID=tmptab[,2],SampleType=tmptab[,3])
                  }
                  if (!is.null(input$Summaryuploadgeneexprfile))
                        Maindata$uploadgeneexpr <- as.matrix(read.table(input$Summaryuploadgeneexprfile$datapath,stringsAsFactors=F,blank.lines.skip=TRUE,row.names=1))
                  geneid <- row.names(Maindata$uploadgeneexpr)
            }
            if (!is.null(Rawdata$genedata)) {
                  Maindata$geneid <- geneid
                  Maindata$genedata <- Rawdata$genedata[Rawdata$genedata[,2] %in% geneid & Rawdata$genedata[,1] %in% input$Selectedgeneset,]
                  Maindata$dim <- length(unique((Maindata$genedata[,1])))
                  Maindata$genesetname <- unique(Maindata$genedata[,1])
            }            
      })
      
      #Render summary of selected dataset
      output$OutputDataSummary <- renderDataTable({
            if (!is.null(Maindata$genedata))
                  if(nrow(Maindata$genedata) != 0) {
                        sumtable <- data.frame(Genesetname=Maindata$genesetname,Originalgene=0,Compendiumgene=0)
                        for (i in 1:nrow(sumtable)) {
                              sumtable[i,2] <- sum(Rawdata$genedata[,1] == sumtable[i,1])
                              sumtable[i,3] <- sum(Maindata$genedata[,1] == sumtable[i,1])
                        }
                        names(sumtable) <- c("Gene Set","Original Number of Genes","Number of Genes in Compendium")
                        sumtable
                  }
      })    
      
      #Missing genedata report
      output$Outputmissinggenesetreport <- renderUI({
            if (!is.null(Maindata$dim))
                  if (Maindata$dim != length(input$Selectedgeneset)) {
                        tagList(
                              br(h4("Following gene sets do not contain any gene in the selected compendium, thus omitted:")),
                              textOutput("Outputmissinggenesetreporttext")
                        )
                  }
      })
      
      output$Outputmissinggenesetreporttext <- renderText({
            do.call("paste",c(as.list(setdiff(input$Selectedgeneset,Maindata$genesetname)),sep=","))
      })
      
      #Geneset breakdown
      
      genebreakdata <- reactiveValues()
      
      output$genesetbreakdownnameui <- renderUI({
            if (!is.null(Maindata$genedata)) {
                  namelist <- list()
                  for (i in Maindata$genesetname) 
                        if (sum(Maindata$genedata$Genesetname==i) > 1)
                              eval(parse(text=paste0("namelist <- c(namelist,`",i,"`=i)")))
                  if (length(namelist) > 0)
                        radioButtons("genesetbreakdownchoose","Choose one of the gene sets with at least two genes",namelist)
            }
      })
      
      output$genesetbreakdowntreenumui <- renderUI({
            if (!is.null(genebreakdata$expr))
                  sliderInput("genesetbreakdowntreenum","Choose number of clusters",min=1,max=min(12,nrow(genebreakdata$expr)),value=1,step=1)
      })
      
      observe({
            if (!is.null(Maindata$genedata) && nrow(Maindata$genedata) > 1) {
                  genebreakdata$geneset <- tmpgeneset <- Maindata$genedata[Maindata$genedata$Genesetname==input$genesetbreakdownchoose,]
                  
                  if (nrow(tmpgeneset) > 1) {
                        if (input$Summarycompmethod=='available') {
                              path <- system.file("extdata",package=paste0("Affy",input$Summarycompselect,"Expr"))
                              tmpgeneexpr <- t(h5read(paste0(path,"/data.h5"),"expr",index=list(NULL,match(tmpgeneset[,2],Maindata$geneid))))/1000
                        } else {
                              tmpgeneexpr <- Maindata$uploadgeneexpr[tmpgeneset[,2],]
                        }        
                        if (input$Summarycompscale) {                              
                              if(input$Summarycompscalemet=="zmuv") {
                                    tmpgeneexpr <- t(apply(tmpgeneexpr,1,scale))   
                              } else if(input$Summarycompscalemet=="zm") {
                                    tmpgeneexpr <- t(apply(tmpgeneexpr,1,scale,scale=F))
                              } else if(input$Summarycompscalemet=="uv") {
                                    tmpgeneexpr <- t(apply(tmpgeneexpr,1,scale,center=F))
                              }
                        }                              
                        tmpgeneexpr <- sweep(tmpgeneexpr,1,tmpgeneset[,3],"*")       
                        rownames(tmpgeneexpr) <- tmpgeneset[,2]
                        genebreakdata$expr <- tmpgeneexpr
                        genebreakdata$hc <- hclust(dist(tmpgeneexpr))
                  }
            }
      })
      
      observe({
            if (!is.null(genebreakdata$hc) & !is.null(input$genesetbreakdowntreenum))
                  genebreakdata$cutree <- cutree(genebreakdata$hc,k=as.numeric(input$genesetbreakdowntreenum))      
      })
      
      output$genesetbreakdownclustplot <- renderPlot({
            if (!is.null(genebreakdata$hc) && !is.null(input$genesetbreakdowntreenum)) {          
                  treenum <- as.numeric(input$genesetbreakdowntreenum)
                  if (treenum == 1) {
                        labelColors <- "black"
                  } else if (treenum == 2) {
                        labelColors <- c("black","red")
                  } else {
                        labelColors <- brewer.pal(treenum,"Paired")
                  }
                  colLab <- function(n) {
                        if (is.leaf(n)) {
                              a <- attributes(n)
                              labCol <- labelColors[genebreakdata$cutree[a$label]]
                              attr(n, "nodePar") <- c(a$nodePar, lab.col=labCol)
                        }
                        n
                  }
                  clusDendro <- dendrapply(as.dendrogram(genebreakdata$hc), colLab)
                  plot(clusDendro)
            }
      })
      
      observe({
            if (input$genesetbreakdownaddbutton > 0)
                  isolate({
                        tmpgeneset <- genebreakdata$geneset  
                        tmpnum <- 1
                        while (paste0(tmpgeneset[1,1],"_breakdown",tmpnum,"_",genebreakdata$cutree[as.character(tmpgeneset[1,2])]) %in% Rawdata$genedata$Genesetname) {
                              tmpnum <- tmpnum + 1
                        }
                        for (i in 1:nrow(tmpgeneset)) {                              
                              tmpgeneset[i,1] <- paste0(tmpgeneset[i,1],"_breakdown",tmpnum,"_",genebreakdata$cutree[as.character(tmpgeneset[i,2])])
                        }
                        Rawdata$genedata <- rbind(Rawdata$genedata,tmpgeneset)
                  })
      })
      
      ######  Mainmethod : GSCA  ######
      
      #When reentering GSCA page, GSCAmethod should be set to GSCAdefault
      observe({
            if (input$Mainmethod == "Input" | input$Mainmethod == "Select" | input$Mainmethod == "Download")
                  updateRadioButtons(session,"GSCAmethod",choices=list("Numeric POI"="GSCAdefault","Interactive POI"="GSCAinteractive"),selected="GSCAdefault")
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
                        message("No significantly enriched biological contexts found.")
                  } else {
                        FIN <- cbind(1:nrow(FIN),FIN)
                        colnames(FIN)[1] <- "Rank"
                        rownames(FIN) <- 1:nrow(FIN)
                  }
            }
            return(FIN)
      }      
      
      #Calculating GSCA Score and hclust dendrogram if there are more than three genesets
      observe({
            if (input$Mainmethod=='GSCA') {
                  isolate({
                        if (!is.null(Maindata$genedata) && nrow(Maindata$genedata) != 0) {
                              scoremat <- matrix(0,nrow=Maindata$dim,ncol=nrow(Maindata$tab))
                              row.names(scoremat) <- rep("tmp",Maindata$dim)
                              for (genesetid in 1:Maindata$dim) {
                                    ###Calculate Score
                                    singlegeneset <- Maindata$genesetname[genesetid]
                                    currentgeneset <- Maindata$genedata[Maindata$genedata[,1] == singlegeneset,]
                                    if (input$Summarycompmethod=='available') {
                                          path <- system.file("extdata",package=paste0("Affy",input$Summarycompselect,"Expr"))
                                          tmpgeneexpr <- t(h5read(paste0(path,"/data.h5"),"expr",index=list(NULL,match(currentgeneset[,2],Maindata$geneid))))/1000
                                    } else {
                                          tmpgeneexpr <- Maindata$uploadgeneexpr[currentgeneset[,2],]
                                    }        
                                    if (input$Summarycompscale) {                              
                                          if(input$Summarycompscalemet=="zmuv") {
                                                tmpgeneexpr <- t(apply(tmpgeneexpr,1,scale))   
                                          } else if(input$Summarycompscalemet=="zm") {
                                                tmpgeneexpr <- t(apply(tmpgeneexpr,1,scale,scale=F))
                                          } else if(input$Summarycompscalemet=="uv") {
                                                tmpgeneexpr <- t(apply(tmpgeneexpr,1,scale,center=F))
                                          }
                                    }  
                                    tmpgeneexpr <- sweep(tmpgeneexpr,1,currentgeneset[,3],"*")            
                                    if (input$Summarygenesetactmethod == "average") {
                                          score <- colMeans(tmpgeneexpr)         
                                    } else {
                                          score <- apply(tmpgeneexpr,2,median)
                                    }
                                    scoremat[genesetid,] <- score
                                    row.names(scoremat)[genesetid] <- singlegeneset
                              }
                              Maindata$GSCAscore <- scoremat
                              if (Maindata$dim >= 3) {
                                    Maindata$GSCAclust <- hclust(dist(t(scoremat)))
                                    Rowv <- rowMeans(Maindata$GSCAscore)
                                    hcr <- hclust(dist(Maindata$GSCAscore))
                                    ddr <- as.dendrogram(hcr)
                                    Maindata$GSCArowclust <- reorder(ddr, Rowv)
                              }
                              if (onesampleslidervalue[1,3] != ncol(Maindata$GSCAscore)) {
                                    onesampleslidervalue[,3] <<- rep(ncol(Maindata$GSCAscore),5)
                                    onesampleslidervalue[,1] <<- rep(0,5)
                                    onesampleslidervalue[,2] <<- rep(ncol(Maindata$GSCAscore),5)
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
      
      #Numeric POI
      output$numericpoiui <- renderUI({
            if (!is.null(Maindata$dim) && input$Mainmethod=='GSCA') {
                  if (input$numericpoimethod == "slider") {
                        lapply(1:Maindata$dim, function(i) {
                              tmpmin <- mean(Maindata$GSCAscore[i,])+2*sd(Maindata$GSCAscore[i,])
                              sliderInput(inputId = paste0("GSCAnumericpoislider",i), label = Maindata$genesetname[i], min = min(Maindata$GSCAscore[i,]), max = max(Maindata$GSCAscore[i,]), value = c(tmpmin, max(Maindata$GSCAscore[i,])))                        
                        })      
                  } else {
                        lapply(1:Maindata$dim, function(i) {
                              tmpmin <- mean(Maindata$GSCAscore[i,])+2*sd(Maindata$GSCAscore[i,])
                              list(textInput(inputId = paste0("GSCAnumericpoitextlower",i), label = paste(Maindata$genesetname[i],"Lower Bound"), tmpmin)
                                   ,textInput(inputId = paste0("GSCAnumericpoitextupper",i), label = paste(Maindata$genesetname[i],"Upper Bound"), max(Maindata$GSCAscore[i,])))
                        })
                  }
            }                  
      })
      
      output$numericpoitext <- renderUI({
            if (!is.null(Maindata$dim)) {
            lapply(1:Maindata$dim, function(i) {
                  helpText(paste0(Maindata$genesetname[i],": ",round(mean(eval(parse(text=paste0("input$GSCAnumericpoislider",i,"[1]")))>Maindata$GSCAscore[i,]),3),"-",round(mean(eval(parse(text=paste0("input$GSCAnumericpoislider",i,"[2]")))>Maindata$GSCAscore[i,]),3)))
            })
            }
      })
      
      observe({
            if (!is.null(Maindata$dim) && !is.null(Maindata$GSCAscore) && input$Mainmethod=='GSCA') {
                  selectsample <- 1:nrow(Maindata$tab)
                  cutoffval <- matrix(0,Maindata$dim,2)
                  if (input$numericpoimethod == "slider") {
                        for (i in 1:Maindata$dim) {
                              if(!is.null(eval(parse(text=paste0("input$GSCAnumericpoislider",i,"[1]"))))) {
                                    eval(parse(text=paste0("cutoffval[",i,",1] <- input$GSCAnumericpoislider",i,"[1]")))
                                    eval(parse(text=paste0("cutoffval[",i,",2] <- input$GSCAnumericpoislider",i,"[2]")))
                                    selectsample <- intersect(selectsample,which(Maindata$GSCAscore[i,] >= cutoffval[i,1] & Maindata$GSCAscore[i,] <= cutoffval[i,2]))
                              }
                        }            
                  } else {
                        for (i in 1:Maindata$dim) {
                              if(!is.null(eval(parse(text=paste0("input$GSCAnumericpoitextlower",i))))) {
                                    eval(parse(text=paste0("cutoffval[",i,",1] <- as.numeric(input$GSCAnumericpoitextlower",i,")")))
                                    eval(parse(text=paste0("cutoffval[",i,",2] <- as.numeric(input$GSCAnumericpoitextupper",i,")")))
                                    selectsample <- intersect(selectsample,which(Maindata$GSCAscore[i,] >= cutoffval[i,1] & Maindata$GSCAscore[i,] <= cutoffval[i,2]))
                              }
                        }   
                  }
                  Maindata$defaultsample <- selectsample
                  Maindata$cutoffval <- cutoffval
            }
      })
      
      output$numericpoimoreopgenesetnameui <- renderUI({
            if (!is.null(Maindata$genesetname))
                  tagList(
                        selectInput("numericpoimoreopselectgene","Select Gene Set",Maindata$genesetname)
                  )
      })
      
      #Numeric POI more options
      observe({
            if (input$numericpoimoreopbutton > 0)
                  isolate({
                        tmpscore <- Maindata$GSCAscore[input$numericpoimoreopselectgene,]
                        id <- which(Maindata$genesetname == input$numericpoimoreopselectgene)
                        tmpinput <- min(as.numeric(input$numericpoimoreopvalue),1)
                        tmpinput <- max(0,tmpinput)
                        if (input$numericpoimoreopcutofftype == "sd") {
                              value <- mean(tmpscore) + as.numeric(input$numericpoimoreopvalue) * sd(tmpscore)
                        } else if (input$numericpoimoreopcutofftype == "Norm") {
                              if (tmpinput == 0) {
                                    value <- min(tmpscore)
                              } else if (tmpinput == 1) {
                                    value <- max(tmpscore)
                              } else {
                                    value <- qnorm(tmpinput,mean(tmpscore),sd(tmpscore))      
                              }
                        } else if (input$numericpoimoreopcutofftype == "Quantile") {
                              value <- as.numeric(quantile(tmpscore,tmpinput))
                        }
                        if (input$numericpoimoreopbound == 'Upper') {
                              if (input$numericpoimethod == "slider") {
                                    eval(parse(text=paste0('updateSliderInput(session,"GSCAnumericpoislider', id ,'",label = Maindata$genesetname[',id,'],value=c(input$GSCAnumericpoislider', id, '[1],value))')))                                    
                              } else {
                                    eval(parse(text=paste0('updateTextInput(session,"GSCAnumericpoitextupper', id ,'",label = paste(Maindata$genesetname[',id,'],"Upper Bound"),value=value)')))                                                            
                              }
                        } else {
                              if (input$numericpoimethod == "slider") {
                                    eval(parse(text=paste0('updateSliderInput(session,"GSCAnumericpoislider', id ,'",label = Maindata$genesetname[',id,'],value=c(value,input$GSCAnumericpoislider', id, '[2]))')))                                    
                              } else {
                                    eval(parse(text=paste0('updateTextInput(session,"GSCAnumericpoitextlower', id ,'",label = paste(Maindata$genesetname[',id,'],"Lower Bound"),value=value)')))
                              }
                        }
                        
                  })
      })
      
      #Establish Maindata$selectsample
      observe({
            if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAinteractive' & !is.null(Maindata$dim)) {
                  if (Maindata$dim == 1) {
                        if(!is.null(input$GSCAinteractiveoneupdate)) {
                              isolate({
                                    selectsample <- NULL
                                    for (i in 1:GSCAoneinfo$sampleslidernum) {
                                          if(!is.null(eval(parse(text=paste0("input$GSCAonesampleslider",i,"[1]"))))) 
                                                eval(parse(text=paste0("selectsample <- union(selectsample,input$GSCAonesampleslider",i,"[1]:input$GSCAonesampleslider",i,"[2])")))
                                    }
                                    Maindata$selectsample <- order(Maindata$GSCAscore)[selectsample]
                              })
                        }
                  } else if (Maindata$dim == 2) {
                        Maindata$selectsample <- inpoly$tf
                  } else {
                        if(!is.null(input$GSCAinteractivethreeupdate)) {
                              isolate({
                                    selectsample <- NULL
                                    for (i in 1:GSCAthreeinfo$sampleslidernum) {
                                          if(!is.null(eval(parse(text=paste0("input$GSCAthreesampleslider",i,"[1]"))))) 
                                                eval(parse(text=paste0("selectsample <- union(selectsample,Maindata$GSCAclust$order[input$GSCAthreesampleslider",i,"[1]:input$GSCAthreesampleslider",i,"[2]])")))
                                    }
                                    Maindata$selectsample <- selectsample
                              })
                        }
                  }
                  Maindata$interactsample <- Maindata$selectsample
            } else {
                  Maindata$selectsample <- Maindata$defaultsample
            }
      })
      
      #GSCA Ranking
      observe({
            if (length(Maindata$selectsample)!=0) {
                  AllRanking <- enrichfunc(Maindata$GSCAscore,Maindata$selectsample,Maindata$tab)
                  SelectRanking <- AllRanking[as.numeric(AllRanking[,"Adj.Pvalue"]) <= as.numeric(input$Inputpvalco) & as.numeric(AllRanking[,"FoldChange"]) >= as.numeric(input$Inputfoldchangeco),]
                  if (is.null(SelectRanking) || nrow(SelectRanking) == 0) {
                        Maindata$Ranking <- NULL
                  } else {
                        Maindata$Ranking <- SelectRanking
                  }
            } else {
                  Maindata$Ranking <- NULL
            }
            if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAdefault') {
                  Maindata$defaultranking <- Maindata$Ranking
            } else if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAinteractive'){
                  Maindata$interactranking <- Maindata$Ranking
            }    
      })
      
      #Specify Context
      output$InputGSCAspecifycontextui <- renderUI({
            if (!is.null(Maindata$Ranking))
                  selectInput("GSCAspecifycontext","",Maindata$Ranking[,"SampleType"],multiple=T)
      })
      
      #Establish GSCA Context
      observe({
            if(!is.null(Maindata$Ranking)) {
                  if (input$Inputcontexttype=='Toprank') {
                        if (!is.null(input$InputN)) {
                              Maindata$GSCAcontext <- Maindata$Ranking[1:as.numeric(input$InputN),"SampleType"]
                        } else {
                              Maindata$GSCAcontext <- NULL
                        }
                  } else {
                        Maindata$GSCAcontext <- input$GSCAspecifycontext 
                  }
            } else {
                  Maindata$GSCAcontext <- NULL
            }
            
            if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAdefault') {
                  Maindata$defaultcontext <- Maindata$GSCAcontext
            } else if (input$Mainmethod == 'GSCA' & input$GSCAmethod == 'GSCAinteractive'){
                  Maindata$interactcontext <- Maindata$GSCAcontext
            }
      })
      
      #GSCA sidebar
      output$InputNslider <- renderUI({
            if (!is.null(Maindata$Ranking)) {
                  maxval <- min(30,nrow(Maindata$Ranking))
                  defaultval <- min(5,nrow(Maindata$Ranking))
                  if (maxval > 0)
                        sliderInput("InputN","Number of top ranked context to display",min=0,max=maxval,value=defaultval,step=1)
            }
      })
      
      output$InputGSCAsidebar <- renderUI({
            if (!is.null(Maindata$genedata) && nrow(Maindata$genedata) != 0) {
                  if (Maindata$dim == 1) {
                        p(actionButton("GSCAinteractiveoneupdate","Update Sample Selection"))
                  } else if (Maindata$dim == 2) {
                        tagList(
                              helpText("Draw polygon by clicking on the scatterplot"),
                              p(actionButton('finishpolygon', 'Finish Drawing Polygon')),
                              p(actionButton('addpolygon', 'Add New Polygon')),
                              p(actionButton('undo', 'Undo Last Operation')),
                              p(actionButton('reset', 'Reset'))
                        )
                  } else if (Maindata$dim >= 3) {             
                        p(actionButton("GSCAinteractivethreeupdate","Update Sample Selection"))
                  }
            }
      })
      
      #GSCA default plot
      output$GSCAdefaultplot <- renderUI({
            if (!is.null(Maindata$genedata) && nrow(Maindata$genedata) != 0) {
                  if (Maindata$dim == 1) {
                        tagList(
                              h4(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                              plotOutput("GSCAdefaultplotoneone",height=300),
                              plotOutput("GSCAdefaultplotonetwo",height=300*length(Maindata$GSCAcontext))
                        )
                  } else if (Maindata$dim == 2) {
                        div(align="center",
                            h4(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                            plotOutput("GSCAdefaultplottwo",width=900,height=900),
                            helpText(paste0("Correlation: ",Maindata$twocorr)),
                            helpText(paste0("Pearson Correlation Test p-value: ",Maindata$twocorrp)),
                            helpText(paste0("Regression Slope: ",Maindata$twoslope)),
                            helpText(paste0("Regression Slope t-test p-value: ",Maindata$twoslopep))
                        )
                  } else {
                        tagList(
                              plotOutput("GSCAdefaultplotthree"),
                              h4(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                              plotOutput("GSCAdefaultplotthreeplus")
                        )
                  }
            }      
      })
      
      output$GSCAdefaultplotoneone <- renderPlot({
            if (Maindata$dim == 1) {
                  hist(Maindata$GSCAscore,xlab="Sample Score",xlim=range(Maindata$GSCAscore),main="All samples in the compendium",cex.main=2)
                  abline(v=Maindata$cutoffval[1,1], lty=2)
                  abline(v=Maindata$cutoffval[1,2], lty=2)
                  GSCAstatus$status <- 0
            }
      })
      
      output$GSCAdefaultplotonetwo <- renderPlot({
            if (Maindata$dim == 1) {
                  par(mfrow=c(length(Maindata$GSCAcontext),1),oma=c(0,0,2,0))
                  if (!is.null(Maindata$GSCAcontext))
                        for(INDEX in Maindata$GSCAcontext) {
                              hist(Maindata$GSCAscore[Maindata$tab$SampleType %in% INDEX],xlim=range(Maindata$GSCAscore),main=substr(INDEX,1,25),xlab="Sample Score",cex.main=2)
                              abline(v=Maindata$cutoffval[1,1], lty=2)
                              abline(v=Maindata$cutoffval[1,2], lty=2)
                        }
            }
      })
      
      #Checkbox whether to display enriched context in all area
      output$plotenrichedareaui <- renderUI({
            if (!is.null(Maindata$dim))
                  if (Maindata$dim == 2)
                        checkboxInput("Inputenrichedareaonly","Show enriched context only in POI")
      })
      
      output$GSCAdefaultplottwo <- renderPlot({  
            if (Maindata$dim == 2) {
                  plot(Maindata$GSCAscore[1,],Maindata$GSCAscore[2,],col="#00000022",pch=20,xlab=Maindata$genesetname[1],cex=0.7,ylab=Maindata$genesetname[2],cex.lab=1.5,ylim=c(min(Maindata$GSCAscore[2,]),1.15*max(Maindata$GSCAscore[2,])))
                  toprankingsample <- NULL
                  if (!is.null(Maindata$GSCAcontext)) {
                        for(INDEX in Maindata$GSCAcontext) {
                              toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                        }
                  }
                  points(Maindata$GSCAscore[1,setdiff(Maindata$selectsample,toprankingsample)],Maindata$GSCAscore[2,setdiff(Maindata$selectsample,toprankingsample)],cex=0.7,pch=20)                  
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
                  lines(c(Maindata$cutoffval[1,1],Maindata$cutoffval[1,2]),c(Maindata$cutoffval[2,1],Maindata$cutoffval[2,1]), lty=2,lwd=4,col="blue")
                  lines(c(Maindata$cutoffval[1,1],Maindata$cutoffval[1,2]),c(Maindata$cutoffval[2,2],Maindata$cutoffval[2,2]), lty=2,lwd=4,col="blue")
                  lines(c(Maindata$cutoffval[1,1],Maindata$cutoffval[1,1]),c(Maindata$cutoffval[2,1],Maindata$cutoffval[2,2]), lty=2,lwd=4,col="blue")
                  lines(c(Maindata$cutoffval[1,2],Maindata$cutoffval[1,2]),c(Maindata$cutoffval[2,1],Maindata$cutoffval[2,2]), lty=2,lwd=4,col="blue")
                  GSCAstatus$status <- 0
            }
      })
      
      output$GSCAdefaultplotthree <- renderPlot({
            if (Maindata$dim >= 3 && !all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                  par(oma=c(0.5,0,0.5,max(nchar(Maindata$genesetname))/1.5))
                  colcolorall <- rep("cyan",ncol(Maindata$GSCAscore))
                  colcolorall[Maindata$defaultsample] <- "blue"
                  tmprowv <- F
                  if (input$heatmapthreerowv)
                        tmprowv <- Maindata$GSCArowclust
                  if (!all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                        heatmap.2(Maindata$threesupscore,col=bluered(100),symbreaks=F,Colv=as.dendrogram(Maindata$GSCAclust),dendrogram="none",trace="none",Rowv=tmprowv,labCol=NA,ColSideColors=colcolorall,main="All Samples",useRaster=T)
                  }
                  legend("bottomleft",legend=c("Selected Samples","Not Selected Samples"),lwd=1,col=c("blue","cyan"))
                  GSCAstatus$status <- 0
            }
      })
      
      output$GSCAdefaultplotthreeplus <- renderPlot({
            if (Maindata$dim >= 3 && length(Maindata$defaultsample) > 0 && !all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                  colcolorselect <- rep("white",ncol(Maindata$GSCAscore))
                  par(oma=c(0.5,0,0.5,max(nchar(Maindata$genesetname))/1.5))
                  tmprowv <- F                  
                  if (input$heatmapthreerowv)
                        tmprowv <- Maindata$GSCArowclust
                  if (!is.null(Maindata$GSCAcontext)) {
                        i <- 1
                        for(INDEX in Maindata$GSCAcontext) {
                              colcolorselect[Maindata$tab$SampleType %in% INDEX] <- COLORS[i]
                              i <- i+1
                        }
                        heatmap.2(Maindata$threesupscore[,Maindata$defaultsample],symbreaks=F,col=bluered,labCol=NA,Rowv=tmprowv,dendrogram="none",trace="none",ColSideColors=colcolorselect[Maindata$defaultsample],main="Selected Samples",useRaster=T,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101))
                        leg.txt <- substr(Maindata$GSCAcontext,1,25)
                        legend("bottomleft",legend=leg.txt,lwd=1,col=COLORS)
                  } else {
                        if (length(Maindata$defaultsample) > 1)
                              heatmap.2(Maindata$threesupscore[,Maindata$defaultsample],symbreaks=F,col=bluered,labCol=NA,Rowv=tmprowv,dendrogram="none",trace="none",main="Selected Samples",useRaster=T,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101))
                  }
            }
      })
      
      #GSCA interactive plot
      output$GSCAinteractiveplot <- renderUI({
            if (!is.null(Maindata$genedata) && nrow(Maindata$genedata) != 0) {
                  if (Maindata$dim == 1) {
                        tagList(
                              helpText("Select sample range using sliders. The range would be the union of multiple sliders. After selection, click 'Update Sample Selection' button on the left sidepanel."),
                              helpText('If the slider and the plot are not aligned, please zoom in the webpage. Windows: Hold "Control" and press "-"; Mac: Hold "Command" and press "-"'),
                              actionButton("GSCAonesampleaddslider","Add Slider"),
                              actionButton("GSCAonesampledeleteslider","Delete Slider"),
                              lapply(1:ifelse(is.null(GSCAoneinfo$sampleslidernum),1,GSCAoneinfo$sampleslidernum) , function(i) {
                                    sliderInput(inputId = paste0("GSCAonesampleslider",i),"",min=0,max=ncol(Maindata$GSCAscore),value=c(onesampleslidervalue[i,1],onesampleslidervalue[i,2]),step=1,width='1000px')
                              }),
                              plotOutput("GSCAinteractiveplotoneone",height=300,width='1000px'),
                              plotOutput("GSCAinteractiveplotonetwo",height=300*length(Maindata$GSCAcontext)),
                              h4(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found",""))
                        )
                  } else if (Maindata$dim == 2) {
                        div(align="center",
                            plotOutput("GSCAinteractiveplottwo",clickId="coords",width=900,height=900),
                            plotOutput("GSCAinteractiveplottwoplus",width=900,height=900),
                            h4(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found",""))
                        )      
                  } else {
                        tagList(          
                              helpText("Select sample range using sliders. The range would be the union set of multiple sliders. After selection, click 'Update Sample Selection' button on the left sidepanel."),
                              helpText('If the slider and the heatmap are not aligned, please zoom in the webpage. Windows: Hold "Control" and press "-"; Mac: Hold "Command" and press "-"'),
                              actionButton("GSCAthreesampleaddslider","Add Slider"),
                              actionButton("GSCAthreesampledeleteslider","Delete Slider"),                              
                              lapply(1:ifelse(is.null(GSCAthreeinfo$sampleslidernum),1,GSCAthreeinfo$sampleslidernum) , function(i) {
                                    sliderInput(inputId = paste0("GSCAthreesampleslider",i),"",min=0,max=ncol(Maindata$GSCAscore),value=c(threesampleslidervalue[i,1],threesampleslidervalue[i,2]),step=1,width='1000px')
                              }),                              
                              plotOutput("GSCAinteractiveplotthreecolbar",height=20,width='1000px'),                              
                              plotOutput("GSCAinteractiveplotthreeheatmap",width='1000px'),
                              #plotOutput("GSCARinteractiveplotthreeheatmaprowlab"),
                              plotOutput("GSCAinteractiveplotthreecolbarunder",height=20,width='1000px')
                        )                        
                  }     
            }      
      })
      
      output$GSCAinteractiveplotthreeplotplusui <- renderUI({
            if (!is.null(Maindata$dim) && Maindata$dim > 2 && input$GSCAmethod=='GSCAinteractive')
                  tagList(
                        h4(ifelse(is.null(Maindata$GSCAcontext),"No significantly enriched biological contexts found","")),
                        plotOutput("GSCAinteractiveplotthreeplus"))
      })
      
      #One geneset case
      
      GSCAoneinfo <- reactiveValues()
      
      #initiate number of sample slider, record current slider value
      observe({
            if (is.null(GSCAoneinfo$sampleslidernum))
                  GSCAoneinfo$sampleslidernum <- 1
            for (i in 1:GSCAoneinfo$sampleslidernum) 
                  if(!is.null(eval(parse(text=paste0("input$GSCAonesampleslider",i,"[1]"))))) {
                        onesampleslidervalue[i,1] <<- eval(parse(text=paste0("input$GSCAonesampleslider",i,"[1]"))) 
                        onesampleslidervalue[i,2] <<- eval(parse(text=paste0("input$GSCAonesampleslider",i,"[2]")))
                  }
      })
      
      #set up action for add slider button in one genedata case
      observe({
            if (!is.null(input$GSCAonesampleaddslider) && input$GSCAonesampleaddslider>0)
                  isolate({
                        if (GSCAoneinfo$sampleslidernum < 5) {
                              GSCAoneinfo$sampleslidernum <- GSCAoneinfo$sampleslidernum + 1
                        }
                  })
      })
      
      #set up action for delete slider button in one genedata case
      observe({
            if (!is.null(input$GSCAonesampledeleteslider) && input$GSCAonesampledeleteslider>0)
                  isolate({
                        if (GSCAoneinfo$sampleslidernum != 1) {
                              GSCAoneinfo$sampleslidernum <- GSCAoneinfo$sampleslidernum - 1
                        }
                  })
      })
      
      output$GSCAinteractiveplotoneone <- renderPlot({
            if (Maindata$dim == 1) {
                  selectsample <- NULL
                  for (i in 1:GSCAoneinfo$sampleslidernum) {
                        if(!is.null(eval(parse(text=paste0("input$GSCAonesampleslider",i,"[1]"))))) 
                              eval(parse(text=paste0("selectsample <- union(selectsample,input$GSCAonesampleslider",i,"[1]:input$GSCAonesampleslider",i,"[2])")))
                  }
                  selectsample <- order(Maindata$GSCAscore)[selectsample]
                  par(oma=c(0,0,0,0),mar=c(5.1,0,0,0))
                  colcolorall <- rep("black",ncol(Maindata$GSCAscore))
                  colcolorall[selectsample] <- "cyan"
                  plot(1:ncol(Maindata$GSCAscore),sort(Maindata$GSCAscore),xaxs="i",yaxs="i",cex.main=2,xlab="",ylab="",pch=16,col=colcolorall[order(Maindata$GSCAscore)])
                  legend("topleft",legend=c("Selected samples","Not selected samples"),pch=c(20,20),pt.bg=c("cyan","black"),col=c("cyan","black"),cex=1.2)    
            }
      })
      
      output$GSCAinteractiveplotonetwo <- renderPlot({
            if (Maindata$dim == 1) {
                  if(!is.null(input$GSCAinteractiveoneupdate) && input$GSCAinteractiveoneupdate>0) {
                        Maindata$Ranking
                        Maindata$GSCAcontext
                        isolate({
                              par(mfrow=c(length(Maindata$GSCAcontext),1),oma=c(0,0,2,0))
                              if (!is.null(Maindata$GSCAcontext)) {
                                    for(INDEX in Maindata$GSCAcontext) {
                                          hist(Maindata$GSCAscore[Maindata$tab$SampleType %in% INDEX],xlim=range(Maindata$GSCAscore),main=substr(INDEX,1,25),xlab="Sample Score",cex.main=2)
                                    }     
                              }
                        })
                  }
            }
      })
      
      #Two genesets case
      inpoly <- reactiveValues()
      
      get.coords <- reactive({
            data.frame(x=input$coords$x, y=input$coords$y)
      })
      
      observe({
            if (!is.null(Maindata$GSCAscore)) {
                  inpoly$tf <- rep(F,ncol(Maindata$GSCAscore))
                  polycord <<- NULL
                  polynum <<- 1
                  actiontaken <<- 1
            }
      })
      
      output$GSCAinteractiveplottwo <- renderPlot({  
            if (Maindata$dim == 2) { 
                  inputcoords <- get.coords()
                  input$GSCAinteractiveloadbutton
                  par(cex=1.5,mar=c(5,5,5,2),xpd=T)
                  plot(Maindata$GSCAscore[1,], Maindata$GSCAscore[2,],pch=20, xlab=Maindata$genesetname[1],ylab=Maindata$genesetname[2],cex.lab=1.5)
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
                                    if (is.matrix(tmpcord) && nrow(tmpcord) > 1) {
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
                                          lines(c(tmpcord[i,1],tmpcord[i+1,1]),c(tmpcord[i,2],tmpcord[i+1,2]),type="l",col="blue",lwd=4)
                                    if (j==polynum && finishpolygonaction == 0) {
                                          lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lty=2,lwd=4)
                                    } else {
                                          lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lwd=4)
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
                        plot(Maindata$GSCAscore[1,], Maindata$GSCAscore[2,],ylim=c(min(Maindata$GSCAscore[2,]),1.15*max(Maindata$GSCAscore[2,])),xlab=Maindata$genesetname[1],ylab=Maindata$genesetname[2],pch=20,cex=0.7,col="#00000022",cex.lab=1.5)            
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
                        for (j in 1:polynum) {
                              tmpcord <- polycord[polycord[,3]==j,]
                              if (is.matrix(tmpcord) && nrow(tmpcord) > 2) {
                                    for (i in 1:(nrow(tmpcord)-1)) 
                                          lines(c(tmpcord[i,1],tmpcord[i+1,1]),c(tmpcord[i,2],tmpcord[i+1,2]),type="l",col="blue",lty=2,lwd=4)
                                    lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lty=2,lwd=4)
                              }
                        }
                  }
            }
      })
      
      #More than two genesets case
      
      GSCAthreeinfo <- reactiveValues()
      
      # Suppress color in heatmap
      output$heatmapcolorsuppressui <- renderUI({
            if (!is.null(Maindata$dim))
                  if (Maindata$dim > 2) {
                        tagList(
                              checkboxInput("heatmapcolorsuppresscheck","Suppress heatmap color range"),
                              conditionalPanel("input.heatmapcolorsuppresscheck == 1",
                                               helpText("Note: this function only impacts heatmap display. GSCA results will not change."),
                                               sliderInput("heatmapcolorsuppressslider","Choose color boundaries",min=min(Maindata$GSCAscore),max=max(Maindata$GSCAscore),value=c(min(Maindata$GSCAscore),max(Maindata$GSCAscore)))
                              )
                        )
                  }
      })
      
      observe({
            if (!is.null(input$heatmapcolorsuppresscheck) && input$heatmapcolorsuppresscheck) {
                  supscore <- Maindata$GSCAscore
                  supscore[supscore<input$heatmapcolorsuppressslider[1]]<-input$heatmapcolorsuppressslider[1]
                  supscore[supscore>input$heatmapcolorsuppressslider[2]]<-input$heatmapcolorsuppressslider[2]
                  Maindata$threesupscore <- supscore      
            } else {
                  Maindata$threesupscore <- Maindata$GSCAscore      
            }
            
      })
      
      #Checkbox whether to cluster on rows for heatmaps
      output$heatmapthreerowvui <- renderUI({
            if (!is.null(Maindata$dim))
                  if (Maindata$dim > 2)
                        checkboxInput("heatmapthreerowv","Cluster on rows")
      })
      
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
      
      #       output$GSCARinteractiveplotthreeheatmaprowlab <- renderText({
      #             par(mar = c(0,0,0,0))
      #             plot(c(1, 1), c(Maindata$dim, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      #             pos <- seq(1.5-1/2/Maindata$dim,Maindata$dim-0.5+1/2/Maindata$dim,length.out=Maindata$dim)
      #             if (input$heatmapthreerowv) {
      #                   Rowv <- rowMeans(Maindata$GSCAscore)
      #                   hcr <- hclust(dist(Maindata$GSCAscore))
      #                   ddr <- as.dendrogram(hcr)
      #                   ddr <- reorder(ddr, Rowv)
      #                   rowInd <- order.dendrogram(ddr)
      #                   names <- rev(row.names(Maindata$GSCAscore)[rowInd])
      #             } else {
      #                   names <- row.names(Maindata$GSCAscore)
      #             }
      #             for (i in 1:Maindata$dim) {
      #                   text(1,pos[i],names[i],cex=1.3)
      #             }
      #       })
      
      output$GSCAinteractiveplotthreecolbar <- renderPlot({
            if (Maindata$dim >= 3) { 
                  selectsample <- NULL
                  for (i in 1:GSCAthreeinfo$sampleslidernum) {
                        if(!is.null(eval(parse(text=paste0("input$GSCAthreesampleslider",i,"[1]"))))) 
                              eval(parse(text=paste0("selectsample <- union(selectsample,Maindata$GSCAclust$order[input$GSCAthreesampleslider",i,"[1]:input$GSCAthreesampleslider",i,"[2]])")))
                  }
                  colcolorall <- rep("cyan",ncol(Maindata$GSCAscore))
                  colcolorall[selectsample] <- "blue"
                  par(mar=c(0,0,0,0))
                  image(cbind(1:ncol(Maindata$GSCAscore)),col=colcolorall[Maindata$GSCAclust$order],axes=F,useRaster=T)
            }
      })
      
      output$GSCAinteractiveplotthreecolbarunder <- renderPlot({
            if (Maindata$dim >= 3) {
                  selectsample <- NULL
                  for (i in 1:GSCAthreeinfo$sampleslidernum) {
                        if(!is.null(eval(parse(text=paste0("input$GSCAthreesampleslider",i,"[1]"))))) 
                              eval(parse(text=paste0("selectsample <- union(selectsample,Maindata$GSCAclust$order[input$GSCAthreesampleslider",i,"[1]:input$GSCAthreesampleslider",i,"[2]])")))
                  }
                  colcolorall <- rep("cyan",ncol(Maindata$GSCAscore))
                  colcolorall[selectsample] <- "blue"
                  par(mar=c(0,0,0,0))
                  image(cbind(1:ncol(Maindata$GSCAscore)),col=colcolorall[Maindata$GSCAclust$order],axes=F,useRaster=T)
            }
      })
      
      output$GSCAinteractiveplotthreeheatmap <- renderPlot({
            if (Maindata$dim >= 3) {
                  par(mar=c(0,0,0,0))
                  if (input$heatmapthreerowv) {
                        rowInd <- order.dendrogram(Maindata$GSCArowclust)
                        image(1:ncol(Maindata$GSCAscore),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order]),col=bluered(100),axes=F,useRaster=T)
                  } else {
                        image(1:ncol(Maindata$GSCAscore),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order]),col=bluered(100),axes=F,useRaster=T)
                  }
                  pos <- 1:Maindata$dim
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
                        text(0.5*ncol(Maindata$GSCAscore),pos[i],names[i],cex=1.5)
                  }
            }
      })
      
      
      output$GSCAinteractiveplotthreezoominallpartsui <- renderUI({
            if (!is.null(Maindata$genedata) && nrow(Maindata$genedata) != 0) {
                  if (Maindata$dim >= 3 && input$GSCAmethod=='GSCAinteractive') {
                        tagList(
                              checkboxInput("GSCAinteractiveplotthreeheatmapzoomincheck","Heatmap zoom in"),
                              conditionalPanel("input.GSCAinteractiveplotthreeheatmapzoomincheck==1",
                                               helpText("Zoom in selected part of the heatmap"),
                                               checkboxInput("GSCAinteractiveplotthreeheatmapzoominrealtime","Real time zoom in",value = T),
                                               conditionalPanel("input.GSCAinteractiveplotthreeheatmapzoominrealtime==0",
                                                                p(actionButton("GSCAinteractiveplotthreeheatmapzoominupdate","Update Sample Selection"))),
                                               uiOutput("GSCAinteractiveplotthreeheatmapzoominplotui")
                              )
                        )
                  }
            }
      })
      
      output$GSCAinteractiveplotthreeheatmapzoominplotui <- renderUI({
            if (GSCAthreeinfo$sampleslidernum == 1) {
                  tags$div(class="row-fluid",
                           tags$div(class="span11",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot1")),
                           tags$div(class="span1",plotOutput("GSCAinteractiveplotthreeheatmapzoominplotlab"))
                  )
            } else if (GSCAthreeinfo$sampleslidernum == 2) {
                  tags$div(class="row-fluid",
                           tags$div(class="span5",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot1")),
                           tags$div(class="span5",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot2")),
                           tags$div(class="span1",plotOutput("GSCAinteractiveplotthreeheatmapzoominplotlab"))
                  )
            } else if (GSCAthreeinfo$sampleslidernum == 3) {
                  tags$div(class="row-fluid",
                           tags$div(class="span3",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot1")),
                           tags$div(class="span3",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot2")),
                           tags$div(class="span3",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot3")),
                           tags$div(class="span1",plotOutput("GSCAinteractiveplotthreeheatmapzoominplotlab"))
                  )
            } else if (GSCAthreeinfo$sampleslidernum == 4) {
                  tags$div(class="row-fluid",
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot1")),
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot2")),
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot3")),
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot4")),
                           tags$div(class="span1",plotOutput("GSCAinteractiveplotthreeheatmapzoominplotlab"))
                  )
            } else if (GSCAthreeinfo$sampleslidernum == 5) {
                  tags$div(class="row-fluid",
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot1")),
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot2")),
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot3")),
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot4")),
                           tags$div(class="span2",plotOutput("GSCAinteractiveplotthreeheatmapzoominplot5")),
                           tags$div(class="span1",plotOutput("GSCAinteractiveplotthreeheatmapzoominplotlab"))
                  )
            }
            
      })
      
      #zoom in heatmaps
      output$GSCAinteractiveplotthreeheatmapzoominplot1 <- renderPlot({
            if (input$GSCAinteractiveplotthreeheatmapzoomincheck && !all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                  par(mar=c(0,0,0,0))
                  if (input$GSCAinteractiveplotthreeheatmapzoominrealtime) {
                        if (input$heatmapthreerowv) {
                              rowInd <- order.dendrogram(Maindata$GSCArowclust)
                              image(1:(input$GSCAthreesampleslider1[2]-input$GSCAthreesampleslider1[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider1[1]:input$GSCAthreesampleslider1[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        } else {
                              image(1:(input$GSCAthreesampleslider1[2]-input$GSCAthreesampleslider1[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider1[1]:input$GSCAthreesampleslider1[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        }
                  } else {
                        if (input$GSCAinteractiveplotthreeheatmapzoominupdate > 0) {
                              if (input$heatmapthreerowv) {
                                    rowInd <- order.dendrogram(Maindata$GSCArowclust)
                                    isolate(image(1:(input$GSCAthreesampleslider1[2]-input$GSCAthreesampleslider1[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider1[1]:input$GSCAthreesampleslider1[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              } else {
                                    isolate(image(1:(input$GSCAthreesampleslider1[2]-input$GSCAthreesampleslider1[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider1[1]:input$GSCAthreesampleslider1[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              }
                        }
                  }
            }
      })
      
      output$GSCAinteractiveplotthreeheatmapzoominplot2 <- renderPlot({
            if (input$GSCAinteractiveplotthreeheatmapzoomincheck & GSCAthreeinfo$sampleslidernum >= 2 && !all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                  par(mar=c(0,0,0,0))
                  if (input$GSCAinteractiveplotthreeheatmapzoominrealtime) {
                        if (input$heatmapthreerowv) {
                              rowInd <- order.dendrogram(Maindata$GSCArowclust)
                              image(1:(input$GSCAthreesampleslider2[2]-input$GSCAthreesampleslider2[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider2[1]:input$GSCAthreesampleslider2[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        } else {
                              image(1:(input$GSCAthreesampleslider2[2]-input$GSCAthreesampleslider2[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider2[1]:input$GSCAthreesampleslider2[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        }
                  } else {
                        if (input$GSCAinteractiveplotthreeheatmapzoominupdate > 0) {
                              if (input$heatmapthreerowv) {
                                    rowInd <- order.dendrogram(Maindata$GSCArowclust)
                                    isolate(image(1:(input$GSCAthreesampleslider2[2]-input$GSCAthreesampleslider2[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider2[1]:input$GSCAthreesampleslider2[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              } else {
                                    isolate(image(1:(input$GSCAthreesampleslider2[2]-input$GSCAthreesampleslider2[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider2[1]:input$GSCAthreesampleslider2[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              }
                        }
                  }
            }
      })
      
      output$GSCAinteractiveplotthreeheatmapzoominplot3 <- renderPlot({
            if (input$GSCAinteractiveplotthreeheatmapzoomincheck & GSCAthreeinfo$sampleslidernum >= 3 && !all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                  par(mar=c(0,0,0,0))
                  if (input$GSCAinteractiveplotthreeheatmapzoominrealtime) {
                        if (input$heatmapthreerowv) {
                              rowInd <- order.dendrogram(Maindata$GSCArowclust)
                              image(1:(input$GSCAthreesampleslider3[2]-input$GSCAthreesampleslider3[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider3[1]:input$GSCAthreesampleslider3[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        } else {
                              image(1:(input$GSCAthreesampleslider3[2]-input$GSCAthreesampleslider3[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider3[1]:input$GSCAthreesampleslider3[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        }
                  } else {
                        if (input$GSCAinteractiveplotthreeheatmapzoominupdate > 0) {
                              if (input$heatmapthreerowv) {
                                    rowInd <- order.dendrogram(Maindata$GSCArowclust)
                                    isolate(image(1:(input$GSCAthreesampleslider3[2]-input$GSCAthreesampleslider3[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider3[1]:input$GSCAthreesampleslider3[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              } else {
                                    isolate(image(1:(input$GSCAthreesampleslider3[2]-input$GSCAthreesampleslider3[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider3[1]:input$GSCAthreesampleslider3[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              }
                        }
                  }
            }
      })
      
      output$GSCAinteractiveplotthreeheatmapzoominplot4 <- renderPlot({
            if (input$GSCAinteractiveplotthreeheatmapzoomincheck & GSCAthreeinfo$sampleslidernum >= 4 && !all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                  par(mar=c(0,0,0,0))
                  if (input$GSCAinteractiveplotthreeheatmapzoominrealtime) {
                        if (input$heatmapthreerowv) {
                              rowInd <- order.dendrogram(Maindata$GSCArowclust)
                              image(1:(input$GSCAthreesampleslider4[2]-input$GSCAthreesampleslider4[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider4[1]:input$GSCAthreesampleslider4[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        } else {
                              image(1:(input$GSCAthreesampleslider4[2]-input$GSCAthreesampleslider4[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider4[1]:input$GSCAthreesampleslider4[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        }
                  } else {
                        if (input$GSCAinteractiveplotthreeheatmapzoominupdate > 0) {
                              if (input$heatmapthreerowv) {
                                    rowInd <- order.dendrogram(Maindata$GSCArowclust)
                                    isolate(image(1:(input$GSCAthreesampleslider4[2]-input$GSCAthreesampleslider4[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider4[1]:input$GSCAthreesampleslider4[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              } else {
                                    isolate(image(1:(input$GSCAthreesampleslider4[2]-input$GSCAthreesampleslider4[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider4[1]:input$GSCAthreesampleslider4[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              }
                        }
                  }
            }
      })
      
      output$GSCAinteractiveplotthreeheatmapzoominplot5 <- renderPlot({
            if (input$GSCAinteractiveplotthreeheatmapzoomincheck & GSCAthreeinfo$sampleslidernum >= 5 && !all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                  par(mar=c(0,0,0,0))
                  if (input$GSCAinteractiveplotthreeheatmapzoominrealtime) {
                        if (input$heatmapthreerowv) {
                              rowInd <- order.dendrogram(Maindata$GSCArowclust)
                              image(1:(input$GSCAthreesampleslider5[2]-input$GSCAthreesampleslider5[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider5[1]:input$GSCAthreesampleslider5[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        } else {
                              image(1:(input$GSCAthreesampleslider5[2]-input$GSCAthreesampleslider5[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider5[1]:input$GSCAthreesampleslider5[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T)
                        }
                  } else {
                        if (input$GSCAinteractiveplotthreeheatmapzoominupdate > 0) {
                              if (input$heatmapthreerowv) {
                                    rowInd <- order.dendrogram(Maindata$GSCArowclust)
                                    isolate(image(1:(input$GSCAthreesampleslider5[2]-input$GSCAthreesampleslider5[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[rowInd,Maindata$GSCAclust$order[input$GSCAthreesampleslider5[1]:input$GSCAthreesampleslider5[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              } else {
                                    isolate(image(1:(input$GSCAthreesampleslider5[2]-input$GSCAthreesampleslider5[1]+1),1:nrow(Maindata$GSCAscore),t(Maindata$threesupscore[nrow(Maindata$GSCAscore):1,Maindata$GSCAclust$order[input$GSCAthreesampleslider5[1]:input$GSCAthreesampleslider5[2]]]),col=bluered(100),axes=F,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101),useRaster=T))
                              }
                        }
                  }
            }
      })
      
      #       output$GSCAinteractiveplotthreeheatmapzoominplotlab <- renderPlot({
      #             par(mar = c(0,0,0,0))
      #             plot(c(1, 1), c(Maindata$dim, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      #             pos <- seq(1.5-1/2/Maindata$dim,Maindata$dim-0.5+1/2/Maindata$dim,length.out=Maindata$dim)
      #             if (input$heatmapthreerowv) {
      #                   Rowv <- rowMeans(Maindata$GSCAscore)
      #                   hcr <- hclust(dist(Maindata$GSCAscore))
      #                   ddr <- as.dendrogram(hcr)
      #                   ddr <- reorder(ddr, Rowv)
      #                   rowInd <- order.dendrogram(ddr)
      #                   names <- row.names(Maindata$GSCAscore)[rowInd]
      #             } else {
      #                   names <- rev(row.names(Maindata$GSCAscore))
      #             }
      #             for (i in 1:Maindata$dim) {
      #                   text(1,pos[i],names[i],cex=1.2)
      #             }
      #       })
      
      output$GSCAinteractiveplotthreeplus <- renderPlot({
            if (Maindata$dim >= 3 && length(Maindata$selectsample) > 0 && !all(Maindata$threesupscore==Maindata$threesupscore[1,1])) {
                  if(!is.null(input$GSCAinteractivethreeupdate) && input$GSCAinteractivethreeupdate > 0) {
                        Maindata$Ranking
                        Maindata$GSCAcontext
                        input$heatmapthreerowv
                        isolate({
                              par(oma=c(0.5,0,0.5,max(nchar(Maindata$genesetname))/1.5))
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
                                    heatmap.2(Maindata$threesupscore[,Maindata$selectsample],symbreaks=F,col=bluered,labCol=NA,Rowv=tmprowv,dendrogram="none",trace="none",ColSideColors=colcolorselect[Maindata$selectsample],main="Selected Samples",useRaster=T,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101))
                                    legend("bottomleft",legend=leg.txt,lwd=1,col=COLORS)
                              } else {
                                    if (length(Maindata$selectsample) > 1)
                                          heatmap.2(Maindata$threesupscore[,Maindata$selectsample],symbreaks=F,col=bluered,labCol=NA,Rowv=tmprowv,dendrogram="none",trace="none",main="Selected Samples",useRaster=T,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101))
                              }
                        })
                  }
            }
      })
      
      output$GSCArankingtable <- renderDataTable(Maindata$Ranking)
      
      #####3D scatterplot
      #       output$RGLplot <- renderWebGL({
      #             if (!is.null(Maindata$dim) && Maindata$dim == 3) {
      #                   dotcolor <- rep("gray",ncol(Maindata$GSCAscore))
      #                   dotcolor[Maindata$selectsample] <- "black"
      #                   i <- 1
      #                   for(INDEX in Maindata$GSCAcontext) {
      #                         dotcolor[Maindata$tab$SampleType %in% INDEX] <- COLORS[i]
      #                         i <- i+1
      #                   }
      #                   points3d(Maindata$GSCAscore[1,],Maindata$GSCAscore[2,],Maindata$GSCAscore[3,],col=dotcolor)
      #                   axes3d()      
      #                   mtext3d(row.names(Maindata$GSCAscore)[1],edge="x",size=2)
      #                   mtext3d(row.names(Maindata$GSCAscore)[2],edge="y",size=2)
      #                   mtext3d(row.names(Maindata$GSCAscore)[3],edge="z",size=2)
      #             }      
      #       })
      
      ##### save and load POI
      
      output$GSCAinteractivesavebutton <- downloadHandler(
            filename = function() { "POIfile.txt" },
            content = function(file) {
                  if (input$GSCAmethod=='GSCAdefault') {
                        write.table(Maindata$cutoffval,file,row.names=F,col.names=F)
                  } else {                  
                        if (Maindata$dim == 1) {
                              tmp <- NULL
                              for (i in 1:GSCAoneinfo$sampleslidernum) {
                                    if(!is.null(eval(parse(text=paste0("input$GSCAonesampleslider",i,"[1]"))))) 
                                          eval(parse(text=paste0("tmp <- rbind(tmp,c(input$GSCAonesampleslider",i,"[1],input$GSCAonesampleslider",i,"[2]))")))
                              }
                              write.table(tmp,file,row.names=F,col.names=F)
                        } else if (Maindata$dim == 2) {
                              write.table(polycord,file,row.names=F,col.names=F)
                        } else {
                              tmp <- NULL
                              for (i in 1:GSCAthreeinfo$sampleslidernum) {
                                    if(!is.null(eval(parse(text=paste0("input$GSCAthreesampleslider",i,"[1]"))))) 
                                          eval(parse(text=paste0("tmp <- rbind(tmp,c(input$GSCAthreesampleslider",i,"[1],input$GSCAthreesampleslider",i,"[2]))")))
                              }
                              write.table(tmp,file,row.names=F,col.names=F)
                        }
                  }
                  
            }
      )
      
      observe({            
            if (input$GSCAinteractiveloadbutton > 0)
                  isolate({
                        if (!is.null(input$GSCAinteractiveload)) { 
                              tmp <- read.table(input$GSCAinteractiveload$datapath)
                              if (input$GSCAmethod=='GSCAdefault') {
                                    if (input$numericpoimethod == "slider") {
                                          for (i in 1:Maindata$dim) {
                                                eval(parse(text=paste0("updateSliderInput(session,'GSCAnumericpoislider",i,"',value=c(tmp[",i,",1],tmp[",i,",2]))")))
                                          }
                                    } else {
                                          for (i in 1:Maindata$dim) {
                                                eval(parse(text=paste0('updateTextInput(session,"GSCAnumericpoitextlower', i,'",value=tmp[',i,',1])')))
                                                eval(parse(text=paste0('updateTextInput(session,"GSCAnumericpoitextupper', i,'",value=tmp[',i,',2])')))
                                          }     
                                    }
                              } else {
                                    if (Maindata$dim == 1) {
                                          GSCAoneinfo$sampleslidernum <- nrow(tmp)
                                          for (i in 1:GSCAoneinfo$sampleslidernum) {
                                                eval(parse(text=paste0("updateSliderInput(session,'GSCAonesampleslider",i,"',value=c(tmp[",i,",1],tmp[",i,",2]))")))
                                          } 
                                    } else if (Maindata$dim == 2) {
                                          polycord <<- as.matrix(tmp)
                                          polynum <<- max(polycord[,3])
                                          actiontaken <<- 1
                                    } else {
                                          GSCAthreeinfo$sampleslidernum <- nrow(tmp)
                                          for (i in 1:GSCAthreeinfo$sampleslidernum) {
                                                eval(parse(text=paste0("updateSliderInput(session,'GSCAthreesampleslider",i,"',value=c(tmp[",i,",1],tmp[",i,",2]))")))
                                          }                                  
                                    }
                              }
                        }
                  })            
      })
      
      #####   Mainmethod : Download   #####
      
      output$Downloadsidebarui <- renderUI({
            if (!is.null(Maindata$genedata)) {
                  if (Maindata$dim == 1) {
                        tagList(
                              wellPanel(
                                    selectInput("Downloadplotonetype","Select File Type",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotonefilename","File Name","GSCA Plot"),
                                    textInput("Downloadplotonefilewidth","Plot Width (inches)",10),
                                    textInput("Downloadplotonefileheight","Plot Height (inches)",3*(as.numeric(input$InputN)+1)),
                                    p(downloadButton("Downloadplotone","Save Plot"))
                              ),
                              wellPanel(
                                    helpText("Change Plotting Details"),
                                    textInput("Downloadmaintitleone","Enter Main Title",""),
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
                                    selectInput("Downloadplottwotype","Select File Type",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplottwofilename","File Name","GSCA Plot"),
                                    textInput("Downloadplottwofilewidth","Plot Width (inches)",15),
                                    textInput("Downloadplottwofileheight","Plot Height (inches)",15),
                                    p(downloadButton("Downloadplottwo","Save Plot"))
                              ),
                              wellPanel(
                                    helpText("Change Plotting Details"),
                                    textInput("Downloadmaintitletwo","Enter Main Title",""),
                                    textInput("Downloadpointsizetwo","Point size",1),
                                    textInput("Downloadcextwo","Font size",2),
                                    textInput("Downloadxlabtwo","Title for X Axis",Maindata$genesetname[1]),
                                    textInput("Downloadylabtwo","Title for Y Axis",Maindata$genesetname[2]),
                                    textInput("Downloadxlimmintwo","Minimum Value of X Axis",min(Maindata$GSCAscore[1,])),
                                    textInput("Downloadxlimmaxtwo","Maximum Value of X Axis",max(Maindata$GSCAscore[1,])),
                                    textInput("Downloadylimmintwo","Minimum Value of Y Axis",min(Maindata$GSCAscore[2,])),
                                    textInput("Downloadylimmaxtwo","Maximum Value of Y Axis",max(Maindata$GSCAscore[2,])),
                                    textInput("Downloadlegcextwo","Legend Size",0.5),
                                    checkboxInput("Downloadlegpostftwo","Change Legend Position",value=F),
                                    conditionalPanel(condition="input.Downloadlegpostftwo=='1'",
                                                     textInput("Downloadlegposxtwo","x-axis position",0),
                                                     textInput("Downloadlegposytwo","y-axis position",0)
                                    ),
                                    checkboxInput("Downloadcorvaluetwo","Show Correlation and p-value")
                              )
                        )
                  } else {
                        tagList(
                              wellPanel(
                                    selectInput("Downloadplotthreeplotonetype","File Type for Heatmap One",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotthreeplotonefilename","File Name","GSCA Heatmap One"),
                                    textInput("Downloadplotthreeplotonefilewidth","Plot Width (inches)",20),
                                    textInput("Downloadplotthreeplotonefileheight","Plot Height (inches)",8),
                                    p(downloadButton("Downloadplotthreeplotone","Save Heatmap One")),
                                    selectInput("Downloadplotthreeplottwotype","File Type for Heatmap Two",choices=c("pdf","png","ps","jpeg","bmp","tiff")),
                                    textInput("Downloadplotthreeplottwofilename","File Name","GSCA Heatmap Two"),
                                    textInput("Downloadplotthreeplottwofilewidth","Plot Width (inches)",20),
                                    textInput("Downloadplotthreeplottwofileheight","Plot Height (inches)",8),
                                    p(downloadButton("Downloadplotthreeplottwo","Save Heatmap Two"))
                              ),
                              wellPanel(
                                    helpText("Change Plotting Details"),
                                    textInput("Downloadmaintitlethreeplotone","Main Title For Heatmap One","All Samples"),
                                    textInput("Downloadmaintitlethreeplottwo","Main Title For Heatmap Two","Selected Samples"),
                                    selectInput("Downloadthreecolplotone","Palette for Heatmap One",choices=c("Bluered","Heat","Rainbow","Terrain","Topo","CM","gray")),
                                    selectInput("Downloadthreecolplottwo","Palette for Heatmap Two",choices=c("Bluered","Heat","Rainbow","Terrain","Topo","CM","gray")),
                                    checkboxInput("Downloadthreecoldendplotone","Display Column Histogram in Heatmap One (Could take some time)"),
                                    checkboxInput("Downloadthreecoldendplottwo","Display Column Histogram in Heatmap Two"),
                                    checkboxInput("Downloadthreerowvplotone","Cluster on Rows in Heatmap One"),
                                    uiOutput("Downloadthreerowdendplotoneui"),
                                    checkboxInput("Downloadthreerowvplottwo","Cluster on Rows in Heatmap Two"),
                                    uiOutput("Downloadthreerowdendplottwoui"),
                                    checkboxInput("Downloadthreerotatelabplotone","Rotate Row Label in Heatmap One"),
                                    checkboxInput("Downloadthreerotatelabplottwo","Rotate Row Label in Heatmap Two"),
                                    textInput("Downloadlegcexthreeone","Legend Size in Heatmap One",1),
                                    textInput("Downloadlegcexthreetwo","Legend Size in Heatmap Two",1)
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
                  if (input$Downloadregionselect == 'Numeric') {
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
                        plotOutput("downloadtworenderplot",width=900,height=900)
                  } else {
                        tagList(
                              h4("Heatmap One"),
                              plotOutput("downloadthreeonerenderplot"),
                              h4("Heatmap Two"),
                              plotOutput("downloadthreetworenderplot")
                        )
                  }
      })
      
      output$downloadonerenderplot <- renderPlot({
            downloadonefunc()
            GSCAstatus$status <- 0
      })
      
      output$downloadtworenderplot <- renderPlot({
            downloadtwofunc()
            GSCAstatus$status <- 0
      })
      
      output$downloadthreeonerenderplot <- renderPlot({
            downloadthreeonefunc()
            GSCAstatus$status <- 0
      })
      
      output$downloadthreetworenderplot <- renderPlot(downloadthreetwofunc())
      
      downloadonefunc <- function() {
            colone <- input$Downloadcolone
            if (colone == "NULL")
                  colone <- NULL
            par(mfrow=c(length(Maindata$downloadcontext)+1,1),oma=c(0,0,2,0),cex.main=2,cex.lab=2,cex.axis=2,mar=c(5,5,4,2))
            if (input$Downloadregionselect == 'Numeric') {
                  hist(Maindata$GSCAscore,xlab=input$Downloadxlabone,ylab=input$Downloadylabone,xlim=as.numeric(c(input$Downloadxlimminone,input$Downloadxlimmaxone)),col=colone,main="All Biological contexts")
                  abline(v=Maindata$cutoffval[1,1], lty=2)
                  abline(v=Maindata$cutoffval[1,2], lty=2)     
            } else {
                  selectsample <- NULL
                  for (i in 1:GSCAoneinfo$sampleslidernum) {
                        if(!is.null(eval(parse(text=paste0("input$GSCAonesampleslider",i,"[1]"))))) 
                              eval(parse(text=paste0("selectsample <- union(selectsample,input$GSCAonesampleslider",i,"[1]:input$GSCAonesampleslider",i,"[2])")))
                  }
                  selectsample <- order(Maindata$GSCAscore)[selectsample]
                  colcolorall <- rep("black",ncol(Maindata$GSCAscore))
                  colcolorall[selectsample] <- "cyan"
                  plot(1:ncol(Maindata$GSCAscore),sort(Maindata$GSCAscore),xaxs="i",yaxs="i",cex.main=2,xlab="",ylab="Gene expression value",pch=16,col=colcolorall[order(Maindata$GSCAscore)])
                  legend("topleft",legend=c("Selected Samples","Not Selected Samples"),pch=c(20,20),pt.bg=c("cyan","black"),col=c("cyan","black"),cex=0.8)
            }
            if (!is.null(Maindata$downloadcontext)) {
                  for(INDEX in Maindata$downloadcontext) {
                        hist(Maindata$GSCAscore[Maindata$tab$SampleType %in% INDEX],xlab=input$Downloadxlabone,ylab=input$Downloadylabone,xlim=as.numeric(c(input$Downloadxlimminone,input$Downloadxlimmaxone)),col=colone,main=substr(INDEX,1,25))
                        if (input$Downloadregionselect == 'Numeric') {
                              abline(v=Maindata$cutoffval[1,1], lty=2)      
                              abline(v=Maindata$cutoffval[1,2], lty=2)      
                        }
                  }     
            }
            title(input$Downloadmaintitleone,outer=T)
      } 
      
      downloadtwofunc <- function() {
            par(cex=input$Downloadcextwo,mar=c(6,6,6,2),xpd=T)
            cortext <- paste0("Correlation: ",round(Maindata$twocorr,3),"; ","Correlation p-value: ",round(Maindata$twocorrp,3),"; ","Slope: ",round(Maindata$twoslope,3),"; ","Slope p-value: ",round(Maindata$twoslopep,3))
            if (input$Downloadregionselect == 'Numeric') {
                  plot(Maindata$GSCAscore[1,],Maindata$GSCAscore[2,],cex=as.numeric(input$Downloadpointsizetwo),col="#00000022",pch=20,xlab=input$Downloadxlabtwo,ylab=input$Downloadylabtwo,xlim = as.numeric(c(input$Downloadxlimmintwo,input$Downloadxlimmaxtwo)), ylim = as.numeric(c(input$Downloadylimmintwo,input$Downloadylimmaxtwo)),main=input$Downloadmaintitletwo)
                  toprankingsample <- NULL
                  if (!is.null(Maindata$downloadcontext)) {
                        for(INDEX in Maindata$downloadcontext) {
                              toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                        }
                  }
                  points(Maindata$GSCAscore[1,setdiff(Maindata$downloadsample,toprankingsample)],Maindata$GSCAscore[2,setdiff(Maindata$downloadsample,toprankingsample)],cex=as.numeric(input$Downloadpointsizetwo),pch=20)                  
                  if (!is.null(Maindata$downloadcontext)) {
                        i <- 1
                        for(INDEX in Maindata$downloadcontext) {
                              if (input$Inputenrichedareaonly == TRUE) {
                                    points(Maindata$GSCAscore[1,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$downloadsample)],Maindata$GSCAscore[2,intersect(which(Maindata$tab$SampleType %in% INDEX), Maindata$downloadsample)],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i],cex=as.numeric(input$Downloadpointsizetwo))
                              } else {
                                    points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX],Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX],
                                           col=COLORS[i],pch=STYLES[i],bg=COLORS[i],cex=as.numeric(input$Downloadpointsizetwo))
                              }
                              i <- i+1
                        }     
                  }
                  leg.txt <- c(substr(Maindata$downloadcontext,1,25),"Selected Samples","Not Selected Samples")
                  if (input$Downloadlegpostftwo) {
                        legtmpx <- as.numeric(input$Downloadlegposxtwo)
                        legtmpy <- as.numeric(input$Downloadlegposytwo)
                  } else {
                        legtmpx <- "topleft"
                        legtmpy <- NULL
                  }
                  if (length(leg.txt)==2) {
                        legend(legtmpx,legtmpy,legend=leg.txt,pch=c(20,20),pt.bg=c("black","#00000022"),col=c("black","#00000022"),cex=as.numeric(input$Downloadlegcextwo))                        
                  } else {
                        legend(legtmpx,legtmpy,legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=as.numeric(input$Downloadlegcextwo))   
                  }
                  if (input$Downloadcorvaluetwo)
                        mtext(cortext)
                  lines(c(Maindata$cutoffval[1,1],Maindata$cutoffval[1,2]),c(Maindata$cutoffval[2,1],Maindata$cutoffval[2,1]), lty=2,lwd=4,col="blue")
                  lines(c(Maindata$cutoffval[1,1],Maindata$cutoffval[1,2]),c(Maindata$cutoffval[2,2],Maindata$cutoffval[2,2]), lty=2,lwd=4,col="blue")
                  lines(c(Maindata$cutoffval[1,1],Maindata$cutoffval[1,1]),c(Maindata$cutoffval[2,1],Maindata$cutoffval[2,2]), lty=2,lwd=4,col="blue")
                  lines(c(Maindata$cutoffval[1,2],Maindata$cutoffval[1,2]),c(Maindata$cutoffval[2,1],Maindata$cutoffval[2,2]), lty=2,lwd=4,col="blue")
            } else {
                  if (sum(inpoly$tf) != 0) {
                        plot(Maindata$GSCAscore[1,], Maindata$GSCAscore[2,], xlim = as.numeric(c(input$Downloadxlimmintwo,input$Downloadxlimmaxtwo)), ylim = as.numeric(c(input$Downloadylimmintwo,input$Downloadylimmaxtwo)),xlab=input$Downloadxlabtwo,ylab=input$Downloadylabtwo,col="#00000022",pch=20,cex=as.numeric(input$Downloadpointsizetwo),main=input$Downloadmaintitletwo)            
                        toprankingsample <- NULL
                        if (!is.null(Maindata$downloadcontext)) {
                              for(INDEX in Maindata$downloadcontext) {
                                    toprankingsample <- union(toprankingsample,which(Maindata$tab$SampleType %in% INDEX))
                              }
                        }
                        points(Maindata$GSCAscore[1,setdiff(which(inpoly$tf),toprankingsample)], Maindata$GSCAscore[2,setdiff(which(inpoly$tf),toprankingsample)],cex=as.numeric(input$Downloadpointsizetwo),pch=20)
                        if (!is.null(Maindata$downloadcontext)) {
                              i <- 1
                              for (INDEX in Maindata$downloadcontext) {    
                                    if (input$Inputenrichedareaonly == TRUE) {
                                          points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0], Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX&inpoly$tf>0],col = COLORS[i], pch = STYLES[i], bg = COLORS[i],cex=as.numeric(input$Downloadpointsizetwo))
                                    } else {
                                          points(Maindata$GSCAscore[1,Maindata$tab$SampleType %in% INDEX], Maindata$GSCAscore[2,Maindata$tab$SampleType %in% INDEX],col = COLORS[i], pch = STYLES[i], bg = COLORS[i],cex=as.numeric(input$Downloadpointsizetwo))
                                    }
                                    i <- i+1
                              }     
                        }
                        leg.txt <- c(substr(Maindata$downloadcontext,1,25),"Selected Samples","Not Selected Samples")
                        if (input$Downloadlegpostftwo) {
                              legtmpx <- as.numeric(input$Downloadlegposxtwo)
                              legtmpy <- as.numeric(input$Downloadlegposytwo)
                        } else {
                              legtmpx <- "topleft"
                              legtmpy <- NULL
                        }
                        if (length(leg.txt)==2) {
                              legend(legtmpx,legtmpy,legend=leg.txt,pch=c(20,20),pt.bg=c("black","#00000022"),col=c("black","#00000022"),cex=as.numeric(input$Downloadlegcextwo))                        
                        } else {
                              legend(legtmpx,legtmpy,legend=leg.txt,pch=c(STYLES[1:(length(leg.txt)-2)],20,20),pt.bg=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),col=c(COLORS[1:(length(leg.txt)-2)],"black","#00000022"),cex=as.numeric(input$Downloadlegcextwo))     
                        }
                        if (input$Downloadcorvaluetwo)
                              mtext(cortext)
                        for (j in 1:polynum) {
                              tmpcord <- polycord[polycord[,3]==j,]
                              if (is.matrix(tmpcord)) {
                                    for (i in 1:(nrow(tmpcord)-1)) 
                                          lines(c(tmpcord[i,1],tmpcord[i+1,1]),c(tmpcord[i,2],tmpcord[i+1,2]),type="l",col="blue",lty=2,lwd=4)
                                    lines(c(tmpcord[1,1],tmpcord[nrow(tmpcord),1]),c(tmpcord[1,2],tmpcord[nrow(tmpcord),2]),type="l",col="blue",lty=2,lwd=4)
                              }
                        }
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
            par(oma=c(0.5,0,0.5,max(nchar(Maindata$genesetname))/ifelse(input$Downloadthreerotatelabplotone,2,1.5)))
            colcolorall <- rep("cyan",ncol(Maindata$GSCAscore))
            colcolorall[Maindata$downloadsample] <- "blue"
            tmprowv <- F
            if (input$Downloadthreerowvplotone)
                  tmprowv <- Maindata$GSCArowclust
            heatmap.2(Maindata$threesupscore,symbreaks=F,col=threepalette(input$Downloadthreecolplotone),Colv=as.dendrogram(Maindata$GSCAclust),srtRow=ifelse(input$Downloadthreerotatelabplotone,-45,0),dendrogram=threeonedendro,trace="none",Rowv=tmprowv,labCol=NA,ColSideColors=colcolorall,main=input$Downloadmaintitlethreeplotone)
            legend("bottomleft",legend=c("Selected Samples","Not Selected Samples"),lwd=1,col=c("blue","cyan"),cex=as.numeric(input$Downloadlegcexthreeone))           
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
                  par(oma=c(0.5,0,0.5,max(nchar(Maindata$genesetname))/ifelse(input$Downloadthreerotatelabplotone,2,1.5)))
                  i <- 1
                  for(INDEX in Maindata$downloadcontext) {
                        colcolorselect[Maindata$tab$SampleType %in% INDEX] <- COLORS[i]
                        i <- i+1
                  }
                  tmprowv <- F
                  if (input$Downloadthreerowvplottwo)
                        tmprowv <- Maindata$GSCArowclust
                  heatmap.2(Maindata$threesupscore[,Maindata$downloadsample],symbreaks=F,col=threepalette(input$Downloadthreecolplottwo),labCol=NA,Rowv=tmprowv,srtRow=ifelse(input$Downloadthreerotatelabplottwo,-45,0),dendrogram=threetwodendro,trace="none",ColSideColors=colcolorselect[Maindata$downloadsample],main=input$Downloadmaintitlethreeplottwo,breaks=seq(min(Maindata$threesupscore),max(Maindata$threesupscore),length.out=101))
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
      
      #####   Mainmethod : Utilities   #####
      
      data("geneIDdata")
      Utidata <- reactiveValues()
      
      observe({
            if (input$Utimethod=="Input") {
                  UtiFileHandle <- input$UtiFile  
                  if (!is.null(UtiFileHandle)) {
                        Utidata$Maindata <- Utidata$rawdata <- read.table(UtiFileHandle$datapath,header=input$Utiheader,sep=input$Utisep,quote=input$Utiquote,stringsAsFactors=F,blank.lines.skip=TRUE)            
                  }
            }
      })
      
      output$utishowdata <- renderDataTable(Utidata$Maindata)
      
      output$Uticonvertselectcolui <- renderUI({
            if (!is.null(Utidata$Maindata))
                  selectInput("Uticonvertselectcol","Select column to be converted",1:ncol(Utidata$Maindata))
      })
      
      observe({
            if (input$Uticonvertbut > 0) {
                  isolate({
                        colid <- as.numeric(input$Uticonvertselectcol)
                        fromname <- paste(input$Utifromspecies,input$Utifromtype,sep="_")
                        toname <- paste(input$Utitospecies,input$Utitotype,sep="_")                        
                        Utidata$Maindata[,colid] <- geneIDdata[match(Utidata$Maindata[,colid], geneIDdata[,fromname]),toname] 
                        Utidata$Maindata <- Utidata$Maindata[complete.cases(Utidata$Maindata),,drop=F]
                  })
            }
      })
      
      observe({
            if (input$Utiresetbut > 0) {
                  isolate({
                        Utidata$Maindata <- Utidata$rawdata
                  })
            }
      })
      
      output$Utidownloadbut <- downloadHandler(
            filename = function() { "ConvertedID.csv" },
            content = function(file) {
                  write.csv(Utidata$Maindata,file,row.names=F)     
            }
      )
      
      
})

