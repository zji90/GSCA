######################################################
##           GSCA:Gene Set Context Analysis         ##
##             Interactive User Interface           ##
##                     UI File                      ##
##           Author:Zhicheng Ji, Hongkai Ji         ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################

shinyUI(pageWithSidebar(
      
      headerPanel('GSCA: Gene Set Context Analysis'),
      
      sidebarPanel(
            
            tags$head(
                  #      tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
                  #      tags$style(type="text/css", "select { max-width: 200px; }"),
                  #      tags$style(type="text/css", "textarea { max-width: 70px; }")
                  #      tags$style(type="text/css", ".jslider { max-width: 200px; }"),
                  #      tags$style(type='text/css', ".well { max-width: 310px; }"),
                  #      tags$style(type='text/css', ".span4 { max-width: 310px; }")
            ),
            
            helpText(a("Youtube short video demo",href="https://www.youtube.com/watch?v=wqv_dmlxdcI",target="_blank")),
            helpText(a("Show User Manual",href="GSCAmanual.pdf",target="_blank")),
            wellPanel(
                  radioButtons("Mainmethod","Main Menu",
                               list("Input Gene Set"="Input",
                                    "Select Gene Set and Compendium"="Select",
                                    "GSCA Analysis"="GSCA",
                                    "Save Results"="Download",
                                    "Utilities"="Utilities",
                                    "About"="About")
                  )
            ),
            
            conditionalPanel(condition="input.Mainmethod=='Input'",
                             wellPanel(
                                   h5("Input Gene Set"),
                                   radioButtons("InputGenesetmethod","",
                                                list("Specify Gene ID"="InputspecifyGeneset",
                                                     "Upload Gene Set File"="InputuploadGeneset"
                                                )
                                   ),   
                                   textInput("InputGenesetname","Gene Set name","Geneset 1"),
                                   conditionalPanel(condition="input.InputGenesetmethod == 'InputspecifyGeneset'",
                                                    helpText("Multiple Entrez GeneID should be seperated by ;"),
                                                    textInput("InputActGeneID","Specify Entrez GeneID for Positive Genes"),
                                                    textInput("InputRepGeneID","Specify Entrez GeneID for Negative Genes")
                                   ),
                                   conditionalPanel(condition="input.InputGenesetmethod == 'InputuploadGeneset'",
                                                    helpText("First column: Entrez GeneID"),
                                                    helpText("Second column: Weight (numeric)"),
                                                    helpText("(Optional) Third column: Gene Set Name"),
                                                    fileInput('InputGenesetFile', 'Choose File'),
                                                    uiOutput("InputGenesetcolnumui"),
                                                    checkboxInput('InputGenesetheader', 'Header', FALSE),
                                                    radioButtons('InputGenesetsep', 'Separator',
                                                                 c('Comma(csv)'=',',
                                                                   'Semicolon'=';',
                                                                   'Tab'='\t'),
                                                                 ','),
                                                    radioButtons('InputGenesetquote', 'Quote',
                                                                 c('None'='',
                                                                   'Double Quote'='"',
                                                                   'Single Quote'="'"),
                                                                 '')
                                   ),
                                   p(actionButton("Inputgenesetadd","Add Gene Set"))
                             ),
                             helpText("Save current gene sets as csv file"),
                             p(downloadButton("Savegenedatafile","Save")),
                             wellPanel(
                                   h5("Delete existing gene set"),
                                   uiOutput("Inputgenesetdeleteui"),
                                   p(actionButton("Inputgenesetdelete","Delete Selected Gene Set")),
                                   p(actionButton("Inputgenesetreset","Reset All Gene Sets"))
                             )
            ),
            
            conditionalPanel(condition="input.Mainmethod=='Select'",
                             wellPanel(
                                   wellPanel(
                                         h5("Select Gene Set"),
                                         uiOutput("Summarydataselect")
                                   ),
                                   checkboxInput("Summaryswitchordertf","Switch gene set order"),
                                   conditionalPanel(condition = "input.Summaryswitchordertf==1",
                                                    helpText("Switch the order of two gene sets"),
                                                    uiOutput("Summaryswitchorderui"),
                                                    p(actionButton("Summaryswitchbut","Switch"))
                                                    ),
                                   wellPanel(
                                         h5("Select Compendium"),
                                         radioButtons("Summarycompmethod","",
                                                      list("Select available GSCA compendium"="available",
                                                           "Upload user's own compendium"="upload")),
                                         conditionalPanel(condition="input.Summarycompmethod=='available'",
                                                          uiOutput("Summarycompselectui"),
                                                          uiOutput("Summarycompinfo")
                                         ),
                                         conditionalPanel(condition="input.Summarycompmethod=='upload'",
                                                          h5("Read the instructions on the right!"),
                                                          fileInput('Summaryuploadgeneexprfile', 'Choose gene expression file'),
                                                          fileInput('Summaryuploadtabfile', 'Choose annotation file')
                                         )
                                   ),
                                   wellPanel(
                                         h5("Scaling Options"),
                                         checkboxInput("Summarycompscale","Scaling and centering expression values across samples"),
                                         conditionalPanel("input.Summarycompscale==1",
                                         radioButtons("Summarycompscalemet","",list("Centering and Scaling"="zmuv","Only Centering"="zm","Only Scaling"="uv"))),
                                         radioButtons("Summarygenesetactmethod","Choose averaging method for gene set activity",
                                                      list("Weighted average"="average",
                                                           "Median"="median"))
                                   )
                             )
            ),
            
            conditionalPanel(condition="input.Mainmethod=='GSCA'",
                             wellPanel(
                                   radioButtons("GSCAmethod","",
                                                list("Numeric POI"="GSCAdefault",
                                                     "Interactive POI"="GSCAinteractive")
                                   ),
                                   wellPanel(
                                         conditionalPanel(condition="input.GSCAmethod=='GSCAdefault'",
                                                          wellPanel(
                                                                h5("Numeric POI"),
                                                                radioButtons("numericpoimethod","",list("Slider Bar"="slider","Exact Number"="number")),
                                                                uiOutput("numericpoiui"),
                                                                helpText("Quantile:"),
                                                                uiOutput("numericpoitext"),
                                                                checkboxInput("numericpoimoreopcheck","More POI cutoff options",value=T),
                                                                conditionalPanel(condition="input.numericpoimoreopcheck==1",
                                                                                 uiOutput("numericpoimoreopgenesetnameui"),
                                                                                 selectInput("numericpoimoreopbound","Choose upper or lower bound",
                                                                                             list("Upper bound" = "Upper",
                                                                                                  "Lower bound" = "Lower")),
                                                                                 selectInput("numericpoimoreopcutofftype","Choose cutoff type",
                                                                                             list("Standard deviation from mean" = "sd",
                                                                                                  "Normal fit quantile" = "Norm",
                                                                                                  "Quantile" = "Quantile")),
                                                                                 textInput("numericpoimoreopvalue","Enter Value","2"),
                                                                                 p(actionButton("numericpoimoreopbutton","Apply New Cutoff"))
                                                                )
                                                          )
                                         ),
                                         conditionalPanel(condition="input.GSCAmethod=='GSCAinteractive'",uiOutput("InputGSCAsidebar")),
                                         uiOutput("plotenrichedareaui"),
                                         uiOutput("heatmapcolorsuppressui"),
                                         uiOutput("heatmapthreerowvui"),
                                         wellPanel(
                                               h5("Specify biological contexts"),
                                               textInput("Inputpvalco","Enrichment adjusted p-value cutoff","0.05"),
                                               textInput("Inputfoldchangeco","Enrichment foldchange cutoff","1.5"),
                                               radioButtons("Inputcontexttype","Choose biological context displaying method",
                                                            list("Display top ranked contexts"="Toprank",
                                                                 "Display specified contexts"="Specify")),
                                               conditionalPanel(condition="input.Inputcontexttype=='Toprank'",
                                                                uiOutput("InputNslider")),
                                               conditionalPanel(condition="input.Inputcontexttype=='Specify'",
                                                                uiOutput("InputGSCAspecifycontextui"))
                                         ),
                                         wellPanel(
                                               h5("Save Current POI"), 
                                               p(downloadButton('GSCAinteractivesavebutton','Save Current POI')),
                                               h5("Load POI"), 
                                               fileInput('GSCAinteractiveload', 'Choose POI file'),
                                               p(actionButton('GSCAinteractiveloadbutton','Load POI'))
                                         )
                                   )
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='Download'",
                             wellPanel(                               
                                   radioButtons("Downloadregionselect","Choose POI type",choices=c("Numeric","Interactive")),
                                   wellPanel(
                                         h5("Download ranking table"),
                                         selectInput("Downloadranktabletype","File Type",choices=c("csv","txt")),
                                         textInput("Downloadranktablefilename","File Name","GSCA Ranking Table"),
                                         p(downloadButton("Downloadranktable","Save Ranking Table"))
                                   ),
                                   wellPanel(
                                         h5("Download plots"),
                                         uiOutput("Downloadsidebarui"))
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='Utilities'",                              
                             h5("ENTREZ ID Conversion Tool"),
                             wellPanel(
                                   radioButtons("Utimethod","",c("Input File"="Input","Convert and Download"="Convert")),
                                   conditionalPanel(condition="input.Utimethod=='Input'",                                                    
                                                    fileInput('UtiFile', 'Choose File'),
                                                    checkboxInput('Utiheader', 'Header', FALSE),
                                                    radioButtons('Utisep', 'Separator',
                                                                 c('Comma(csv)'=',',
                                                                   'Semicolon'=';',
                                                                   'Tab'='\t'),
                                                                 ','),
                                                    radioButtons('Utiquote', 'Quote',
                                                                 c('None'='',
                                                                   'Double Quote'='"',
                                                                   'Single Quote'="'"),
                                                                 '')                  
                                   ),
                                   conditionalPanel(condition="input.Utimethod=='Convert'",
                                                    uiOutput("Uticonvertselectcolui"),
                                                    selectInput("Utifromspecies","Select original species",list("Human"="human","Mouse"="mouse")),
                                                    selectInput("Utifromtype","Select original data type",list("ENTREZ ID"="ENTREZ","Gene Name"="genename")),
                                                    selectInput("Utitospecies","Select target species",list("Human"="human","Mouse"="mouse")),
                                                    selectInput("Utitotype","Select target data type",list("ENTREZ ID"="ENTREZ","Gene Name"="genename")),
                                                    p(actionButton("Uticonvertbut","Convert"),actionButton("Utiresetbut","Reset")),
                                                    downloadButton("Utidownloadbut","Download")
                                   )                                   
                             )
                             
                             
                             
                             
                             
            )
            ,width=3),      
      
      
      mainPanel(
            uiOutput("GSCAstatusui"),
            conditionalPanel(condition="input.Mainmethod=='Input'",
                             tabsetPanel(
                                   tabPanel("Input Gene Set",uiOutput("OutputCurrentGenedatainstui"),h4("Current gene set"),dataTableOutput("OutputCurrentGenedata"),br(h4("All gene sets in GSCA")),textOutput("OutputGenedataname")),
                                   tabPanel("Gene Set Summary", dataTableOutput("OutputAllGenedata")), 
                                   tabPanel("Gene Set Details", uiOutput("Indigenesetnameui"), dataTableOutput("Indigeneset"))
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='Select'",
                             tabsetPanel(
                                   tabPanel("Gene Set Summary",
                                            conditionalPanel(condition="input.Summarycompmethod=='upload'",
                                                             checkboxInput("Summarycompuploadinfo","Hide instructions"),
                                                             conditionalPanel(condition="input.Summarycompuploadinfo==0",
                                                                              h4("Important! Read instructions before preparing files"),
                                                                              p('GSCA requires rigorous file format if users want to upload their own gene expression data and annotation files'),
                                                                              p('For gene expression data file: the file should be in txt format and should be separated by space. Each row corresponds to the expression of one single gene. The first column should be gene ENTREZ ID and all other columns are the gene expressions in all samples. All fields should be numeric. DO NOT include header in the file.'),
                                                                              p('For annotation file: the file should be in txt format and should be separated by space. Each row corresponds to one single sample. First column: sample ID; Second column: experiment ID; Third column: sampleType. NO space is allowed in any of the entries and they SHOULD be replaced by other separators like "_". DO NOT include header in the file.'),
                                                                              p('Each column in the gene expression file (except the first geneID column) SHOULD correspond to each row in the annotation file in order.'),
                                                                              p('Files unable to meet the requirements could fail to be read in or lead to unpredictable error.'),
                                                                              p('Example for gene expression data file:'),
                                                                              p(br('10000 -0.315 0.457 0.658 0.685 0.651 0.677'),br('10001 0.166 0.009 0.098 1.108 1.183 1.446'),br('10002 0.303 -0.39 -0.149 0.686 1.068 0.066')),
                                                                              p('Example for annotation file:'),
                                                                              p(br('GSM132917 GSE5681 skidlcl_cells:normal'),br('GSM132918 GSE5681 skidlcl_cells:normal'),br('GSM132920 GSE5681 skidlcl_cells:normal'),br('GSM148748 GSE6475 skin:normal'),br('GSM148763 GSE6475 skin:normal'),br('GSM148765 GSE6475 skin:normal')),
                                                                              p('The expression of gene 10001 in sample GSM132918 is thus 0.009')
                                                             )
                                            ),
                                            dataTableOutput("OutputDataSummary"),
                                            uiOutput("Outputmissinggenesetreport")
                                   ),
                                   tabPanel("Gene Set Breakdown",
                                            uiOutput("genesetbreakdownnameui"),
                                            uiOutput("genesetbreakdowntreenumui"),
                                            p(actionButton("genesetbreakdownaddbutton","Add sub gene sets")),
                                            plotOutput("genesetbreakdownclustplot")
                                   )
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='GSCA'",
                             tabsetPanel(
                                   tabPanel("Plot",
                                            conditionalPanel(condition="input.GSCAmethod=='GSCAdefault'",uiOutput("GSCAdefaultplot")),  
                                            conditionalPanel(condition="input.GSCAmethod=='GSCAinteractive'",uiOutput("GSCAinteractiveplot")),
                                            uiOutput("GSCAinteractiveplotthreezoominallpartsui"),
                                            uiOutput("GSCAinteractiveplotthreeplotplusui")
                                   ),
                                   tabPanel("Ranking Table",dataTableOutput("GSCArankingtable"))
                                   #tabPanel("3D scatterplot",helpText("Should have X11 installed on your computer; Only available with three genesets"),webGLOutput("RGLplot",height="800px"))
                             )                
            ),
            conditionalPanel(condition="input.Mainmethod=='Download'",
                             tabsetPanel(
                                   tabPanel("Plot",uiOutput("Downloadshowplotui")),
                                   tabPanel("Ranking Table",dataTableOutput("Downloadshowrankingtable"))
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='Utilities'",
                              dataTableOutput("utishowdata")                
            ),
            conditionalPanel(condition="input.Mainmethod=='About'",
                             p('GSCA: Gene Set Context Analysis'),
                             p('Current Version: 1.5.0'),
                             p('Release Date: 2015-1-6'),
                             p('Author: Zhicheng Ji,Hongkai Ji'),
                             p('Maintainer: Zhicheng Ji <zji4@jhu.edu>'),
                             p(a("GSCA Github home page",href="https://github.com/zji90/GSCA",target="_blank")),
                             p(a("Visit my home page",href="http://www.biostat.jhsph.edu/~zji4/",target="_blank")),
                             p(a("Visit web page of our lab",href="http://www.biostat.jhsph.edu/~hji/",target="_blank"))
            )
            
      )
)
)