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
            
            helpText(a("Youtube short video demo",href="https://www.youtube.com/watch?v=1OeZ1PAUMhw",target="_blank")),
            helpText(a("Show User Manual",href="GSCAmanual.pdf",target="_blank")),
            wellPanel(
                  radioButtons("Mainmethod","Main Menu",
                               list("Input geneset"="Input",
                                    "Select geneset and compendium"="Select",
                                    "GSCA analysis"="GSCA",
                                    "Save results"="Download",
                                    "Utilities"="Utilities",
                                    "About"="About")
                  )
            ),
            
            conditionalPanel(condition="input.Mainmethod=='Input'",
                             wellPanel(
                                   h5("Input Geneset"),
                                   radioButtons("InputGenesetmethod","",
                                                list("Specify gene ID"="InputspecifyGeneset",
                                                     "Upload geneset file"="InputuploadGeneset"
                                                )
                                   ),   
                                   textInput("InputGenesetname","Input geneset name","Geneset 1"),
                                   conditionalPanel(condition="input.InputGenesetmethod == 'InputspecifyGeneset'",
                                                    helpText("Multiple Entrez GeneID should be seperated by ;"),
                                                    textInput("InputActGeneID","Specify activated Entrez GeneID"),
                                                    textInput("InputRepGeneID","Specify repressed Entrez GeneID")
                                   ),
                                   conditionalPanel(condition="input.InputGenesetmethod == 'InputuploadGeneset'",
                                                    helpText("First column: Entrez GeneID"),
                                                    helpText("Second column: Weight (numeric)"),
                                                    helpText("(Optional) Third column: Geneset Name"),
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
                                                                 '"')
                                   ),
                                   p(actionButton("Inputgenesetadd","Add geneset"))
                             ),
                             helpText("Save current genesets as csv file"),
                             p(downloadButton("Savegenedatafile","Save")),
                             wellPanel(
                                   h5("Delete existing geneset"),
                                   uiOutput("Inputgenesetdeleteui"),
                                   p(actionButton("Inputgenesetdelete","Delete selected geneset")),
                                   p(actionButton("Inputgenesetreset","Reset all genesets"))
                             )
            ),
            
            conditionalPanel(condition="input.Mainmethod=='Select'",
                             wellPanel(
                                   wellPanel(
                                         h5("Select Geneset"),
                                         uiOutput("Summarydataselect")
                                   ),
                                   wellPanel(
                                         h5("Select Compendium"),
                                         radioButtons("Summarycompmethod","",
                                                      list("Select available GSCA compendium"="available",
                                                           "Upload your own compendium"="upload")),
                                         conditionalPanel(condition="input.Summarycompmethod=='available'",
                                                          uiOutput("Summarycompselectui"),
                                                          uiOutput("Summarycompinfo")
                                         ),
                                         conditionalPanel(condition="input.Summarycompmethod=='upload'",
                                                          h5("See instructions on the right!"),
                                                          fileInput('Summaryuploadgeneexprfile', 'Choose gene expression file'),
                                                          fileInput('Summaryuploadtabfile', 'Choose annotation file')
                                         )
                                   ),
                                   wellPanel(
                                         h5("Other options"),
                                         checkboxInput("Summarycompscale","Scale expression values across samples"),
                                         radioButtons("Summarygenesetactmethod","Choose method of defining geneset activity",
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
                                                                uiOutput("numericpoisliderui"),
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
                                                                                 textInput("numericpoimoreopvalue","Enter value","2"),
                                                                                 p(actionButton("numericpoimoreopbutton","Apply new cutoff"))
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
                                               h5("Save/Load POI"), 
                                               p(downloadButton('GSCAinteractivesavebutton','Save current POI')),
                                               fileInput('GSCAinteractiveload', 'Load exact POI file'),
                                               p(actionButton('GSCAinteractiveloadbutton','Load POI'))
                                         )
                                   )
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='Download'",
                             wellPanel(
                                   h5("Download GSCA outputs"),
                                   radioButtons("Downloadregionselect","Choose POI type",choices=c("Numeric","Interactive")),
                                   wellPanel(
                                         h5("Download ranking table"),
                                         selectInput("Downloadranktabletype","Select File Type",choices=c("csv","txt")),
                                         textInput("Downloadranktablefilename","Enter File Name","GSCA Ranking Table"),
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
                                                                 '"')                  
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
                                   tabPanel("Current Geneset Data",uiOutput("OutputCurrentGenedatawarnui"),dataTableOutput("OutputCurrentGenedata"),br(h4("All genedata in GSCA:")),textOutput("OutputGenedataname")),
                                   tabPanel("Input Summary", h4("Geneset Summary"), dataTableOutput("OutputAllGenedata"), h4("Precise Pattern Summary"), dataTableOutput("OutputAllPattern")), 
                                   tabPanel("Individual Geneset", uiOutput("Indigenesetnameui"), dataTableOutput("Indigeneset"))
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='Select'",
                             tabsetPanel(
                                   tabPanel("Genesets summary",
                                            dataTableOutput("OutputDataSummary"),
                                            uiOutput("Outputmissinggenesetreport"),
                                            conditionalPanel(condition="input.Summarycompmethod=='upload'",
                                                             checkboxInput("Summarycompuploadinfo","Hide instruction for upload"),
                                                             conditionalPanel(condition="input.Summarycompuploadinfo==0",
                                                                              h4("Important! See instructions before preparing files"),
                                                                              p('GSCA requires rigorous file format if you want to upload your own gene expression data and annotation files'),
                                                                              p('For gene expression data file: the file should be in txt format and entries should be separated by space. Each row stands for expression for one gene. Row order SHOULD correspond to the column order in the gene expression data file (For example, first row should annotate the first column in gene expression data). The first column should be gene ENTREZ ID and all other columns stands for expression in samples. Data types in all entries should be numeric. DO NOT include header in the file.'),
                                                                              p('For annotation file: the file should be in txt format and entries should be separated by space. Each row stands for a sample. First column: sample ID; Second column: experiment ID; Third column: sampleType. NO space is allowed in any of the entries and they SHOULD be replaced by other separators like "_". DO NOT include header in the file.'),
                                                                              p('Files unable to meet the requirements could fail to be read in or lead to unpredictable error.'),
                                                                              p('Example for gene expression data file:'),
                                                                              p(br('10000 -0.315 -0.457 -0.658 -0.685 -0.651 -0.677'),br('10001 0.166 0.009 0.098 -1.108 -1.183 -1.446'),br('10002 -0.303 -0.39 -0.149 -0.686 -1.068 0.066')),
                                                                              p('Example for annotation file:'),
                                                                              p(br('GSM132917 GSE5681 skidlcl_cells:normal'),br('GSM132918 GSE5681 skidlcl_cells:normal'),br('GSM132920 GSE5681 skidlcl_cells:normal'),br('GSM148748 GSE6475 skin:normal'),br('GSM148763 GSE6475 skin:normal'),br('GSM148765 GSE6475 skin:normal'))
                                                             )
                                            )
                                   ),
                                   tabPanel("Geneset breakdown",
                                            uiOutput("genesetbreakdownnameui"),
                                            uiOutput("genesetbreakdowntreenumui"),
                                            p(actionButton("genesetbreakdownaddbutton","Add sub genesets")),
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
                             p('Current Version: 1.0.0'),
                             p('Release Date: 2014-3-20'),
                             p('Author: Zhicheng Ji,Hongkai Ji'),
                             p('Maintainer: Zhicheng Ji <zji4@jhu.edu>'),
                             p(a("Visit web page of our lab",href="http://www.biostat.jhsph.edu/~hji/",target="_blank"))
            )
            
      )
)
)