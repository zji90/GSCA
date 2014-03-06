######################################################
##           GSCA:Gene Set Context Analysis         ##
##             Interactive User Interface           ##
##                     UI File                      ##
##           Author:Zhicheng Ji, Hongkai Ji         ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################

library(shiny)

sidebarPanel3 <- function (...) 
{
      div(class = "span3", tags$form(class = "well", ...))
}

shinyUI(pageWithSidebar(
      
      headerPanel('GSCA: Gene Set Context Analysis'),
      
      sidebarPanel3(
            
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
                               list("Input Geneset Data"="Input",
                                    "Select Geneset and Compendium"="Select",
                                    "GSCA"="GSCA",
                                    "Download"="Download",
                                    "About"="About")
                  )
            ),
            
            conditionalPanel(condition="input.Mainmethod=='Input'",
                             wellPanel(
                                   h4("Input Geneset"),
                                   radioButtons("InputGenesetmethod","",
                                                list("Specify Gene ID"="InputspecifyGeneset",
                                                     "Upload Geneset File"="InputuploadGeneset"
                                                )
                                   ),   
                                   textInput("InputGenesetname","Input Geneset Name","Geneset 1"),
                                   conditionalPanel(condition="input.InputGenesetmethod == 'InputspecifyGeneset'",
                                                    helpText("Multiple Entrez GeneID should be seperated by ;"),
                                                    textInput("InputActGeneID","Specify Activated Entrez GeneID"),
                                                    textInput("InputRepGeneID","Specify Repressed Entrez GeneID")
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
                                                                 'Comma(csv)'),
                                                    radioButtons('InputGenesetquote', 'Quote',
                                                                 c(None='',
                                                                   'Double Quote'='"',
                                                                   'Single Quote'="'"),
                                                                 'Double Quote')
                                   ),
                                   h4("Input Precise Pattern"),
                                   helpText("Multiple uploaded genesets share the same pattern"),
                                   selectInput("Inputgenesetpatternactivity","Choose Pattern",
                                               list("High activity" = "High",
                                                    "Low activity" = "Low")),
                                   selectInput("Inputgenesetpatterncotype","Cutoff Type",
                                               list("Norm" = "Norm",
                                                    "Quantile" = "Quantile",
                                                    "Exprs" = "Exprs")),
                                   textInput("Inputgenesetpatternco","Cutoff Value","0.1"),
                                   ####### Incubation: more exhaustive options for specifying pattern #######
                                   #                                    radioButtons("InputPatternmethod","",
                                   #                                                 list("Specify Pattern"="InputspecifyPattern",
                                   #                                                      "Upload Pattern File"="InputuploadPattern"
                                   #                                                 )
                                   #                                    ),                                   
                                   #                                    conditionalPanel(condition="input.InputPatternmethod == 'InputspecifyPattern'",
                                   #                                                     helpText("If multiple genesets are uploaded simultaneously, they will share the same pattern"),
                                   #                                                     selectInput("Inputgenesetpatternactivity","Geneset Regulatory Pattern",
                                   #                                                                 list("High activity" = "High",
                                   #                                                                      "Low activity" = "Low")),
                                   #                                                     selectInput("Inputgenesetpatterncotype","Cutoff Type",
                                   #                                                                 list("Norm" = "Norm",
                                   #                                                                      "Quantile" = "Quantile",
                                   #                                                                      "Exprs" = "Exprs")),
                                   #                                                     textInput("Inputgenesetpatternco","Cutoff Value","0.1")
                                   #                                    ),
                                   #                                    conditionalPanel(condition="input.InputPatternmethod == 'InputuploadPattern'",
                                   #                                                     helpText("First column: Gene Entrez ID"),
                                   #                                                     helpText("Second column:1 for activated gene, -1 for repressed gene"),
                                   #                                                     fileInput('InputGenesetFile', 'Choose File'),
                                   #                                                     checkboxInput('InputGenesetheader', 'Header', FALSE),
                                   #                                                     radioButtons('InputGenesetsep', 'Separator',
                                   #                                                                  c('Comma(csv)'=',',
                                   #                                                                    'Semicolon'=';',
                                   #                                                                    'Tab'='\t'),
                                   #                                                                  'Comma(csv)'),
                                   #                                                     radioButtons('InputGenesetquote', 'Quote',
                                   #                                                                  c(None='',
                                   #                                                                    'Double Quote'='"',
                                   #                                                                    'Single Quote'="'"),
                                   #                                                                  'Double Quote')
                                   #                                    ),
                                   p(actionButton("Inputgenesetadd","Add Genedata Information"))
                             ),
                             wellPanel(
                                   h4("Delete Existing Genedata"),
                                   uiOutput("Inputgenesetdeleteui"),
                                   p(actionButton("Inputgenesetdelete","Delete selected genedata")),
                                   p(actionButton("Inputgenesetreset","Reset all genedata"))
                             ),
                             p(downloadButton("Savegenedatafile","Save current genedata as csv file"))
                             
            ),
            
            conditionalPanel(condition="input.Mainmethod=='Select'",
                             wellPanel(
                                    wellPanel(
                                         h4("Select Geneset"),
                                         uiOutput("Summarydataselect")
                                         ),
                                    wellPanel(
                                         h4("Select Compendium"),
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
                                    checkboxInput("Summarycompscale","Scale expression values across samples"),
                                    radioButtons("Summarygenesetactmethod","Choose method of defining geneset activity",
                                                list("Weighted average"="average",
                                                     "Medium"="medium"))
                             )
            ),
            
            conditionalPanel(condition="input.Mainmethod=='GSCA'",
                             wellPanel(
                                   radioButtons("GSCAmethod","",
                                                list("Precise pattern selection"="GSCAdefault",
                                                     "Interactive pattern selection"="GSCAinteractive")
                                   ),
                                   wellPanel(
                                         conditionalPanel(condition="input.GSCAmethod=='GSCAinteractive'",
                                                          checkboxInput("GSCAinteractivesaveload","save/load interactive POI"),
                                                          conditionalPanel(condition="input.GSCAinteractivesaveload==1",
                                                                           p(downloadButton('GSCAinteractivesavebutton','Save current POI')),
                                                                           fileInput('GSCAinteractiveload', 'Load exact POI file'),
                                                                           p(actionButton('GSCAinteractiveloadbutton','Load POI'))
                                                            ),
                                                          uiOutput("InputGSCAsidebar")),
                                         uiOutput("plotenrichedareaui"),
                                         ######## Incubation: Suppress color in heatmap ############
                                         #uiOutput("heatmapcolorsuppressui"),
                                         uiOutput("heatmapthreerowvui"),
                                         helpText("Context Cutoff and Display options"),
                                         wellPanel(
                                               textInput("Inputpvalco","Enrichment P-value cutoff","0.05"),
                                               textInput("Inputfoldchangeco","Enrichment Foldchange cutoff","1.5"),
                                               radioButtons("Inputcontexttype","Choose Biological Context displaying Method",
                                                            list("Display top ranked contexts"="Toprank",
                                                                 "Display specified contexts"="Specify")),
                                               conditionalPanel(condition="input.Inputcontexttype=='Toprank'",
                                                                uiOutput("InputNslider")),
                                               conditionalPanel(condition="input.Inputcontexttype=='Specify'",
                                                                uiOutput("InputGSCAspecifycontextui"))
                                         )
                                   )
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='Download'",
                             wellPanel(
                                   h4("Download GSCA outputs"),
                                   radioButtons("Downloadregionselect","Choose pattern of interest type",choices=c("Precise","Interactive")),
                                   wellPanel(
                                         selectInput("Downloadranktabletype","Select File Type",choices=c("csv","txt")),
                                         textInput("Downloadranktablefilename","Enter File Name","GSCA Ranking Table"),
                                         p(downloadButton("Downloadranktable","Save Ranking Table"))
                                   ),
                                   uiOutput("Downloadsidebarui")
                             )
            )               
      ),
      
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
                             tabPanel("Analysis data summary",
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
                                            uiOutput("GSCAinteractiveplotthreezoominallpartsui")
                                   ),
                                   tabPanel("Ranking Table",dataTableOutput("GSCArankingtable")),
                                   tabPanel("3D scatterplot",helpText("Should have X11 installed on your computer; Only available with three genesets"),webGLOutput("RGLplot",height="800px"))
                             )                
            ),
            conditionalPanel(condition="input.Mainmethod=='Download'",
                             tabsetPanel(
                                   tabPanel("Plot",uiOutput("Downloadshowplotui")),
                                   tabPanel("Ranking Table",dataTableOutput("Downloadshowrankingtable"))
                             )
            ),
            conditionalPanel(condition="input.Mainmethod=='About'",
                             p('GSCA: Gene Set Context Analysis'),
                             p('Current Version: 0.99.1'),
                             p('Release Date: 2014-2-15'),
                             p('Author: Zhicheng Ji,Hongkai Ji'),
                             p('Maintainer: Zhicheng Ji <zji4@jhu.edu>')
            )
            
      )
)
)