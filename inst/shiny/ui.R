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
      div(class = "span2", tags$form(class = "well", ...))
}

shinyUI(pageWithSidebar(
      
      headerPanel('GSCA: Gene Set Context Analysis'),
      
      sidebarPanel3(
            
            tags$head(
            #      tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
            #      tags$style(type="text/css", "select { max-width: 200px; }"),
            #      tags$style(type="text/css", "textarea { max-width: 70px; }")
            #      tags$style(type="text/css", ".jslider { max-width: 200px; }"),
            #     tags$style(type='text/css', ".well { max-width: 310px; }"),
            #      tags$style(type='text/css', ".span4 { max-width: 310px; }")
            ),
            
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
                                                    helpText("Second column: 1 for activated gene, -1 for repressed gene"),
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
                                   h4("Input Activity Pattern"),
                                   helpText("Multiple uploaded genesets share the same pattern"),
                                   selectInput("Inputgenesetpatternactivity","Geneset Regulatory Pattern",
                                               list("High activity" = "High",
                                                    "Low activity" = "Low")),
                                   selectInput("Inputgenesetpatterncotype","Cutoff Type",
                                               list("Norm" = "Norm",
                                                    "Quantile" = "Quantile",
                                                    "Exprs" = "Exprs")),
                                   textInput("Inputgenesetpatternco","Cutoff Value","0.1"),
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
                              )
            ),
            
            conditionalPanel(condition="input.Mainmethod=='Select'",
                             wellPanel(
                                    uiOutput("Summarydataselect"),
                                    uiOutput("Summarycompselectui")
                             )
            ),
             
            conditionalPanel(condition="input.Mainmethod=='GSCA'",
                             wellPanel(
                                    radioButtons("GSCAmethod","",
                                                list("Default Enrichment Region Selection"="GSCAdefault",
                                                     "Interactive Enrichment Region Selection"="GSCAinteractive")
                                    ),
                                    wellPanel(
                                          conditionalPanel(condition="input.GSCAmethod=='GSCAinteractive'",
                                                           helpText("Choose your interested region"),
                                                           wellPanel(uiOutput("InputGSCAsidebar"))),
                                          uiOutput("plotenrichedareaui"),
                                          uiOutput("heatmapthreerowvui"),
                                          helpText("Context Cutoff and Display options"),
                                          wellPanel(
                                                textInput("Inputpvalco","Enrichment P-value cutoff","0.05"),
                                                textInput("Inputfoldchangeco","Enrichment Foldchange cutoff","1.5"),
                                                radioButtons("Inputcontexttype","Choose Biological Context displaying Method",
                                                       list("Display top ranked Contexts"="Toprank",
                                                            "Specify Contexts"="Specify")),
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
                                    radioButtons("Downloadregionselect","Select Genedata Acitivity Region",choices=c("GSCAdefault","GSCAinteractive")),
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
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$h5("GSCA status: Calculating...Please Wait...")),
            conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                             tags$h5("GSCA status: idle")),
            conditionalPanel(condition="input.Mainmethod=='Input'",
                              tabsetPanel(
                                   tabPanel("Current Geneset Data",uiOutput("OutputCurrentGenedatawarnui"),dataTableOutput("OutputCurrentGenedata"),br(h4("All genedata in GSCA:")),textOutput("OutputGenedataname")),
                                   tabPanel("All Geneset Data", dataTableOutput("OutputAllGenedata")), 
                                   tabPanel("All Geneset Pattern", dataTableOutput("OutputAllPattern"))
                              )
            ),
            conditionalPanel(condition="input.Mainmethod=='Select'",
                              dataTableOutput("OutputDataSummary"),
                              uiOutput("Outputmissinggenesetreport")
            ),
            conditionalPanel(condition="input.Mainmethod=='GSCA'",
                             tabsetPanel(
                                    tabPanel("Plot",
                                             conditionalPanel(condition="input.GSCAmethod=='GSCAdefault'",uiOutput("GSCAdefaultplot")),  
                                             conditionalPanel(condition="input.GSCAmethod=='GSCAinteractive'",uiOutput("GSCAinteractiveplot"))
                                    ),
                                    tabPanel("Ranking Table",dataTableOutput("GSCArankingtable"))
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
                             p('Current Version: 0.99.1.'),
                             p('Release Date: 2014-2-15'),
                             p('Author: Zhicheng Ji,Hongkai Ji'),
                             p('Maintainer: Zhicheng Ji <zji4@jhu.edu>')
                             )
            
            )
      )
)
