library(shiny)
library(ggplot2)
#library(htmlwidgets)
library(plotly)
library(shinycustomloader)
library(markdown)

library(DT)
library(tumorcomparer)

# Common code goes in global.R

shinyUI(
  navbarPage("TumorComparer",
             header = list(tags$head(includeScript("www/js/google-analytics.js"))),
             tabPanel("Pre-Computed",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          selectInput("preComputedType", "Cancer Type", choices=tcgaTypes)
                        ),
                        mainPanel(
                          h3("Results Plot"),
                          div(align="left", plotlyOutput("preComputedPlot", height=600, width=600)), 
                          h3("Results Table"),
                          downloadLink("preComputedDownload", "Download Table as Tab-Delimited File"),
                          DT::dataTableOutput("preComputedTable")
                        )
                      )
             ),
             tabPanel("User Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          fileInput('datasetFile', 'Choose Dataset File',
                                    accept=c('application/zip', '.zip')),
                          helpText("Download: ", a(href="ovarian_tcga_cclp.zip", 
                                     target="_blank", download="ovarian_tcga_cclp.zip", "Sample Ovarian Dataset (.zip)")),
                        ),
                        mainPanel(
                          h3("Results Plot"),
                          div(align="left", 
                              withLoader(plotlyOutput("userPlot", height=600, width=600), type="html", loader="loader3")
                          ), 
                          h3("Results Table"),
                          withLoader(DT::dataTableOutput("userTable"), type="html", loader="loader3")
                        )
                      )
             ),
             tabPanel("About",
                      includeMarkdown("www/files/about.md")
             )
  )
)
