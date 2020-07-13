library(shiny)
library(ggplot2)
library(markdown)
library(htmlwidgets)
library(plotly)
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
                          div(align="center", plotlyOutput("preComputedPlot", height=600, width=600)), 
                          h3("Data Table"),
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
                                     target="_blank", download="ovarian_tcga_cclp.zip", "Sample Ovarian Dataset")),
                        ),
                        mainPanel(
                          div(align="center", plotlyOutput("userPlot", height=600, width=600)), 
                          DT::dataTableOutput("userTable")
                        )
                      )
             ),
             tabPanel("About",
                      includeMarkdown("www/files/about.md")
             )
  )
)
