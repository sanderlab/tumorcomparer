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
                          selectInput("preComputedType", "Cancer Type", choices=tcgaTypes),
                          selectizeInput("preComputedDisplayCategory", "Display Match Categories", 
                                         choices=c("All", categorizations), 
                                         multiple=TRUE, 
                                         selected = "All")
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
                          selectInput("userType", "Cancer Type", choices=tcgaTypes), 
                          fileInput('mutFile', 'Choose Mutation File',
                                    accept=c('text/plain', '.txt')),
                          fileInput('cnaFile', 'Choose Copy Number File',
                                    accept=c('text/plain', '.txt')),
                          helpText("Example: ", a(href="single_cell_line_example/cell_line_MUT.txt", 
                                                  target="_blank", download="cell_line_MUT.txt", "Mutation File")),
                          helpText("Example: ", a(href="single_cell_line_example/cell_line_CNA.txt", 
                                                  target="_blank", download="cell_line_CNA.txt", "Copy Number File"))
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
