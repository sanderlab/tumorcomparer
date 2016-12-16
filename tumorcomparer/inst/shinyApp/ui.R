library(shiny)
library(ggplot2)
library(markdown)
library(metabologram)
library(htmlwidgets)
library(plotly)
library(DT)

# Common code goes in global.R

shinyUI(
  navbarPage("Pan-Cancer Metabolism Data Explorer",
             tabPanel("Data Overview",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          selectInput("volcanoPlotStudy", "Cancer Type", choices=unique(fc[,"Study"])), 
                          selectInput("volcanoPlotType", "Query Type", choices=c("Druggability"="db", "HMDB"="hmdb")),
                          uiOutput("overview")
                        ),
                        mainPanel(
                          div(align="center", plotlyOutput("volcanoPlot")),
                          DT::dataTableOutput("ptsTable")
                        )
                      )
             ),
             tabPanel("Univariate Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          selectInput("uniType", "Plot Data", choices=c("Fold Change", "All Data")), 
                          uiOutput("uniUi")
                        ),
                        mainPanel(
                          div(align="center", plotOutput("uniPlot", height=600, width=600))
                        )
                      )
             ),
             tabPanel("Bivariate Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          selectInput("biType", "Plot Data", choices=c("Fold Change", "All Data")),
                          uiOutput("biUi")
                        ),
                        mainPanel(
                          div(align="center", plotOutput("biPlot", height=600, width=600))
                        )
                      )
             ),
             tabPanel("Clinical Variables",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          selectInput("clinMetab", "Metabolite", choices=colnames(metdata)[1:(ncol(metdata)-2)], selected = "kynurenine"),
                          selectInput("clinFea", "Clinical Feature", choices=c("Stage", "Grade")),
                          selectInput("clinStudies", "Studies", choices=unique(metdata[,"Study"]))
                        ),
                        mainPanel(
                          div(align="center", plotOutput("clinPlot", height=600, width=600))
                          #div(align="center", plotlyOutput("clinPlot"))
                        )
                      )
             ),
             # tabPanel("Metabologram",
             #          sidebarLayout(
             #            sidebarPanel(
             #              width=3,
             #              selectInput("pathway", "Pathway Name", choices=names(allmgrams_KIRC),selected = names(allmgrams_KIRC)[1] ),
             #              br(),
             #              selectInput("study", "Study", 
             #                          choices=c('BLCA','BRCA','BRCATang','COAD','KIRC','OV','PAAD','PAADHussein1','PAADHussein2','PRAD','PRADLODA','STAD'), 
             #                          selected="KIRC")  
             #            ),
             #            mainPanel(
             #              metabologramOutput("metabologram", height=600, width=600)
             #            )
             #          )
             # ),
             tabPanel("About",
                      includeMarkdown("www/files/about.md")
             ),
             header = list(
               tags$head(tags$style(".shiny-plot-output{height:90vh !important;}"))
             )
  )
)
