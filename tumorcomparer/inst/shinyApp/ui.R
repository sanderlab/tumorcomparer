library(shiny)
library(ggplot2)
library(markdown)
library(htmlwidgets)
library(plotly)
library(DT)

# Common code goes in global.R

shinyUI(
  navbarPage("TumorComparer",
             tabPanel("Pre-Computed",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          selectInput("preComputedType", "Cancer Type", choices=c("Ovarian", "Lung", "Bladder"))
                        ),
                        mainPanel(
                          div(align="center", plotOutput("preComputedPlot", height=600, width=600)), 
                          DT::dataTableOutput("preComputedTable")
                        )
                      )
             ),
             tabPanel("User Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          selectInput("userType", "Cancer Type", choices=c("Ovarian", "Lung", "Bladder"))
                        ),
                        mainPanel(
                          div(align="center", plotOutput("userPlot", height=600, width=600)), 
                          DT::dataTableOutput("userTable")
                        )
                      )
             ),
             tabPanel("About",
                      includeMarkdown("www/files/about.md")
             )
  )
)
