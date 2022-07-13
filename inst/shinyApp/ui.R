library(shiny)
library(ggplot2)
#library(htmlwidgets)
library(plotly)
library(shinycustomloader)
library(markdown)

#library(cyjShiny)

library(DT)
library(tumorcomparer)

# Common code goes in global.R

shinyUI(
  navbarPage("TumorComparer",
             header = list(
               tags$head(includeScript("www/js/google_analytics.js"), 
               tags$head(includeScript("www/js/nav_append_static.js")),
               tags$meta(name="description", content="Compare experimental model systems (e.g., cell lines) to patient samples by various -omic profiles (e.g., expression, mutation, copy number)."))
             ),
             tabPanel("Pre-Computed Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          selectInput("preComputedType", "Cancer Type", choices=tcgaTypes),
                          selectInput("gene_set", "Select Gene Set", choices=genesets, selected = "Most Variable Genes")
                        ),
                        mainPanel(
                          p("TumorComparer uses tumor and cell line data to guide experimental model choice; 24 cancers are assessed using NCI TCGA and Sanger COSMIC data. Please cite: ", a(href="https://doi.org/10.1016/j.crmeth.2021.100039", "Sinha R, et al., 2021.")),
                          h3("Similarity Plot"),
                          div(align="left", plotlyOutput("preComputedPlot", height=600, width=800)), 
                          h3("Similarity Table"),
                          downloadLink("preComputedDownload", "Download Table as Tab-Delimited File"),
                          DT::dataTableOutput("preComputedTable")
                        )
                      )
             ),
             tabPanel("User Data Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          fileInput('datasetFile', 'Choose Dataset File',
                                    accept=c('application/zip', '.zip')),
                          numericInput("default_weight", "Default (Background) Weight:", 0.01, min = 0, max = 100, step = 0.01),
                          helpText(br(), a(href="read_data_for_running_tc.zip", 
                                           target="_blank", download="read_data_for_running_tc.zip", "Sample Rectum Adenocarcinoma (READ) Dataset (.zip)")),
                          helpText(br(), a(href="https://zenodo.org/record/4627644", 
                                           target="_blank", "Datasets for multiple cancers from TumorComparer publication (PMID: 35475239)"))
                        ),
                        mainPanel(
                          # Results showing in Tabs (can use navlistPanel to show on left)
                          tabsetPanel(
                            tabPanel(
                              "Ranked Results",
                              h3("Similarity Plot"),
                              div(align="left", 
                                  withLoader(plotlyOutput("userPlot", height=600, width=600), type="html", loader="loader3")
                              ), 
                              h3("Similarity Table"),
                              downloadLink("userDownload", "Download Table as Tab-Delimited File"),
                              withLoader(DT::dataTableOutput("userTable"), type="html", loader="loader3")
                            ),
                            tabPanel(
                              "MDS Plot",
                              p("Tumors: Small blue points; Cell Lines: Labeled points; See documentation for more details"),
                              textOutput("userStress"), 
                              withLoader(plotlyOutput("userMdsPlot", height=600, width=600), type="html", loader="loader3")
                            )
                            #tabPanel(
                            #  "Similarity Network",
                            #  div(style="display: inline-block;vertical-align:top; width: 200px; margin-top: 1px;", selectizeInput("selected_dist_mat", label = "Select Distance Matrix", choices = c("Combined"), selected = "Combined")),
                            #  div(style="display: inline-block;vertical-align:top; width: 200px; margin-top: 1px;", numericInput("corr_threshold", "Similarity Threshold", value = 0.85, min = 0, max = 1, step = 0.01)),
                            #  cyjShinyOutput('corr_network_out', width = "100%", height = "800px")
                            #  # visNetworkOutput("corr_network_out",  height = "800px")
                            #)
                          )
                        )
                      )
             ),
             tabPanel("About",
                      includeMarkdown("www/files/about.md"),
                      h1("Version"),
                      p(paste0("TumorComparer: ", packageVersion("tumorcomparer")))
             )
  )
)
