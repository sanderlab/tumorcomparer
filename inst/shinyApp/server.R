library(shiny)
library(ggplot2)
library(markdown)
library(plotly)

library(magrittr)
library(dplyr)
library(DT)

library(tumorcomparer)

# Common code goes in global.R

shinyServer(function(input, output, session) {
  # PRECOMPUTED ----
  preComputedDat <- reactive({
    tcgaType <- input$preComputedType
    df <- mtc_dataset
    
    # Filter
    idx <- which(df$Tumor_Cancer_Type == tcgaType)
    df <- df[idx, ]
    
    # Round numeric columns 
    is.num <- sapply(df, is.numeric)
    df[is.num] <- lapply(df[is.num], round, 3)
    
    return(df)
  })
  
  output$preComputedPlot <- renderPlotly({
    df <- preComputedDat()
    tcgaType <- input$preComputedType
    
    validate(need(nrow(df) > 0, "ERROR: Cancer type does not have pre-computed results"))
    
    p <- plot_balloon_plot(df, tcgaType)

    #print(p)
    #cat("A", colnames(df))
    ggplotly(p, tooltip=c("toolTip"), source="preComputedPlot") %>% 
      # Remove buttons from: https://community.plotly.com/t/how-to-remove-mode-bar-buttons/34002
      config(cloud=FALSE, displaylogo=FALSE, modeBarButtonsToRemove = c(
        "select2d", "sendDataToCloud", "pan2d", "resetScale2d",
        "hoverClosestCartesian", "hoverCompareCartesian",
        "lasso2d", "zoomIn2d", "zoomOut2d", "toggleSpikelines"
      ))
  })
    
  output$preComputedTable <- DT::renderDataTable({
    df <- preComputedDat()
    
    # Filter dataset
    df <- df[, mtc_selected_columns]
    colnames(df) <- names(mtc_selected_columns)
    
    DT::datatable(df, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE)
  })
  
  output$preComputedDownload <- downloadHandler(
    filename = function() {
      paste0("table_", input$preComputedType, ".txt")
    },
    content = function(file) {
      df <- preComputedDat()
      df <- df[df$Sample_Type == "Cell_Line",]
      
      write.table(df, file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    }
  )
    
  # USER INPUT ----
  userDat <- reactive({
    dataset_file <- input$datasetFile
    dataset_path <- dataset_file$datapath
    
    dataset_path <- system.file("shinyApp", "www", "ovarian_tcga_cclp.zip", package="tumorcomparer")
    tmp_file_names <- unzip(dataset_path, list=TRUE)
    tmp_dir <- tempdir()
    
    files <- c("cell_line_cna.txt", "cell_line_exp.txt", "cell_line_mut.txt", 
               "default_weights_for_known_cancer_genes_cna.txt", "default_weights_for_known_cancer_genes_exp.txt", "default_weights_for_known_cancer_genes_mut.txt",
               "genes_and_weights_cna.txt", "genes_and_weights_exp.txt", "genes_and_weights_mut.txt", 
               "tumor_cna.txt", "tumor_exp.txt", "tumor_mut.txt")
    
    # Check files are present
    are_files_present <- all(sapply(files, function(file) { any(grepl(file, tmp_file_names$Name)) }, USE.NAMES = FALSE))
    validate(
      need(are_files_present, 
           paste0("ERROR: Missing file. Case-sensitive. Following ", length(files), " files needed: ", paste(files, collapse=", ")))
    )
    
    unzip(dataset_path, junkpaths = TRUE, exdir = tmp_dir)
    
    # NOTE: Example files are embedded in the package and are accessible with system.file()
    tumor_mut_file <- file.path(tmp_dir, "tumor_mut.txt")
    tumor_cna_file <- file.path(tmp_dir, "tumor_cna.txt")
    tumor_exp_file <- file.path(tmp_dir, "tumor_exp.txt")
    
    cell_line_mut_file <- file.path(tmp_dir, "cell_line_mut.txt")
    cell_line_cna_file <- file.path(tmp_dir, "cell_line_cna.txt")
    cell_line_exp_file <- file.path(tmp_dir, "cell_line_exp.txt")
    
    known_cancer_gene_weights_mut_file <- file.path(tmp_dir, "default_weights_for_known_cancer_genes_mut.txt")
    known_cancer_gene_weights_cna_file <- file.path(tmp_dir, "default_weights_for_known_cancer_genes_cna.txt")
    known_cancer_gene_weights_exp_file <- file.path(tmp_dir, "default_weights_for_known_cancer_genes_exp.txt")
    
    cancer_specific_gene_weights_mut_file <- file.path(tmp_dir, "Genes_and_weights_mut.txt")
    cancer_specific_gene_weights_cna_file <- file.path(tmp_dir, "Genes_and_weights_cna.txt")
    cancer_specific_gene_weights_exp_file <- file.path(tmp_dir, "Genes_and_weights_exp.txt")
    
    comparison_result <- run_comparison(
      available_data_types=c("mut", "cna", "exp"), 
      mut_data_type_weight = 1/3,
      cna_data_type_weight = 1/3,
      exp_data_type_weight = 1/3,
      cna_default_weight=0.01, 
      mut_default_weight=0.01,
      exp_default_weight=0.01,
      cna_known_cancer_gene_weight=0.1, 
      mut_known_cancer_gene_weight=0.1, 
      exp_known_cancer_gene_weight=0.1, 
      tumor_mut_file=tumor_mut_file, 
      tumor_cna_file=tumor_cna_file, 
      tumor_exp_file=tumor_exp_file, 
      cell_line_mut_file=cell_line_mut_file, 
      cell_line_cna_file=cell_line_cna_file, 
      cell_line_exp_file=cell_line_exp_file, 
      known_cancer_gene_weights_mut_file=known_cancer_gene_weights_mut_file, 
      known_cancer_gene_weights_cna_file=known_cancer_gene_weights_cna_file, 
      known_cancer_gene_weights_exp_file=known_cancer_gene_weights_exp_file, 
      cancer_specific_gene_weights_mut_file=cancer_specific_gene_weights_mut_file, 
      cancer_specific_gene_weights_cna_file=cancer_specific_gene_weights_cna_file, 
      cancer_specific_gene_weights_exp_file=cancer_specific_gene_weights_exp_file,
      distance_similarity_measures=c("generalized_jaccard", 
                                     "generalized_jaccard", 
                                     "weighted_correlation"))
    
    # See returned outputs
    names(comparison_result)
    
    tmpDf
  })
    
  output$userPlot <- renderPlotly({
    validate(need(input$mutFile, 'Upload mutation file'))
    validate(need(input$cnaFile, 'Upload copy number file'))
    
    df <- userDat()
    
    # DEBUG
    cat("START userPlot")
    str(df)
    cat("L: ", class(df))
    cat("N: ", names(df))
    cat("END userPlot")
    
    toolTip <- paste(df$Sample_ID, '</br>Coordinate 1:', df$points.1, '</br>Coordinate 2: ', df$points.2, '</br>Category: ', df$Category)
    df$toolTip <- toolTip
    
    p1 <- ggplot() + 
      geom_point(data=df, 
                 aes(x=points.1, y=points.2, color=Sample_Type, text=toolTip),
                 size = 2, 
                 alpha=0.4) + 
      theme_minimal() +
      scale_color_manual(name="", values=c("orange", "blue")) +
      xlab("Coordinate 1") + ylab("Coordinate 2")
    
    ggplotly(p1, tooltip=c("toolTip"), source="userPlot") %>% 
      config(cloud=FALSE, displaylogo=FALSE)
    
    #print(p1)
  })
    
  output$userTable <- DT::renderDataTable({
    df <- userDat()
    df <- df[df$Sample_Type == "Cell_Line",]
    
    # Only the first is the user data
    df <- df[1,]
    
    DT::datatable(df, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE)
  })
})
