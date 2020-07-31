library(shiny)
library(ggplot2)
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
    
    plot_dat <- make_balloon_plot_data_from_mtc(df, tcgaType)
    p <- plot_balloon_plot(plot_dat, paste0(plot_title_prefix, ": ", tcgaType))

    #print(p)
    #cat("A", colnames(df))
    ggplotly(p, tooltip=c("toolTip"), source="preComputedPlot") %>% 
      # Remove buttons from: https://community.plotly.com/t/how-to-remove-mode-bar-buttons/34002
      config(cloud=FALSE, displaylogo=FALSE, modeBarButtonsToRemove = plotlyModeBarButtonsToRemove)
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
    
    validate(
      need(nchar(dataset_path) > 1,"ERROR: Please upload dataset as ZIP file.")
    )
    
    cat("DP: ", dataset_path, "\n")
    
    #dataset_path <- system.file("shinyApp", "www", "ovarian_tcga_cclp.zip", package="tumorcomparer")
    tmp_file_names <- unzip(dataset_path, list=TRUE)
    cat("FILES: ", paste(tmp_file_names, collapse = "; "), "\n")
    
    tmp_dir <- tempdir()
    
    # files <- c("cell_line_cna.txt", "cell_line_exp.txt", "cell_line_mut.txt", 
    #            "default_weights_for_known_cancer_genes_cna.txt", "default_weights_for_known_cancer_genes_exp.txt", "default_weights_for_known_cancer_genes_mut.txt",
    #            "genes_and_weights_cna.txt", "genes_and_weights_exp.txt", "genes_and_weights_mut.txt", 
    #            "tumor_cna.txt", "tumor_exp.txt", "tumor_mut.txt")
    # 
    # # Check files are present
    # are_files_present <- all(sapply(files, function(file) { any(grepl(file, tmp_file_names$Name)) }, USE.NAMES = FALSE))
    # validate(
    #   need(are_files_present, 
    #        paste0("ERROR: Missing file. Case-sensitive. Following ", length(files), " files needed: ", paste(files, collapse=", ")))
    # )

    unzip(dataset_path, junkpaths = TRUE, exdir = tmp_dir)
    cat("DIR: ", paste(dir(tmp_dir), collapse = "; "))
    
    # tmp_dir <- dir("inst/extdata/ovarian_tcga_cclp/")
    has_mut <- any(grepl("mut", dir(tmp_dir)))
    has_cna <- any(grepl("cna", dir(tmp_dir)))
    has_exp <- any(grepl("exp", dir(tmp_dir)))
    
    validate(
      need(any(c(has_mut, has_cna, has_exp)), 
           paste0("ERROR: Missing files. File names are all lowercase. Please consult documentation."))
    )
    
    available_data_types <- NULL 
    distance_similarity_measures <- NULL 
    
    if(has_mut) {
      tumor_mut_file <- file.path(tmp_dir, "tumor_mut.txt")
      cell_line_mut_file <- file.path(tmp_dir, "cell_line_mut.txt")
      known_cancer_gene_weights_mut_file <- file.path(tmp_dir, "default_weights_for_known_cancer_genes_mut.txt")
      cancer_specific_gene_weights_mut_file <- file.path(tmp_dir, "genes_and_weights_mut.txt")
      available_data_types <- c(available_data_types, "mut")
      distance_similarity_measures <- c(distance_similarity_measures, "generalized_jaccard")
    }

    if(has_cna) {
      tumor_cna_file <- file.path(tmp_dir, "tumor_cna.txt")
      cell_line_cna_file <- file.path(tmp_dir, "cell_line_cna.txt")
      known_cancer_gene_weights_cna_file <- file.path(tmp_dir, "default_weights_for_known_cancer_genes_cna.txt")
      cancer_specific_gene_weights_cna_file <- file.path(tmp_dir, "genes_and_weights_cna.txt")
      available_data_types <- c(available_data_types, "cna")
      distance_similarity_measures <- c(distance_similarity_measures, "generalized_jaccard")
    }
    
    if(has_exp) {
      tumor_exp_file <- file.path(tmp_dir, "tumor_exp.txt")
      cell_line_exp_file <- file.path(tmp_dir, "cell_line_exp.txt")
      known_cancer_gene_weights_exp_file <- file.path(tmp_dir, "default_weights_for_known_cancer_genes_exp.txt")
      cancer_specific_gene_weights_exp_file <- file.path(tmp_dir, "genes_and_weights_exp.txt")
      available_data_types <- c(available_data_types, "exp")
      distance_similarity_measures <- c(distance_similarity_measures, "weighted_correlation")
    }
    
    mut_data_type_weight <- ifelse(has_mut, 1/length(available_data_types), NULL)
    cna_data_type_weight <- ifelse(has_cna, 1/length(available_data_types), NULL)
    exp_data_type_weight <- ifelse(has_exp, 1/length(available_data_types), NULL)
    
    default_weight <- 0.01
    known_cancer_gene_weight <- 0.1
    
    # Start timer
    tic() 
    
    comparison_result <- run_comparison(
      available_data_types = available_data_types, 
      mut_data_type_weight = mut_data_type_weight,
      cna_data_type_weight = cna_data_type_weight,
      exp_data_type_weight = exp_data_type_weight,
      cna_default_weight = default_weight, 
      mut_default_weight = default_weight,
      exp_default_weight = default_weight,
      cna_known_cancer_gene_weight = known_cancer_gene_weight, 
      mut_known_cancer_gene_weight = known_cancer_gene_weight, 
      exp_known_cancer_gene_weight = known_cancer_gene_weight, 
      tumor_mut_file = tumor_mut_file, 
      tumor_cna_file = tumor_cna_file, 
      tumor_exp_file = tumor_exp_file, 
      cell_line_mut_file = cell_line_mut_file, 
      cell_line_cna_file = cell_line_cna_file, 
      cell_line_exp_file = cell_line_exp_file, 
      known_cancer_gene_weights_mut_file = known_cancer_gene_weights_mut_file, 
      known_cancer_gene_weights_cna_file = known_cancer_gene_weights_cna_file, 
      known_cancer_gene_weights_exp_file = known_cancer_gene_weights_exp_file, 
      cancer_specific_gene_weights_mut_file = cancer_specific_gene_weights_mut_file, 
      cancer_specific_gene_weights_cna_file = cancer_specific_gene_weights_cna_file, 
      cancer_specific_gene_weights_exp_file = cancer_specific_gene_weights_exp_file,
      distance_similarity_measures = distance_similarity_measures)
    
    # See returned outputs
    names(comparison_result)
    
    # Stop timer
    toc()
    
    cat("DONE")
    
    comparison_result
  })
    
  output$userPlot <- renderPlotly({
    comparison_result <- userDat()
    
    # DEBUG
    cat("START userPlot")
    str(comparison_result)
    cat("L: ", class(comparison_result))
    cat("N: ", names(comparison_result))
    cat("END userPlot")
    
    # toolTip <- paste(df$Sample_ID, '</br>Coordinate 1:', df$points.1, '</br>Coordinate 2: ', df$points.2, '</br>Category: ', df$Category)
    # df$toolTip <- toolTip
    # 
    # p1 <- ggplot() + 
    #   geom_point(data=df, 
    #              aes(x=points.1, y=points.2, color=Sample_Type, text=toolTip),
    #              size = 2, 
    #              alpha=0.4) + 
    #   theme_minimal() +
    #   scale_color_manual(name="", values=c("orange", "blue")) +
    #   xlab("Coordinate 1") + ylab("Coordinate 2")
    # 
    # ggplotly(p1, tooltip=c("toolTip"), source="userPlot") %>% 
    #   config(cloud=FALSE, displaylogo=FALSE)
    
    # Example plot
    # p <- ggplot(mtcars,aes(x=wt,y=mpg)) + 
    #   geom_point() + 
    #   xlab('Weight (x 1000lbs)') + 
    #   ylab('Miles per Gallon')
    
    plot_dat <- make_balloon_plot_data_from_comparison_result(comparison_result)
    p <- plot_balloon_plot(plot_dat, plot_title_prefix)
    
    #print(p)
    #cat("A", colnames(df))
    ggplotly(p, tooltip=c("toolTip"), source="preComputedPlot") %>% 
      # Remove buttons from: https://community.plotly.com/t/how-to-remove-mode-bar-buttons/34002
      config(cloud=FALSE, displaylogo=FALSE, modeBarButtonsToRemove = plotlyModeBarButtonsToRemove)
    
    ggplotly(p, source="userPlot") %>% 
      config(cloud=FALSE, displaylogo=FALSE)
    
    #print(p1)
  })
    
  output$userTable <- DT::renderDataTable({
    comparison_result <- userDat()
    df <- make_balloon_plot_data_from_comparison_result(comparison_result)
    
    DT::datatable(df, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE)
  })
  
  output$userDownload <- downloadHandler(
    filename = function() {
      paste0("table_results.txt")
    },
    content = function(file) {
      comparison_result <- userDat()
      df <- make_balloon_plot_data_from_comparison_result(comparison_result)

      write.table(df, file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    }
  )
})
