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
  # PRECOMPUTED
  preComputedDat <- reactive({
    tcgaType <- input$preComputedType
    displayCategory <- input$preComputedDisplayCategory
    #displayCategory <- "Great"
    #tcgaType <- "OV"
    
    # Read pre-computed results
    mdsfitDf <- read.table(system.file("extdata", "isomdsfit", paste0("isoMDSfit_", tcgaType, "_coordinates.txt"), package = "tumorcomparer"), header = TRUE, sep="\t")
    categorizationDf <- read.table(system.file("extdata", "cell_line_categorization", paste0("categorization_", tcgaType, ".txt"), package = "tumorcomparer"), header = TRUE, sep="\t")
    
    tmpDf <- merge(mdsfitDf, categorizationDf, all.x=TRUE)
    
    # Filter
    if("All" %in% displayCategory) {
      categoryFilter <- c(NA, categorizations)
    } else {
      categoryFilter <- c(NA, categorizations[categorizations %in% displayCategory])
    }
    
    idx <- which(tmpDf$Category %in% categoryFilter)
    tmpDf <- tmpDf[idx, ]
    
    # Round numeric columns 
    is.num <- sapply(tmpDf, is.numeric)
    tmpDf[is.num] <- lapply(tmpDf[is.num], round, 3)
    
    return(tmpDf)
  })
  
  output$preComputedPlot <- renderPlotly({
    df <- preComputedDat()
    
    # Initialize columns
    df$color <- "blue"
    df$size <- 1
    
    idx <- which(df$Sample_Type == "Cell_Line")
    df$color[idx] <- "orange"
    df$size[idx] <- 3
    
    toolTip <- paste(df$Sample_ID, '</br>Coordinate 1: ', df$Coordinate1, '</br>Coordinate 2: ', df$Coordinate2, '</br>Category: ', df$Category)
    df$toolTip <- toolTip

    p1 <- ggplot() + 
      geom_point(data=df, 
                 aes(x=Coordinate1, y=Coordinate2, color=Sample_Type, text=toolTip),
                 size = 2, 
                 alpha=0.4) + 
      theme_minimal() +
      scale_color_manual(name="", values=c("orange", "blue")) + 
      xlab("Coordinate 1") + ylab("Coordinate 2")

    #print(p1)
    #cat("A", colnames(df))
    ggplotly(p1, tooltip=c("toolTip"), source="preComputedPlot") %>% 
      config(cloud=FALSE, collaborate=FALSE, displaylogo=FALSE)
  })
    
  output$preComputedTable <- DT::renderDataTable({
    df <- preComputedDat()
    
    DT::datatable(df, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE)
  })
    
  # USER 
  userDat <- reactive({
    mutFile <- input$mutFile
    cnaFile <- input$cnaFile
    userType <- input$userType
    
    cell_line_mut_file <- mutFile$datapath
    cell_line_cna_file <- cnaFile$datapath
    
    if (is.null(cell_line_mut_file) || is.null(cell_line_cna_file))
      return(NULL)
    
    pancancer_gene_weights_file <- system.file("extdata", "weights", "default_weights_for_known_cancer_genes.txt", package="tumorcomparer")
    cancer_specific_gene_weights_file <- system.file("extdata", "weights", paste0("Genes_and_weights_TCGA_", userType, "_based.txt"), package="tumorcomparer")
    tumor_mut_file <- system.file("extdata", "tumor_alterations", paste0("TCGA_", userType, "_MUT.txt"), package="tumorcomparer")
    tumor_cna_file <- system.file("extdata", "tumor_alterations", paste0("TCGA_", userType, "_CNA.txt"), package="tumorcomparer")
    cell_line_mut_file <- cell_line_mut_file
    cell_line_cna_file <- cell_line_cna_file
    # Do not save the output
    output_composite_alteration_matrix_file <- NULL
    # Use weighted correlation
    distance_similarity_measure <- "weighted_correlation"

    # DEBUG TEST     
    a <- read.table(pancancer_gene_weights_file, sep="\t")
    
    comparison <- run_comparison(pancancer_gene_weights_file=pancancer_gene_weights_file, 
                                 cancer_specific_gene_weights_file=cancer_specific_gene_weights_file, 
                                 tumor_mut_file=tumor_mut_file, 
                                 tumor_cna_file=tumor_cna_file,
                                 cell_line_mut_file=cell_line_mut_file,
                                 cell_line_cna_file=cell_line_cna_file, 
                                 output_composite_alteration_matrix_file=output_composite_alteration_matrix_file,
                                 distance_similarity_measure=distance_similarity_measure)
    
    categorization <- categorize_cell_lines(fraction_of_tumors_for_comparison=0.1, 
                                            dist_mat=comparison$dist_mat,
                                            composite_mat=comparison$composite_mat,
                                            cell_lines_with_both_MUT_and_CNA=comparison$cell_lines_with_both_MUT_and_CNA, 
                                            tumors_with_both_MUT_and_CNA=comparison$tumors_with_both_MUT_and_CNA,
                                            trim_cell_line_names=TRUE)
  
    tmpDf <- merge(comparison$isomdsfit, categorization$categorization, all.x=TRUE)
    tmpDf$Sample_Type <- c(rep("Cell_Line", 1), rep("Tumor", nrow(tmpDf)-1))
    
    # DEBUG
    cat("START userDat")
    cat(str(comparison$isomdsfit))
    cat(str(categorization$categorization))
    cat(str(tmpDf))
    cat("END userDat")
    
    # Round numeric columns 
    is.num <- sapply(tmpDf, is.numeric)
    tmpDf[is.num] <- lapply(tmpDf[is.num], round, 3)
    
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
      config(cloud=FALSE, collaborate=FALSE, displaylogo=FALSE)
    
    #print(p1)
  })
    
  output$userTable <- DT::renderDataTable({
    df <- userDat()
    
    # Only the first is the user data
    df <- df[1,]
    
    DT::datatable(df, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE)
  })
    
  # # IGNORE
  # # Get reactive data
  # volcanoPlotDat <- reactive({
  #   #study <- "BRCA"
  #   #queryClass <- "Druggable"
  #   #queryType <- "db"
  #   study <- input$volcanoPlotStudy
  #   queryType <- input$volcanoPlotType
  #   
  #   fcCol <- paste0("FC.", study)
  #   padjCol <- paste0("Padj.", study)
  #   
  #   differentialAbundanceSummary <- filter(differentialAbundanceSummary, abs(differentialAbundanceSummary[[fcCol]]) > 0)
  # 
  #   tmpMerge <- merge(differentialAbundanceSummary, drugbankHmdb, by.x="X", by.y="finalName")
  #   
  #   tmpFc <- round(tmpMerge[[fcCol]], digits=2)
  #   tmpNegLogPadj <- round(-log10(tmpMerge[[padjCol]]), digits=2)
  #   
  #   tmpDat <- data.frame(
  #     name=tmpMerge$X, 
  #     FC=tmpFc, 
  #     Padj=tmpMerge[[padjCol]], 
  #     negLogPadj=tmpNegLogPadj, 
  #     db=tmpMerge$drugbankDrug, 
  #     dbGene=tmpMerge$drugbankGene, 
  #     hmdb=tmpMerge$hmdbClass, 
  #     hmdbId=tmpMerge$hmdbId
  #   )
  #   
  #   tmpDat <- tmpDat[with(tmpDat, order(Padj)), ]
  # 
  #   toolTip <- paste(tmpDat$name, '</br>', tmpDat$hmdb, '</br>Druggable (Drugbank): ', tmpDat$db)
  #   tmpDat$toolTip <- toolTip
  #   
  #   #tmpDat <- as_data_frame(tmpDat)
  # 
  #   if(queryType == "db") {
  #     queryClass <- "Druggable"
  #     
  #     tmpDat[, "color"] <- rep(paste0("Not ", queryClass), nrow(tmpDat))
  #     tmpDat[which(!is.na(tmpDat[,queryType])), "color"] <- queryClass
  #   } else {
  #     queryClass <- input$volcanoPlotClass
  #     
  #     tmpDat[, "color"] <- rep(paste0("Not ", queryClass), nrow(tmpDat))
  #     tmpDat[which(tmpDat[,queryType] == queryClass), "color"] <- queryClass
  #   }
  #   
  #   #write.table(tmpDat, "~/Downloads/tmpDat.txt", sep="\t", quote=FALSE, row.names=FALSE)
  #   
  #   if(length(unique(tmpDat$color)) == 2) {
  #     levels <- c(queryClass, paste0("Not ", queryClass))
  #     tmpDat$color <- factor(tmpDat$color, levels=levels)
  #   }
  #   
  #   return(tmpDat)
  # })
  # 
  # output$volcanoPlot <- renderPlotly({
  #   tmpDat <- volcanoPlotDat()
  #   
  #   #write.table(tmpDat, "~/Downloads/tmpPlot.txt", sep="\t", quote=FALSE, row.names=FALSE)
  #   
  #   ## Construct the plot object
  #   # NOTE: Cannot change labels
  #   ptColor <- c("red", "black")
  #   
  #   # SET COLORS
  #   ## Check if only one class
  #   if(length(unique(tmpDat$color)) == 1) {
  #     ptColor <- "black"
  #   }
  #   
  #   ## Check if the legend colors need to be reversed because of the legend labels are alphabetical
  #   #if(length(input$volcanoPlotClass) > 0 && input$volcanoPlotClass > "Not") {
  #   #  ptColor <- c("black", "red")
  #   #} 
  #   
  #   ## Check to make sure the colors are always in the correct order for the drug information
  #   #if(input$volcanoPlotType == "db") {
  #   #  ptColor <- c("red", "black")
  #   #}
  #   
  #   # CREATE PLOT
  #   p <- ggplot() +
  #     geom_point(data=tmpDat, aes(x=FC, y=negLogPadj, color=color, text=toolTip), alpha=0.4, size=2) +
  #     theme_bw() +
  #     xlab("log2 (Tumor/Normal)") + ylab("-log10 P-Value") +
  #     scale_color_manual(name="", values=ptColor)
  #   
  #   # Interactive
  #   ggplotly(p, width=900, tooltip=c("toolTip", "FC", "negLogPadj"), source="volcanoPlot") %>% config(cloud=FALSE, collaborate=FALSE, displaylogo=FALSE) %>% layout(dragmode = "select")
  # })
  # 
  # selectedPoints <- reactive({
  #   #DEBUG
  #   #str(reactType$type)
  #   
  #   # Check if user if last action was on plot or menu; force to NULL any existing selection data (i.e. eventData) if the user made a menu change
  #   if(length(reactType$type) > 0 && reactType$type == "lastUpdatedValue") {
  #     #cat("LASTUPDATEDVALUE")
  #     eventData <- NULL       
  #   } else {
  #     #cat("EVENTDATA")
  #     eventData <- reactValues$eventData 
  #     pts <- as.numeric(eventData$pointNumber)
  #     
  #     #DEBUG
  #     cat("eDa: ", str(eventData), "\n")
  #   }
  #   
  #   tmpDat <- volcanoPlotDat()
  #   
  #   results <- tmpDat
  #   
  #   #DEBUG 
  #   #cat("tDat: ", str(tmpDat), "\n")
  #   
  #   # is.data.frame(x) check accounts for empty selections
  #   if(!is.null(eventData) && is.data.frame(eventData)) {
  #     idx <- which(!grepl("^Not", tmpDat$color))
  #     notIdx <- which(grepl("^Not ", tmpDat$color))
  #     
  #     # Check the case where there are no entries of the queryClass
  #     if(length(notIdx) == nrow(tmpDat)) {
  #       idx <- notIdx
  #     }
  #     
  #     # NOTE: EXTREMELY IMPORTANT FOR DEBUGGING GGPLOTLY AUTOMATICALLY ALPHABETIZES CLASSES 
  #     # RESULTING IN THE NEED TO ACCOUNT FOR FLIPPING
  #     # Flip entries post-Not
  #     #if(input$volcanoPlotClass > "Not") {
  #     #  t1 <- idx 
  #     #  idx <- notIdx
  #     #  notIdx <- t1 
  #     #}
  #     
  #     # DEBUG
  #     #cat("L1: ", str(idx), "\n")
  #     #cat("L2: ", str(notIdx), "\n")
  #     
  #     tmpB1 <- tmpDat[idx, ]
  #     tmpB2 <- tmpDat[notIdx, ]
  #     
  #     class1 <- tmpB1[subset(eventData, curveNumber == 0)$pointNumber + 1,]
  #     class2 <- tmpB2[subset(eventData, curveNumber == 1)$pointNumber + 1,]
  #     
  #     #cat("C1: ", str(class1), "\n")
  #     #cat("C2: ", str(class2), "\n")
  # 
  #     plotSubset <- rbind(class1, class2)
  #     #cat("PS1: ", str(plotSubset), "\n")
  #     
  #     plotSubset <- plotSubset[!is.na(as.character(plotSubset$name)),] 
  #     #cat("PS2: ", str(plotSubset), "\n")
  #     
  #     results <- plotSubset
  #     
  #     # DEBUG
  #     cat("dat: ", str(plotSubset), "\n")
  #     cat("evt: ", eventData[["pointNumber"]], "\n")
  #     cat("pts: ", pts, "\n")
  #   }
  #   
  #   #write.table(tmpDat, "~/Downloads/tmpClick.txt", sep="\t", quote=FALSE, row.names = FALSE)
  # 
  #   return(results) 
  # })
  
  # START: CHECK LAST USER INTERACTION ----
  # NOTE: Starting point http://stackoverflow.com/questions/41022310/detect-specific-ui-change-in-r-shiny
  # The below code sets the value of eventData and records the last menu value as reactive variables
  # A second reactive variable records whether the last user interaction was on the menu (i.e. lastUpdatedValue) or the plot was modified (i.e. eventData)
  # reactValues <- reactiveValues(lastUpdatedValue=NULL, eventData=NULL)
  # reactType <- reactiveValues(type=NULL)
  # 
  # observe({
  #   vars <- c("volcanoPlotStudy", "volcanoPlotType", "volcanoPlotClass")
  #   lapply(vars, function(x) {
  #     observe({
  #       input[[x]]
  #       reactValues$lastUpdatedValue <- input[[x]]
  #     })
  #   })
  # })
  # 
  # observe({
  #   eventData <- event_data("plotly_selected", source="volcanoPlot")
  #   reactValues$eventData <- eventData
  # })
  # 
  # observe({
  #   vars <- c("lastUpdatedValue", "eventData")
  #   lapply(vars, function(x) {
  #     observe({
  #       reactValues[[x]]
  #       reactType$type <- x
  #     })
  #   })
  # })
  # END: CHECK LAST USER INTERACTION ----
})
