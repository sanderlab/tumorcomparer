library(shiny)
library(ggplot2)
library(markdown)
library(plotly)

library(magrittr)
library(dplyr)
library(DT)

# Common code goes in global.R

shinyServer(function(input, output, session) {
    output$preComputedPlot <- renderPlot({
      p1 <- ggplot(data=iris, aes(x=Sepal.Length, y=Sepal.Width)) + geom_point()
      print(p1)
    })
    
    output$preComputedTable <- DT::renderDataTable({
      DT::datatable(iris, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE)
    })
    
    output$userPlot <- renderPlot({
      p1 <- ggplot(data=iris, aes(x=Sepal.Length, y=Sepal.Width)) + geom_point()
      print(p1)
    })
    
    output$userTable <- DT::renderDataTable({
      DT::datatable(iris, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE)
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
