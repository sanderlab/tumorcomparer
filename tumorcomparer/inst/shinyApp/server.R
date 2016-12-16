library(shiny)
library(ggplot2)
library(markdown)
library(plotly)

library(magrittr)
library(dplyr)
#library(tibble)
library(DT)

# Common code goes in global.R

shinyServer(function(input, output, session) {
    # UI
    output$overview <- renderUI({
      queryClass <- NULL
      
      if(input$volcanoPlotType == "hmdb") {
        queryClass <- selectInput(
          "volcanoPlotClass",
          "Query Class",
          choices = sort(unique(drugbankHmdb$hmdbClass)),
          selected = "Fatty Acyls"
        )
      }

      return(queryClass)
    })
  
    output$uniUi <- renderUI({
        uniMetab <- switch(
            input$uniType,
            "Fold Change" = selectInput(
                "uniMetab",
                "Metabolite",
                choices = colnames(fc)[1:(ncol(fc)-1)],
                selected = "kynurenine"
            ),
            "All Data" = selectInput(
                "uniMetab",
                "Metabolite",
                choices = colnames(metdata)[1:(ncol(metdata)-2)],
                selected = "kynurenine"
            )
        )
        
        uniStudies <- switch(
            input$uniType,
            "Fold Change" = selectizeInput(
                "uniStudies", 
                "Studies", 
                choices=c("All Studies", unique(fc[,"Study"])), 
                selected="All Studies", 
                multiple=TRUE
            ),
            "All Data" = selectizeInput(
                "uniStudies", 
                "Studies", 
                choices=c("All Studies", unique(metdata[,"Study"])), 
                selected="All Studies", 
                multiple=TRUE
            )
        )
        
        return(list(uniMetab, uniStudies))
    })
    
    output$biUi <- renderUI({
        biX <- switch(
            input$biType,
            "Fold Change" = selectInput(
                "biX",
                "Metabolite (X-Axis)",
                choices = colnames(fc)[1:(ncol(fc)-1)],
                selected = "glucose"
            ),
            "All Data" = selectInput(
                "biX",
                "Metabolite (X-Axis)",
                choices = colnames(metdata)[1:(ncol(metdata)-2)],
                selected = "glucose"
            )
        )
        
        biY <- switch(
            input$biType,
            "Fold Change" = selectInput(
                "biY",
                "Metabolite (Y-Axis)",
                choices = colnames(fc)[1:(ncol(fc)-1)],
                selected = "fructose"
            ),
            "All Data" = selectInput(
                "biY",
                "Metabolite (Y-Axis)",
                choices = colnames(metdata)[1:(ncol(metdata)-2)],
                selected = "fructose"
            )
        )
        
        biStudies <- switch(
            input$biType,
            "Fold Change" = selectizeInput(
                "biStudies", 
                "Studies", 
                choices=c("All Studies", unique(fc[,"Study"])), 
                selected="All Studies", 
                multiple=TRUE
            ),
            "All Data" = selectizeInput(
                "biStudies", 
                "Studies", 
                choices=c("All Studies", unique(metdata[,"Study"])), 
                selected="All Studies", 
                multiple=TRUE
            )
        )
        
        return(list(biX, biY, biStudies))
    })
    
    # PLOTS
    output$uniPlot <- renderPlot({
        type <- input$uniType
        metabolite <- input$uniMetab
        studies <- input$uniStudies
        
        # Hide error on load
        validate(need(metabolite != "", ""))
        
        if(!("All Studies" %in% studies) && length(studies) > 0) {
            metdata <- metdata[which(metdata$Study %in% studies),]
            fc <- fc[which(fc$Study %in% studies),]
        }
        
        if (type == 'Fold Change') {
            fcdata <- data.frame(x = fc[, metabolite], y = fc$Study)
            colnames(fcdata) <- c('X', 'Study')
            
            # Find studies which have notNA data
            notNA = unique(fcdata[which(!is.na(fcdata[, 1])), 2])
            keepidx = which(fcdata[, 2] %in% notNA)
            
            plotdata = fcdata[keepidx, ]
            plotdata$Study = factor(plotdata$Study)
            
            # Plot
            p1 = ggplot(plotdata, aes(X, fill = Study)) + geom_histogram(aes(y = ..density..)) + theme_classic() +
                geom_density(alpha = 0.5) + ylab('Density') + xlab('Log2 Ratio, Tumor:Normal') +
                geom_vline(xintercept = 0, linetype = 'longdash') +
                facet_grid(Study ~ .)
        }
        
        if (type == 'All Data') {
            alldata <-
                data.frame(x = metdata[, metabolite],
                           y = metdata$Study,
                           z = metdata$Type)
            colnames(alldata) <- c('X', 'Study', 'Type')
            
            # Find studies which have notNA data
            notNA <-
                unique(alldata[which(!is.na(alldata[, 1])), 2])
            keepidx <- which(alldata[, 2] %in% notNA)
            
            plotdata <- alldata[keepidx, ]
            plotdata$Study <- factor(plotdata$Study)
            
            # Plot
            p1 <-
                ggplot(plotdata, aes(X, fill = Type)) + geom_histogram(aes(y = ..density..)) + theme_classic() +
                geom_density(alpha = 0.5) + ylab('Density') + xlab('Log2 Ratio, Tumor:Normal') +
                facet_grid(Study ~ ., scales = 'free')
        }
        
        # Hide message: stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
        suppressMessages(print(p1))
    })
    
    output$biPlot <- renderPlot({
        type <- input$biType
        mx <- input$biX
        my <- input$biY
        studies <- input$biStudies
        
        if(!("All Studies" %in% studies) && length(studies) > 0) {
            metdata <- metdata[which(metdata$Study %in% studies),]
            fc <- fc[which(fc$Study %in% studies),]
        }
        
        # Hide error on load
        validate(need(mx != "" && my != "", ""))
        
        if (type == 'Fold Change') {
            fcdata <-
                data.frame(x = fc[, mx],
                           y = fc[, my],
                           z = fc$Study)
            colnames(fcdata) <- c('X', 'Y', 'Study')
            
            # Find studies which have notNA data
            notNA = unique(fcdata[which(!is.na(fcdata[, 1]) &
                                            !is.na(fcdata[, 2])), 2])
            keepidx = which(fcdata[, 2] %in% notNA)
            
            plotdata = fcdata[keepidx, ]
            plotdata$Study = factor(plotdata$Study)
            
            # Plot
            p1 <- ggplot(plotdata, aes(X, Y), color = Study) + geom_point(aes(color = Study)) +
                theme_classic() + facet_wrap( ~ Study) +
                ylab(input$biY) + xlab(input$biX)
            
        }
        
        if (type == 'All Data') {
            fcdata <-
                data.frame(
                    x = metdata[, mx],
                    y = metdata[, my],
                    z = metdata$Study,
                    t = metdata$Type
                )
            colnames(fcdata) <- c('X', 'Y', 'Study', 'Type')
            
            # Find studies which have notNA data
            notNA = unique(fcdata[which(!is.na(fcdata[, 1]) &
                                            !is.na(fcdata[, 2])), 2])
            keepidx = which(fcdata[, 2] %in% notNA)
            
            plotdata = fcdata[keepidx, ]
            plotdata$Study = factor(plotdata$Study)
            
            # Plot
            p1 <- ggplot(plotdata, aes(X, Y), color = Type) + geom_point(aes(color = Type)) +
                theme_classic() + facet_wrap( ~ Study) +
                ylab(input$biY) + xlab(input$biX)
            
        }
        
        # Hide message: stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
        suppressMessages(print(p1))
    })
    
    output$clinPlot <- renderPlot({
        studies <- input$clinStudies
        feature <- input$clinFea
        metabolite <- input$clinMetab
        
        # DEBUGGING (WORKS)
#         studies <- "KIRC"
#         metabolite <- "kynurenine"
#         feature <- "grade"
        
        # DEBUGGING
        #studies <- "PRADLODA"
        #metabolite <- "1-stearoylglycerophosphoethanolamine"
        #feature <- "grade"
        
        # Hide error on load
        validate(need(metabolite != "", ""))
        
        study <- unlist(lapply(strsplit(rownames(metdata), ":"), "[[", 1))
        id <- unlist(lapply(strsplit(rownames(metdata), ":"), "[[", 2))
        status <- unlist(lapply(strsplit(rownames(metdata), ":"), "[[", 3))
        metab <- metdata[, metabolite]
        
        tmpDf <- data.frame(study=study, id=id, metab=metab, stringsAsFactors = FALSE)
        
        metDf <- tmpDf[tmpDf$study %in% studies, ]
        clinDf <- clinFeaDat[clinFeaDat$study %in% studies, ]
    
        mergedDf <- merge(metDf, clinDf, by.x=c("study", "id"), by.y=c("study", "sample"))
        
        feature <- tolower(trimws(feature))
        
#         cat(colnames(mergedDf))
#         cat("A ", feature)
#         cat(" ", studies)
#         cat(" ", metabolite)
        
        # NOTE: Grade and stage are lower case
        uqclin <- sort(unique(mergedDf[, feature]))
        uqtype <- sort(unique(mergedDf$type))

        mergedDf$clinxtype <- paste(mergedDf[,feature], mergedDf$type,sep = ':')
        clinxtype_levels <- as.vector(outer(uqclin, uqtype, paste, sep=":"))
        mergedDf$clinxtype <- factor(mergedDf$clinxtype, levels=clinxtype_levels)
        
        #mergedDf$grade = factor(mergedDf$grade)
        #mergedDf$stage = factor(mergedDf$stage)
        
        # Plot
        if(feature == "stage") {
            tmpDf <- mergedDf[!is.na(mergedDf$stage),]
            
            p1 <- ggplot(tmpDf, aes(x=clinxtype, y=metab, color=type, label = id)) + geom_boxplot() + 
              ylab(paste('Log2', metabolite)) + xlab("Stage") + theme_bw() + 
              geom_point(position = 'jitter') + 
              scale_color_discrete(name = 'Tissue\n Type') + 
              theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_text(vjust=3))
        } else {
            tmpDf <- mergedDf[!is.na(mergedDf$grade),]
            print( tmpDf )
            p1 <- ggplot(tmpDf, aes(x=clinxtype, y=metab, color=type, label = id)) + geom_boxplot() +
              ylab(paste('Log2', metabolite)) + xlab("Grade") + theme_bw() + 
              geom_point(position = 'jitter') + 
              scale_color_discrete(name = 'Tissue\n Type') + 
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }

        validate(need(nrow(tmpDf) > 0, "No valid entries for the selection."))
        
        # Hide message: stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
        suppressMessages(print(p1))  
    })
    
    # VOLCANO PLOT SECTION
    volcanoPlotDat <- reactive({
      #study <- "BRCA"
      #queryClass <- "Druggable"
      #queryType <- "db"
      study <- input$volcanoPlotStudy
      queryType <- input$volcanoPlotType
      
      fcCol <- paste0("FC.", study)
      padjCol <- paste0("Padj.", study)
      
      differentialAbundanceSummary <- filter(differentialAbundanceSummary, abs(differentialAbundanceSummary[[fcCol]]) > 0)

      tmpMerge <- merge(differentialAbundanceSummary, drugbankHmdb, by.x="X", by.y="finalName")
      
      tmpFc <- round(tmpMerge[[fcCol]], digits=2)
      tmpNegLogPadj <- round(-log10(tmpMerge[[padjCol]]), digits=2)
      
      tmpDat <- data.frame(
        name=tmpMerge$X, 
        FC=tmpFc, 
        Padj=tmpMerge[[padjCol]], 
        negLogPadj=tmpNegLogPadj, 
        db=tmpMerge$drugbankDrug, 
        dbGene=tmpMerge$drugbankGene, 
        hmdb=tmpMerge$hmdbClass, 
        hmdbId=tmpMerge$hmdbId
      )
      
      tmpDat <- tmpDat[with(tmpDat, order(Padj)), ]

      toolTip <- paste(tmpDat$name, '</br>', tmpDat$hmdb, '</br>Druggable (Drugbank): ', tmpDat$db)
      tmpDat$toolTip <- toolTip
      
      #tmpDat <- as_data_frame(tmpDat)

      if(queryType == "db") {
        queryClass <- "Druggable"
        
        tmpDat[, "color"] <- rep(paste0("Not ", queryClass), nrow(tmpDat))
        tmpDat[which(!is.na(tmpDat[,queryType])), "color"] <- queryClass
      } else {
        queryClass <- input$volcanoPlotClass
        
        tmpDat[, "color"] <- rep(paste0("Not ", queryClass), nrow(tmpDat))
        tmpDat[which(tmpDat[,queryType] == queryClass), "color"] <- queryClass
      }
      
      #write.table(tmpDat, "~/Downloads/tmpDat.txt", sep="\t", quote=FALSE, row.names=FALSE)
      
      if(length(unique(tmpDat$color)) == 2) {
        levels <- c(queryClass, paste0("Not ", queryClass))
        tmpDat$color <- factor(tmpDat$color, levels=levels)
      }
      
      return(tmpDat)
    })
    
    output$volcanoPlot <- renderPlotly({
      tmpDat <- volcanoPlotDat()
      
      #write.table(tmpDat, "~/Downloads/tmpPlot.txt", sep="\t", quote=FALSE, row.names=FALSE)
      
      ## Construct the plot object
      # NOTE: Cannot change labels
      ptColor <- c("red", "black")
      
      # SET COLORS
      ## Check if only one class
      if(length(unique(tmpDat$color)) == 1) {
        ptColor <- "black"
      }
      
      ## Check if the legend colors need to be reversed because of the legend labels are alphabetical
      #if(length(input$volcanoPlotClass) > 0 && input$volcanoPlotClass > "Not") {
      #  ptColor <- c("black", "red")
      #} 
      
      ## Check to make sure the colors are always in the correct order for the drug information
      #if(input$volcanoPlotType == "db") {
      #  ptColor <- c("red", "black")
      #}
      
      # CREATE PLOT
      p <- ggplot() +
        geom_point(data=tmpDat, aes(x=FC, y=negLogPadj, color=color, text=toolTip), alpha=0.4, size=2) +
        theme_bw() +
        xlab("log2 (Tumor/Normal)") + ylab("-log10 P-Value") +
        scale_color_manual(name="", values=ptColor)
      
      # Interactive
      ggplotly(p, width=900, tooltip=c("toolTip", "FC", "negLogPadj"), source="volcanoPlot") %>% config(cloud=FALSE, collaborate=FALSE, displaylogo=FALSE) %>% layout(dragmode = "select")
    })
    
    selectedPoints <- reactive({
      #DEBUG
      #str(reactType$type)
      
      # Check if user if last action was on plot or menu; force to NULL any existing selection data (i.e. eventData) if the user made a menu change
      if(length(reactType$type) > 0 && reactType$type == "lastUpdatedValue") {
        #cat("LASTUPDATEDVALUE")
        eventData <- NULL       
      } else {
        #cat("EVENTDATA")
        eventData <- reactValues$eventData 
        pts <- as.numeric(eventData$pointNumber)
        
        #DEBUG
        cat("eDa: ", str(eventData), "\n")
      }
      
      tmpDat <- volcanoPlotDat()
      
      results <- tmpDat
      
      #DEBUG 
      #cat("tDat: ", str(tmpDat), "\n")
      
      # is.data.frame(x) check accounts for empty selections
      if(!is.null(eventData) && is.data.frame(eventData)) {
        idx <- which(!grepl("^Not", tmpDat$color))
        notIdx <- which(grepl("^Not ", tmpDat$color))
        
        # Check the case where there are no entries of the queryClass
        if(length(notIdx) == nrow(tmpDat)) {
          idx <- notIdx
        }
        
        # NOTE: EXTREMELY IMPORTANT FOR DEBUGGING GGPLOTLY AUTOMATICALLY ALPHABETIZES CLASSES 
        # RESULTING IN THE NEED TO ACCOUNT FOR FLIPPING
        # Flip entries post-Not
        #if(input$volcanoPlotClass > "Not") {
        #  t1 <- idx 
        #  idx <- notIdx
        #  notIdx <- t1 
        #}
        
        # DEBUG
        #cat("L1: ", str(idx), "\n")
        #cat("L2: ", str(notIdx), "\n")
        
        tmpB1 <- tmpDat[idx, ]
        tmpB2 <- tmpDat[notIdx, ]
        
        class1 <- tmpB1[subset(eventData, curveNumber == 0)$pointNumber + 1,]
        class2 <- tmpB2[subset(eventData, curveNumber == 1)$pointNumber + 1,]
        
        #cat("C1: ", str(class1), "\n")
        #cat("C2: ", str(class2), "\n")

        plotSubset <- rbind(class1, class2)
        #cat("PS1: ", str(plotSubset), "\n")
        
        plotSubset <- plotSubset[!is.na(as.character(plotSubset$name)),] 
        #cat("PS2: ", str(plotSubset), "\n")
        
        results <- plotSubset
        
        # DEBUG
        cat("dat: ", str(plotSubset), "\n")
        cat("evt: ", eventData[["pointNumber"]], "\n")
        cat("pts: ", pts, "\n")
      }
      
      #write.table(tmpDat, "~/Downloads/tmpClick.txt", sep="\t", quote=FALSE, row.names = FALSE)

      return(results) 
    })
    
    output$ptsTable <- DT::renderDataTable({
      tmpDat <- selectedPoints()

      #DEBUGGING
      #cat("str:", str(tmpDat), "\n")
      
      if(ncol(tmpDat) > 3) {
        selectedColumns <- c("hmdbId", "name", "FC", "negLogPadj", "db", "dbGene", "hmdb")
        tmpDat$hmdbId <- paste0('<a target="_blank" href="http://www.hmdb.ca/metabolites/', tmpDat$hmdbId, '">', tmpDat$hmdbId, '</a>')
        
        tmpDat[grepl("metabolites/NA", tmpDat$hmdbId), "hmdbId"] <- NA
        
        tmpDat <- tmpDat[, selectedColumns]

        colnames(tmpDat) <- c("HMDB", "Name", "log2 (Tumor/Normal)", "-log10 (Padj)", "DrugBank Drug", "DrugBank Gene", "HMDB Classification")
      }
      
      DT::datatable(tmpDat, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE)
    })
    
    # START: CHECK LAST USER INTERACTION ----
    # NOTE: Starting point http://stackoverflow.com/questions/41022310/detect-specific-ui-change-in-r-shiny
    # The below code sets the value of eventData and records the last menu value as reactive variables
    # A second reactive variable records whether the last user interaction was on the menu (i.e. lastUpdatedValue) or the plot was modified (i.e. eventData)
    reactValues <- reactiveValues(lastUpdatedValue=NULL, eventData=NULL)
    reactType <- reactiveValues(type=NULL)
    
    observe({
      vars <- c("volcanoPlotStudy", "volcanoPlotType", "volcanoPlotClass")
      lapply(vars, function(x) {
        observe({
          input[[x]]
          reactValues$lastUpdatedValue <- input[[x]]
        })
      })
    })
    
    observe({
      eventData <- event_data("plotly_selected", source="volcanoPlot")
      reactValues$eventData <- eventData
    })
    
    observe({
      vars <- c("lastUpdatedValue", "eventData")
      lapply(vars, function(x) {
        observe({
          reactValues[[x]]
          reactType$type <- x
        })
      })
    })
    # END: CHECK LAST USER INTERACTION ----
    
    # 	output$metabologram <- renderMetabologram({
    #       topdata <- get(paste('allmgrams_',input$study,sep=''))
    #       sampleMetabologramData <- topdata[[input$pathway]]
    #
    # 	  sampleMetabologramBreaks <- c(
    # 	    -3.00, -2.75, -2.50, -2.25, -2.00, -1.75, -1.50, -1.25, -1.00,
    # 	    -0.75, -0.50, -0.25, 0.00,  0.25,  0.50,  0.75,  1.00,
    # 	    1.25,  1.50,  1.75,  2.00,  2.25,  2.50,  2.75, 3.00
    # 	  )
    #
    # 	  sampleMetabologramColors <- c(
    # 	    "#0000FF", "#1616FF", "#2C2CFF", "#4242FF", "#5858FF", "#6E6EFF", "#8585FF",
    # 	    "#9B9BFF", "#B1B1FF", "#C7C7FF", "#DDDDFF", "#F3F3FF", "#FFF3F3", "#FFDDDD",
    # 	    "#FFC7C7", "#FFB1B1", "#FF9B9B", "#FF8585", "#FF6E6E", "#FF5858", "#FF4242",
    # 	    "#FF2C2C", "#FF1616", "#FF0000"
    # 	  )
    #
    # 	  print(metabologram(
    # 	    sampleMetabologramData,
    # 	    width=1200,
    # 	    height=1200,
    # 	    main=input$pathway,
    # 	    showLegend=TRUE,
    # 	    fontSize=10,
    # 	    legendBreaks=sampleMetabologramBreaks,
    # 	    legendColors=sampleMetabologramColors,
    # 	    legendText="Log2 Ratio Tumor Abundance:Normal Abundance"
    # 	  ))
    # 	})
    
    # 	output$downloadNormalData <- downloadHandler(
    # 	    filename="normalData.tsv",
    # 	    content = function(file) {
    # 	        write.table(normal, file, quote=FALSE, sep="\t", col.names=NA)
    # 	    }
    # 	)
})
