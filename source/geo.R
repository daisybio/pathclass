# Download and process data from GEO
getGEO_data <- eventReactive(input$geo_dl_button, {
    # Create a Progress object
    progress <- shiny::Progress$new(style="old")

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Downloading data from GEO", value = 0)

    # check if it is a proper GEO id
    if(!grepl("GSE[0-9]*$", input$geo_id)){
        shinysky::showshinyalert(session, "geo_error", "Only GEO GSE ids are supported, e.g. GSE45827." ,"danger")
        return(NULL)
    }

    # load series and platform data from GEO
    gset <- tryCatch({
      getGEO(input$geo_id,
             GSEMatrix = TRUE,
             AnnotGPL=TRUE)[[1]]

    }, error = function(e) {
      shinysky::showshinyalert(session, "geo_error", paste("Download failed: ", e$message, sep="") ,"danger")
      return(NULL)})

    if(is.null(gset)) return(NULL)

    progress$set(message = "Final processing of GEO data", value = 0.8)

    # make proper column names to match toptable
    fvarLabels(gset) <- make.names(fvarLabels(gset))

    return(gset)
})

# Only when we have results do we want to show the respective result tabs...
startGEO <- observeEvent(input$geo_button, {
    getGEO_subtypes()

    shinyjs::show(selector = "#geo_nav li a[data-value='Predictions Plot']")
    shinyjs::show(selector = "#geo_nav li a[data-value='Predictions Table']")
    shinyjs::show(selector = "#geo_nav li a[data-value='Performance (subtypes)']")
    shinyjs::show(selector = "#geo_nav li a[data-value='Performance (predictors)']")
    shinyjs::show(selector = "#geo_nav li a[data-value='Confusion matrices']")
})

# Reactive element that checks if a file has been uploaded already...
output$geo_file_uploaded <- reactive({
    return(!is.null(getGEO_data()))
})

# which we need to know on the front end. We can sneak in this variable in this way.
outputOptions(output, 'geo_file_uploaded', suspendWhenHidden=FALSE)

# Get the column names in the GEO meta data
getGEO_phenoColumns <- observeEvent(input$geo_dl_button,{
    gset <- getGEO_data()

    if(is.null(gset)) return(NULL)

    pheno_cols <- gset@phenoData@varMetadata
    updateSelectInput(session,
                      "geo_pheno_columns",
                      "Select a class label column (optional)",
                      rownames(pheno_cols))
})

# Assign class labels from the selected meta data column in GEO
getGEO_classLabels <- observeEvent(input$geo_pheno_columns, {
    gset <- getGEO_data()

    labels <- unique(as.character(gset@phenoData@data[,input$geo_pheno_columns]))

    for(subtype in brca_subtypes){
        updateSelectInput(session, paste("subtype", subtype, sep="_"), subtype, labels)
    }
})

# Map the class labels selected above to the ones we expect.
getGEO_mappedClassLabels <- reactive({
    lumA <- input$subtype_LumA
    lumB <- input$subtype_LumB
    basal <- input$subtype_Basal
    her2 <- input$subtype_Her2

    gset <- getGEO_data()
    labels <- as.character(gset@phenoData@data[,input$geo_pheno_columns])

    for(subtype in brca_subtypes){
        labels[which(labels == input[[paste("subtype", subtype, sep = "_")]])] <- subtype
    }

    labels[-which(labels %in% brca_subtypes)] <- "Other"

    return(labels)
})

# Where the magic happens...
getGEO_subtypes <- eventReactive(input$geo_button, {

    # Create a Progress object
    progress <- shiny::Progress$new(style="old")

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    # load series and platform data from GEO
    gset <- getGEO_data()

    progress$set(message = "Processing GEO data", value = 0.1)

    # transform array signal
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }

    #scaling
    ex <- t(scale(t(ex), center = TRUE, scale = FALSE))

    # get gene ids
    entrez_ids <- featureData(gset)@data$Gene.ID
    entrez_ids <- cbind(entrez_ids = as.character(entrez_ids), probe = rownames(ex))

    #if there are several entrez ids for one probe. we make one row for each
    progress$set(message = "Processing GEO data", value = 0.15,
                 detail = "Expanding data for multiple gene matches")

    entrez_ids <- tidyr::separate_rows(as.data.frame(entrez_ids), entrez_ids, sep = "///")

    ex <- ex[entrez_ids$probe,]
    rownames(ex) <- entrez_ids$entrez_ids

    # map gene ids. when several probes match the same gene use the one
    # that is most variable
    progress$set(message = "Processing GEO data", value = 0.2,
                 detail = "Selecting most variable probes for each gene")

    ex <- fn_filter_most_variable(ex, entrez_ids$entrez_ids)

    progress$set(message = "Performing pathway enrichment analysis", value = 0.2)

    # Use ssGSEA and the existing random forest predictor to get class labels
    result <- fn_ssgsea(ex,
                        progress = progress,
                        selected.pathways = selected.pathways,
                        all.rf = all.rf,
                        pathway.sources = input$selected_predictors)

    result <- spread(result, "pathways", "subtype")

    progress$set(message = "Predicting subtype labels", value = 0.9)

    # get gene expression based predictions, e.g. PAM50
    gfu_result <- geneFu(t(ex))

    # get the reference class labels
    annotated_result <- getGEO_mappedClassLabels()

    #put it all together
    final_result <- cbind(result, gfu_result[result$Sample,],
                          `GEO annotation` = annotated_result)
    final_result <- as.data.frame(final_result)

    progress$set(message = "Returning results", value = 1.0)

    #select 1st result tab
    updateTabsetPanel(session, "geo_nav",
        selected = "Predictions Table"
    )

    return(final_result)
})
