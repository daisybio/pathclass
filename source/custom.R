getCUSTOM_data <- reactive({

    inFile <- input$custom_file

    if (is.null(inFile))
        return(NULL)

    tryCatch(read.table(inFile$datapath,
               header = TRUE,
               sep = input$sep),
             error = function(e) {
                 shinysky::showshinyalert(session, "custom_error", paste("Upload failed: ", e$message, sep="") ,"danger")
                 return(NULL)
             })
})

output$custom_file_uploaded <- reactive({
    return(!is.null(getCUSTOM_data()))
})

outputOptions(output, 'custom_file_uploaded', suspendWhenHidden=FALSE)

getCUSTOM_classLabels <- reactive({

    inFile <- input$custom_class_labels
    if(is.null(inFile)) return(NULL)
    labels <- scan(inFile$datapath, what="character", sep="\n")

    if(length(labels) != (ncol(getCUSTOM_data()) - 1)) stop("Number of labels does not fit number of samples.")

    return(labels)
})

observeEvent(input$custom_file, {
    if(!is.null(getCUSTOM_data())){
        shinyjs::show(selector = "#custom_nav li a[data-value='Data Table']")
        shinyjs::show(selector = "#custom_nav li a[data-value='Gene ID mapping']")
    }
})

output$custom_labels_uploaded <- reactive({
    return(!is.null(getCUSTOM_classLabels()))
})

outputOptions(output, 'custom_labels_uploaded', suspendWhenHidden=FALSE)


updateCUSTOM_classLabels <- observeEvent(input$custom_class_labels,{
    labels <- getCUSTOM_classLabels()
    for(subtype in brca_subtypes){
        updateSelectInput(session, paste("custom", "subtype", subtype, sep="_"),
                          subtype,
                          unique(labels))
    }
})


getCUSTOM_mappedClassLabels <- reactive({
    lumA <- input$custom_subtype_LumA
    lumB <- input$custom_subtype_LumB
    basal <- input$custom_subtype_Basal
    her2 <- input$custom_subtype_Her2

    gset <- getCUSTOM_data()

    labels <- getCUSTOM_classLabels()

    for(subtype in brca_subtypes){
        labels[which(labels == input[[paste("subtype", subtype, sep = "_")]])] <- subtype
    }

    labels[-which(labels %in% brca_subtypes)] <- "Other"

    return(labels)
})

getCUSTOM_data_exprs <- reactive({
    exprs <- as.matrix(getCUSTOM_data()[,-1])
    rownames(exprs) <- getCUSTOM_data()[,1]
    return(exprs)
})

getCUSTOM_data_genes <- reactive({
    as.character(getCUSTOM_data()[,1])
})

getCUSTOM_data_mapped_genes <- reactive({
    AnnotationDbi::mapIds(org.Hs.eg.db, getCUSTOM_data_genes(),
                          column = "ENTREZID",
                          keytype = input$custom_gene_id_type)
})


getCUSTOM_data_mapped_genes_table <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Mapping gene ids", value = 0.1)

    AnnotationDbi::select(org.Hs.eg.db, getCUSTOM_data_genes(),
                          column = "ENTREZID",
                          keytype = input$custom_gene_id_type)
})

startCUSTOM <- observeEvent(input$custom_button, {
    getCUSTOM_subtypes()

    show(selector = "#custom_nav li a[data-value='Predictions Plot']")
    show(selector = "#custom_nav li a[data-value='Predictions Table']")
    show(selector = "#custom_nav li a[data-value='Performance (subtypes)']")
    show(selector = "#custom_nav li a[data-value='Performance (predictors)']")
    show(selector = "#custom_nav li a[data-value='Confusion matrices']")
})

getCUSTOM_subtypes <- reactive({

    # Create a Progress object
    progress <- shiny::Progress$new()

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    # load series and platform data from GEO
    ex <- getCUSTOM_data_exprs()

    progress$set(message = "Processing uploaded data", value = 0.1)

    # transform array signal
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }

    # get gene ids
    entrez_ids <- getCUSTOM_data_mapped_genes()
    rownames(ex) <- entrez_ids

    # map gene ids. when several probes match the same gene use the one
    # that is most variable
    progress$set(message = "Processing uploaded data", value = 0.2,
                 detail = "Selecting most variable entry for each gene")

    ex <- fn_filter_most_variable(ex, unlist(entrez_ids))

    progress$set(message = "Performing pathway enrichment analysis", value = 0.2)

    result <- fn_ssgsea(ex,
                        progress = progress,
                        selected.pathways = selected.pathways,
                        all.rf = all.rf,
                        pathway.sources = pathway.sources)

    result <- spread(result, "pathways", "subtype")

    progress$set(message = "Predicting subtype labels", value = 0.9)

    gfu_result <- geneFu(t(ex))

    annotated_result <- getCUSTOM_mappedClassLabels()

    if(length(annotated_result) > 0){
        final_result <- cbind(result, gfu_result[result$Sample,],
                              CUSTOM = annotated_result)
    }else {
        final_result <- cbind(result, gfu_result[result$Sample,])
    }

    final_result <- as.data.frame(final_result)

    progress$set(message = "Returning results", value = 1.0)

    return(final_result)
})
