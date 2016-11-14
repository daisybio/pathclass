# Collect class labels for TCGA samples
getTCGA_subtypes <- reactive({

    # Create a Progress object
    progress <- shiny::Progress$new(style="old")

    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Processing TCGA data", value = 0.1)

    result <- foreach(pathway.source = pathway.sources,
                      .combine = cbind,
                      .final = function(x){
                          colnames(x) <- pathway.sources
                          rownames(x) <- names(all.rf$KPM_HTRIdb$predicted)
                          return(x)
                      }) %do%{
                          as.character(all.rf[[pathway.source]]$predicted)
                      }
    gfu_result <- geneFu(brca.gene.expr)

    final_result <- cbind(result, gfu_result[rownames(result),])
    final_result <- as.data.frame(final_result) %>%
        tibble::rownames_to_column("Sample")
    return(final_result)
})
