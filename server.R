# Dependencies

library(shiny)
library(shinyjs)
library(GSVA)
library(randomForest)
library(foreach)
library(tidyr)
library(Biobase)
library(GEOquery)
library(ggplot2)
library(genefu)
library(matrixStats)
library(org.Hs.eg.db)
library(RColorBrewer)

source("source/gsva_helper_methods.R")
source("source/shiny_helper_methods.R")
source("source/plot_methods.R")

# Load the final models and gene expression dataset
load("data/final_rf_models_BRCA_TCGA_EXPR.RData")

# Shiny options. Allow more than the default file size for upload
options(shiny.maxRequestSize = 100*1024^2)

shinyServer(function(input, output, session) {

    # Reactive elements to better describe the pathways used as features in the
    # various predictors
    pathway_names <- reactive({
        names(selected.pathways[[input$selp_dataset]])
    })

    pathway_members <- reactive({
        genes <- selected.pathways[[input$selp_dataset]][[input$selp_pathway]]
        AnnotationDbi::select(org.Hs.eg.db,
                              keys = genes,
                              keytype = "ENTREZID",
                              columns = c("ENTREZID", "SYMBOL"))
    })

    pathway_numbers <- reactive({

        foreach(pathway_source = names(selected.pathways), .combine = rbind) %do%
        {
            foreach(pathway = names(selected.pathways[[pathway_source]]), .combine = rbind) %do% {
                data.frame(predictor = pathway_source, pathway = pathway,
                           number_of_genes = length(selected.pathways[[pathway_source]][[pathway]]))
            }
        }
    })

    observeEvent(input$selp_dataset, {
        updateSelectInput(session, "selp_pathway", choices = pathway_names())
    })

    output$selected_pathway_members <- renderDataTable({
        pathway_members()
    })

    output$pathways_overview <- renderPlot({
        ggplot(pathway_numbers(), aes(x = predictor, y = number_of_genes)) +
            geom_violin(fill = "skyblue") + theme_bw() + ylab("Number of genes")
    })

    # reactive elements and observers for the different functionalities
    source("source/geo.R", local = TRUE)
    source("source/tcga.R", local = TRUE)
    source("source/custom.R", local = TRUE)

    # output elements such as data tables and plots for the different functionalities
    source("output/geo_output.R", local = TRUE)
    source("output/tcga_output.R", local = TRUE)
    source("output/custom_output.R", local = TRUE)
})
