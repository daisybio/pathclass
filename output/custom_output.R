output$raw_data_id_map <- renderDataTable({
    getCUSTOM_data_mapped_genes_table()
})

output$raw_data_custom <- renderDataTable({
    getCUSTOM_data()
})

output$gssea_result_custom <- renderDataTable({
    getCUSTOM_subtypes()
})

output$gssea_result_custom_heatmap <- renderPlot({
    plot_heatmap(getCUSTOM_subtypes(), input$custom_predictors)
})

output$class_error_subt_plot_custom <- renderPlot({
    plot_class_errors_subtypes(getCUSTOM_subtypes(), input$custom_predictors)
})

output$class_error_pred_plot_custom <- renderPlot({
    plot_class_errors_predictors(getCUSTOM_subtypes(), input$custom_predictors)
})

output$class_error_cm_plot_custom <- renderPlot({
    plot_confusion_matrices(getCUSTOM_subtypes(), input$custom_predictors)
})

output$download_custom_subtypes <- downloadHandler(
    filename = function() { "predicted_subtypes.txt" },
    content = function(file) {
        data <- getCUSTOM_subtypes()
        write.table(data, file, row.names=F, sep="\t", quote=F)
    }
)

observeEvent(input$custom_class_labels, {
    updateSelectInput(session, "custom_predictors",
                      choices = c(all.sources, "CUSTOM"),
                      selected = "CUSTOM")
})
