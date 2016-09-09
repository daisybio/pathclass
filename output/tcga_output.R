output$class_error_subt_plot_tcga_array <- renderPlot({
    plot_class_errors_subtypes(getTCGA_subtypes(), input$tcga_predictors)
})

output$class_error_pred_plot_tcga_array <- renderPlot({
    plot_class_errors_predictors(getTCGA_subtypes(), input$tcga_predictors)
})

output$class_error_cm_plot_tcga_array <- renderPlot({
    plot_confusion_matrices(getTCGA_subtypes(), input$tcga_predictors)
})

output$gssea_result_tcga <- renderDataTable({
    getTCGA_subtypes()
})

output$gssea_result_tcga_heatmap <- renderPlot({
    plot_heatmap(getTCGA_subtypes(), input$tcga_predictors)
})

output$download_tcga_subtypes <- downloadHandler(
    filename = function() { "TCGA_predicted_subtypes.txt" },
    content = function(file) {
        data <- getTCGA_subtypes()
        write.table(data, file, row.names=F, sep="\t", quote=F)
    }
)
