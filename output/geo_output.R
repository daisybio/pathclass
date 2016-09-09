output$gssea_result_geo <- renderDataTable({
    getGEO_subtypes()
})

output$gssea_result_geo_heatmap <- renderPlot({
    plot_heatmap(getGEO_subtypes(), input$geo_predictors)
})

output$class_error_subt_plot_geo <- renderPlot({
    plot_class_errors_subtypes(getGEO_subtypes(), input$geo_predictors)
})

output$class_error_pred_plot_geo <- renderPlot({
    plot_class_errors_predictors(getGEO_subtypes(), input$geo_predictors)
})

output$class_error_cm_plot_geo <- renderPlot({
    plot_confusion_matrices(getGEO_subtypes(), input$geo_predictors)
})

output$download_geo_subtypes <- downloadHandler(
    filename = function() { paste(input$geo_id, "_predicted_subtypes.txt", sep="") },
    content = function(file) {
        data <- getGEO_subtypes()
        write.table(data, file, row.names=F, sep="\t", quote=F)
    }
)
