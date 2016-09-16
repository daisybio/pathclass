# A heatmap-like plot in which each subtype is represented by one color. 
# Predictors are on the x axis and samples on y
# Samples are sorted according to the selected reference
plot_heatmap <- function(subtypes, orderBy){
  subtypes <- subtypes %>% dplyr::select(-scmgene, -scmod1, -scmod2)
  subtypes <- transform(subtypes, Sample=reorder(Sample, as.integer(as.factor(subtypes[[orderBy]]))))
  plot.data <- tidyr::gather(subtypes, "predictor", "subtype", -Sample)
  ggplot(plot.data, aes(y = Sample, x = as.factor(predictor), fill = as.factor(subtype))) +
    geom_tile() +
    theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) +
    xlab("Predictor") + labs(fill = "Breast cancer subtype") +
    scale_fill_brewer(palette="Accent")
}

# This is a plot of F measure, precision and recall with focus on different subtypes
plot_class_errors_subtypes <- function(subtypes, reference){
  class_errors <- calculate_class_errors(subtypes, reference)
  class_errors <- tidyr::gather(class_errors, key = "measure", value = "value", -subtype, -predictor, -reference)
  
  ggplot(class_errors, aes(x = predictor, y = value, fill = predictor)) +
    facet_grid(measure~subtype) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

# This is a plot of F measure, precision and recall with focus on different predictors
# This is essentially the same plot as the one above but grouping results by predictor
plot_class_errors_predictors <- function(subtypes, reference){
  class_errors <- calculate_class_errors(subtypes, reference)
  class_errors <- tidyr::gather(class_errors, key = "measure", value = "value", -subtype, -predictor, -reference)
  
  ggplot(class_errors, aes(x = subtype, y = value, fill = subtype)) +
    facet_grid(measure~predictor) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

# This is the ggplot2 version of a confusion matrix
plot_confusion_matrices <- function(subtypes, reference){
  plot.data <- calculate_confusion_matrices(subtypes, reference)
  
  ggplot(plot.data, aes(x = Predicted, y = Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "orange") +
    facet_wrap(~predictor) +
    theme_bw()
}