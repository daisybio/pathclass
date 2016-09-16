plot_heatmap <- function(subtypes, orderBy){
    subtypes <- subtypes %>% dplyr::select(-scmgene, -scmod1, -scmod2, -intClust)
    subtypes <- transform(subtypes, Sample=reorder(Sample, as.integer(as.factor(subtypes[[orderBy]]))))
    plot.data <- tidyr::gather(subtypes, "predictor", "subtype", -Sample)
    ggplot(plot.data, aes(y = Sample, x = as.factor(predictor), fill = as.factor(subtype))) +
        geom_tile() +
        theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) +
        xlab("Predictor") +
        scale_fill_discrete(name = "Breast cancer subtype") +
        scale_fill_brewer(palette="Accent")
}

fn_filter_most_variable <- function(ex, entrez_ids){
    row.variance <- rowVars(ex)
    names(row.variance) <- rownames(ex)

    most.variable <- foreach(id = unique(entrez_ids), .combine = c) %do% {
        probes.to.test <- row.variance[id]
        names(probes.to.test)[which.max(probes.to.test)]
    }

    ex[most.variable,]
}

calculate_class_errors <- function(subtypes, reference){
    subtypes <- subtypes %>% dplyr::select(-scmgene, -scmod1, -scmod2, -intClust)

    subtype_reference <- subtypes[[reference]]
    subtype_reference[!(subtype_reference %in% brca_subtypes)] <- NA
    subtype_reference <- factor(as.character(subtype_reference), levels = brca_subtypes)

    foreach(predictor = colnames(subtypes), .combine = rbind) %do% {

        if(predictor %in% c(reference, "Sample")) return(NULL)
        #browser()
        predicted <- subtypes[[predictor]]
        predicted[!(predicted %in% brca_subtypes)] <- NA
        predicted <- factor(as.character(predicted), levels = brca_subtypes)

        # create the confusion matrix
        cm = as.matrix(table(Reference = subtype_reference,
                             Predicted = predicted))

        n = sum(cm) # number of instances
        nc = nrow(cm) # number of classes
        diag = diag(cm) # number of correctly classified instances per class
        rowsums = apply(cm, 1, sum) # number of instances per class
        colsums = apply(cm, 2, sum) # number of predictions per class
        p = rowsums / n # distribution of instances over the actual classes
        q = colsums / n # distribution of instances over the predicted classes

        accuracy = sum(diag) / n
        precision = diag / (colsums + 1e-64)
        recall = diag / (rowsums + 1e-64)
        f1 = 2 * precision * recall / (precision + recall + 1e-64)

        data.frame(predictor = predictor, reference = reference, precision, recall, `F-measure` = f1) %>%
            tibble::rownames_to_column("subtype")
    }
}

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

calculate_confusion_matrices <- function(subtypes, reference){
    subtypes <- subtypes %>% dplyr::select(-scmgene, -scmod1, -scmod2, -intClust)

    subtype_reference <- subtypes[[reference]]
    subtype_reference[!(subtype_reference %in% brca_subtypes)] <- NA
    subtype_reference <- factor(as.character(subtype_reference), levels = rev(brca_subtypes))

    foreach(predictor = colnames(subtypes), .combine = rbind) %do% {

        if(predictor %in% c(reference, "Sample")) return(NULL)
        #browser()
        predicted <- subtypes[[predictor]]
        predicted[!(predicted %in% brca_subtypes)] <- NA
        predicted <- factor(as.character(predicted), levels = brca_subtypes)

        if(predictor %in% c(reference, "Sample")) return(NULL)

        # create the confusion matrix
        cm = as.matrix(table(Reference = subtype_reference,
                             Predicted = predicted))

        cbind(predictor = predictor, reference = reference, as.data.frame(cm))
    }
}

plot_confusion_matrices <- function(subtypes, reference){
    plot.data <- calculate_confusion_matrices(subtypes, reference)

    ggplot(plot.data, aes(x = Predicted, y = Reference, fill = Freq)) +
        geom_tile() +
        geom_text(aes(label = Freq), color = "orange") +
        facet_wrap(~predictor) +
        theme_bw()
}

geneFu <- function(expression_data){
    sbt.models = c("scmgene", "scmod1", "scmod2",
                   "pam50", "ssp2006", "ssp2003")
    annot <- data.frame(EntrezGene.ID = as.integer(colnames(expression_data)))
    rownames(annot) <- colnames(expression_data)
    annot$probe <- colnames(expression_data)

    foreach(sbt.model = sbt.models,
            .combine = cbind,
            .final = function(x){
                colnames(x) <- sbt.models
                rownames(x) <- rownames(expression_data)
                return(x)
            }) %do%{

                as.character(molecular.subtyping(sbt.model = sbt.model,
                                                 data = expression_data,
                                                 annot = annot, do.mapping = TRUE)$subtype)
            }
}

fn_ssgsea <- function(exprs, progress, pathway.sources, selected.pathways, all.rf, sel.samples = NULL){

    if(is.null(sel.samples)) sel.samples <- rownames(exprs)
    current_progress = 0.1

    foreach(pathway.source = pathway.sources, .combine = rbind) %do% {

        current_progress <- current_progress + 0.1

        progress$set(message = "Performing pathway enrichment analysis",
                     value = current_progress,
                     detail = paste("for", pathway.source))

        pathways <- selected.pathways[[pathway.source]]

        # For ssGSEA, we require all genes in the pathways to have
        # expression values. Since this is not the case, we create
        # "default" values for the missing genes in the expression
        # in this case we set them to the mean of the gene expression
        # for that sample.

        unique.ids <- unique(unlist(pathways))

        not.in.expr <- unique.ids[!unique.ids %in% rownames(exprs)]
        df.values <- matrix(rep(colMeans(exprs, na.rm=T), length(not.in.expr)), ncol = ncol(exprs))
        rownames(df.values) <- not.in.expr

        exprs.plus <- rbind(exprs, df.values)

        # We create the pathway scores for the sample
        # (this can also be done for multiple samples in parallel)
        ssgsea_result <- run.gsva(exprs.plus, pathways, cores = 2, method = "ssgsea")

        # We get the random forest model we built for MSIG features
        # and predict the sample label
        rf.obj <- all.rf[[pathway.source]]
        data.frame(Sample = colnames(exprs.plus),
                   pathways = pathway.source,
                   subtype = as.character(predict(rf.obj, newdata = t(ssgsea_result)))
        )
    }

}
