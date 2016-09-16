# when we have several probes per gene we use the most variable one
fn_filter_most_variable <- function(ex, entrez_ids){
    row.variance <- rowVars(ex)
    names(row.variance) <- rownames(ex)

    most.variable <- foreach(id = unique(entrez_ids), .combine = c) %do% {
        probes.to.test <- row.variance[id]
        names(probes.to.test)[which.max(probes.to.test)]
    }

    ex[most.variable,]
}

# we calculate the prediction error using a reference
# includes F measure, precision and recall
calculate_class_errors <- function(subtypes, reference){
    subtypes <- subtypes %>% dplyr::select(-scmgene, -scmod1, -scmod2)

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

# Confusion matrix
calculate_confusion_matrices <- function(subtypes, reference){
    subtypes <- subtypes %>% dplyr::select(-scmgene, -scmod1, -scmod2)

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

# Uses the bioconductor geneFu package to get predictions of classical gene expression signatures
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

# Performs for each pathway predictor selected single sample gene set enrichment analysis
# The enrichment scores are used for predicting subtypes in random forest, i.e. the features.
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
        ssgsea_result <- run.gsva(exprs.plus, pathways, cores = max(parallel::detectCores() - 1, 1), method = "ssgsea")

        # We get the random forest model we built for MSIG features
        # and predict the sample label
        rf.obj <- all.rf[[pathway.source]]
        data.frame(Sample = colnames(exprs.plus),
                   pathways = pathway.source,
                   subtype = as.character(predict(rf.obj, newdata = t(ssgsea_result)))
        )
    }

}
