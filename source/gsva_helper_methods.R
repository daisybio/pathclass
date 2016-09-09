run.gsva <- function(expression, geneSets, method = "gsva", cores) {
  require(GSVA)
  require(parallel)
  gsva.object <- gsva(as.matrix(expression), geneSets, method = method, parallel.sz = cores)#$es.obs
  return(gsva.object)
}

### Features are rows ###
add.unmapped.nodes <- function(graph.nodes, data, method = "mean") {
  not.mapped <- setdiff(graph.nodes, rownames(data))
  na.matrix <- data.frame(matrix(NA, nrow = length(not.mapped), ncol=ncol(data)))
  colnames(na.matrix) <- colnames(data)
  rownames(na.matrix) <- not.mapped
  data.plus <- rbind(data, na.matrix)
  if (method == "mean") {  
    for (i in which(sapply(data.plus, is.numeric))) {
      data.plus[is.na(data.plus[, i]), i] <- mean(data.plus[, i],  na.rm = TRUE)
    }
  } else if (method == "median") {
    for (i in which(sapply(data.plus, is.numeric))) {
      data.plus[is.na(data.plus[, i]), i] <- median(data.plus[, i],  na.rm = TRUE)
    }
  }
  return(data.plus)
}
# t(apply(X, 1, function(x) (x-mean(x))/sd(x)))

## rows must be features to aggregate
aggregate.mean <- function(data, genesets, abs = TRUE, cores = 1) {
  require(parallel)
  require(foreach)
  require(doMC)
  
  if (abs) {
    data <- abs(data)
  }
  registerDoMC(cores)
  data.mean <- foreach(i = 1:length(genesets), .combine = rbind) %dopar% {
    gene.rows <- genesets[[i]]
    ret <- apply(data[gene.rows,], 2, function(x) mean)
    return(ret)
  }
  rownames(data.mean) <- names(genesets)
  colnames(data.mean) <- colnames(data)
  return(data.mean)
}

aggregate.median <- function(data, genesets, abs = TRUE, cores = 1) {
  require(parallel)
  require(foreach)
  require(doMC)
  
  if (abs) {
    data <- abs(data)
  }
  registerDoMC(cores)
  data.mean <- foreach(i = 1:length(genesets), .combine = rbind) %dopar% {
    gene.rows <- genesets[[i]]
    ret <- apply(data[gene.rows,], 2, function(x) median)
    return(ret)
  }
  rownames(data.mean) <- names(genesets)
  colnames(data.mean) <- colnames(data)
  return(data.mean)
}

aggregate.zscore <- function(data, genesets, abs = TRUE, cores = 1) {
  require(parallel)
  require(foreach)
  require(doMC)
  
  
  
  if (abs) {
    data <- abs(data)
  }
  registerDoMC(cores)
  data.mean <- foreach(i = 1:length(genesets), .combine = rbind) %dopar% {
    gene.rows <- genesets[[i]]
    ret <- apply(data[gene.rows,], 2, function(x) median)
    return(ret)
  }
  rownames(data.mean) <- names(genesets)
  colnames(data.mean) <- colnames(data)
  return(data.mean)
}