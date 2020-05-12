### installation ###

#CRAN packages
install.packages(c("shiny", "shinydashboard", "shinyjs", "dplyr", "ggplot2", 
                   "gplots", "scales", "gridExtra", "lazyeval",
                   "randomForest", "foreach", "tidyr", "stringr", "RColorBrewer",
                   "iterators", "htmlwidgets", "matrixStats",
                   "networkD3", "devtools", "XML", "R.utils", "BiocManager"))

#bioconductor packages
BiocManager::install(c("Biobase", "GEOquery", "genefu", "org.Hs.eg.db", "GSVA"))

devtools::install_github("AnalytixWare/ShinySky")

