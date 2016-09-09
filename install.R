### installation ###

#CRAN packages
install.packages(c("shiny", "shinyjs", "shinysky", "dplyr", "ggplot2", 
                   "gplots", "scales", "gridExtra", "lazyeval",
                   "randomForest", "foreach", "tidyr", "stringr", "RColorBrewer",
                   "iterators", "htmlwidgets", "matrixStats",
                   "networkD3", "devtools", "XML", "R.utils"))

#bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite("Biobase", ask=F)
biocLite("GEOquery", ask=F)
biocLite("genefu", ask=F)
biocLite("org.Hs.eg.db", ask=F)
biocLite("GSVA", ask=F)
