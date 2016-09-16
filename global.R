pathway.sources <- c("KPM_HTRIdb", "KPM_HNET", "KPM_HPRD", "KPM_I2D", "CPDB", "MSIG")
all.sources <- c(pathway.sources, c("scmgene", "scmod1", "scmod2","pam50", "ssp2006", "ssp2003"))
brca_subtypes <- c("LumA", "LumB", "Basal", "Her2")

# load about page
aboutFileName <- "www/about.html"
aboutText <- readChar(aboutFileName, file.info(aboutFileName)$size)
