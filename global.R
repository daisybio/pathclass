pathway.sources <- c("KPM_HTRIdb", "KPM_HNET", "KPM_HPRD", "KPM_I2D", "CPDB", "MSIG")
all.sources <- c(pathway.sources, "pam50")
brca_subtypes <- c("LumA", "LumB", "Basal", "Her2")

# load about page
aboutFileName <- "www/about.html"
aboutText <- readChar(aboutFileName, file.info(aboutFileName)$size)

# load user guide page
guideFileName <- "www/guide.html"
guideText <- readChar(guideFileName, file.info(guideFileName)$size)
