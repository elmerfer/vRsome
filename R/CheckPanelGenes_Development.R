rm(list=ls())
library(data.table)
library(openxlsx)
source("R/vRsome.R", echo=FALSE)



xlsx.file <- rstudioapi::selectFile(path = "~/DATA/NGS/Paneles", filter = ".xlsx")
print(xlsx.file)
if(!is.null(xlsx.file)){

  CheckGeneOnPanels(xlsx.file)

}

