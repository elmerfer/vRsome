library(data.table)
library(openxlsx)

hg19_genes <- readRDS("Data/hg19_genes.RDS")

xlsx.file <- rstudioapi::selectFile(, filter = ".xlsx")

if(!is.null(xlsx.file)){
  panel <- read.xlsx(xlsxFile = xlsx.file)
  estas<-panel$GeneSymbol %in% hg19_genes$GeneName
  if(all(estas)){
    cat(paste0("Todos los genes del panel han sido encontrados"))
  }else{
    panel$NotFound <- NA
    panel$NotFound[!estas] <- "REVISAR"
    write.xlsx(panel, file = xlsx.file)
    cat("los siguientes genes no fueron encontrados\n")
    cat(paste0(panel$GeneSymbol[!estas],"\n"))
    cat("En el archivo original se agrego una columna -NotFound- con la identificacion de REVISAR estos genes")
  }
}

