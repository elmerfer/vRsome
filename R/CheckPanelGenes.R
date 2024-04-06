library(data.table)
library(openxlsx)

hg19_genes <- readRDS("Data/hg19_genes.RDS")
hg19_snpEFF_genes <- readRDS("Data/hg19_snpEFF_genes.RDS")

xlsx.file <- rstudioapi::selectFile(path = "/home/elmer/Elmer/FLENI/Paneles/", filter = ".xlsx")
print(xlsx.file)
if(!is.null(xlsx.file)){
  panel <- read.xlsx(xlsxFile = xlsx.file)
  estas<-panel$GeneSymbol %in% hg19_genes$GeneName
  estas_en_snpEFF <- panel$GeneSymbol %in% hg19_snpEFF_genes$geneName
  
  if(all(c(estas,estas_en_snpEFF))){
    cat("Todos los genes del panel han sido encontrados tanto en hg19 como en la base de datos de snpEFF")
  }else{
    panel$NotInHg19 <- NA
    panel$NotInHg19snpEff <-NA
    panel$NotInHg19[!estas] <- "REVISAR"
    panel$NotInHg19snpEff[!estas_en_snpEFF] <- "REVISAR"
    write.xlsx(panel, file = xlsx.file)
    if(!all(estas)){
      cat("los siguientes genes NO fueron encontrados en hg19\n")
      cat(paste0(panel$GeneSymbol[!estas],"\n"))  
    }
    if(!all(estas_en_snpEFF)){
      cat("los siguientes genes NO fueron encontrados en hg19 snpEFF\n")
      cat(paste0(panel$GeneSymbol[!estas],"\n"))
    }
  
    cat("En el archivo original se agrego una columna -NotInHg19 and NotInHg19snpEFF- con la identificacion de REVISAR estos genes")
    browseURL(xlsx.file)
    
  }
}

