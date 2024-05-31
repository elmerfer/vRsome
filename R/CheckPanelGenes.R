library(data.table)
library(openxlsx)
## Actualizado : 18/04/24

hg19_genes <- readRDS("Data/hg19_genes.RDS")
hg19_snpEFF_genes <- readRDS("Data/hg19_snpEFF_genes.RDS")

xlsx.file <- rstudioapi::selectFile(path = "~/DATA/NGS/Paneles", filter = ".xlsx")
print(xlsx.file)
if(!is.null(xlsx.file)){
  
  wb <- openxlsx::loadWorkbook(file = xlsx.file)
  panel <- openxlsx::read.xlsx(wb, detectDates = TRUE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  
  num_cols <- ncol(panel)
  
  col_id <- which(stringr::str_detect(colnames(panel), "VCF|Cobertura"))
  
  estas<-panel$GeneSymbol %in% hg19_genes$GeneName
  estas_en_snpEFF <- panel$GeneSymbol %in% hg19_snpEFF_genes$geneName
  
  if(all(c(estas,estas_en_snpEFF))){
    cat("Todos los genes del panel han sido encontrados tanto en hg19 como en la base de datos de snpEFF")
  }else{
    
    panel$NotInHg19 <- NA
    panel$NotInHg19snpEff <-NA
    
    panel$NotInHg19[!estas] <- "REVISAR"
    panel$NotInHg19snpEff[!estas_en_snpEFF] <- "REVISAR"
    
    if(length(col_id)<1){
      cat("Agrega las nuevas columnas\n")
      openxlsx::writeData(wb, sheet =1, startCol = num_cols+1, startRow = 1, x= c("Atencion VCF", panel$NotInHg19))
      openxlsx::writeData(wb, sheet =1, startCol = num_cols+2, startRow = 1, x= c("Atencion Cobertura", panel$NotInHg19snpEff))
    }else{
      cat("Reescribe las nuevas columnas\n")
      openxlsx::writeData(wb, sheet =1, startCol = col_id[1], startRow = 2, x= panel$NotInHg19)
      openxlsx::writeData(wb, sheet =1, startCol = col_id[2], startRow = 2, x= panel$NotInHg19snpEff)
    }
    
    openxlsx::saveWorkbook(wb, xlsx.file, overwrite = TRUE)
    if(all(estas)==FALSE){
      cat("los siguientes genes NO fueron encontrados en hg19\n")
      cat(paste0(panel$GeneSymbol[!estas],"\n"))  
    }
    if(all(estas_en_snpEFF)==FALSE){
      cat("los siguientes genes NO fueron encontrados en hg19 snpEFF\n")
      cat(paste0(panel$GeneSymbol[!estas_en_snpEFF],"\n"))
    }
  
    cat("En el archivo original se agrego una columna -Atencion VCF y Atencion Cobertura- con la identificacion de REVISAR estos genes")
    browseURL(xlsx.file)
    
  }
}

