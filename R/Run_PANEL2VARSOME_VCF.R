source("~/Dropbox/Fleni/RStudio/libraries/vRsome/R/vRsome.R")
pat_path <- "/home/biomolecular/DATA/NGS/Pacientes/"


library(rChoiceDialogs)

library(stringr)
select <- rstudioapi::selectDirectory(
  caption = "Select Directory",
  label = "Select",
  path = "/home/biomolecular/DATA/NGS/RNAseq/Pacientes"
)

if(length(select)>0){
  if(stringr::str_detect(select,".xlsx")){
    stop("CUIDADO: DEBE SER UN ARCHIVO xlsx")
  }
  BuildVarsomeVCF(openxlsx::read.xlsx(select))

}else{
  cat("Not selected")
}
