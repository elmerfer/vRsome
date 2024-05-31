source("R/vRsome.R")
pat_path <- "/home/biomolecular/DATA/NGS/Pacientes/"




library(stringr)
select <- rstudioapi::selectFile(
  caption = "Select Directory",
  filter = "Excel Files (*.xlsx)",
  path = pat_path
)

if(length(select)>0){
  if(stringr::str_detect(select,".xlsx")==FALSE){
    stop("CUIDADO: DEBE SER UN ARCHIVO xlsx")
  }
  BuildVarsomeVCF(select)

}else{
  cat("Not selected")
}
