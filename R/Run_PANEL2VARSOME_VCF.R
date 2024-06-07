rm(list = ls())
source("R/vRsome.R")
if(Sys.info()["nodename"] == "elmer-pc"){
  pat_path <- "/home/elmer/Elmer/FLENI/EXOMAS/"
}else{
  pat_path <- "/home/biomolecular/DATA/NGS/Pacientes/"  
}





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
  # BuildVarsomeVCF(select)
  BuildVarsomeVCF_FULL(select)

}else{
  cat("Not selected")
}
