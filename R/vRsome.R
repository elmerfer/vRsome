#' BuildVarsomeVCF
#' It buils the appropriate VFC file to feed Varsome Clinical
#' from an MODApy excel panel
#' @param xlsxFile the full file path
#' @return it creates a VCF with the same name as the original file 
#' @export
#' 
BuildVarsomeVCF <- function(xlsxFile){
  if(!file.exists(xlsxFile)){
    stop("file does not exists")
  }
  panel <- openxlsx::read.xlsx(xlsxFile,sheet=1)
  ## prepare the VCF header
  
  panel$QUAL="."
  panel$FILTER="PASS"
  panel$INFO="."
  vcf.header <- c("##fileformat=VCFv4.2",
                  paste0("##fileDate=",Sys.time()),
                  "##source=vRsome",
                  paste0("##reference=",xlsxFile),
                  paste0(c("#CHROM", "POS", "ID","REF","ALT","QUAL","FILTER","INFO"),collapse = "\t"))
  
  vcf.header <- c(vcf.header,apply(panel[,c("CHROM","POS","RSID","REF","ALT","QUAL","FILTER","INFO")],MARGIN=1, FUN=function(x) paste0(x,collapse = "\t")))
  writeLines(vcf.header,stringr::str_replace_all(xlsxFile,".xlsx","_vRsome.vcf"))
  if(file.exists(stringr::str_replace_all(xlsxFile,".xlsx","_vRsome.vcf"))){
    cat(paste0("\nSaved as: ",stringr::str_replace_all(xlsxFile,".xlsx","_vRsome.vcf")))
  }else{
    stop("failed")
  }
}

