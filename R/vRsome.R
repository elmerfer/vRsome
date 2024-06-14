#' BuildVarsomeVCF
#' It buils the appropriate VFC file to feed Varsome Clinical
#' from an MODApy excel panel
#' @param xlsxFile the full file path
#' @return it creates a VCF with only those variants selected by NON EMPTY VARSOME excel column
#' @export
#' 
BuildVarsomeVCF <- function(xlsxFile){
  if(!file.exists(xlsxFile)){
    stop("file does not exists")
  }
  panel <- openxlsx::read.xlsx(xlsxFile,sheet=1)
  selected_variants <- which(!is.na(panel$VARSOME))
  if(length(selected_variants)<1){
    stop("candidate variants NOT selected")
  }
  ## prepare the VCF header
  
  panel$QUAL="."
  panel$FILTER="PASS"
  panel$INFO="."
  vcf.header <- c("##fileformat=VCFv4.2",
                  paste0("##fileDate=",Sys.time()),
                  "##source=vRsome",
                  paste0("##reference=",xlsxFile),
                  paste0(c("#CHROM", "POS", "ID","REF","ALT","QUAL","FILTER","INFO"),collapse = "\t"))
  
  vcf.header <- c(vcf.header,apply(panel[selected_variants,c("CHROM","POS","RSID","REF","ALT","QUAL","FILTER","INFO")],MARGIN=1, FUN=function(x) paste0(x,collapse = "\t")))
  writeLines(vcf.header,stringr::str_replace_all(xlsxFile,".xlsx","_vRsome.vcf"))
  if(file.exists(stringr::str_replace_all(xlsxFile,".xlsx","_vRsome.vcf"))){
    cat(paste0("\nSaved as: ",stringr::str_replace_all(xlsxFile,".xlsx","_vRsome.vcf")))
  }else{
    stop("failed")
  }
}



#' BuildVarsomeVCF_FULL
#' It builds the appropriate VFC file to feed Varsome Clinical
#' from an MODApy excel panel, FULL implies that all the varian information is used in the VCF
#' @param xlsxFile the full file path
#' @return it creates a VCF with only those variants selected by NON EMPTY VARSOME excel column
#' @export
#' 
BuildVarsomeVCF_FULL <- function(xlsxFile){
  if(!file.exists(xlsxFile)){
    stop("file does not exists")
  }
  panel <- openxlsx::read.xlsx(xlsxFile,sheet=1)
  selected_variants <- which(!is.na(panel$VARSOME))
  if(length(selected_variants)<1){
    stop("candidate variants NOT selected")
  }
  if(!any(stringr::str_detect(colnames(panel),"GENE_NAME")==TRUE)){
    stop("The Excel file PANEL should contain a columna with the gene symbola named GENE_NAME")
  }
  genes <- unique(panel$GENE_NAME[selected_variants])
  
  vcf_file <- list.files(dirname(xlsxFile), pattern = "final.vcf", full.names = T)
  patient_ID <- unlist(stringr::str_split(basename(xlsxFile),"_MODApy"))[1]
  vcf_file <- vcf_file[stringr::str_detect(vcf_file,patient_ID)]
  if(length(vcf_file)<1){
    stop(paste("VCF final NOT FOUND for patient ", patient_ID))
  }
  if(!file.exists(vcf_file)){
    stop(paste("File NOT FOUND ", vcf_file))
  }
  
  # header_vcf <- vcfR::read.vcfR(vcf_file,nrows = 1)
  header_vcf <- data.table::fread(cmd= paste("grep", "\\##", vcf_file),
                                  sep="",sep2="",header=F,data.table=F)[,,drop=T]
  lines_to_skip <- length(header_vcf)
  ##elimino la informacion de las lineas de comando
  header_vcf<-header_vcf[-which(stringr::str_detect(header_vcf,"Cmd"))]
  header_vcf <- header_vcf[-which(stringr::str_detect(header_vcf,"CommandLine"))]
  header_vcf <- header_vcf[-which(stringr::str_detect(header_vcf,"contig"))]
  
  vcf_col_names <- data.table::fread(input=vcf_file, skip=lines_to_skip,nrow=1)
  
  ##busco los genes seleccionados del panel
  
  variants_DT <- plyr::ldply(genes, function(gen){
    var <- data.table::fread(cmd=paste("grep","-e",gen, vcf_file))
  })
  colnames(variants_DT) <- colnames(vcf_col_names)
  
  ## busco solo las variantes seleccionadas basandome en la HGSV.C que siempre deberia estar
  vDT<- plyr::ldply(panel$HGVS.C[selected_variants],function(hgvs.c){
    variants_DT[stringr::str_detect(variants_DT$INFO,hgvs.c),]
  })
  
  if(nrow(vDT)==length(selected_variants)){
    print("EXITO")
  }
  ##creo el objeto para guardar
  vcf2save <- new("vcfR", meta = header_vcf, fix = as.matrix(vDT[,-c(ncol(vDT)-1,ncol(vDT))]), gt = as.matrix(vDT[,c(ncol(vDT)-1,ncol(vDT))]))
  file_to_save <- stringr::str_replace_all(xlsxFile,".xlsx","_vRsome_FULL.vcf.gz")
  write.vcf(x= vcf2save, file =file_to_save)
  gunzip(file_to_save)

}


vcf_header <- function(vcf_file){
  header_vcf <- data.table::fread(cmd= paste("grep", "\\##", vcf_file),
                                  sep="",sep2="",header=F,data.table=F)[,,drop=T]
}

#' CheckGeneOnPanels
#' Check that the gene in the panel are present in the hg19 and snpEFF reference gene list and database respectively
#' @param xlsx.file file path to the exel file containing the gene panel
#' @return a data frame with those genes to be revisited (not found in the databses )
CheckGeneOnPanels <- function(xlsx.file, verbose =TRUE){
  if(!is.null(xlsx.file)){
    hg19_genes <- readRDS("Data/hg19_genes.RDS")
    hg19_snpEFF_genes <- readRDS("Data/hg19_snpEFF_genes.RDS")
    
    wb <- openxlsx::loadWorkbook(xlsxFile = xlsx.file)
    panel <- openxlsx::readWorkbook(wb,skipEmptyRows = FALSE,
                                    skipEmptyCols = FALSE,detectDates=TRUE)
    if(all(c("Atencion.Cobertura","Atencion.VCF") %in% colnames(panel))==TRUE){
      n_cols<-which(colnames(panel)  %in% c("Atencion.Cobertura","Atencion.VCF"))
    }else{
      n_cols <- c(1:2)+ncol(panel)
    }
    
    estas<-panel$GeneSymbol %in% hg19_genes$GeneName
    estas_en_snpEFF <- panel$GeneSymbol %in% hg19_snpEFF_genes$geneName
    
    if(all(c(estas,estas_en_snpEFF))){
      cat("Todos los genes del panel han sido encontrados tanto en hg19 como en la base de datos de snpEFF")
    }else{
      panel$Atencion.Cobertura <- NA
      panel$Atencion.VCF <-NA
      panel$Atencion.Cobertura[!estas] <- "REVISAR"
      panel$Atencion.VCF[!estas_en_snpEFF] <- "REVISAR"
      panel$GeneSymbol[!estas_en_snpEFF]
     
    }
    if(verbose==TRUE){
      cat("\n---------------------------------------------\n")
      cat(paste0("The Panel is : ", basename(xlsx.file),"\n"))
      cat("\n---------------------------------------------\n")
      print(panel[!estas_en_snpEFF | !estas_en_snpEFF ,c("GeneSymbol","Atencion.Cobertura","Atencion.VCF")])
      cat("\n---------------------------------------------\n")      
    }
     invisible(return(panel[!estas_en_snpEFF | !estas_en_snpEFF ,c("GeneSymbol","Atencion.Cobertura","Atencion.VCF")]))
  }  
  return(invisible(NA))
  
}
