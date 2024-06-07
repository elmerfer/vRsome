##Author: Elmer A. Fernández
##last update: 22/04/2024
# Changes : trim galore, BreathCoverage
#This file is stored at source('~/Dropbox/Fleni/R_Pipelines/NGS.Utils.R')
##Funcionalidades utiles
library(stringr)##manipulacion de strings
library(tools)
library(tidyxl)##manipulacion de archivos excel
library(crayon)## colores de letras
library(Rsubread)
library(openxlsx)
library(Rsamtools)
library(BiocParallel)
library(parallel)

##
library(rstudioapi)

#Configuration stuff ####
warn <- magenta $ underline $ bold
error <- red $ bold
msg <- blue $ bold
#------------------------------------------------------------------
LoadOrderFromXlsx <- function(path){
  #It reads the excel file sent by macrogen with the order information
  #Args:
  #     path  : the path to the excel file. 
  #Returns :
  #     a data frame with columns HTTP, Patient, MD5. HTTP: the URL address for downloading data, 
  #           Patient: a character with the patient identifier and MD5 number for file checking     
  
  extension <- rev(unlist(str_split(basename(path),"\\.")))[1]
  if( !str_detect(extension, "xlsx")) {
    stop(red$bold("Please convert the file in a XLSX format"))
  }
  ex <- xlsx_cells( path )
  num.of.orders <- max( subset( ex, sheet != "Download_Address" )$row )-1
  
  hyperlinks <- do.call(rbind,lapply(subset(ex, sheet == "Download_Address")$formula, 
                                     function(x) {
                                       hp <- str_detect(x, "HYPERLINK")
                                       hp <- if(is.na(hp)) {FALSE} else hp
                                       tar <- str_detect(x, ".tar")
                                       tar <- if(is.na(tar)) {FALSE} else tar
                                       
                                       if( hp ) {
                                         if( tar ) {
                                           hyp <- unlist(str_split(str_replace_all(unlist(str_split(str_replace_all(x, "\"", ""), "\\("))[2],"[)]",""),","))[1]
                                           return(c(hyp, unlist(str_split(basename(hyp),"\\."))[1]))
                                           
                                         }else return(c(NA,NA))
                                       }else return(c(NA,NA))
                                       
                                     }))
  
  
  hyperlinks <- as.data.frame(na.omit(hyperlinks))
  colnames(hyperlinks) <- c("HTTP","Patient")
  ## Buscar la posicion de md5sum
  md5 <- unlist(as.data.frame(subset(ex, sheet == "Download_Address" & character == "md5sum")[,c("row","col")]))
  rows <- md5[1]+c(1,nrow(hyperlinks))
  hyperlinks$MD5 <- subset(ex, sheet == "Download_Address" & row %in% as.character(rows) & col == as.character(md5[2]) )$character
  
  return(hyperlinks)
  
}
#--------------------------------------------------
DownLoadMacrogenFilesFromOrder <- function(path, patientDirectory){
  hyp.urls <- LoadOrderFromXlsx(path)
  for( f in nrow(hyp.urls)){
    DownLoadMacrogenFiles( as.character(hyp.urls$HTTP[ f ]) , patientDirectory, as.character(hyp.urls$MD5[ f ]))
  }
}

#--------------------------------------------------
DownLoadMacrogenFiles <- function(url, patientDirectory, md5 = NULL){
  patient <- unlist(str_split( basename(url) , "\\."))[1]
  patient.dir <- paste(patientDirectory,patient, sep="")
  if(dir.exists(patient.dir)) {
    stop(red$bold("ERROR: Patient already downloaded"))
  }
  dir.create(patient.dir)
  if (dir.exists(patient.dir)) {
    cat(paste("\nCreated Directory :", patient.dir,"\n" ))
  }
  
  cat(blue(paste("\nStarting download...")))
  destination.file <- paste(patient.dir,"/",basename(url),sep="") ## esto esta mal
  cat(paste("Saving file :", destination.file," ... \n" ))
  comm <- paste("-P ",patient.dir,"/ ", url, sep="")
  cat("\n Command passed to console\n wget ")
  cat(comm)
  system2(command = "wget", args = comm, wait = FALSE, stdout = FALSE)
  #if( download.file(url, destfile = destination.file, method = "wget", extra = "-c -o /reporte.txt ") !=0) stop(red$bold("FAILED download"))
  # if( download.file(url, destfile = destination.file, method = "wget") !=0) stop(red$bold("FAILED download"))
  # 
  # if( is.null( md5 ) ) {
  #   cat(red$bold(paste("\nPlease Check MD5SUM :" , md5sum(destination.file), "\n")))
  # }else {
  #   Md5SumCheckAndMsg(fileName = destination.file, md5Num = md5, 
  #                     stopMsg = paste("Uncompres",destination.file,  "file Process aborted"))
  # }
  # 
  # UncompressMacrogenTarFile( patientName = patient, directory = patientDirectory)  
}
#------------------------------------------------------------------------------------------
UncompressMacrogenTarFile <- function(patientName, directory = "/mnt/A8FCD3E3FCD3AA38/NGS/Pacientes/"){
# UncompresMacrogenFile: Each macrogen work provides a file named XXXX.tar. This file contains all the 
# data (raw sequencing data as well as results data)
#This function helps you to decompress requiered files for Variant Analysis
# Args :
#   patientName : character given the patient name. This is found in the "Download Address" sheet
#                 of the provided excel file from Macrogen
#   diretcory : character providing the rooot patient directory. 
#               The full patient directory is given as directory + patientname 
# example :
#     UncompressMAcrogenTarFile("PUFRA")
  if(!str_detect(patientName,".tar")) {
    if( !dir.exists(paste(directory,patientName,sep="") ) ) {
        stop(paste(directory,patientName,sep="") %+% red$bold("Directory do not exists"))
    }
    setwd(paste(directory,patientName,sep=""))#assume directory existence.
    tar.file.name<- paste(patientName,".tar",sep="")
  }else {
    if( !dir.exists(paste(directory,unlist(str_split(patientName,"\\."))[1],sep="") ) ){
        stop(red$bold("Directory do not exists"))
    }
    setwd(paste(directory,unlist(str_split(patientName,"\\."))[1],sep=""))
    tar.file.name <- patientName
  }
  
  if(!file.exists(tar.file.name)) stop(paste(getwd(),tar.file.name, "-----------> FILE not exists"))
  cat(paste("\n Uncompressing ", tar.file.name, " \n"))
  system(paste("tar -xvf ", tar.file.name))
  ## aca descomprimio el .tar
  
  
  
  ##
  fastq1 <- str_replace(tar.file.name,".tar","_1.fastq.gz")
  if (!file.exists(fastq1)) {
    stop(paste(fastq1,"---> FILE not exists"))
  } 
  cat(paste("\n Uncompressing ", fastq1, " \n"))
  system(paste("gunzip ", fastq1, sep=""))
  if ( !file.exists(str_replace(fastq1, "fastq.gz","fastq")) ) {
    stop(paste(fastq1,"---> Uncompression Failed"))
  }
    
  fastq2 <- str_replace(tar.file.name,".tar","_2.fastq.gz")
  
  if (!file.exists(fastq2)) {
    stop(paste(fastq2,"---> FILE not exists"))
  } 
  cat(paste("\n Uncompressing ", fastq2, " \n"))
  system(paste("gunzip ", fastq2, sep="" ))
  if (!file.exists(str_replace(fastq2, "fastq.gz","fastq")) ) {
    stop(paste(fastq2,"---> Uncompression Failed"))
  }
    
  cat(paste("\n## Uncompressed files at :" , getwd(),"\n"))
  cat(paste(list.files(), collapse = "\n"))
  
  ## Verificamos los md5
  md5.info <- read.table(paste(patientName,".md5sum.file.txt",sep=""), colClasses = "character")
  if (nrow(md5.info) < 1) stop(red$bold("Faile to open md5 information"))
  colnames(md5.info) <- c("Md5","Files")
  
    md5verify <- md5sum(md5.info$Files)
    if ( all(md5verify == md5.info$Md5, na.rm = TRUE)) {
      cat(blue("\n md5sum OK \n"))
    } else {
      cat(red$bold("\n md5sum FAILED \n"))
      id.false <- which(md5verify != md5.info$Md5)
      for( j in id.false){
        cat( md5verify[j] %+% " : " %+% red$bold("FAILED") )
      }
   }
}


## Auxiliary functions

Md5SumCheckAndMsg <- function(fileName, md5Num, stopMsg){
  #check for md5sum and report error mesgs
  #if the md5 check fails, the program stops
  #Args :
  #     fileName  : a character or string with the file name to be checked
  #     md5Num    : the md5sum check
  #     stopMsg   : the message to be displayed on fail.
  # Returns:
  #       the program stop on failure or display ok messages.
  if( md5sum(fileName) == md5Num ) {
    cat( blue$bold(fileName) %+% " : " %+% md5Num %+% " : " %+% green$bold("OK"))
  }else {
    cat( blue$bold(fileName) %+% " : " %+% md5Num %+% " : " %+% red$bold("FAILED"))
    stop(red$bold(stopMsg))
  }
}
# General puporpose functions ------------------
# 

BuildPatientFiles <- function(mainPath, subj, gzip = TRUE){
  sbjDir <- file.path(mainPath,subj)
  if(gzip){
    files <- c(file1 = file.path(mainPath,subj,paste(subj,"_1.fastq.gz",sep="")), 
               file2= file.path(mainPath,subj,paste(subj,"_2.fastq.gz",sep="")))  
  }else{
    files <- c(file1 = file.path(mainPath,subj,paste(subj,"_1.fastq",sep="")), 
               file2= file.path(mainPath,subj,paste(subj,"_2.fastq",sep="")))
  }
  
  
  all.fastq.files <- list.files(sbjDir, full.names = TRUE)
  if(gzip){
    all.fastq.files <- all.fastq.files[str_detect(all.fastq.files,"fastq.gz")]
    
  }else{
    all.fastq.files <- all.fastq.files[str_detect(all.fastq.files,"fastq")]
  }

  if(all(file.exists( files ))){
    cat("\nfiles detected\n")
    return(list(subjectDir=sbjDir,files=files, fastq.files = all.fastq.files))
  }else{
    warn <- magenta $ underline $ bold
    error <- red $ bold
    cat(error("\n...................................................\n"))
    cat(warn(paste("Subject:", subj, "Not Found")))
    cat("\n")
    cat(warn("The following files where not found :\n"))
    cat(warn(paste(files,"\n")))
    cat(error("\n...................................................\n"))
    cat(error("Please check the file names in Patient directory\n"))
    cat(error("\n...................................................\n"))
    return(NA)  
  }
  
}

######### TEST 
BuildPatientFiles2 <- function(mainPath, subj){
  sbjDir <- file.path(mainPath,subj)
  
  if(!dir.exists(sbjDir)) stop(error("\nsubject not found\n"))
  
  all.fastq.files <- list.files(sbjDir, full.names = TRUE)
  ##
  fastq <- c(which(str_detect(all.fastq.files,"fastq")),
             which(str_detect(all.fastq.files,"fq")) )
  
  #check macrogen WES aligned files. "recal.bam
  macrogen.recal.bam.bai <- all.fastq.files[str_detect(all.fastq.files, "recal.bam.bai")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files, "recal.bam.bai")]
  macrogen.recal.bam <- all.fastq.files[str_detect(all.fastq.files, "recal.bam")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files, "recal.bam")]
  
  #Check for Rsubread allign DNA - WES
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files, "ralign_PE.bam.summary")]
  ralign.indel.vcf <- all.fastq.files[str_detect(all.fastq.files, "ralign_PE.bam.indel.vcf")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files, "ralign_PE.bam.indel.vcf")]
  
  ralign.breakpoints.vcf <- all.fastq.files[str_detect(all.fastq.files, "ralign_PE.bam.breakpoints.vcf")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files, "ralign_PE.bam.breakpoints.vcf")]
  
  ralign.bam.bai <- all.fastq.files[str_detect(all.fastq.files, "ralign_PE.bam.bai")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files, "ralign_PE.bam.bai")]
  
  ralign.bam <- all.fastq.files[str_detect(all.fastq.files, "ralign_PE.bam")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files, "ralign_PE.bam")]
  
  all.fastq.files <- all.fastq.files[fastq]
  ## remove summry files
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files,".txt")]
  ##check for aligned subjunc files (RNA)
  subjunc.bam.summary <- all.fastq.files[str_detect(all.fastq.files,".subjunc.BAM.summary")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files,".subjunc.BAM.summary")]
  
  subjunc.bam.indels.vcf <- all.fastq.files[str_detect(all.fastq.files,".subjunc.BAM.indel.vcf")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files,".subjunc.BAM.indel.vcf")]
  
  subjunc.bam.junction.bed <- all.fastq.files[str_detect(all.fastq.files,".subjunc.BAM.junction.bed")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files,".subjunc.BAM.junction.bed")]
  
  subjunc.bam.breackpoints.vcf <- all.fastq.files[str_detect(all.fastq.files,".subjunc.BAM.breakpoints.vcf")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files,".subjunc.BAM.breakpoints.vcf")]

  subjunc.bam <- all.fastq.files[str_detect(all.fastq.files,".subjunc.BAM")]
  all.fastq.files <- all.fastq.files[!str_detect(all.fastq.files,".subjunc.BAM")]
  
  #Check for Rsubread allign DNA - WES
  
  #check macrogen WES aligned files. "recal.bam
  
  
  ## remove trimming_report
  # cat(msg("\n verifiyin the existence of the following files:\n"))
  # cat(msg(paste(" * ",all.fastq.files,"\n",sep="")))
  # 
  if(all(file.exists( all.fastq.files ))){
    cat("\nfiles detected\n")
    return(list(subjectDir=sbjDir, subject = subj, 
                fastq.files = all.fastq.files, 
                fastq.files.raw = all.fastq.files[str_detect(all.fastq.files,".fastq")],
                fastq.trimmed.files = all.fastq.files[str_detect(all.fastq.files,"val")],
                subjunc.bam = subjunc.bam,
                subjunc.bam.breackpoints.vcf = subjunc.bam.breackpoints.vcf,
                subjunc.bam.junction.bed = subjunc.bam.junction.bed,
                subjunc.bam.indels.vcf = subjunc.bam.indels.vcf,
                subjunc.bam.summary = subjunc.bam.summary,
                ralign.bam = ralign.bam,
                ralign.breakpoints.vcf = ralign.breakpoints.vcf,
                ralign.indel.vcf = ralign.indel.vcf,
                macrogen.recal.bam = macrogen.recal.bam,
                macrogen.recal.bam.bai = macrogen.recal.bam.bai))
  }else{

    cat(error("\n...................................................\n"))
    cat(warn(paste("Subject:", subj, "Not Found")))
    cat("\n")
    cat(warn("The following files where not found :\n"))
    cat(warn(paste(files,"\n")))
    cat(error("\n...................................................\n"))
    cat(error("Please check the file names in Patient directory\n"))
    cat(error("\n...................................................\n"))
    return(NA)  
  }
  
}

######
# > fastq.files <- list.files(sbj$subjectDir)
# > fastq.files <- fastq.files[str_detect(fastq.files,"fastq.gz")]
# > fastq.files
#######################################
### Control de Calidad ####
#######################################
RUN_FastQC <- function(subjPath, subj, fastQcPath = fastqc.path, allFiles = FALSE , gzip = TRUE){
  sbj <- BuildPatientFiles(subjPath,subj, gzip = gzip)
  
  if(any(is.na(sbj))) {
    stop("Error")
  }
  
  if(!dir.exists(file.path(sbj$subjectDir,"FASTQC"))){
    dir.create(file.path(sbj$subjectDir,"FASTQC"))  
    if(dir.exists(file.path(sbj$subjectDir,"FASTQC"))) {
      cat(paste("FASTQC directory created: Results will be saved in",file.path(sbj$subjectDir,"FASTQC"),"\n"))
    }else{
      cat(paste("FASTQC directory creation FAILED: ",file.path(sbj$subjectDir,"FASTQC"),"\n"))
      return()
    }
  }else{
    cat(paste("FASTQC directory already created : ",file.path(sbj$subjectDir,"FASTQC"),"\n"))
  }
  
  t<-Sys.time()
  # print(c(sbj$files,"-t 2", paste("-o ",file.path(sbj$subjectDir,"FASTQC"),sep="")))
  if(allFiles == TRUE){
    cat("\n All fastq files\n")
    system2(fastqc.path, args = c(sbj$fastq.files,"-t 2", paste("-o ",file.path(sbj$subjectDir,"FASTQC"),sep="")))
  }else{
    cat("\n Only raw files\n")
    system2(fastqc.path, args = c(sbj$files,"-t 2", paste("-o ",file.path(sbj$subjectDir,"FASTQC"),sep="")))
  }
    
  tiempo <- Sys.time() - t
  cat("\n###########################\n")
  cat(tiempo)
  cat("\n###########################\n")
      
}

RUN_FastQC2 <- function(subjPath, subj, fastQcPath = fastqc.path, allFiles = c("raw","trim","all"), wait = FALSE){
  sbj <- BuildPatientFiles2(subjPath,subj)
  allFiles <- match.arg(allFiles, c("raw","trim","all"))
  if(any(is.na(sbj))) {
    stop("Error")
  }
  
  if(!dir.exists(file.path(sbj$subjectDir,"FASTQC"))){
    dir.create(file.path(sbj$subjectDir,"FASTQC"))  
    if(dir.exists(file.path(sbj$subjectDir,"FASTQC"))) {
      cat(msg(paste("FASTQC directory created: Results will be saved in",file.path(sbj$subjectDir,"FASTQC"),"\n")))
    }else{
      cat(error(paste("FASTQC directory creation FAILED: ",file.path(sbj$subjectDir,"FASTQC"),"\n")))
      return()
    }
  }else{
    cat(msg(paste("FASTQC directory already created : ",file.path(sbj$subjectDir,"FASTQC"),"\n")))
  }
 
  if(allFiles == "all"){
    arguments <- c(sbj$fastq.file,
                   paste("-t ",length(sbj$fastq.file),sep=""), 
                   paste("-o ",file.path(sbj$subjectDir,"FASTQC"),sep=""))
    cat("\n All fastq files\n")
  }
  if(allFiles == "trim"){
     files <- sbj$fastq.files[str_detect(sbj$fastq.files,"val")]
     if(length(files)==0) stop(error("trimmed files NOT found"))
      arguments <- c(files,
                   paste("-t ",length(files),sep=""), 
                   paste("-o ",file.path(sbj$subjectDir,"FASTQC"),sep=""))
      cat("\n trimmed files\n")
  }
  if(allFiles == "raw"){
    files <- sbj$fastq.files[!str_detect(sbj$fastq.files,"val")]
    arguments <- c(files,
                   paste("-t ",length(files),sep=""), 
                   paste("-o ",file.path(sbj$subjectDir,"FASTQC"),sep=""))
    
  }
  # arguments <- c(sbj$fastq.file,
  #                paste("-t ",length(sbj$fastq.file),sep=""), 
  #                paste("-o ",file.path(sbj$subjectDir,"FASTQC"),sep=""))
   print(arguments)
  
    
  t<-Sys.time()
      system2(fastqc.path, args = arguments, wait = wait)
  tiempo <- Sys.time() - t
  cat("\n###########################\n")
  cat(paste("Elapsed time (mins): ",tiempo))
  cat("\n###########################\n")
  
}
###-----------------------------------
RUN_TrimGalore <- function(subjPath,subject, trimgalorePath,trimVersion = "0.6.4"){
  sbj <- BuilPatientFiles(subjPath,subject)
  trimVersion <- str_replace_all(trimVersion,"\\.","_")
  comm <- file.path(trimgalorePath,paste("trim_galore_",trimVersion,sep=""))
  arguments <- c("--paired",sbj$files, "--cores 2 " ,"--path_to_cutadap ~/.local/bin/cutadapt ", paste("-o ",sbj$subjectDir,sep=""))
  print(c(comm,arguments))
  
  t<-Sys.time()
  
  system2(command = comm, args = arguments)  
  
  tiempo <- Sys.time() - t
  cat("\n###########################\n")
  cat(paste("Elapsed time (mins): ",tiempo))
  cat("\n###########################\n")
  
  # ./trim_galore_0_6_4 
  # --paired 
  # --fastqc 
  # /home/elmer/FLENI/Pacientes/34708/34708_1.fastq.gz /home/elmer/FLENI/Pacientes/34708/34708_2.fastq.gz
}

TrimGalore <- function(pathFiles, baseNameFile, trimgalorePath,trimVersion = "0.6.4"){
  # sbj <- BuilPatientFiles(subjPath,subject)
  sbj <- vector(mode = "list", length = 2)
  names(sbj) <- c("files", "subjectDir")
  sbj$files <- c(file.path(pathFiles, paste(baseNameFile,"_1.fastq",sep="")),
                file.path(pathFiles, paste(baseNameFile,"_2.fastq",sep="")))
  sbj$subjectDir <- pathFiles
  trimVersion <- str_replace_all(trimVersion,"\\.","_")
  comm <- file.path(trimgalorePath,paste("trim_galore_",trimVersion,sep=""))
  arguments <- c("--paired",sbj$files, "--path_to_cutadap ~/.local/bin/cutadapt", paste("-o ",sbj$subjectDir,sep=""))
  print(c(comm,arguments))
  cat("\n All fastq files\n")
  t<-proc.time()
  system2(command = comm, args = arguments)
  tiempo <- proc.time() - t
  cat("\n###########################\n")
  cat(tiempo)
  cat("\n###########################\n")
}

RNAseq_allign <- function(subjPath,subject, pathToIndex, useTrimmed = TRUE){
  sbj <- BuildPatientFiles2(subjPath,subject)
  if(useTrimmed == TRUE){
    res <-subjunc(index = pathToIndex,
             readfile1 = sbj$fastq.trimmed.files[str_detect(sbj$fastq.trimmed.files,"val_1")],
             readfile2 = sbj$fastq.trimmed.files[str_detect(sbj$fastq.trimmed.files,"val_2")],
             nthreads = 6,
             reportAllJunctions = TRUE,
             )   
  }
  return(res)
}

DNAseq_allign <- function(subjPath,subject, pathToIndex, useTrimmed = TRUE){
  sbj <- BuildPatientFiles2(subjPath,subject)
  if(useTrimmed == TRUE){
    res <-align(index = pathToIndex,
                  readfile1 = sbj$fastq.trimmed.files[str_detect(sbj$fastq.trimmed.files,"val_1")],
                  readfile2 = sbj$fastq.trimmed.files[str_detect(sbj$fastq.trimmed.files,"val_2")],
                  type = "dna",
                  detectSV = TRUE,
                  annot.inbuilt = "hg19",
                  nthreads = 6
                )   
  }else{
    res <-align(index = pathToIndex,
                readfile1 = sbj$fastq.files.raw[str_detect(sbj$fastq.files.raw,"_1.")],
                readfile2 = sbj$fastq.files.raw[str_detect(sbj$fastq.files.raw,"_2.")],
                type = "dna",
                detectSV = TRUE,
                annot.inbuilt = "hg19",
                nthreads = 6
    )   
  }
  return(res)
}
  
## Auxiliary functions ####
##Only internal functions
GetBAMfileFromAlligner <- function(sbj, alligner = c("subjunc", "allign","bwa","wes.macrogen")){
  alligner <- match.arg(alligner,c("subjunc", "ralign","bwa","wes.macrogen"))
  bam.file <- switch( alligner,
          "subjunc" = sbj$subjunc.bam,
          "ralign" = sbj$ralign.bam,
          "bwa" = sbj$bwa.bam,
          "wes.macrogen" = sbj$macrogen.recal.bam)
}

CountsToTPM <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
####
RNAseq_counts <- function(subjPath,subject, featureType = c("gene","exon"),
                          genome = c("hg38","hg19"),
                          alligner = c("subjunc", "rallign","bwa","wes.macrogen"),
                          ncores = 4){
  alligner <- match.arg(alligner,c("subjunc", "ralign","bwa","wes.macrogen"))
  featureType <- match.arg(featureType, c("gene","exon"))
  genome <- match.arg(genome , c("hg38","hg19") )
  
  ncores <- abs(ncores)
  if(ncores==0) ncores <- 1
  
  ns <- length(subject)
  if(ns == 1){
    sbj <- BuildPatientFiles2(subjPath,subject)
    bam.files <- GetBAMfileFromAlligner(sbj,alligner) 
  }else{
    bam.files <- unlist(lapply(subject, function(x){
      sbj <- BuildPatientFiles2(subjPath,x)
      return(GetBAMfileFromAlligner(sbj,alligner))
    }))
  }
  
  fcounts <- featureCounts(files = bam.files,
                           annot.inbuilt = genome,
                           useMetaFeatures = featureType == "gene",
                           isPairedEnd = TRUE,
                           requireBothEndsMapped = TRUE,
                           countChimericFragments = FALSE, #un chimerico no es un gen
                           nthreads = ncores)
  fcounts$TPM <- CountsToTPM(fcounts$counts, fcounts$annotation$Length )
  if(ns==1){
    count.file <- file.path(sbj$subjectDir, paste(sbj$subject,"_Counts_",featureType,"_",genome,".xlsx",sep=""))
    write.xlsx(fcounts, file = count.file)
    cat(paste("\n Counts saved at: ",count.file,"\n") )
  }else{
    count.file <- paste(subjPath,"Counts_",featureType,"_",genome,".xlsx",sep="")
    write.xlsx(fcounts, file = count.file)
    cat(paste("\n Counts saved at: ",count.file,"\n") )
  }
  return(invisible(fcounts))
}


###Sequence Statistics Functions ####
#BreathRegionInPanel: This function estimates the covered % of the sequencing on the specified 
#panel
sel.canonico <- function(bed){
  transcripts <- unique(bed$NR)
  if( length(transcripts) == 1){
    return(bed)
  }
  
  lt <- unlist(lapply(transcripts, function(x) sum(bed$NR==x)))
  mp <- which.max(lt)[1]
  return(subset(bed, NR == transcripts[mp]))
}

fung <- function(gene.search, bed, bam.file, maxDepth = 1000, nucleotide = FALSE, offset = TRUE){
  if(any(stringr::str_detect(bed$GeneSymbol, gene.search) )==FALSE){##arreglado any 22/04/24
    return(c(rep(NA,10)))
  }
  ##el ofset no me queda muy claro, pero pareceria que empiezan en -1
  off <- ifelse(offset,1,0)
  
  print(gene.search)
  gene.bed <- subset(bed, GeneSymbol == gene.search)
  if(any(dim(gene.bed)==0)){#added 10/03/2024
    return(c(rep(NA,10)))
  }
  ## viene del hg19 que esta en la referencia
  if(any(stringr::str_detect(colnames(gene.bed),"NR"))==TRUE){
    gene.bed <- sel.canonico(gene.bed)  
  }
  
  chrs <- which(gene.bed$Chr %in% paste0("chr",c(as.character(rep(1:22)), "X", "Y")))
  gene.bed <- gene.bed[chrs,]
  gene.bed$Chr <- as.character(gene.bed$Chr)
  
  
  gene.gr <- GRanges(seqname = gene.bed$Chr, IRanges(start = gene.bed$start+off,end = gene.bed$end))
  
  sbp <- ScanBamParam(which=gene.gr)
  
  p_param <- PileupParam(distinguish_strands=FALSE,
                         distinguish_nucleotides = nucleotide,
                         include_insertions = TRUE,
                         min_nucleotide_depth=1,
                         max_depth = maxDepth)
  
  pl <- pileup(bam.file, scanBamParam = sbp, pileupParam = p_param)
  
  # return(pl)
  stats.summary <- summary(pl$count)[-1]
  # names(stats.summary) <- c("1Q","2Q","Mean","3Q","Max")
  Length <- sum(width(gene.gr),na.rm=T)
  cLength <- sum(pl$count>0,na.rm=T)
  cLength_10x <- sum(pl$count>=10,na.rm=T)
  per <- round(100*cLength/Length,2)
  per_10x <- round(100*cLength_10x/Length,2)
  res <- c(Start = min(start(gene.gr)),
           End = max(end(gene.gr)), 
           Length = Length,
           cLength = cLength,
           cLength_10x = cLength_10x,
           "physCov" = per,
           "physCov_10x" = per_10x,
           round(stats.summary,2)
  )
  
  return(res)
}

#' BreathRegionInPanel
#' This function calculatethe Phicical coverage region (the true length that was sequenced for each exon)
#' @param patDir string whith the patient directory main root
#' @param subjectID string with the patientID (folder name where the bam file is stored)
#' @param panel string, the panel file name
#' @param bed the bed file or the bed data.frame
#' @param depth the minimal depth to consider in the analysis
#' @BioinfProvider string c("MODApy","Macrogen"), to indicate whoch provider performed the alignement. (default MODApy)
#' @return an excel file with the covarage statistics. The file name is SubjectID_MODApy_CoverageStat.xlsx
#' in addition, the resulting is invisible returned 
#' 
BreathRegionInPanel <- function(patDir, subjectID, panel, bed, depth, BioinfProvider = c("MODApy","Macrogen")){
  BioinfProvider <- match.arg(BioinfProvider, choices = c("MODApy","Macrogen"))
  bam.file <- normalizePath(list.files(file.path(patDir,paste0(subjectID,"_MODApy")),pattern = ".bam", full.names = T))
   if(BioinfProvider == "MODApy"){
     # bam.file <- file.path(patDir,paste0(subjectID,"_MODApy"),paste(subjectID,"_MODApy_realigned_reads_recal.bam",sep=""))
     bam.file <- normalizePath(list.files(file.path(patDir,paste0(subjectID,"_MODApy")),pattern = ".bam", full.names = T))
   }else{
     bam.file <- file.path(patDir,subjectID,paste(subjectID,".recal.bam",sep=""))  
   }
  
 cat(bam.file)
  
  if(file.exists(bam.file)==FALSE){
    stop(paste("BAM FILE not found:",bam.file))
  }
  if(file.exists(panel)== FALSE){
    stop(error(paste("Panel file not found :", panel)))
  }
  # 
  panel.gsymbol <- try(read.xlsx(panel, sheet = "GeneList")$GeneSymbol)
  if(class(panel.gsymbol) == "try-error" ){
    stop(error(paste("GeneList sheet or GeneSymbol column not found ", basename(panel))))
  }
  panel.gsymbol <- na.omit(panel.gsymbol)
  mc.cores <- max(1,detectCores()-2)
  
  cat(msg(paste("\n... Using ", mc.cores," CPUs...\n")))
  cat(msg(paste("\nProcessing ,", length(panel.gsymbol),"genes \n")))


  pileup.stats <- data.frame(do.call(rbind,bplapply(panel.gsymbol, fung, bed=bed,
                                         bam.file=bam.file, maxDepth = 4*depth,
                                         BPPARAM = MulticoreParam(workers =  mc.cores))
                                     )
                             )

  
  # fung(gene.search = "MAPT", bed=bed, bam.file=bam.file, maxDepth = 4*depth)
  
                               
  #
  ss <- summary(pileup.stats$physCov)
  rownames(pileup.stats) <- panel.gsymbol  
        cat("\n.... Panel coverage Statistics ....")
  cat(paste("\nSubject             : ", subjectID))
  cat(paste("\nPanel file name     : ", basename(panel)))
  cat(paste("\nGenes in panel      : ", nrow(pileup.stats)))
        cat("\nPhysical coverage   = #of gen positions with reads / gene length")
  for(n in names(ss)) {
    nn <- ss[n]
    if(n == "NA's") {
      n ="# not found genes"
    }
    l <- 18-str_length(n)
    cat(paste("\n",n,paste(rep(" ",l),collapse = ""), ": ",round(nn,2),sep=""))
    }
  
  cat("\n....................................................\n")
  if(any(is.na(pileup.stats$physCov))){
    cat("\n.... Gene Symbols not detected in BED file ....\n")
    for(i in which(is.na(pileup.stats$physCov)))
      cat(error(paste(rownames(pileup.stats)[i],"\n")))   
    cat("Please check Gene Symbols in panel file!!!")
  }
  
  if(BioinfProvider == "MODApy"){
    res.file <- file.path(patDir,paste0(subjectID,"_MODApy"),paste(subjectID,"_",str_remove(basename(panel),".xlsx"),"_CoverageStat.xlsx",sep=""))
  }else{
    res.file <- file.path(patDir,subjectID,paste(subjectID,"_",str_remove(basename(panel),".xlsx"),"_CoverageStat.xlsx",sep=""))
  }
  pileup.stats<-rbind(pileup.stats,"---"=rep("---", ncol(pileup.stats)-1), 
        Total = round(apply(pileup.stats,2,mean),2), 
        "Total>10" = round(apply(subset(pileup.stats, Median >=10),2,mean),2))
  write.xlsx(list("Panel Coverage Statistics" = pileup.stats), 
             file = res.file, rowNames = TRUE, overwrite = TRUE)
  cat(paste("\nResults saved at ", res.file,"\n"))
  cat("\\.................................................................../\n")
   return(invisible(pileup.stats))
}

PileUpGenePlot <- function(gene, bed, bam.file, maxDepth = 1000, offset = TRUE){
  if(any(str_detect(bed$GeneSymbol, gene) )==FALSE){
    return(c(rep(NA,10)))
  }
  ##el ofset no me queda muy claro, pero pareceria que empiezan en -1
  off <- ifelse(offset,1,0)
  
  print(gene)
  gene.bed <- subset(bed, GeneSymbol == gene)
  gene.bed <- sel.canonico(gene.bed)
  gene.gr <- GRanges(seqname = gene.bed$Chr, IRanges(start = gene.bed$start+off,end = gene.bed$end))
  
  sbp <- ScanBamParam(which=gene.gr)
  
  p_param <- PileupParam(distinguish_strands=FALSE,
                         distinguish_nucleotides = FALSE,
                         include_insertions = TRUE,
                         min_nucleotide_depth=1,
                         max_depth = maxDepth)
  
  pl <- pileup(bam.file, scanBamParam = sbp, pileupParam = p_param)
  gene.bed$exon <- stringr::str_remove_all(gene.bed$exon, "exon_")
  gene.bed$exonID <- paste(as.character(gene.bed$Chr[1]),":",gene.bed$start+off,"-",gene.bed$end,sep="")
  pl <- merge(pl,gene.bed[,c("exonID","exon")],by.x="which_label",by.y = "exonID",sort=FALSE)
  ggplot(pl, aes(x=pos, y = count)) + geom_line() + facet_wrap(~exon, ncol = ifelse(length(unique(pl$exon))>20,10,5), scale = "free")
}
