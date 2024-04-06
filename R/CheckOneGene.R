##checkOneGene
library(data.table)
library(openxlsx)


### --------------------------
### MODIFICAR AQUI
gene_to_check <- "JUAN" #Poner el nombre entre los parentesis
#############################



hg19_genes <- readRDS("Data/hg19_genes.RDS")
hg19_snpEFF_genes <- readRDS("Data/hg19_snpEFF_genes.RDS")

cat("Evaluando existencia contra HG19 bed\n")
if(length(which(hg19_genes$GeneName == gene_to_check))<1){
  cat(paste0("El Gen dado :", gene_to_check, " NO ESTA en Hg19"))
}else{
  cat(paste0((hg19_genes[which(hg19_genes$GeneName == gene_to_check),])))  
}
cat("\n------------------------\n")
cat("Evaluando existencia contra HG19 del snpEFF - Anotador de variantes\n")
if(length(which(hg19_snpEFF_genes$geneName == gene_to_check))<1){
  cat(paste0("El Gen dado :", gene_to_check, " NO ESTA en Hg19 del anotador snpEFF"))
}else{
  cat(paste0((hg19_snpEFF_genes[which(hg19_snpEFF_genes$geneName == gene_to_check),])))
}
