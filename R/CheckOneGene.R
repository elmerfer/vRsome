##checkOneGene
library(data.table)
library(openxlsx)

hg19_genes <- readRDS("Data/hg19_genes.RDS")

gene_to_check <- "ADAM33" #Poner el nombre entre los parentesis

print(hg19_genes[which(hg19_genes$GeneName == gene_to_check),])
