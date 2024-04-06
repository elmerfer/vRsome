library(rtracklayer)
library(data.table)

bed<-fread("/home/elmer/Elmer/FLENI/PanelVCF_VARSOM/hg19-from-snpeff/hg19-canon-txt.bed")
bed2 <-fread("/home/elmer/Elmer/FLENI/PanelVCF_VARSOM/hg19-from-snpeff/hg19-txt.bed")
colnames(bed)
head(bed)
subset(bed, geneName == "C2-AS1")
length(unique(bed2$geneName))

bed_genes <- tidytable::distinct(bed, geneName, .keep_all = T)
dim(bed_genes)

dim(hg19_genes)
colnames(hg19_genes)
colnames(bed_genes)
sum(hg19_genes$GeneName %in% bed_genes$geneName)
sum(bed_genes$geneName %in% hg19_genes$GeneName)

panel <- openxlsx::read.xlsx("/home/elmer/Elmer/FLENI/Paneles/GeneList_Rasopatias.xlsx")
View(panel)

estas<-panel$GeneSymbol %in% hg19_genes$GeneName
estas
saveRDS(bed_genes,"Data/hg19_snpEFF_genes.RDS")
panel$GeneSymbol %in% bed_genes$geneName
