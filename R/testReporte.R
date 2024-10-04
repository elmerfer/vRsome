file.copy(from = "/home/biomolecular/DATA/NGS/ONCO/PACIENTES/49644_D_MODApy/49644_D_MODApy.final.vcf", 
          to = "/home/efernandez/49644_D_MODApy.final.vcf", overwrite = FALSE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
BuildVarsomeVCF_FULL("/home/efernandez/49644_D_MODApy_GeneList_Neuroonco.annotated_TEST.xlsx")


tabla_reporte <- as.data.table(openxlsx::read.xlsx("/home/efernandez/BM23-48847_Alzheimer-exported-variants.xlsx"))
require(data.table) ## v1.9.6+

View(tabla_reporte)
tabla_reporte[, HGVS:= stringr::str_replace_all(HGVS," ",":") ]
tabla_reporte[, HGVS:= stringr::str_replace_all(HGVS,"&gt;",">") ]
tabla_reporte[, HGVS:= stringr::str_remove_all(HGVS,"\\)|\\(") ]
tabla_reporte$HGVS
tabla_reporte[, c("NM", "Var","Prot") := tstrsplit(HGVS, ":", fixed=TRUE)]

d2t <- data.table(GEN = tabla_reporte$Genes,
                  VARIANTE = tabla_reporte$Var,
                  NM = tabla_reporte$NM,
                  "PREDICCION PROTEICA"= tabla_reporte$Prot,
                  CIGOSIDAD = toupper(stringr::str_sub(tabla_reporte$Zygosity,end=3)),
                  "FRECUENCIA ALELICA" = tabla_reporte$Frequency,
                  ACMG = tabla_reporte$Germline.Class)
tinytex::install_tinytex()
