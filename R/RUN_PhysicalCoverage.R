#Autor: Elmer Fernández
#Date: 06/03/2024
#source('~/NGS/R_Pipelines/NGS.Utils.R')
# source('~/Dropbox/Fleni/R_Pipelines/NGS.Utils.R')
main.path <- "/home/biomolecular/DATA"
file.exists('~/DATA/Dropbox/Fleni/R_Pipelines/NGS.Utils.R')
source(file.path(main.path,'Dropbox/Fleni/R_Pipelines/NGS.Utils.R'))
df <- readRDS("~/DATA/NGS/References/hg19/hg19.exon.bed.RDS")

# pats.path <- "/home/elmer/FLENI/Pacientes"
# pats.path <- "/home/biomolecular/DATA/NGS/Pacientes"

panel.dir <- "~/DATA/NGS/Paneles"
#######################################
## DONDE ESTA EL PACIENTE ?
####################################### 
pats.path <- "/home/biomolecular/DATA/NGS/Pacientes"
#######################################
## ES ONCOLOGICO? PONER TRUE O FALSE SI NO LO ES
####################################### 
isONCO <- TRUE #FALSE
#######################################
## QUIEN ES EL PROVEEDOR DE BIOINFORMATICA
## MODApy o Macrogen ? NO OLVIDE CAMBIAR LA LINEA QUE SIGUE
####################################### 
BioinformaticProvider <- "MODApy" 
#######################################
## SOLO MODIFICAR PACIENTE ID
#######################################
## Seleccionar paciente 
PacienteID <- "BM24-49914"
PanelName <- "GeneList_TESTEO_ELA_NO_USAR.xlsx"
#######################################
st <-BreathRegionInPanel(patDir = pats.path, 
                    subjectID =  PacienteID,
                    panel = file.path(panel.dir,PanelName), 
                    depth = ifelse(isONCO, 500,100),
                    bed = df,
                    BioinfProvider = BioinformaticProvider)
