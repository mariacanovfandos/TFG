# SCRIPT 1: PREPROCESAMIENTO datos de bases de datos
#==================================================

# 1. Cargar librerías necesarias
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("here")) install.packages("here")

library(tidyverse)
library(here)

ruta_descargables <- here("../DESCARGABLES")
ruta_resultados   <- here("../RESULTADOS")

# Crear la carpeta de resultados si no existe
if (!dir.exists(ruta_resultados)) dir.create(ruta_resultados)

# ==========================================
# PASO 1: Cargar y limpiar DrugBank
# ==========================================

# 1. Leer y ver el archivo de DrugBank
drugbank_inicial <- read_csv(file.path(ruta_descargables, "DRUGBANK_Vocabulary.csv"))
View(drugbank_inicial)

# 2. Seleccionar columnas de interés (estandarización)
drugbank_interes <- drugbank_inicial %>%
  # Mantener solo las columnas de interés
  select(`DrugBank ID`, `Common name`, `Standard InChI Key`) %>%
  # Cambiar el nombre para que no tengan espacios
  rename(
    ID_DrugBank = `DrugBank ID`,
    Nombre_DrugBank = `Common name`,
    InChIKey = `Standard InChI Key`
  ) %>%
  # Descartar filas vacías: Si no tiene InChIKey, no sirve para cruzar datos
  drop_na(InChIKey)

View(drugbank_interes)
# La lista drugbank_interes contiene los químicos que en la lista inicial
# tenían InChIKey, los químicos que no lo tenían se han descartado


# ==========================================
# PASO 2: Cargar y limpiar CTD
# ==========================================

# 1. Leer el archivo ignorando los comentarios (#)
ctd_inicial <- read_tsv(file.path(ruta_descargables, "CTD_Chemicals.tsv"), 
                        comment = "#", 
                        col_names = FALSE)

# 2. Asignar nombres a las columnas
colnames(ctd_inicial) <- c("ChemicalName", "ChemicalID", "CasRN", "PubChemCID", "PubChemSID", "DTXSID", "InChIKey", "Definition", "ParentIDs", "TreeNumbers", "ParentTreeNumbers", "MESHSynonyms", "CTDCuratedSynonyms")

View(ctd_inicial)

# 3. Limpiar y seleccionar columnas de interés
ctd_interes <- ctd_inicial %>%
  # Seleccionar y mantener las tres variables de interés (InChIKey se queda igual) (estandarización)
  select(ChemicalID, ChemicalName, InChIKey) %>%
  rename(
    ID_MeSH = ChemicalID,
    Nombre_CTD = ChemicalName
  ) %>%
  # Quitar las filas que no tengan InChIKey
  drop_na(InChIKey)

View(ctd_interes)


# ==========================================
# PASO 3: Cargar y limpiar ChEBI
# ==========================================

# 1. Leer los dos archivos (en ChEBI son 2 archivos, uno de estrucutra y otro de nombres)
chebi_compounds_inicial <- read_tsv(file.path(ruta_descargables, "CHEBI_Compounds.tsv"))
chebi_structures_inicial <- read_tsv(file.path(ruta_descargables, "CHEBI_Structures.tsv"))
View(chebi_compounds_inicial)
View(chebi_structures_inicial)

# 2. Limpiar tabla estructuras (tiene InChIKey) (estandarización)
chebi_structures_interes <- chebi_structures_inicial %>%
  select(compound_id, standard_inchi_key) %>%
  rename(
    ID_ChEBI = compound_id,
    InChIKey = standard_inchi_key
  ) %>%
  # Quitar los que no tengan InChIKey
  drop_na(InChIKey)

# 3. Limpiar la tabla de compunds (no tiene InChIKey) (estandarización)
chebi_compounds_interes <- chebi_compounds_inicial %>%
  select(name, chebi_accession) %>%
  rename(
    Nombre_ChEBI = name,
    ID_ChEBI = chebi_accession
  )

view(chebi_structures_interes)
view(chebi_compounds_interes)


# 4. Corregir listas para poder hacer la unión
# Quitar "CHEBI:", convertir a texto y quitar espacios para poder unir
chebi_compounds_interes$ID_ChEBI <- trimws(gsub("CHEBI:", "", as.character(chebi_compounds_interes$ID_ChEBI)))
View(chebi_compounds_interes)

chebi_structures_interes$ID_ChEBI <- trimws(gsub("CHEBI:", "", as.character(chebi_structures_interes$ID_ChEBI)))
View(chebi_structures_interes)

# 5. Unir tablas ChEBI en una sola por ID
chebi_unificado_interes <- inner_join(chebi_structures_interes, chebi_compounds_interes, by = "ID_ChEBI")
View(chebi_unificado_interes)


# ==========================================
# PASO 4: Crear tabla con las 3 bases de datos juntas
# ==========================================

quimicos_interes_todos <- ctd_interes %>%
  # 1. Pegar DrugBank usando el InChIKey
  full_join(drugbank_interes, by = "InChIKey") %>%
  # 2. Pegar ChEBI usando también el InChIKey
  full_join(chebi_unificado_interes, by = "InChIKey")

view(quimicos_interes_todos)


# ==========================================
# PASO 5: Guardar datos
# ==========================================

write_csv(quimicos_interes_todos, file.path(ruta_resultados, "TFG_quimicos_interes.csv"))
saveRDS(quimicos_interes_todos, file.path(ruta_resultados, "TFG_quimicos_interes.rds"))

message("Proceso finalizado. Archivos guardados en: ", ruta_resultados)
