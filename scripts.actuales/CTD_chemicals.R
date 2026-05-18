# ------------- LIMPIAR Y SERIALIZAR CTD_chemicals -------------

# ------------- LIMPIEZA --------------

# 1. Cargar librerías
library(readr)
library(dplyr)
library(data.table)
library(tidyverse)
library(rdflib)

# 2. Leer archivo
datos_quimicos <- fread("https://ctdbase.org/reports/CTD_chemicals.csv.gz", 
                        skip = 27, stringsAsFactors = FALSE)

# 3. Limpieza archivo
datos_quimicos <- datos_quimicos %>%
  set_names(c("ChemicalName", "ChemicalID", "CasRN", "PubChemCID", "PubChemSID", 
              "DTXSID", "InChIKey", "Definition", "ParentIDs", 
              "TreeNumbers", "ParentTreeNumbers", "MESHSynonyms", "CTDCuratedSynonyms")) %>%
  # Seleccionar columnas de interés
  select(ChemicalName, ChemicalID, InChIKey, Definition) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~na_if(., ""))) %>%
  # MEJORA: Filtrar compuestos sin InChIKey ANTES del bucle para evitar colapsar la memoria de R
  filter(!is.na(InChIKey)) %>%
  distinct()

# ----------- SERIALIZACIÓN ------------

# 1. Iniciar grafo vacío
grafo_quimicos <- rdf()

# 2. Definir los prefijos estandarizados
sio <- "http://semanticscience.org/resource/"
rdfs <- "http://www.w3.org/2000/01/rdf-schema#"
rdf_base <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
mesh_prefix <- "https://id.nlm.nih.gov/mesh/" 
schema <- "https://schema.org/"
owl <- "http://www.w3.org/2002/07/owl#"
inchikey_prefix <- "https://identifiers.org/inchikey/"
skos <- "http://www.w3.org/2004/02/skos/core#"

# 3. Bucle para procesar toda la tabla optimizada
for (i in 1:nrow(datos_quimicos)) {
  
  ChemicalName <- datos_quimicos$ChemicalName[i]
  ID_inicial <- datos_quimicos$ChemicalID[i]  
  InChIKey <- datos_quimicos$InChIKey[i]
  Definicion <- datos_quimicos$Definition[i] 
  
  ID_solo_codigo <- sub("MESH:", "", ID_inicial)
  uri_quimico <- paste0(mesh_prefix, ID_solo_codigo)
  
  # --- TRIPLETAS ---
  
  # A. Tipo (Químico) - Usando el prefijo base unificado
  rdf_add(grafo_quimicos, 
          subject = uri_quimico, 
          predicate = paste0(rdf_base, "type"), 
          object = paste0(sio, "SIO_010004"))
  
  # B. Label
  rdf_add(grafo_quimicos, 
          subject = uri_quimico, 
          predicate = paste0(rdfs, "label"), 
          object = ChemicalName)
  
  # C. Identificador
  rdf_add(grafo_quimicos, 
          subject = uri_quimico, 
          predicate = paste0(schema, "identifier"), 
          object = ID_solo_codigo)
  
  # D. sameAs (InChIKey) con blindaje de espacios
  InChIKey_limpio <- trimws(InChIKey)
  uri_inchi <- paste0(inchikey_prefix, InChIKey_limpio)
  
  rdf_add(grafo_quimicos, 
          subject = uri_quimico, 
          predicate = paste0(owl, "sameAs"), 
          object = uri_inchi)
  
  # E. Definición (SKOS)
  if (!is.na(Definicion)) {
    rdf_add(grafo_quimicos, 
            subject = uri_quimico, 
            predicate = paste0(skos, "definition"), 
            object = Definicion)
  }
  
  # Avance en la consola de R
  if (i %% 10000 == 0) {
    message(paste("Procesados", i, "químicos de CTD"))
  }
}

# 4. Guardar el archivo
ruta_guardado_quimicos <- "/Users/mersmac/Desktop/TFG/RESULTADOS/Chemicals_ctd.ttl"
rdf_serialize(grafo_quimicos, doc = ruta_guardado_quimicos, format = "turtle")
