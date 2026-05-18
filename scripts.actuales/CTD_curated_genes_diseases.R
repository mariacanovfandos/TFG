# ------------- LIMPIAR Y SERIALIZAR CTD_curated_genes_diseases -------------

# ------------- LIMPIEZA --------------
# 1. Cargar librerías
library(dplyr)
library(data.table)
library(tidyverse)

# 2. Cargar el archivo saltando las líneas
datos_gen_enfermedad <- fread("https://ctdbase.org/reports/CTD_curated_genes_diseases.csv.gz", 
                              skip = 27, stringsAsFactors = FALSE)

# 3. Limpieza archivo
datos_gen_enfermedad <- datos_gen_enfermedad %>%
  set_names(c("GeneSymbol", "GeneID", "DiseaseName", "DiseaseID", "DirectEvidence", 
              "OmimIDs", "PubMedIDs")) %>%
  # Seleccionar columnas de interés
  select(GeneSymbol, DiseaseID, DirectEvidence, PubMedIDs) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~na_if(., ""))) %>%
  # Comprobar que hay GeneSymbol y DiseaseID
  filter(!is.na(GeneSymbol) & !is.na(DiseaseID)) %>%
  distinct()


# ----------- SERIALIZACIÓN ------------
grafo_gen_enfermedad <- rdf()

# 2. Definir prefijos
rdf <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
biolink <- "https://w3id.org/biolink/vocab/"
sio <- "http://semanticscience.org/resource/"
schema <- "https://schema.org/"
mesh_prefix <- "https://id.nlm.nih.gov/mesh/"
omim_prefix <- "http://purl.bioontology.org/ontology/OMIM/"
obo <- "http://purl.obolibrary.org/obo/"
pubmed_base <- "https://pubmed.ncbi.nlm.nih.gov/"
assoc_base <- "http://tfg.org/gene_disease_association/"
gene_base <- "http://rdf.biogateway.eu/gene/9606/"

# 3. Bucle
for (i in 1:nrow(datos_gen_enfermedad)) {
  
  # --- EXTRAER LOS DATOS ---
  gen_symbol <- datos_gen_enfermedad$GeneSymbol[i] 
  enfermedad_id <- datos_gen_enfermedad$DiseaseID[i] 
  evidencia <- datos_gen_enfermedad$DirectEvidence[i] 
  pubmed <- datos_gen_enfermedad$PubMedIDs[i] 
  
  # --- CORRECCIÓN: MANEJO INTELIGENTE DE MESH Y OMIM ---
  if (grepl("MESH:", enfermedad_id)) {
    id_enf_limpio <- sub("MESH:", "", enfermedad_id)
    uri_enfermedad <- paste0("<", mesh_prefix, id_enf_limpio, ">")
  } else if (grepl("OMIM:", enfermedad_id)) {
    id_enf_limpio <- sub("OMIM:", "", enfermedad_id)
    uri_enfermedad <- paste0("<", omim_prefix, id_enf_limpio, ">")
  } else {
    # Por si CTD introduce un formato inesperado
    id_enf_limpio <- enfermedad_id
    uri_enfermedad <- paste0("<", mesh_prefix, id_enf_limpio, ">")
  }
  
  uri_asociacion <- paste0("<", assoc_base, gen_symbol, "--", id_enf_limpio, ">")
  uri_gen <- paste0("<", gene_base, gen_symbol, ">")
  
  # Vector temporal para las líneas
  lineas <- c()
  
  # --- CONSTRUIR TRIPLETAS ---
  # Type
  lineas <- c(lineas, paste(uri_asociacion, 
                            paste0("<", rdf, "type>"), 
                            paste0("<", sio, "SIO_000983> .")))
  
  # A. El Sujeto es el Gen (Tu corrección aplicada)
  lineas <- c(lineas, paste(uri_asociacion, 
                            paste0("<", rdf, "subject>"), 
                            uri_gen, "."))
  
  # B. El Objeto es la Enfermedad (Tu corrección aplicada)
  lineas <- c(lineas, paste(uri_asociacion, 
                            paste0("<", rdf, "object>"), 
                            uri_enfermedad, "."))
  
  # Relación directa gen --> enfermedad
  lineas <- c(lineas, paste(uri_gen, 
                            paste0("<", obo, "RO_0002331>"), 
                            uri_enfermedad, "."))
  
  # C. Evidencia Directa (si existe)
  if (!is.na(evidencia)) {
    lista_evidencias <- strsplit(evidencia, "\\|")[[1]]
    for (ev in lista_evidencias) {
      lineas <- c(lineas, paste(uri_asociacion, 
                                paste0("<", biolink, "has_evidence>"), 
                                paste0("\"", ev, "\" .")))
    }
  }
  
  # D. PubMed
  if (!is.na(pubmed)) {
    lista_pubmeds <- strsplit(pubmed, "\\|")[[1]]
    for (pmid in lista_pubmeds) {
      uri_articulo <- paste0("<", pubmed_base, pmid, ">")
      
      # Conectar asociación con artículo
      lineas <- c(lineas, paste(uri_asociacion, 
                                paste0("<", sio, "SIO_000772>"), 
                                uri_articulo, "."))
      
      # Definir el tipo de artículo
      lineas <- c(lineas, paste(uri_articulo, 
                                paste0("<", rdf, "type>"), 
                                paste0("<", sio, "SIO_000154> .")))
      
      # Identificador del artículo
      lineas <- c(lineas, paste(uri_articulo, 
                                paste0("<", schema, "identifier>"), 
                                paste0("\"", pmid, "\" .")))
    }
  }
  # Avance
  if (i %% 10000 == 0) {
    message(paste("Procesadas", i, "filas de CTD Curated Genes-Diseases"))
  }
}

ruta_guardado <- "/Users/mersmac/Desktop/TFG/RESULTADOS/Curated_genes_diseases.ttl"
rdf_serialize(grafo_gen_enfermedad, doc = ruta_guardado, format = "turtle")
