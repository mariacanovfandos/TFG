# ------------- LIMPIAR Y SERIALIZAR CTD_chem_gene_interactions -------------

# ------------- LIMPIEZA --------------
# 1. Cargar librerías
library(dplyr)
library(readr)

# 2. Leer archivo
datos_quimico_gen <- read.csv("TFG/DB/CTD_chem_gene_interactions.csv", skip = 27)

# 3. Limpiar archivo
enfermedad_ruta_limpio <- datos_enfermedad_ruta %>%
  rename(DiseaseName = 1) %>%
  # Seleccionar columnas de interés
  select(DiseaseName, DiseaseID, PathwayName, PathwayID) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~na_if(., ""))) %>%
  slice(-1) %>%
  # Filtro de seguridad
  filter(!is.na(DiseaseID) & !is.na(PathwayID))


# ----------- SERIALIZACIÓN ------------

# 1. Guardar grafo
ruta_guardado_cg <- "RESULTADOS/chem_gene_association_ctd.ttl"
cat("", file = ruta_guardado_cg)

# 2. Prefijos
rdf_ns <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
rdfs <- "http://www.w3.org/2000/01/rdf-schema#"
sio <- "http://semanticscience.org/resource/"
schema <- "https://schema.org/"
bao <- "http://www.bioassayontology.org/bao#"

mesh_prefix <- "http://id.nlm.nih.gov/mesh/"
gene_base <- "http://rdf.biogateway.eu/gene/9606/" 
assoc_base <- "http://tfg.org/asociacion_cg/" 
pubmed_base <- "https://pubmed.ncbi.nlm.nih.gov/"

# 3. Bucle
for (i in 1:nrow(quimico_gen_limpio)) {
  
  # Extraer datos
  chem_id <- quimico_gen_limpio$ChemicalID[i]
  gene_symbol <- quimico_gen_limpio$GeneSymbol[i]
  interaccion_val <- quimico_gen_limpio$InteractionActions[i]
  pubmed <- quimico_gen_limpio$PubMedIDs[i]
  
  # Crear URIs
  uri_asociacion <- paste0("<", assoc_base, i, ">")
  uri_quimico <- paste0("<", mesh_prefix, sub("MESH:", "", chem_id), ">")
  uri_gen <- paste0("<", gene_base, gene_symbol, ">")
  
  # Vector temporal para guardar las líneas de texto de esta iteración
  lineas <- c()
  
  # Tripletas Estructurales
  lineas <- c(lineas, paste(uri_asociacion, 
                            paste0("<", rdf_ns, "type>"), 
                            paste0("<", sio, "SIO_001257> .")))
  
  lineas <- c(lineas, paste(uri_asociacion, 
                            paste0("<", rdfs, "subject>"), 
                            uri_quimico, "."))
  lineas <- c(lineas, paste(uri_quimico, 
                            paste0("<", rdf_ns, "type>"), 
                            paste0("<", sio, "SIO_010004> .")))
  
  lineas <- c(lineas, paste(uri_asociacion, 
                            paste0("<", rdfs, "object>"), 
                            uri_gen, "."))
  lineas <- c(lineas, paste(uri_gen, 
                            paste0("<", rdf_ns, "type>"), 
                            paste0("<", sio, "SIO_010035> .")))
  
  # Tripletas de InteractionActions (BAO)
  if (!is.na(interaccion_val) && interaccion_val != "") {
    lista_acciones <- strsplit(interaccion_val, "\\|")[[1]]
    for (accion in lista_acciones) {
      # Para literales (texto), se ponen comillas dobles en lugar de los brackets de las URIs
      lineas <- c(lineas, paste(uri_asociacion, 
                                paste0("<", bao, "BAO_0000196>"), 
                                paste0("\"", accion, "\" .")))
    }
  }
  
  # Tripletas de PubMed
  if (!is.na(pubmed) && pubmed != "") {
    lista_pubmeds <- strsplit(pubmed, "\\|")[[1]]
    for (pmid in lista_pubmeds) {
      uri_articulo <- paste0("<", pubmed_base, pmid, ">")
      
      lineas <- c(lineas, paste(uri_asociacion, 
                                paste0("<", sio, "SIO_000772>"), 
                                uri_articulo, "."))
      lineas <- c(lineas, paste(uri_articulo, 
                                paste0("<", rdf_ns, "type>"), 
                                paste0("<", sio, "SIO_000154> .")))
      lineas <- c(lineas, paste(uri_articulo, 
                                paste0("<", schema, "identifier>"), 
                                paste0("\"", pmid, "\" .")))
    }
  }
  
  # Escribir bloque de texto en el archivo
  cat(paste(lineas, collapse = "\n"), "\n", file = ruta_guardado_cg, append = TRUE)
  
  # Avance
  if (i %% 10000 == 0) {
    message(paste("Procesadas", i, "filas de", nrow(quimico_gen_limpio)))
  }
}
