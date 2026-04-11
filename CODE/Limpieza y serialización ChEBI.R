# ------------- LIMPIAR Y SERIALIZAR CHEBI -------------

#------ LIMPIAR -------

library(dplyr)
library(data.table)
library(tidyverse)
library(rdflib)

# 1. Leer el archivo de estructuras comprimido online
datos_chebi <- fread("https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/structures.tsv.gz", 
                     stringsAsFactors = FALSE)

# 2. Limpieza y Filtrado
datos_chebi <- datos_chebi %>%
  # Seleccionar columnas de interés
  select(compound_id, standard_inchi_key) %>%
  rename(ChEBI_ID = compound_id,
         InChIKey = standard_inchi_key) %>%
  # Convertir todo a texto 
  mutate(across(everything(), as.character)) %>%
  # Convertir espacios en blanco a NA reales
  mutate(across(everything(), ~na_if(., ""))) %>%
  # Que no haya filas vacías
  filter(!is.na(ChEBI_ID) & !is.na(InChIKey))



#--------- SERIALIZAR ----------

# 1. Inicializar grafo vacío
grafo_chebi <- rdf()

# 2. Prefijos Oficiales
chebi_prefix <- "http://purl.obolibrary.org/obo/CHEBI_"
inchikey_prefix <- "https://identifiers.org/inchikey/"
owl <- "http://www.w3.org/2002/07/owl#"
rdf_type <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
sio <- "http://semanticscience.org/resource/"

# 3. Bucle
for (i in 1:nrow(datos_chebi)) {
  
  # Extraer datos de la fila
  chebi_id <- datos_chebi$ChEBI_ID[i]
  inchi <- datos_chebi$InChIKey[i]
  
  # Escudo protector: Si falta algún dato vital, pasamos a la siguiente fila
  if (length(chebi_id) == 0 || is.na(chebi_id) || length(inchi) == 0 || is.na(inchi)) {
    next 
  }
  
  # URIs
  # Nota: En ChEBI, el COMPOUND_ID viene como un número (ej. 16236). 
  # Le concatenamos "CHEBI_" delante para que sea la URI oficial.
  uri_chebi <- paste0(chebi_prefix, chebi_id)
  
  # Limpiamos el InChIKey por si acaso el texto viniera con "InChIKey=" delante
  inchi_limpio <- sub("InChIKey=", "", inchi)
  uri_inchi <- paste0(inchikey_prefix, inchi_limpio)
  
  # --- TRIPLETAS ---
  
  # A. Crear el nodo de ChEBI (Es un Químico)
  rdf_add(grafo_chebi, 
          subject = uri_chebi, 
          predicate = rdf_type, 
          object = paste0(sio, "SIO_010004"))
  
  # B. La Integración (owl:sameAs apuntando al InChIKey universal)
  rdf_add(grafo_chebi, 
          subject = uri_chebi, 
          predicate = paste0(owl, "sameAs"), 
          object = uri_inchi)
  
  # Mensaje de progreso
  if (i %% 20000 == 0) {
    message(paste("Procesadas", i, "filas de ChEBI"))
  }
}

# 5. Guardar el archivo final
ruta_guardado <- "/Users/mersmac/Desktop/TFG/RESULTADOS/chebi_structures.ttl"
rdf_serialize(grafo_chebi, doc = ruta_guardado, format = "turtle")
