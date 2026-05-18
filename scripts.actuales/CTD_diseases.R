# ------------- LIMPIAR Y SERIALIZAR CTD_diseases -------------

# ------------- LIMPIEZA --------------
# 1. Cargar paquetes
library(rdflib)
library(dplyr)
library(data.table)
library(tidyverse)

# 2. Leer archivo
datos_enfermedades <- fread("https://ctdbase.org/reports/CTD_diseases.csv.gz", 
                            skip = 27, stringsAsFactors = FALSE)

# 3. Limpiar archivo
datos_enfermedades <- datos_enfermedades %>%
  set_names(c("DiseaseName", "DiseaseID", "AltDiseaseIDs", "Definition", "ParentIDs",
              "TreeNumbers", "ParentTreeNumbers", "Synonyms", "SlimMappings")) %>%
  select(DiseaseName, DiseaseID, Definition) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~na_if(., ""))) %>%
  filter(!is.na(DiseaseID)) %>%
  distinct()


# ----------- SERIALIZACIÓN ------------

# 1. Iniciar grafo vacío
mi_grafo <- rdf()

# 2. Definir los prefijos estandarizados
sio <- "http://semanticscience.org/resource/"
rdfs <- "http://www.w3.org/2000/01/rdf-schema#"
rdf_base <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
obo <- "http://purl.obolibrary.org/obo/"
schema <- "https://schema.org/"
mesh_prefix <- "https://id.nlm.nih.gov/mesh/" 
omim_prefix <- "http://purl.bioontology.org/ontology/OMIM/"  # MANTENIDO COMO PEDISTE
skos <- "http://www.w3.org/2004/02/skos/core#"

# 3. Bucle
for (i in 1:nrow(datos_enfermedades)) {
  
  # Extraer datos de fila actual
  DiseaseName <- datos_enfermedades$DiseaseName[i]
  ID_inicial <- datos_enfermedades$DiseaseID[i]
  Definicion <- datos_enfermedades$Definition[i]
  
  # Crear URI según la base de datos de origen y limpiar ID
  if (grepl("^MESH:", ID_inicial)) {
    id_limpio <- sub("MESH:", "", ID_inicial)
    uri_enfermedad <- paste0(mesh_prefix, id_limpio)
  } else if (grepl("^OMIM:", ID_inicial)) {
    id_limpio <- sub("OMIM:", "", ID_inicial)
    uri_enfermedad <- paste0(omim_prefix, id_limpio)
  } else {
    id_limpio <- ID_inicial # Fallback de seguridad
    uri_enfermedad <- paste0("http://unknown.org/disease/", id_limpio)
  }
  
  # --- TRIPLETAS ---
  
  # Tripleta A: nodo es de tipo "Enfermedad" (SIO_010299)
  rdf_add(mi_grafo, 
          subject = uri_enfermedad, 
          predicate = paste0(rdf_base, "type"), 
          object = paste0(sio, "SIO_010299"))
  
  # Tripleta B: asignar su nombre en texto normal
  rdf_add(mi_grafo, 
          subject = uri_enfermedad, 
          predicate = paste0(rdfs, "label"), 
          object = DiseaseName)
  
  # Tripleta C: asignar su identificador (usando el ID limpio)
  rdf_add(mi_grafo, 
          subject = uri_enfermedad, 
          predicate = paste0(schema, "identifier"), 
          object = id_limpio)
  
  # Tripleta D: Definition (SKOS)
  if (!is.na(Definicion)) {
    rdf_add(mi_grafo, 
            subject = uri_enfermedad, 
            predicate = paste0(skos, "definition"), 
            object = Definicion)
  }
  
  # Avance
  if (i %% 10000 == 0) {
    message(paste("Procesadas", i, "filas de enfermedades de CTD"))
  }
}

# 5. Guardar el grafo definitivo en un archivo
ruta_guardado <- "/Users/mersmac/Desktop/TFG/RESULTADOS/Diseases_ctd.ttl"
rdf_serialize(mi_grafo, doc = ruta_guardado, format = "turtle")
