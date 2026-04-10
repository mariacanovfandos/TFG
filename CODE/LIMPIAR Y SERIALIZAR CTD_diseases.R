# ------------- LIMPIAR Y SERIALIZAR CTD_diseases -------------

# ------------- LIMPIEZA --------------

# 1. Instalar paquetes y cargarlos
install.packages(c("rdflib", "dplyr"))
library(rdflib)
library(dplyr)
library(data.table)
library(tidyverse)

datos_enfermedades <- fread("https://ctdbase.org/reports/CTD_diseases.csv.gz", skip = 27, stringsAsFactors = FALSE)

# 2. Limpiar archivo
datos_enfermedades <- datos_enfermedades %>%
  set_names(c("DiseaseName", "DiseaseID", "AltDiseaseIDs", "Definition", "ParentIDs",
              "TreeNumbers", "ParentTreeNumbers", "Synonyms", "SlimMappings")) %>%
  select(DiseaseName, DiseaseID) %>%
  filter(DiseaseName != "#") %>%          # Elimina la fila del símbolo #
  filter(DiseaseID != "MESH:C")           # Elimina la categoría general "Diseases"


# ----------- SERIALIZACIÓN ------------

# 1. Iniciar grafo vacío
mi_grafo <- rdf()

# 2. Definir los prefijos
sio <- "http://semanticscience.org/resource/"
rdfs <- "http://www.w3.org/2000/01/rdf-schema#"
rdf_type <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
obo <- "http://purl.obolibrary.org/obo/"
schema = "https://schema.org/"

# 3. Bucle
for (i in 1:nrow(datos_enfermedades)) {
  
  # Extraer datos de fila actual
  DiseaseName <- datos_enfermedades$DiseaseName[i]
  ID_inicial <- datos_enfermedades$DiseaseID[i] # Ej: "MESH:C538288"
  
  # Sustituir dos puntos por guion bajo para que sea una URL válida
  ID_formateado <- sub(":", "_", ID_inicial)
  
  # Crear la URI final de esta enfermedad específica
  uri_enfermedad <- paste0(obo, ID_formateado)
  
  # --- TRIPLETAS ---
  
  # Tripleta A: Declarar que este nodo es de tipo "Enfermedad" (SIO_010299)
  rdf_add(mi_grafo, 
          subject = uri_enfermedad, 
          predicate = rdf_type, 
          object = paste0(sio, "SIO_010299"))
  
  # Tripleta B: Asignar su nombre en texto normal
  rdf_add(mi_grafo, 
          subject = uri_enfermedad, 
          predicate = paste0(rdfs, "label"), 
          object = DiseaseName)
  
  # Tripleta C: schema:identifier -> DiseaseID (asignar su identificador)
  rdf_add(mi_grafo, 
          subject = uri_enfermedad, 
          predicate = paste0(schema, "identifier"), 
          object = ID_inicial)
}


# Guardar grafo en la carpeta "RESULTADOS"
rdf_serialize(mi_grafo, doc = "TFG/RESULTADOS/diseases_ctd.ttl", format = "turtle")
