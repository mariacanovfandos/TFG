# ------------- LIMPIAR Y SERIALIZAR CTD_diseases -------------

# ------------- LIMPIEZA --------------

# 1. Instalar paquetes y cargarlos
install.packages(c("rdflib", "dplyr"))
library(rdflib)
library(dplyr)

datos_quimicos <- read.csv("TFG/DB/CTD_chemicals.csv", skip = 27, stringsAsFactors = FALSE)

quimicos_limpio <- datos_quimicos %>%
  rename(ChemicalName = 1) %>%
  # Seleccionar columnas de interés
  select(ChemicalName, ChemicalID, InChIKey, Definition) %>% 
  # Quitar la fila de separadores y el posible nodo raíz "MESH:D" (Chemicals)
  filter(ChemicalName != "#") %>% 
  filter(ChemicalID != "MESH:D") 

# ----------- SERIALIZACIÓN ------------

# 1. Iniciar grafo vacío
grafo_quimicos <- rdf()

# 2. Definir los prefijos
sio <- "http://semanticscience.org/resource/"
rdfs <- "http://www.w3.org/2000/01/rdf-schema#"
rdf_type <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
mesh_prefix <- "http://id.nlm.nih.gov/mesh/" 
schema <- "https://schema.org/"
owl <- "http://www.w3.org/2002/07/owl#"
inchikey_prefix <- "https://identifiers.org/inchikey/"
skos <- "http://www.w3.org/2004/02/skos/core#"

# 3. Bucle para procesar toda la tabla
for (i in 1:nrow(quimicos_limpio)) {
  
  ChemicalName <- quimicos_limpio$ChemicalName[i]
  ID_inicial <- quimicos_limpio$ChemicalID[i]  
  InChIKey <- quimicos_limpio$InChIKey[i]
  Definicion <- quimicos_limpio$Definition[i] 
  
  ID_solo_codigo <- sub("MESH:", "", ID_inicial)
  uri_quimico <- paste0(mesh_prefix, ID_solo_codigo)
  
  # --- TRIPLETAS ---
  
  # A. Tipo (Químico)
  rdf_add(grafo_quimicos, 
          subject = uri_quimico, 
          predicate = rdf_type, 
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
  
  # D. sameAs (InChIKey)
  if (!is.na(InChIKey)) {
    uri_inchi <- paste0(inchikey_prefix, InChIKey)
    rdf_add(grafo_quimicos, 
            subject = uri_quimico, 
            predicate = paste0(owl, "sameAs"), 
            object = uri_inchi)
  }
  
  # E. Definición (SKOS). Solo añadir tripleta si la molécula tiene una definición (no es NA)
  if (!is.na(Definicion)) {
    rdf_add(grafo_quimicos, 
            subject = uri_quimico, 
            predicate = paste0(skos, "definition"), 
            object = Definicion)
  }
}


# Guardar grafo en la carpeta "RESULTADOS" de tu proyecto
rdf_serialize(grafo_quimicos, doc = "TFG/RESULTADOS/chemicals_ctd.ttl", format = "turtle")
