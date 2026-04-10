# ------------- LIMPIAR Y SERIALIZAR CTD_diseases -------------

# ------------- LIMPIEZA --------------

# 1. Instalar paquetes y cargarlos
install.packages(c("rdflib", "dplyr"))
library(rdflib)
library(dplyr)

datos_enfermedades <- read.csv("TFG/DB/CTD_diseases.csv", skip = 27, stringsAsFactors = FALSE)

# 2. Corregir nombre de 1ª columna "X...DiseaseName" por "DiseaseName"
colnames(datos_enfermedades)[1] <- "DiseaseName"

# 3. Crear nueva tabla con columnas de interés
enfermedades_limpio <- datos_enfermedades %>%
  select(DiseaseName, DiseaseID)

# 4. Excluir fila con # y el nodo raíz "MESH:C"
enfermedades_limpio <- enfermedades_limpio %>%
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

# 3. Hacer un bucle para toda la tabla
for (i in 1:nrow(enfermedades_limpio)) {
  
  # Extraer datos de fila actual
  DiseaseName <- enfermedades_limpio$DiseaseName[i]
  ID_inicial <- enfermedades_limpio$DiseaseID[i] # Ej: "MESH:C538288"
  
  # Sustituir dos puntos por guion bajo para que sea una URL válida
  ID_formateado <- sub(":", "_", ID_inicial)
  
  # Creamos la URI final de esta enfermedad específica
  uri_enfermedad <- paste0(obo, ID_formateado)
  
  # --- AÑADIR LOS TRIPLES AL GRAFO ---
  
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


# Guardar grafo en la carpeta "RESULTADOS" de tu proyecto
rdf_serialize(mi_grafo, doc = "RESULTADOS/diseases_ctd.ttl", format = "turtle")
