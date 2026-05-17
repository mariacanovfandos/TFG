# ------------- CHEBI -------------

# -------- LIMPIAR -------
# 1. Cargar librerías y archivos
library(dplyr)
library(data.table)
library(tidyverse)
library(rdflib)

datos_compuestos <- fread("https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/compounds.tsv.gz", 
                          stringsAsFactors = FALSE)
datos_relaciones <- fread("https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/relation.tsv.gz",
                          stringsAsFactors = FALSE)
datos_estructuras <- fread("https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/structures.tsv.gz",
                           stringsAsFactors = FALSE)

# 2. Procesar archivos
nombres_chebi <- datos_compuestos %>%
  select(id, name, stars) %>%
  rename(ChEBI_ID = id, Chebi_Name = name)

roles_asignados <- datos_relaciones %>%
  filter(relation_type_id == 4) %>%
  select(init_id, final_id) %>%
  rename(ChEBI_ID = init_id, Role_ID = final_id)

nombres_roles <- datos_compuestos %>%
  select(id, name) %>%
  rename(Role_ID = id, Role_Name = name)

chebi_roles_final <- roles_asignados %>%
  inner_join(nombres_roles, by = "Role_ID") %>%
  group_by(ChEBI_ID) %>%
  summarise(Chemical_Role = paste(Role_Name, collapse = " | "))

datos_chebi <- datos_estructuras %>%
  select(compound_id, standard_inchi_key) %>%
  rename(ChEBI_ID = compound_id, InChIKey = standard_inchi_key) %>%
  inner_join(nombres_chebi, by = "ChEBI_ID") %>%
  left_join(chebi_roles_final, by = "ChEBI_ID") %>%
  mutate(across(everything(), as.character)) %>%
  filter(!is.na(InChIKey)) %>%
  distinct()


# --------- SERIALIZAR ----------

grafo_chebi <- rdf()

# Prefijos
chebi_prefix <- "http://purl.obolibrary.org/obo/CHEBI_"
inchikey_prefix <- "https://identifiers.org/inchikey/"
sio <- "http://semanticscience.org/resource/"
rdfs <- "http://www.w3.org/2000/01/rdf-schema#"
owl <- "http://www.w3.org/2002/07/owl#"
rdf_type <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"


for (i in 1:nrow(datos_chebi)) {
  
  uri_chebi <- paste0(chebi_prefix, 
                      datos_chebi$ChEBI_ID[i])
  uri_inchi <- paste0(inchikey_prefix, 
                      sub("InChIKey=", 
                          "", 
                          datos_chebi$InChIKey[i]))
  
  # A. Tipo y Equivalencia
  rdf_add(grafo_chebi, 
          subject = uri_chebi, 
          predicate = rdf_type, 
          object = paste0(sio, "SIO_010004"))
  
  rdf_add(grafo_chebi, 
          subject = uri_chebi, 
          predicate = paste0(owl, "sameAs"), 
          object = uri_inchi)
  
  # B. Nombre de la molécula
  rdf_add(grafo_chebi, 
          subject = uri_chebi, 
          predicate = paste0(rdfs, "label"), 
          object = datos_chebi$Chebi_Name[i])
  
  # C. Añadir los Roles como literales independientes
  if (!is.na(datos_chebi$Chemical_Role[i])) {
    roles_separados <- unlist(strsplit(datos_chebi$Chemical_Role[i], " \\| "))
    for (rol in roles_separados) {
      rdf_add(grafo_chebi, 
              subject = uri_chebi, 
              predicate = paste0(sio, "SIO_000008"), # Propiedad 'has attribute'
              object = trimws(rol))
    }
  }
  
  # Avance
  if (i %% 10000 == 0) {
    message(paste("Procesadas", i, "filas de ChEBI"))
  }
}

# 5. Guardar el archivo final en la carpeta de resultados
rdf_serialize(grafo_chebi, doc = "../RESULTADOS/ChEBI_DB.ttl", format = "turtle")
