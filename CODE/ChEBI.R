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

# 2. Seleccionar variables de interés de compuestos
nombres_chebi <- datos_compuestos %>%
  select(id, name, stars) %>%
  rename(ChEBI_ID = id, Chebi_Name = name)

# 3. Sustituir el ID del rol por su nombre
nombres_roles <- datos_compuestos %>%
  select(id, name) %>%
  rename(Role_ID = id, Role_Name = name)


# 4. Seleccionar variables de interés de relaciones y filtrar por "has_role"
roles_asignados <- datos_relaciones %>%
  filter(relation_type_id == 4) %>%
  select(init_id, final_id) %>%
  rename(ChEBI_ID = init_id, Role_ID = final_id)


# Unir 3 y 4
chebi_roles_final <- roles_asignados %>%
  inner_join(nombres_roles, by = "Role_ID") %>%
  group_by(ChEBI_ID) %>%
  summarise(Chemical_Role = paste(Role_Name, collapse = " | "))

# 3. Cruzar archivos limpios
datos_chebi <- datos_estructuras %>%
  select(compound_id, standard_inchi_key) %>%
  rename(ChEBI_ID = compound_id, InChIKey = standard_inchi_key) %>%
  inner_join(nombres_chebi, by = "ChEBI_ID") %>%
  left_join(chebi_roles_final, by = "ChEBI_ID") %>%
  mutate(across(everything(), as.character)) %>%
  filter(!is.na(InChIKey)) %>%
  distinct()


#--------- SERIALIZAR ----------

# 1. Inicializar grafo vacío
grafo_chebi <- rdf()

# 2. Prefijos Oficiales
chebi_prefix <- "http://purl.obolibrary.org/obo/CHEBI_"
inchikey_prefix <- "https://identifiers.org/inchikey/"
owl <- "http://www.w3.org/2002/07/owl#"
rdf_type <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
sio <- "http://semanticscience.org/resource/"

for (i in 1:nrow(datos_chebi)) {
  
  uri_chebi <- paste0(chebi_prefix, datos_chebi$ChEBI_ID[i])
  uri_inchi <- paste0(inchikey_prefix, sub("InChIKey=", "", datos_chebi$InChIKey[i]))
  
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
          oredicate = paste0(rdfs, "label"), 
          object = datos_chebi$Chebi_Name[i])
  
  # C. Añadir los Roles como literales independientes
  if (!is.na(datos_chebi$Chemical_Role[i])) {
    
    # 1. Separar texto usando la barra como delimitador
    roles_separados <- unlist(strsplit(datos_chebi$Chemical_Role[i], " \\| "))
    
    # 2. Bucle para crear una tripleta por cada rol
    for (rol in roles_separados) {
      rdf_add(grafo_chebi, 
              subject = uri_chebi, 
              predicate = paste0(sio, "SIO_000008"), # Propiedad 'has attribute'
              object = trimws(rol))
    }
  }
}
  
  # Mensaje de progreso
  if (i %% 10000 == 0) {
    message(paste("Procesadas", i, "filas de ChEBI"))
  }
}

# 5. Guardar el archivo final en la carpeta de resultados
rdf_serialize(grafo_chebi, doc = "../RESULTADOS/chebi_DB.ttl", format = "turtle")
