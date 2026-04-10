# ------------- LIMPIAR Y SERIALIZAR CTD_curated_genes_diseases -------------

# ------------- LIMPIEZA --------------
# 1. Cargar librería
library(dplyr)
library(readr)
library(rdflib)

# 2. Cargar el archivo saltando las líneas
datos_gen_enfermedad <- read.csv("TFG/DB/CTD_curated_genes_diseases.csv", skip = 27)

gen_enfermedad_limpio <- datos_gen_enfermedad %>%
  rename(GeneSymbol = 1) %>%
  # Seleccionar columnas de interés
  select(GeneSymbol, DiseaseID, DirectEvidence, PubMedIDs) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~na_if(., ""))) %>%
  slice(-1) %>%
  # Comprobar que hay GeneSymbol
  filter(!is.na(GeneSymbol) & !is.na(DiseaseID))

head(gen_enfermedad_limpio)


# ----------- SERIALIZACIÓN ------------

# 1. Inicializar el grafo de asociaciones
grafo_gen_enfermedad <- rdf()

# 2. Definir prefijos exactos
rdf_ns <- "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
biolink <- "https://w3id.org/biolink/vocab/"
sio <- "http://semanticscience.org/resource/"
schema <- "https://schema.org/"
mesh_prefix <- "http://id.nlm.nih.gov/mesh/" # URI Oficial de Enfermedades

# Prefijos para crear URIs nuevas
pubmed_base <- "https://pubmed.ncbi.nlm.nih.gov/"
assoc_base <- "http://tfg.org/asociacion_ge/"
gene_base <- "http://rdf.biogateway.eu/gene/9606/"

# 3. Bucle
for (i in 1:nrow(gen_enfermedad_limpio)) {
  
  # --- EXTRAEMOS LOS DATOS ---
  gen_symbol <- gen_enfermedad_limpio$GeneSymbol[i] 
  enfermedad_id <- gen_enfermedad_limpio$DiseaseID[i] 
  evidencia <- gen_enfermedad_limpio$DirectEvidence[i] 
  pubmed <- gen_enfermedad_limpio$PubMedIDs[i] 
  
  # --- CREAMOS LAS URIs DE ESTA FILA ---
  uri_asociacion <- paste0(assoc_base, i)
  uri_gen <- paste0(gene_base, gen_symbol)
  
  # Limpiamos MESH:
  id_enf_limpio <- sub("MESH:", "", enfermedad_id)
  uri_enfermedad <- paste0(mesh_prefix, id_enf_limpio)
  
  # --- CONSTRUIR TRIPLETAS ---
  
  # A. El Sujeto es el Gen
  rdf_add(grafo_gen_enfermedad, 
          subject = uri_asociacion, 
          predicate = paste0(rdf_ns, "subject"), 
          object = uri_gen)
  
  # B. El Objeto es la Enfermedad
  rdf_add(grafo_gen_enfermedad, 
          subject = uri_asociacion, 
          predicate = paste0(rdf_ns, "object"), 
          object = uri_enfermedad)
  
  # C. Evidencia Directa (si existe)
  if (!is.na(evidencia) && evidencia != "") {
    rdf_add(grafo_gen_enfermedad, 
            subject = uri_asociacion, 
            predicate = paste0(biolink, "has_evidence"), 
            object = evidencia)
  }
  
  # D. PubMed (Con separación automática de la barra vertical '|')
  if (!is.na(pubmed) && pubmed != "") {
    
    lista_pubmeds <- strsplit(pubmed, "\\|")[[1]]
    
    for (pmid in lista_pubmeds) {
      uri_articulo <- paste0(pubmed_base, pmid)
      
      # Conectar asociación con artículo
      rdf_add(grafo_gen_enfermedad, 
              subject = uri_asociacion, 
              predicate = paste0(sio, "SIO_000772"), 
              object = uri_articulo)
      
      # Definir el tipo de artículo
      rdf_add(grafo_gen_enfermedad, 
              subject = uri_articulo, 
              predicate = paste0(rdf_ns, "type"), 
              object = paste0(sio, "SIO_000154"))
      
      # Identificador del artículo
      rdf_add(grafo_gen_enfermedad, 
              subject = uri_articulo, 
              predicate = paste0(schema, "identifier"), 
              object = pmid)
    }
  }
}

# 4. Guardar el grafo definitivo en un archivo
rdf_serialize(grafo_gen_enfermedad, doc = "RESULTADOS/gene_disease_association_ctd.ttl", format = "turtle")
