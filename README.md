# Grafo de Conocimiento Semántico: Impacto Químico en la Salud

Este repositorio contiene el código fuente, los scripts de procesamiento y las consultas SPARQL asociados al Trabajo de Fin de Grado (TFG) para la construcción de un Grafo de Conocimiento Semántico (Knowledge Graph). 

El proyecto extrae, limpia y serializa datos de la *Comparative Toxicogenomics Database* (CTD) para modelar computacionalmente los efectos de diversos compuestos químicos en la salud, poniendo especial énfasis en el análisis de fármacos quimioterápicos y agentes dopantes.

## Resumen y Objetivos
El objetivo principal de este trabajo es transformar datos tabulares masivos en una red semántica interconectada, utilizando ontologías biomédicas estándar (OMIM, SIO, BAO, Biolink, OBO...). 
El sistema diseñado permite:
- **Modelado Toxicológico (Bottom-Up):** Rastrear cómo quimioterápicos y sustancias dopantes alteran la expresión genética y afectan a rutas biológicas concretas.
- **Farmacología Inversa (Top-Down):** Identificar compuestos químicos (xenobióticos) ocultos a partir de fenotipos o alteraciones clínicas específicas.
- **Inferencia Diagnóstica (CDSS):** Emular el razonamiento clínico identificando posibles patologías a partir del peso genético de un conjunto de síntomas (midiendo el efecto pleiotrópico).

## Estructura del Repositorio

- `/CODE/`: Scripts en R para la extracción, limpieza (`dplyr`, `data.table`) y serialización a tripletas RDF (`rdflib`). Contiene los módulos secuenciales para procesar químicos, enfermedades, genes y sus interacciones.
- `/SPARQL/`: Código de las queries listas para ejecutar en GraphDB.
- `/RESULTADOS/`: Directorio de salida donde se alojan los archivos `.ttl` generados, listos para su importación al motor de grafos.
- `/DIAGRAMA CONCEPTUAL/`: Documentación visual y técnica de la red.
    - `Diagrama Conceptual.png`: Representación gráfica de la arquitectura del grafo.
    - `Prefijos Oontologias.xlsx`: Diccionario de prefijos y namespaces utilizados (MESH, OMIM, SIO, etc.).
    - `README.md`: Explicación detallada de la lógica de relaciones y la jerarquía de las ontologías.

## Guía de Ejecución

### 1. Pre-requisitos
El preprocesamiento de datos está programado en **R**. Se requieren las siguientes librerías:
``R install.packages(c("dplyr", "readr", "data.table", "tidyverse", "rdflib"))´´

### 2. Generación del Grafo
Los scripts descargan automáticamente los datos de CTD. Al ejecutar los scripts de /CODE/, se realiza el filtrado de duplicados y la estandarización de URIs, generando el grafo en formato Turtle.

### 3. Exploración en GraphDB
Cargar los archivos .ttl en un repositorio de Ontotext GraphDB y utilizar las archivos de la carpeta /SPARQL/ para interrogar al grafo y reproducir las consultas del trabajo.
















