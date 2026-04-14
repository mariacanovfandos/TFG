Los prefijos de las ontologías usadas son los siguientes:

sio --> http://semanticscience.org/resource/
bao--> http://www.bioassayontology.org/bao#
biolink --> https://w3id.org/biolink/vocab/
obo--> http://purl.obolibrary.org/obo/
rdfs--> http://www.w3.org/2000/01/rdf-schema#
owl--> http://www.w3.org/2002/07/owl#
inchikey --> https://identifiers.org/inchikey:
skos --> http://www.w3.org/2004/02/skos/core#

Explicación del diagrama
La figura muestra el diagrama conceptual final, que organiza la información en cinco dominios principales (“chemical entity”, “pathway”, “gene”, “protein” y “disease”), representados como clases en color verde. Cada uno incluye atributos (en violeta) y se conecta mediante propiedades (en negro) y nodos intermedios (en naranja). Este modelo permite representar información clave sobre las sustancias químicas, como sus dianas moleculares, sus efectos sobre ellas y su relación con distintas enfermedades.
Delante de cada clase y propiedad se observan los prefijos, componentes esenciales de las URIs. Estos prefijos simplifican la escritura y lectura de las propiedades y clases en el diagrama, asegurando su claridad, facilidad de lectura y legibilidad. Así, las URIs extensas como http://semanticscience.org/resource/SIO_010004, que representa la clase de “chemical entity”, se pueden abreviar usando el prefijo sio para representar http://semanticscience.org/resource/ y SIO_010004 para la clase de “chemical entity”, de forma que el nodo principal se representa como sio: SIO_010004.
Para la creación de este diagrama se han reutilizado términos procedentes de distintas ontologías, destacando especialmente OBO y SIO. SIO representa una ontología concebida para describir objetos, atributos y procesos en diversos dominios científicos, facilitando una representación de los datos científicos y biomédicos consistente y ampliable. Por otro lado, OBO constituye un conjunto de ontologías orientadas al ámbito biomédico, dentro del cual se incluyen la Gene Ontology (GO) y la Disease Ontology (DO). Este sistema ofrece un marco estructurado y normalizado que permite la anotación y el intercambio de datos biológicos, favoreciendo así la interoperabilidad y la integración de información de múltiples fuentes.

