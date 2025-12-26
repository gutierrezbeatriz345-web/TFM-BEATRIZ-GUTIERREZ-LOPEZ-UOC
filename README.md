# TFM – Análisis transcriptómico del nervio ciático

Este repositorio contiene todos los scripts que he utilizado para el análisis transcriptómico de nervio ciático lesionado, integrando datos de secuenciación unicelular (scRNA-seq) y RNA-seq bulk.

El objetivo del trabajo es caracterizar cómo el tipo de lesión (aplastamiento *crush* vs sección *cut*) y la edad (jóvenes vs envejecidos) afectan la respuesta inmunitaria y glial en el nervio periférico, así como estudiar el papel de las células T reguladoras (Treg) en la inflamación y regeneración.

He organizado los códigos en carpetas según la función de cada análisis, para que sea más fácil de seguir y reproducir.

## Estructura del repositorio

- **data_preparation/**  
  Scripts para preparación de datos brutos (descarga de datos SRA y cuantificación con Salmon).

- **scRNAseq/**  
  Scripts de análisis de datos single-cell (scRNA-seq). Incluye el procesamiento de datos celulares individuales, generación de UMAPs, anotación de tipos celulares y renombrado de clusters.

- **bulk_RNAseq/**  
  Scripts de análisis de datos bulk RNA-seq. Contiene el pipeline de análisis diferencial con DESeq2, visualizaciones (volcanos, heatmaps, PCA), identificación de genes housekeeping y análisis funcional por módulos, así como un script separado para enriquecimiento funcional (Gene Ontology y KEGG).

## Descripción de los scripts

- **data_preparation/descarga_y_cuantificacion.ps1**  
  *(PowerShell)*  
  Descarga los datos crudos de secuenciación desde el repositorio (archivos SRA), los convierte a formato FASTQ usando SRA Toolkit (`prefetch` y `fasterq-dump`) y ejecuta Salmon (en entorno WSL Ubuntu) para cuantificar la expresión génica. Automatiza la obtención de los archivos `quant.sf`.

- **scRNAseq/01_umap_por_dia.R**  
  *(R)*  
  Carga objetos de expresión unicelular previamente procesados (`.rds`) para un experimento de lesión por aplastamiento en ratones jóvenes. Aplica control de calidad, normalización (SCTransform), reducción de dimensionalidad (PCA), clustering y genera UMAPs por día post-lesión (0, 3 y 7).

- **scRNAseq/02_anotacion_y_marcadores.R**  
  *(R)*  
  Anota automáticamente los clusters celulares usando SingleR y extrae marcadores genéticos por cluster. Genera UMAPs coloreados por tipo celular, tablas de marcadores (incluyendo top 5 por cluster) y un ejemplo de análisis diferencial entre clusters con volcán.

- **scRNAseq/03_renombrar_clusters_y_umaps.R**  
  *(R)*  
  Renombra manualmente los clusters con etiquetas biológicas descriptivas y genera UMAPs finales por día post-lesión, así como un panel combinado.

- **bulk_RNAseq/01_analisis_bulk_RNAseq.R**  
  *(R)*  
  Script principal para el análisis de RNA-seq bulk. Realiza control de calidad, normalización y análisis de expresión diferencial con DESeq2. Incluye PCA, heatmaps, identificación automática de genes housekeeping, definición de módulos funcionales (inmunidad, senescencia, Treg, rutas de señalización) y cálculo de scores por muestra. Genera tablas de resultados y un archivo Excel consolidado.

- **bulk_RNAseq/02_enriquecimiento_GO_KEGG.R**  
  *(R)*  
  Realiza análisis de enriquecimiento funcional GO (Biological Process) y KEGG a partir de los genes diferencialmente expresados de cada módulo funcional y contraste experimental. Guarda tablas y gráficos de los términos enriquecidos.

## Notas importantes

- Los datos brutos y resultados voluminosos no están incluidos en este repositorio.  
- Los scripts asumen que los datos están disponibles localmente y que las rutas pueden necesitar ajuste.  
- Cada script incluye comentarios y mensajes para facilitar su comprensión.

Este repositorio recoge el código utilizado para mi Trabajo Final de Máster y está organizado para facilitar su revisión y reproducibilidad.
