# Describir los géneros en las muestras clínicas.
# No se repiten gráficas exploratorias por PM10 ni análisis inferenciales.

suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(readr); library(tibble)})
raiz_codigo <- Sys.getenv("TFM_CODE_ROOT")
source(file.path(raiz_codigo, "01_analisis/modulos/taxonomia.R"))
ps <- readRDS("outputs/analisis_microbiota/objetos_phyloseq/ps_saliva_206_sin_filtrado.rds")
conteos <- agregar_rango(ps, "Genus")
res_genus <- resumen_taxones(conteos)
carpeta <- "outputs/analisis_microbiota/composicion_taxonomica_v3/tablas"
dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)
write_csv(res_genus, file.path(carpeta, "03_resumen_taxonomico_genus.csv"), na = "")
prev_generos <- res_genus %>% mutate(
  presente_10pct = prevalencia_prop >= 0.10,
  presente_25pct = prevalencia_prop >= 0.25,
  presente_50pct = prevalencia_prop >= 0.50,
  presente_75pct = prevalencia_prop >= 0.75)
write_csv(prev_generos, file.path(carpeta, "13_prevalencia_generos_saliva206.csv"), na = "")
