# Guardar la tabla analítica utilizada por las etapas siguientes.

local({
  f <- "outputs/analisis_microbiota/reanalisis_completo_MITECO/tablas/01_dataset_metricas_MITECO_analitico.csv"
  if (!file.exists(f)) stop("Falta la salida del script 33: ", f)
  data <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
  needed <- c("id_muestra", "age", "sex", "island", "pool", "antibiotic", "cavities", "gum",
              "exac6m", "fecha", "year", "sin_doy", "cos_doy", "cluster_isla_fecha",
              "Observed", "Shannon", "Simpson", "aitchison_centro_isla", "pm10_media_10d")
  missing <- setdiff(needed, names(data))
  if (length(missing)) stop("No se crea maestra; faltan columnas: ", paste(missing, collapse = ", "))
  if (anyDuplicated(data$id_muestra)) stop("Hay identificadores duplicados.")
  out <- "outputs/datos_procesados_finales"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(out, "TABLA_MAESTRA_PACIENTES_205_COMPLETA.csv")
  if (file.exists(dest)) stop("La tabla maestra ya existe. No se sobrescribe: ", dest)
  write.csv(data, dest, row.names = FALSE, na = "")
  clinical_dest <- file.path(out, "TABLA_MAESTRA_COHORTE_206.csv")
  if (!file.exists(clinical_dest)) {
    clinical_source <- "outputs/01_cohorte_microbiota_para_exposicion_ambiental.csv"
    if (!file.exists(clinical_source)) stop("Falta la salida de preparación de la cohorte.")
    clinical <- read.csv(clinical_source, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("id_muestra", "edad", "sexo", "pool", "caries", "chicle_30min", "antibiotico_2m", "exacerbacion_6m")
    if (!all(required %in% names(clinical))) stop("La tabla clínica no tiene las columnas requeridas por DESeq2.")
    write.csv(clinical, clinical_dest, row.names = FALSE, na = "")
  }
  message("Tabla maestra creada SIN modificar columnas o estimaciones: ", dest)
})
