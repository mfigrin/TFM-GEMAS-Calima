# Ejecución solo ambiental. 
local({
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "TFM-GEMAS-Calima.Rproj"))) stop("Ejecutar desde la raíz del repositorio.")
  work <- file.path(root, "environment_workspace")
  if (dir.exists(work)) stop("environment_workspace ya existe. Renombrarlo antes de repetir; no se sobrescribe.")
  source(file.path(root, "01_analisis/utilidades/ejecucion.R"), local = TRUE)
  dir.create(work)
  writeLines("SOLO SERIES AMBIENTALES; SIN DATOS DE MUESTRAS", file.path(work, "ENVIRONMENT_ONLY.txt"))
  tfm_stage_environment(root, work)
  tfm_run_isolated(root, work, c("04_datos_ambientales/procesar/01_normalizar_red_canaria.R",
    "04_datos_ambientales/procesar/02_extraer_miteco.R"))
  generated <- read.csv2(file.path(work, "data_raw/CALIMA/auditoria_estaciones/03_datos_diarios_normalizados.csv"), stringsAsFactors = FALSE)
  reference <- read.csv2(file.path(root, "04_datos_ambientales/datos/processed/red_canaria/03_datos_diarios_normalizados.csv"), stringsAsFactors = FALSE)
  generated <- generated[order(generated$estacion, generated$fecha), ]
  reference <- reference[order(reference$estacion, reference$fecha), ]
  rownames(generated) <- rownames(reference) <- NULL
  if (!identical(generated$estacion, reference$estacion) || !identical(generated$fecha, reference$fecha)) stop("Discrepancia de claves respecto a la auditoría independiente.")
  numeric_cols <- c("n_horas_pm10", "pm10_media_diaria", "pm10_max_horario", "n_horas_pm25", "pm25_media_diaria", "pm25_max_horario", "n_horas_no2", "no2_media_diaria", "n_horas_o3", "o3_media_diaria")
  ok <- vapply(numeric_cols, function(n) isTRUE(all.equal(generated[[n]], reference[[n]], tolerance = 1e-9, check.attributes = FALSE)), logical(1))
  write.csv(data.frame(variable = numeric_cols, coincide = ok), file.path(work, "logs/comparacion_referencia.csv"), row.names = FALSE)
  if (!all(ok)) stop("Algún valor no coincide con la auditoría independiente. Revisar logs/comparacion_referencia.csv.")
  capture.output(sessionInfo(), file = file.path(work, "logs/sessionInfo_control_ambiental.txt"))
  message("Control ambiental superado: originales, fechas y referencia diaria. No se han ajustado modelos clínicos.")
})
