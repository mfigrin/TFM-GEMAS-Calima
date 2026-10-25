# Utilidades de ejecución: solo rutas explícitas, no busca datos reales recursivamente.
tfm_stage_environment <- function(root, work) {
  destination <- file.path(work, "data_raw", "CALIMA")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  src <- file.path(root, "04_datos_ambientales", "datos", "raw", "red_canaria")
  for (year in 2016:2020) {
    source <- file.path(src, paste0("Datos ", year, ".xlsx"))
    if (!file.exists(source)) stop("Falta original ambiental: ", source)
    if (!file.copy(source, file.path(destination, basename(source)), overwrite = FALSE)) stop("No se copió: ", basename(source))
  }
  md <- file.path(destination, "CALIMA MITECO 2016-2020")
  dir.create(md)
  books <- list.files(file.path(root, "04_datos_ambientales/datos/raw/miteco"), pattern = "[.]xls[x]?$", full.names = TRUE)
  if (length(books) != 5L) stop("Se requieren los cinco libros MITECO.")
  if (!all(file.copy(books, md))) stop("No se copiaron todos los libros MITECO.")
  for (name in c("00_MITECO_Canarias_calendario_africanos_2016_2020.csv", "00b_MITECO_Canarias_rangos_oficiales_2016_2020.csv")) {
    if (!file.copy(file.path(root, "04_datos_ambientales/datos/processed/miteco", name), md)) stop("No se copió: ", name)
  }
  invisible(destination)
}
tfm_run_isolated <- function(root, work, scripts) {
  original_wd <- getwd(); on.exit(setwd(original_wd), add = TRUE)
  old <- Sys.getenv("TFM_PROJECT_DIR", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("TFM_PROJECT_DIR") else Sys.setenv(TFM_PROJECT_DIR = old), add = TRUE)
  Sys.setenv(TFM_PROJECT_DIR = work)
  setwd(work)
  dir.create("logs", showWarnings = FALSE)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- data.frame(script = character(), estado = character())
  for (f in scripts) {
    path <- file.path(root, f)
    if (!file.exists(path)) stop("Falta script: ", f)
    log <- file.path(work, "logs", paste0(basename(f), ".log"))
    message("Ejecutando: ", basename(f))
    code <- system2(rscript, c("--vanilla", shQuote(path)), stdout = log, stderr = log)
    status <- rbind(status, data.frame(script = f, estado = if (code == 0) "finalizado" else "error"))
    write.csv(status, file.path(work, "logs/estado_ejecucion.csv"), row.names = FALSE)
    if (code != 0) stop("Etapa fallida. Revisar ", log, ". No se continuará.")
  }
  message("Etapas terminadas. Revisar los registros internos; el código de salida no evalúa la validez científica.")
  invisible(status)
}
