# Construir el objeto phyloseq a partir de las tablas artificiales incluidas.
# La segunda parte deriva las mismas entradas que esperan los scripts del estudio.
local({
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (!file.exists("TFM-GEMAS-Calima.Rproj")) stop("Abrir el proyecto y ejecutar desde su raíz.")
  if (!requireNamespace("phyloseq", quietly = TRUE) || !requireNamespace("ape", quietly = TRUE))
    stop("Se necesitan phyloseq y ape.")
  datos <- file.path(root, "02_ejemplo_sintetico/datos")
  work <- file.path(root, "02_ejemplo_sintetico/trabajo")
  if (dir.exists(work)) stop("02_ejemplo_sintetico/trabajo ya existe; no se sobrescriben sus resultados.")
  dir.create(work, recursive = TRUE)
  writeLines("SOLO DATOS SINTETICOS. No contiene participantes reales.",
             file.path(work, "SYNTHETIC_ONLY.txt"))
  for (d in c("data_raw", "inputs_synthetic", "outputs/datos_procesados_finales",
              "outputs/analisis_microbiota/miteco_integrado/tablas", "logs"))
    dir.create(file.path(work, d), recursive = TRUE, showWarnings = FALSE)
  meta <- read.csv(file.path(datos, "metadatos.csv"), row.names = 1,
                   check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  meta$Date <- as.Date(meta$Date)
  meta$Pool <- as.character(meta$Pool)
  meta$synthetic <- tolower(as.character(meta$synthetic)) %in% c("true", "t", "1")
  if (!all(meta$synthetic) || !all(grepl("^SYN_", rownames(meta)))) stop("Entrada no sintética.")
  counts <- as.matrix(read.csv(file.path(datos, "conteos.csv"), row.names = 1, check.names = FALSE))
  storage.mode(counts) <- "integer"
  taxonomy <- as.matrix(read.csv(file.path(datos, "taxonomia.csv"), row.names = 1, check.names = FALSE))
  tree <- ape::read.tree(file.path(datos, "arbol.nwk"))
  stopifnot(identical(colnames(counts), rownames(meta)),
            identical(rownames(counts), rownames(taxonomy)))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
                          phyloseq::tax_table(taxonomy), phyloseq::sample_data(meta), tree)
  n_saliva <- 206L
  stopifnot(methods::validObject(ps),
            all(phyloseq::sample_sums(ps)[seq_len(n_saliva)] >= 10000))
  saveRDS(ps, file.path(work, "phyloseq_sintetico.rds"), compress = "xz")
  saveRDS(ps, file.path(work, "data_raw/data_phyloseq_salivas_GEMAS_Canarias.rds"), compress = "xz")
  daily <- read.csv(file.path(datos, "ambiente_diario.csv"), stringsAsFactors = FALSE,
                    na.strings = c("", "NA"))
  daily$fecha <- as.Date(daily$fecha)
  if (!all(tolower(as.character(daily$synthetic)) %in% c("true", "t", "1"))) stop("El ambiente de demostración debe ser sintético.")
  clinical <- meta[seq_len(n_saliva), , drop = FALSE]
  clinical$id_muestra <- rownames(clinical)
  env_ids <- which(!is.na(clinical$Date))
  trace <- lapply(env_ids, function(j) {
    ds <- daily[daily$isla == clinical$Island[j], ]
    idx <- match(clinical$Date[j] - seq_len(180), ds$fecha)
    stopifnot(!anyNA(idx))
    x <- ds[idx, , drop = FALSE]
    x$id_muestra <- clinical$id_muestra[j]
    x$fecha_muestra <- clinical$Date[j]
    x$lag <- seq_len(180)
    x
  })
  trace <- do.call(rbind, trace)
  one_summary <- function(id) {
    x <- trace[trace$id_muestra == id, ]
    out <- list(id_muestra = id, fecha_muestra = as.character(x$fecha_muestra[1]),
                isla = x$isla[1], synthetic = TRUE)
    for (window in c(10L, 30L, 90L, 180L)) {
      z <- x[x$lag <= window, ]
      present <- z$miteco_africano == 1
      out[[paste0("miteco_any_", window, "d")]] <- as.integer(any(present))
      out[[paste0("miteco_days_", window, "d")]] <- sum(present)
      out[[paste0("miteco_fraction_", window, "d")]] <- mean(present)
      out[[paste0("miteco_episodes_", window, "d")]] <- length(unique(z$episodio_miteco_id[present]))
      out[[paste0("miteco_recency_", window, "d")]] <- if (any(present)) min(z$lag[present]) else NA_integer_
      out[[paste0("naf_mean_calendar_", window, "d")]] <- mean(z$aporte_naf_usable)
      out[[paste0("naf_sum_calendar_", window, "d")]] <- sum(z$aporte_naf_usable)
      out[[paste0("naf_max_", window, "d")]] <- max(z$aporte_naf_usable)
      out[[paste0("naf_mean_event_", window, "d")]] <- if (any(present)) mean(z$aporte_naf_usable[present]) else 0
      out[[paste0("naf_event_coverage_", window, "d")]] <- 1
    }
    as.data.frame(out, stringsAsFactors = FALSE)
  }
  exposed <- do.call(rbind, lapply(unique(trace$id_muestra), one_summary))
  pm_summary <- lapply(unique(trace$id_muestra), function(id) {
    x <- trace[trace$id_muestra == id, ]
    out <- list(id_muestra = id, fecha_muestra = as.character(x$fecha_muestra[1]), adulto = TRUE)
    for (window in c(10L, 180L)) {
      y <- x$PM10[x$lag <= window]
      out[[paste0("pm10_media_", window, "d")]] <- mean(y)
      out[[paste0("pm10_maximo_", window, "d")]] <- max(y)
      out[[paste0("pm10_acumulado_", window, "d")]] <- sum(y)
      out[[paste0("pm10_auc50_", window, "d")]] <- sum(pmax(y-50, 0))
      ndays <- sum(y >= 50)
      out[[paste0("n_dias_sdi_", window, "d")]] <- ndays
      out[[paste0("sdi_exposicion_", window, "d")]] <- ndays > 0
      out[[paste0("grado_sdi_", window, "d")]] <- if (ndays == 0) "Sin_exposicion" else if (ndays < 4) "Baja_1_3_dias" else "Alta_ge4_dias"
    }
    y <- x$PM10[x$lag <= 10]
    for (lag in seq_len(10)) out[[sprintf("pm10_d%02d", lag)]] <- y[lag]
    out$sdi_temprana_10d <- any(y[1:3] >= 50)
    out$sdi_tardia_10d <- any(y[4:10] >= 50)
    out$patron_temporal_sdi_10d <- if (out$sdi_temprana_10d && out$sdi_tardia_10d) "Temprana_y_tardia" else if (out$sdi_temprana_10d) "Solo_temprana_D1_D3" else if (out$sdi_tardia_10d) "Solo_tardia_D4_D10" else "Sin_SDI"
    month <- as.integer(format(x$fecha_muestra[1], "%m"))
    out$tipo_sdi_estacional_10d <- if (month %in% c(10,11,12,1,2,3)) "Invernal" else "Calido"
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  pm <- do.call(rbind, pm_summary)
  # Prevenir columnas duplicadas al unir PM10 (en phyloseq) con MITECO (externo).
  p <- file.path(work, "outputs/analisis_microbiota/miteco_integrado/tablas")
  write.csv(exposed, file.path(p, "04_exposicion_MITECO_PM10_por_muestra.csv"), row.names = FALSE, na = "")
  write.csv(trace, file.path(p, "02_trazabilidad_MITECO_PM10_D1_D180.csv"), row.names = FALSE, na = "")
  write.csv2(pm, file.path(work, "outputs/17_macrotabla_clinica_ambiental_por_muestra.csv"), row.names = FALSE, na = "")
  cohort <- data.frame(id_muestra = clinical$id_muestra, edad = clinical$Age,
    sexo = clinical$Sex.letter, pool = clinical$Pool, caries = clinical$Cavities,
    chicle_30min = clinical$Chewing.gum.30min, antibiotico_2m = clinical$Antibiotic.2mo,
    exacerbacion_6m = clinical$Exacerbations.6M, synthetic = TRUE)
  write.csv(cohort, file.path(work, "outputs/datos_procesados_finales/TABLA_MAESTRA_COHORTE_206.csv"), row.names = FALSE)
  writeLines(c("fixture=datos_portables_v03", "n_saliva=206", "n_ambiental=205", "n_controles=12", "n_ASV=320",
    "Tablas generadas con NumPy PCG64, semilla 20261024; ver 02_ejemplo_sintetico/generar_tablas.py.",
    "RDS construido con phyloseq::phyloseq; no se ha leido ningun objeto real."),
    file.path(work, "inputs_synthetic/PROVENANCE.txt"))
  capture.output(sessionInfo(), file = file.path(work, "logs/sessionInfo_generador.txt"))
  message("Creado: ", file.path(work, "phyloseq_sintetico.rds"))
  message("Objeto sintético preparado. No contiene datos de GEMAS.")
})
