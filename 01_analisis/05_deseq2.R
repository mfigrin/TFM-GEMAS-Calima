# DESeq2 clínico y ambiental como análisis complementario.

PROJECT_DIR <- normalizePath(Sys.getenv("TFM_PROJECT_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
BASE_DIR <- PROJECT_DIR

if (!dir.exists(PROJECT_DIR)) {
  stop("No se encuentra el directorio del proyecto: ", PROJECT_DIR)
}

# Resultados de este analisis.
OUT_DIR <- file.path(PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Objeto de entrada. La demo usa datos sintéticos bajo el nombre de archivo histórico.
PHYLOSEQ_RDS <- file.path(PROJECT_DIR, "data_raw", "data_phyloseq_salivas_GEMAS_Canarias.rds")

if (!file.exists(PHYLOSEQ_RDS)) {
  stop("No se encuentra el objeto phyloseq de entrada en: ", PHYLOSEQ_RDS)
}

# Covariables del modelo clinico conservado del analisis original del TFM.
COVARS_REPLICA <- c("age", "sex", "pool", "cavities", "gum")

# Ajuste clínico, demográfico y técnico. La estacionalidad se evalúa por separado.
COVARS_CALIMA_BASE <- c("age", "sex", "island", "pool", "antibiotic", "cavities", "gum")
EJECUTAR_SENSIBILIDAD_ESTACIONAL <- TRUE

# Filtros CMM conservados del análisis clínico: prevalencia y media incluidos los ceros.
CMM_PREVALENCIA <- 0.10  # estrictamente >10%
CMM_MEDIA_LECTURAS <- 40 # estrictamente >40 lecturas/muestra
FDR_ALPHA <- 0.05

# 1. Comprobar y cargar los paquetes.
required_pkgs <- c("phyloseq", "DESeq2", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Faltan paquetes: ", paste(missing_pkgs, collapse = ", "),
    "\nInstalarlos antes de ejecutar. Para Bioconductor: BiocManager::install(c('phyloseq','DESeq2'))."
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(ggplot2)
})

# 2. Funciones de preparación y ajuste.
msg <- function(...) cat(sprintf(...), "\n")

find_one_file <- function(pattern, base_dir = BASE_DIR, required = TRUE) {
  x <- list.files(base_dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(x) == 0) {
    if (required) stop("No se encontro archivo con patron: ", pattern)
    return(NA_character_)
  }
  if (length(x) > 1) {
    msg("AVISO: varios archivos coinciden con '%s'. Se usara: %s", pattern, x[1])
  }
  x[1]
}

load_phyloseq_safe <- function() {
# 1) Ruta indicada expresamente
  if (!is.na(PHYLOSEQ_RDS) && file.exists(PHYLOSEQ_RDS)) {
    msg("Cargando objeto phyloseq original desde: %s", PHYLOSEQ_RDS)
    obj <- readRDS(PHYLOSEQ_RDS)
    if (!inherits(obj, "phyloseq")) {
      stop("El archivo indicado en PHYLOSEQ_RDS existe, pero no contiene un objeto phyloseq.")
    }
    msg("Objeto phyloseq cargado correctamente: %d muestras, %d taxa.", nsamples(obj), ntaxa(obj))
    return(obj)
  }

# 2) Si ya hay objetos phyloseq cargados en el entorno, reutilizar el mas grande.
  obj_names <- ls(envir = .GlobalEnv)
  ps_names <- obj_names[vapply(obj_names, function(nm) {
    inherits(get(nm, envir = .GlobalEnv), "phyloseq")
  }, logical(1))]
  if (length(ps_names) > 0) {
    sizes <- vapply(ps_names, function(nm) nsamples(get(nm, envir = .GlobalEnv)), numeric(1))
    chosen <- ps_names[which.max(sizes)]
    msg("Se reutilizara el objeto phyloseq cargado '%s' (%d muestras).", chosen, max(sizes))
    return(get(chosen, envir = .GlobalEnv))
  }

# 3) Buscar candidatos RDS por nombre
  cand <- list.files(BASE_DIR, pattern = "(phyloseq|gemas|(^|_)ps(_|\\.|$)).*\\.rds$",
                     recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(cand) == 0) {
    stop(
      "No se localizo el objeto phyloseq. Guarda el objeto como .rds o indica su ruta en PHYLOSEQ_RDS.\n",
      "Ejemplo: PHYLOSEQ_RDS <- 'datos/phyloseq_GEMAS.rds'"
    )
  }
  if (length(cand) > 1) {
    msg("Se localizaron varios RDS candidatos. Se probara el primero: %s", cand[1])
  }
  obj <- readRDS(cand[1])
  if (!inherits(obj, "phyloseq")) {
    stop("El RDS localizado no es phyloseq. Indique manualmente PHYLOSEQ_RDS.")
  }
  obj
}

read_csv_flexible <- function(path) {
# read.csv funciona con nuestros CSV; check.names=FALSE conserva nombres.
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("", "NA", "NaN", "null", "NULL"))
}

standardize_id <- function(x) trimws(as.character(x))

get_count_matrix <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) m <- t(m)
  storage.mode(m) <- "numeric"
  m
}

subset_ps_ids <- function(ps, ids) {
  ids <- intersect(sample_names(ps), standardize_id(ids))
  prune_samples(ids, ps)
}

prepare_genus_counts <- function(ps) {
  if (!("Genus" %in% rank_names(ps))) {
    stop("El objeto phyloseq no contiene el rango taxonomico 'Genus'. Rangos: ",
         paste(rank_names(ps), collapse = ", "))
  }

  ps_g <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)
  counts <- get_count_matrix(ps_g)
  tax <- as.data.frame(tax_table(ps_g), stringsAsFactors = FALSE)
  tax$taxon_id <- rownames(tax)

  genus <- as.character(tax$Genus)
  genus[is.na(genus) | trimws(genus) == ""] <- paste0("Unclassified_", tax$taxon_id[is.na(genus) | trimws(genus) == ""])
  tax$Genus_label <- genus

# DESeq2 exige enteros. Los conteos de 16S deberian ser enteros ya.
  if (any(abs(counts - round(counts)) > 1e-8, na.rm = TRUE)) {
    warning("Se detectaron conteos no enteros. Se redondearan para DESeq2; revise el objeto fuente.")
  }
  counts <- round(counts)

  list(ps_genus = ps_g, counts = counts, tax = tax)
}

filter_cmm <- function(counts, tax) {
  prevalence <- rowMeans(counts > 0)
  mean_reads <- rowMeans(counts)
  keep <- prevalence > CMM_PREVALENCIA & mean_reads > CMM_MEDIA_LECTURAS

  qc <- data.frame(
    taxon_id = rownames(counts),
    Genus = tax[rownames(counts), "Genus_label"],
    prevalence = prevalence,
    mean_reads = mean_reads,
    CMM = keep,
    stringsAsFactors = FALSE
  )
  list(counts = counts[keep, , drop = FALSE], qc = qc)
}

make_dds <- function(counts, meta, design_formula) {
  common <- intersect(colnames(counts), rownames(meta))
  counts <- counts[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]

  vars <- all.vars(design_formula)
  complete <- complete.cases(meta[, vars, drop = FALSE])
  meta <- meta[complete, , drop = FALSE]
  counts <- counts[, rownames(meta), drop = FALSE]

# Eliminar taxones totalmente ausentes tras complete-case.
  keep_tax <- rowSums(counts) > 0
  counts <- counts[keep_tax, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  dds
}

extract_results <- function(dds, result_name = NULL, contrast = NULL,
                            tax_lookup = NULL, alpha = FDR_ALPHA) {
  if (!is.null(result_name)) {
    rr <- results(dds, name = result_name, alpha = alpha)
  } else if (!is.null(contrast)) {
    rr <- results(dds, contrast = contrast, alpha = alpha)
  } else {
    stop("Debe especificarse result_name o contrast.")
  }

  out <- as.data.frame(rr)
  out$taxon_id <- rownames(out)
  rownames(out) <- NULL
  names(out)[names(out) == "log2FoldChange"] <- "log2FC"
  names(out)[names(out) == "pvalue"] <- "P.Value"
  names(out)[names(out) == "padj"] <- "FDR"

  if (!is.null(tax_lookup)) {
    tax_tmp <- tax_lookup
    tax_tmp$taxon_id <- rownames(tax_tmp)
    keep_cols <- intersect(c("taxon_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Genus_label"), names(tax_tmp))
    out <- merge(tax_tmp[, keep_cols, drop = FALSE], out, by = "taxon_id", all.y = TRUE, sort = FALSE)
  }
  if ("Genus_label" %in% names(out) && !("Genus" %in% names(out))) out$Genus <- out$Genus_label
  if ("Genus_label" %in% names(out)) {
    out$Genus <- ifelse(is.na(out$Genus) | out$Genus == "", out$Genus_label, out$Genus)
  }

  out$significativo_FDR <- !is.na(out$FDR) & out$FDR < alpha
  out <- out[order(out$FDR, out$P.Value, na.last = TRUE), , drop = FALSE]
  out
}

find_coef <- function(dds, must_contain) {
  rn <- resultsNames(dds)
  hit <- rn[vapply(rn, function(x) all(vapply(must_contain, grepl, logical(1), x = x, fixed = TRUE)), logical(1))]
  if (length(hit) == 0) {
    stop("No se encontro coeficiente con: ", paste(must_contain, collapse = " + "),
         "\nCoeficientes disponibles: ", paste(rn, collapse = ", "))
  }
  if (length(hit) > 1) {
    msg("AVISO: varios coeficientes coinciden; se usara: %s", hit[1])
  }
  hit[1]
}

save_volcano <- function(res, path, title) {
  d <- res
  d$minuslog10FDR <- -log10(pmax(d$FDR, .Machine$double.xmin))
  d$sig <- !is.na(d$FDR) & d$FDR < FDR_ALPHA
  p <- ggplot(d, aes(x = log2FC, y = minuslog10FDR)) +
    geom_point(aes(shape = sig), alpha = 0.8) +
    geom_hline(yintercept = -log10(FDR_ALPHA), linetype = 2) +
    labs(x = "log2 fold change", y = "-log10(FDR)", shape = "FDR<0,05", title = title) +
    theme_minimal(base_size = 11)
  ggsave(path, p, width = 7, height = 5, dpi = 300)
}

run_deseq_analysis <- function(counts, tax, meta, design_formula,
                               result_name = NULL, contrast = NULL,
                               analysis_id, cmm_filter = TRUE) {
  if (cmm_filter) {
    fc <- filter_cmm(counts, tax)
    counts_use <- fc$counts
    write.csv(fc$qc, file.path(OUT_DIR, paste0("QC_CMM_", analysis_id, ".csv")), row.names = FALSE)
  } else {
    counts_use <- counts
  }

  dds <- make_dds(counts_use, meta, design_formula)
  res <- extract_results(dds, result_name = result_name, contrast = contrast, tax_lookup = tax)
  res$analysis_id <- analysis_id
  res$n_modelo <- ncol(dds)
  res$n_taxa_modelo <- nrow(dds)

  write.csv(res, file.path(OUT_DIR, paste0("DESEQ2_", analysis_id, ".csv")), row.names = FALSE)
  save_volcano(res, file.path(OUT_DIR, paste0("VOLCANO_", analysis_id, ".png")), analysis_id)

  list(dds = dds, results = res)
}

# 3. archivos.
# Priorizar las tablas maestras finales ya consolidadas.
DATA_FINAL_DIR <- file.path(PROJECT_DIR, "outputs", "datos_procesados_finales")

file_cohort206_exact <- file.path(DATA_FINAL_DIR, "TABLA_MAESTRA_COHORTE_206.csv")
file_master205_exact <- file.path(DATA_FINAL_DIR, "TABLA_MAESTRA_PACIENTES_205_COMPLETA.csv")
file_master_results_exact <- file.path(DATA_FINAL_DIR, "TABLA_MAESTRA_RESULTADOS_ANALISIS.csv")

file_cohort206 <- if (file.exists(file_cohort206_exact)) {
  file_cohort206_exact
} else {
  find_one_file("TABLA_MAESTRA_COHORTE_206.*\\.csv$")
}

file_master205 <- if (file.exists(file_master205_exact)) {
  file_master205_exact
} else {
  find_one_file("TABLA_MAESTRA_PACIENTES_205_COMPLETA.*\\.csv$")
}

file_master_results <- if (file.exists(file_master_results_exact)) {
  file_master_results_exact
} else {
  find_one_file("TABLA_MAESTRA_RESULTADOS_ANALISIS.*\\.csv$", required = FALSE)
}

msg("Archivos de entrada:")
msg("  phyloseq original: %s", PHYLOSEQ_RDS)
msg("  cohorte 206:       %s", file_cohort206)
msg("  tabla 205:         %s", file_master205)
if (!is.na(file_master_results)) {
  msg("  resultados ANCOM:  %s", file_master_results)
}

cohort206 <- read_csv_flexible(file_cohort206)
master205 <- read_csv_flexible(file_master205)
ps0 <- load_phyloseq_safe()

cohort206$id_muestra <- standardize_id(cohort206$id_muestra)
master205$id_muestra <- standardize_id(master205$id_muestra)
sample_names(ps0) <- standardize_id(sample_names(ps0))

msg("Objeto phyloseq: %d muestras, %d taxa", nsamples(ps0), ntaxa(ps0))
msg("Cohorte 206: %d filas", nrow(cohort206))
msg("Tabla ambiental 205: %d filas", nrow(master205))

# 4. subconjuntos.
# Usar explicitamente los IDs de las tablas maestras para evitar incluir controles.
ids206 <- intersect(cohort206$id_muestra, sample_names(ps0))
ids205 <- intersect(master205$id_muestra, sample_names(ps0))

if (length(ids206) < 190) {
  stop("Solo se emparejaron ", length(ids206), " de los 206 IDs con phyloseq. Revisar identificadores.")
}
if (length(ids205) < 190) {
  stop("Solo se emparejaron ", length(ids205), " de los 205 IDs ambientales con phyloseq. Revisar identificadores.")
}

ps206 <- subset_ps_ids(ps0, ids206)
ps205 <- subset_ps_ids(ps0, ids205)

# Agregar a genero por separado porque el CMM depende de la cohorte analizada.
g206 <- prepare_genus_counts(ps206)
g205 <- prepare_genus_counts(ps205)

# 5. metadatos replica gemas.
meta206 <- cohort206[match(sample_names(ps206), cohort206$id_muestra), , drop = FALSE]
rownames(meta206) <- meta206$id_muestra

meta206$age <- as.numeric(meta206$edad)
meta206$sex <- factor(as.character(meta206$sexo), levels = c("Female", "Male"))
meta206$pool <- factor(meta206$pool)
meta206$cavities <- factor(meta206$caries, levels = c(0, 1), labels = c("No", "Si"))
meta206$gum <- factor(meta206$chicle_30min, levels = c(0, 1), labels = c("No", "Si"))
meta206$antibiotic <- factor(meta206$antibiotico_2m, levels = c(0, 1), labels = c("No", "Si"))
meta206$exac <- factor(meta206$exacerbacion_6m, levels = c(0, 1),
                       labels = c("No_exacerbador", "Exacerbador"))

# QC de codificacion
qc_groups <- as.data.frame(table(meta206$exac, useNA = "ifany"))
names(qc_groups) <- c("grupo", "n")
write.csv(qc_groups, file.path(OUT_DIR, "QC_GRUPOS_EXACERBACION_206.csv"), row.names = FALSE)

# 6. replica deseq2 exacerbacion.
# Modelo comparable con gemas: edad + sexo + covariables salivales + exac..
form_replica <- as.formula(paste("~", paste(c(COVARS_REPLICA, "exac"), collapse = " + ")))

replica_cmm <- run_deseq_analysis(
  counts = g206$counts,
  tax = g206$tax,
  meta = meta206,
  design_formula = form_replica,
  contrast = c("exac", "Exacerbador", "No_exacerbador"),
  analysis_id = "EXAC6M_REPLICA_CMM",
  cmm_filter = TRUE
)

# Sensibilidad anadiendo antibioticos recientes..
form_replica_abx <- as.formula(paste("~", paste(c(COVARS_REPLICA, "antibiotic", "exac"), collapse = " + ")))
replica_abx <- run_deseq_analysis(
  counts = g206$counts,
  tax = g206$tax,
  meta = meta206,
  design_formula = form_replica_abx,
  contrast = c("exac", "Exacerbador", "No_exacerbador"),
  analysis_id = "EXAC6M_REPLICA_CMM_AJUSTE_ANTIBIOTICO",
  cmm_filter = TRUE
)

# 8. metadatos calima 205.
meta205 <- master205[match(sample_names(ps205), master205$id_muestra), , drop = FALSE]
rownames(meta205) <- meta205$id_muestra

meta205$age <- as.numeric(meta205$age)
meta205$sex <- factor(as.character(meta205$sex), levels = c("Female", "Male"))
meta205$island <- factor(as.character(meta205$island), levels = c("Tenerife", "La Palma"))
meta205$pool <- factor(meta205$pool)
meta205$antibiotic <- factor(meta205$antibiotic, levels = c(0, 1), labels = c("No", "Si"))
meta205$cavities <- factor(meta205$cavities, levels = c(0, 1), labels = c("No", "Si"))
meta205$gum <- factor(meta205$gum, levels = c(0, 1), labels = c("No", "Si"))
meta205$exac <- factor(meta205$exac6m, levels = c(0, 1), labels = c("No_exacerbador", "Exacerbador"))
meta205$miteco <- factor(meta205$miteco10, levels = c(0, 1), labels = c("No_MITECO", "Si_MITECO"))
meta205$miteco_days10_z <- as.numeric(meta205$miteco_days10_z)
meta205$naf10_z <- as.numeric(meta205$naf10_z)
meta205$sin_doy <- as.numeric(meta205$sin_doy)
meta205$cos_doy <- as.numeric(meta205$cos_doy)

# El CMM de calima se calcula una sola vez sobre los 205 para que todos los
# contrastes de calima partan del mismo conjunto taxonomico.
cmm205 <- filter_cmm(g205$counts, g205$tax)
counts205_cmm <- cmm205$counts
write.csv(cmm205$qc, file.path(OUT_DIR, "QC_CMM_CALIMA_205.csv"), row.names = FALSE)

# 9. deseq2 calima: modelo base.
# Miteco10 binario.
form_miteco <- as.formula(paste("~", paste(c(COVARS_CALIMA_BASE, "miteco"), collapse = " + ")))
res_miteco <- run_deseq_analysis(
  counts205_cmm, g205$tax, meta205, form_miteco,
  contrast = c("miteco", "Si_MITECO", "No_MITECO"),
  analysis_id = "MITECO10",
  cmm_filter = FALSE
)

# Numero de dias MITECO d1-d10, estandarizado.
form_days <- as.formula(paste("~", paste(c(COVARS_CALIMA_BASE, "miteco_days10_z"), collapse = " + ")))
res_days <- run_deseq_analysis(
  counts205_cmm, g205$tax, meta205, form_days,
  result_name = "miteco_days10_z",
  analysis_id = "MITECO_DIAS10_Z",
  cmm_filter = FALSE
)

# Naf / aporte natural 10d estandarizado.
form_naf <- as.formula(paste("~", paste(c(COVARS_CALIMA_BASE, "naf10_z"), collapse = " + ")))
res_naf <- run_deseq_analysis(
  counts205_cmm, g205$tax, meta205, form_naf,
  result_name = "naf10_z",
  analysis_id = "NAF10_Z",
  cmm_filter = FALSE
)

# Interaccion MITECO x sexo.
cov_no_sex <- setdiff(COVARS_CALIMA_BASE, "sex")
form_int_sex <- as.formula(paste("~", paste(c(cov_no_sex, "sex", "miteco", "sex:miteco"), collapse = " + ")))
dds_int_sex <- make_dds(counts205_cmm, meta205, form_int_sex)
coef_int_sex <- find_coef(dds_int_sex, c("sexMale", "mitecoSi_MITECO"))
out_int_sex <- extract_results(dds_int_sex, result_name = coef_int_sex, tax_lookup = g205$tax)
out_int_sex$analysis_id <- "MITECO10_X_SEXO"
out_int_sex$n_modelo <- ncol(dds_int_sex)
out_int_sex$n_taxa_modelo <- nrow(dds_int_sex)
write.csv(out_int_sex, file.path(OUT_DIR, "DESEQ2_MITECO10_X_SEXO.csv"), row.names = FALSE)
save_volcano(out_int_sex, file.path(OUT_DIR, "VOLCANO_MITECO10_X_SEXO.png"), "MITECO10 x sexo")

# Interaccion MITECO x exacerbacion.
form_int_exac <- as.formula(paste("~", paste(c(COVARS_CALIMA_BASE, "exac", "miteco", "exac:miteco"), collapse = " + ")))
dds_int_exac <- make_dds(counts205_cmm, meta205, form_int_exac)
coef_int_exac <- find_coef(dds_int_exac, c("exacExacerbador", "mitecoSi_MITECO"))
out_int_exac <- extract_results(dds_int_exac, result_name = coef_int_exac, tax_lookup = g205$tax)
out_int_exac$analysis_id <- "MITECO10_X_EXAC6M"
out_int_exac$n_modelo <- ncol(dds_int_exac)
out_int_exac$n_taxa_modelo <- nrow(dds_int_exac)
write.csv(out_int_exac, file.path(OUT_DIR, "DESEQ2_MITECO10_X_EXAC6M.csv"), row.names = FALSE)
save_volcano(out_int_exac, file.path(OUT_DIR, "VOLCANO_MITECO10_X_EXAC6M.png"), "MITECO10 x exacerbacion 6m")

# 10. sensibilidad estacional.
if (EJECUTAR_SENSIBILIDAD_ESTACIONAL) {
  cov_season <- c(COVARS_CALIMA_BASE, "sin_doy", "cos_doy")
  form_miteco_season <- as.formula(paste("~", paste(c(cov_season, "miteco"), collapse = " + ")))
  run_deseq_analysis(
    counts205_cmm, g205$tax, meta205, form_miteco_season,
    contrast = c("miteco", "Si_MITECO", "No_MITECO"),
    analysis_id = "MITECO10_AJUSTE_ESTACIONAL",
    cmm_filter = FALSE
  )
}

# 11. resumen deseq2 calima.
all_calima <- list(
  MITECO10 = res_miteco$results,
  MITECO_dias10 = res_days$results,
  NAF10 = res_naf$results,
  MITECO10_x_sexo = out_int_sex,
  MITECO10_x_exac6m = out_int_exac
)

summary_calima <- do.call(rbind, lapply(names(all_calima), function(nm) {
  x <- all_calima[[nm]]
  data.frame(
    analisis = nm,
    n = if ("n_modelo" %in% names(x)) unique(x$n_modelo)[1] else NA_integer_,
    n_taxa = nrow(x),
    n_p005 = sum(!is.na(x$P.Value) & x$P.Value < 0.05),
    n_FDR005 = sum(!is.na(x$FDR) & x$FDR < FDR_ALPHA),
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_calima, file.path(OUT_DIR, "RESUMEN_DESEQ2_CALIMA.csv"), row.names = FALSE)

# 12. comparacion resumida ancom-bc2.
if (!is.na(file_master_results) && file.exists(file_master_results)) {
  mr <- read_csv_flexible(file_master_results)
  ancom <- mr[grepl("ANCOM-BC2", mr$bloque) & mr$nivel_taxonomico == "Genus", , drop = FALSE]

  map <- data.frame(
    analisis = c("MITECO10", "MITECO10_x_exac6m", "MITECO10_x_sexo", "MITECO_dias10", "NAF10"),
    exposicion_ancom = c("MITECO10", "MITECO10_x_exac6m", "MITECO10_x_sexo", "MITECO_dias10", "NAF10"),
    stringsAsFactors = FALSE
  )
  comp_methods <- merge(summary_calima, map, by = "analisis", all.x = TRUE)
  ancom_small <- ancom[, intersect(c("exposicion", "n", "n_taxa", "n_q005", "n_robustos", "estado_inferencia"), names(ancom)), drop = FALSE]
  names(ancom_small)[names(ancom_small) == "exposicion"] <- "exposicion_ancom"
  names(ancom_small)[names(ancom_small) == "n"] <- "n_ANCOMBC2"
  names(ancom_small)[names(ancom_small) == "n_taxa"] <- "n_taxa_ANCOMBC2"
  names(ancom_small)[names(ancom_small) == "n_q005"] <- "n_q005_ANCOMBC2"
  names(ancom_small)[names(ancom_small) == "n_robustos"] <- "n_robustos_ANCOMBC2"

  comp_methods <- merge(comp_methods, ancom_small, by = "exposicion_ancom", all.x = TRUE)
  write.csv(comp_methods, file.path(OUT_DIR, "COMPARACION_RESUMIDA_DESEQ2_VS_ANCOMBC2.csv"), row.names = FALSE)
}

# 13. informe qc en texto.
sink(file.path(OUT_DIR, "00_INFORME_QC_DESEQ2.txt"))
cat("DESEQ2 CLINICO Y AMBIENTAL / TFM\n")
cat("Fecha de ejecucion: ", as.character(Sys.time()), "\n", sep = "")
cat("Muestras phyloseq inicial: ", nsamples(ps0), "\n", sep = "")
cat("Taxa phyloseq inicial: ", ntaxa(ps0), "\n", sep = "")
cat("IDs emparejados cohorte 206: ", length(ids206), "\n", sep = "")
cat("IDs emparejados cohorte 205 ambiental: ", length(ids205), "\n", sep = "")
cat("Distribucion exacerbacion 206:\n")
print(table(meta206$exac, useNA = "ifany"))
cat("\nCMM replica: prevalencia >", CMM_PREVALENCIA * 100, "% y media >", CMM_MEDIA_LECTURAS, " lecturas/muestra\n", sep = "")
cat("Generos CMM subcohorte 206: ", nrow(filter_cmm(g206$counts, g206$tax)$counts), "\n", sep = "")
cat("\nResumen DESeq2 calima:\n")
print(summary_calima)
cat("\nNOTA DE INTERPRETACION:\n")
cat("- En calima, DESeq2 se considera analisis de sensibilidad frente al principal ANCOM-BC2.\n")
sink()

# 14. session info.
capture.output(sessionInfo(), file = file.path(OUT_DIR, "SESSION_INFO.txt"))

msg("ANALISIS FINALIZADO")
msg("Resultados guardados en: %s", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
msg("Revise primero:")
msg("  1) 00_INFORME_QC_DESEQ2.txt")
msg("  5) RESUMEN_DESEQ2_CALIMA.csv")
msg("  6) COMPARACION_RESUMIDA_DESEQ2_VS_ANCOMBC2.csv (si existe)")
