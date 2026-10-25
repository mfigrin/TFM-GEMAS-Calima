# PERMANOVA principal, dispersión y sensibilidad isla-fecha.
# 9.999 permutaciones y semilla 20260831 desde la primera ejecución.
# La tabla de covariables reutiliza el ajuste completo de Bray-Curtis.

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260831)

PROJECT_DIR <- normalizePath(Sys.getenv("TFM_PROJECT_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
N_PERM <- 9999
SEED <- 20260831

if (!dir.exists(PROJECT_DIR)) {
  stop("No se encuentra PROJECT_DIR: ", PROJECT_DIR)
}

setwd(PROJECT_DIR)

OUT_DIR <- file.path(PROJECT_DIR, "outputs", "23_PERMANOVA_FINAL_ESTABLE")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("phyloseq", "vegan", "dplyr", "readr", "tibble")
faltan <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan) > 0) {
  stop("Faltan paquetes: ", paste(faltan, collapse = ", "))
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(readr)
  library(tibble)
})

msg <- function(...) cat(sprintf(...), "\n")

ruta_existente <- function(rutas, etiqueta) {
  ok <- rutas[file.exists(rutas)]
  if (length(ok) == 0) {
    stop("No se encontró ", etiqueta, ".\nRutas probadas:\n", paste(rutas, collapse = "\n"))
  }
  normalizePath(ok[1], winslash = "/", mustWork = TRUE)
}

dist_labels <- function(d) {
  z <- attr(d, "Labels")
  if (is.null(z) || length(z) == 0) z <- rownames(as.matrix(d))
  as.character(z)
}

subdist <- function(d, ids) {
  labs <- dist_labels(d)
  ids <- intersect(as.character(ids), labs)
  m <- as.matrix(d)
  stats::as.dist(m[ids, ids, drop = FALSE])
}

n_niveles <- function(x) length(unique(x[!is.na(x)]))

covariables_activas <- function(data, covars) {
  covars[
    vapply(covars, function(v) v %in% names(data) && n_niveles(data[[v]]) > 1, logical(1))
  ]
}

primero_no_na <- function(x) {
  z <- x[!is.na(x)]
  if (length(z) == 0) return(NA)
  z[1]
}

media_segura <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# Archivos.

ruta_master <- ruta_existente(
  c(
    file.path(PROJECT_DIR, "tablas_maestras", "TABLA_MAESTRA_PACIENTES_205_COMPLETA.csv"),
    file.path(PROJECT_DIR, "outputs", "datos_procesados_finales",
              "TABLA_MAESTRA_PACIENTES_205_COMPLETA.csv")
  ),
  "tabla maestra 205"
)

ruta_ps <- ruta_existente(
  c(
    file.path(PROJECT_DIR, "outputs", "analisis_microbiota",
              "reanalisis_completo_MITECO", "objetos",
              "ps_ambiental_205_filtrado_con_MITECO.rds"),
    file.path(PROJECT_DIR, "outputs", "analisis_microbiota",
              "objetos_phyloseq",
              "ps_ambiental_205_filtrado_prev2_total10.rds")
  ),
  "phyloseq ambiental"
)

ruta_wu <- ruta_existente(
  c(
    file.path(PROJECT_DIR, "outputs", "analisis_microbiota",
              "unifrac_sensibilidad", "objetos",
              "dist_weighted_unifrac_ambiental_205.rds")
  ),
  "weighted UniFrac"
)

ruta_uu <- ruta_existente(
  c(
    file.path(PROJECT_DIR, "outputs", "analisis_microbiota",
              "unifrac_sensibilidad", "objetos",
              "dist_unweighted_unifrac_ambiental_205_rarefactado_10000.rds")
  ),
  "unweighted UniFrac"
)

# Datos.

meta <- read_csv(ruta_master, show_col_types = FALSE) %>%
  mutate(
    id_muestra = as.character(id_muestra),
    age = as.numeric(age),
    sex = factor(as.character(sex), levels = c("Female", "Male")),
    island = factor(as.character(island), levels = c("La Palma", "Tenerife")),
    pool = factor(pool),
    antibiotic = factor(as.numeric(antibiotic), levels = c(0, 1)),
    cavities = factor(as.numeric(cavities), levels = c(0, 1)),
    gum = factor(as.numeric(gum), levels = c(0, 1)),
    exac6m = factor(as.numeric(exac6m), levels = c(0, 1)),
    year = factor(year),
    miteco10 = as.numeric(miteco10),
    cluster_isla_fecha = factor(cluster_isla_fecha)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

rownames(meta) <- meta$id_muestra

ps <- readRDS(ruta_ps)

ids <- intersect(sample_names(ps), meta$id_muestra)
if (length(ids) != 205) {
  stop("Se esperaban 205 IDs comunes y se encontraron ", length(ids))
}

ps <- prune_samples(ids, ps)
meta <- meta[match(sample_names(ps), meta$id_muestra), , drop = FALSE]
rownames(meta) <- meta$id_muestra

# Bray-Curtis sobre abundancia relativa
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
dist_bray <- phyloseq::distance(ps_rel, method = "bray")

distancias <- list(
  Bray_Curtis = dist_bray,
  Weighted_UniFrac = readRDS(ruta_wu),
  Unweighted_UniFrac_rarefactado = readRDS(ruta_uu)
)

# Modelo principal.

covars <- c(
  "age",
  "sex",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m",
  "year",
  "sin_doy",
  "cos_doy"
)

fits_principales <- list()

run_permanova <- function(d, nombre) {
  labs <- dist_labels(d)
  ids <- intersect(labs, meta$id_muestra)

  dd <- meta[match(ids, meta$id_muestra), , drop = FALSE]

  vars <- c("miteco10", covars)
  cc <- complete.cases(dd[, vars, drop = FALSE])
  dd <- droplevels(dd[cc, , drop = FALSE])

  d2 <- subdist(d, dd$id_muestra)
  labs2 <- dist_labels(d2)
  dd <- dd[match(labs2, dd$id_muestra), , drop = FALSE]
  rownames(dd) <- dd$id_muestra

  cov2 <- covariables_activas(dd, covars)

  formula_txt <- paste("d2 ~", paste(c("miteco10", cov2), collapse = " + "))

  set.seed(SEED)

  fit <- adonis2(as.formula(formula_txt), data = dd, permutations = N_PERM, by = "margin")

  fits_principales[[nombre]] <<- fit

  tab <- as.data.frame(fit)
  tab$term <- rownames(tab)

  fila <- tab[tab$term == "miteco10", , drop = FALSE]

  tibble(
    distancia = nombre,
    n = nrow(dd),
    permutaciones = N_PERM,
    R2 = fila$R2,
    pseudo_F = fila$F,
    p = fila$`Pr(>F)`,
    formula = formula_txt
  )
}

msg("Ejecutando PERMANOVA estable con %d permutaciones...", N_PERM)

beta_final <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      msg("  %s", nm)
      run_permanova(distancias[[nm]], nm)
    }
  )
) %>%
  mutate(q_BH_3_distancias = p.adjust(p, method = "BH"))

write_csv(beta_final, file.path(OUT_DIR, "01_PERMANOVA_FINAL_9999.csv"))

# PERMANOVA completa de Bray-Curtis.

fit_bray <- fits_principales[["Bray_Curtis"]]

bray_full <- as.data.frame(fit_bray)
bray_full$term <- rownames(bray_full)
bray_full <- as_tibble(bray_full) %>%
  relocate(term)

idx <- !bray_full$term %in% c("Residual", "Total") &
  !is.na(bray_full$`Pr(>F)`)

bray_full$q_BH_modelo <- NA_real_
bray_full$q_BH_modelo[idx] <- p.adjust(bray_full$`Pr(>F)`[idx], method = "BH")

write_csv(bray_full, file.path(OUT_DIR, "02_PERMANOVA_BRAY_COVARIABLES_9999.csv"))

# Permdisp.

ids_disp <- intersect(dist_labels(dist_bray), meta$id_muestra)
grupo <- meta[match(ids_disp, meta$id_muestra), "miteco10"]
grupo <- factor(ifelse(grupo == 1, "Con intrusión oficial", "Sin intrusión oficial"))

d_disp <- subdist(dist_bray, ids_disp)
grupo <- grupo[match(dist_labels(d_disp), ids_disp)]

bd <- betadisper(d_disp, grupo)

set.seed(SEED)

bd_perm <- permutest(bd, permutations = N_PERM)

permdisp <- tibble(
  n = length(grupo),
  permutaciones = N_PERM,
  F = bd_perm$tab[1, "F"],
  p = bd_perm$tab[1, "Pr(>F)"]
)

write_csv(permdisp, file.path(OUT_DIR, "03_PERMDISP_BRAY_9999.csv"))

# Sensibilidad agregada por isla-fecha.

otu_rel <- as(otu_table(ps_rel), "matrix")
if (!taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)

clusters <- as.character(meta$cluster_isla_fecha)
names(clusters) <- meta$id_muestra

ids_otu <- colnames(otu_rel)
clusters_otu <- clusters[ids_otu]
niveles_cluster <- unique(clusters_otu)

perfil_cluster <- sapply(
  niveles_cluster,
  function(cl) {
    ii <- ids_otu[clusters_otu == cl]
    if (length(ii) == 1) {
      otu_rel[, ii, drop = TRUE]
    } else {
      rowMeans(otu_rel[, ii, drop = FALSE])
    }
  }
)

colnames(perfil_cluster) <- niveles_cluster

d_cluster <- vegdist(t(perfil_cluster), method = "bray")

meta_cluster <- meta %>%
  group_by(cluster_isla_fecha) %>%
  summarise(
    n_muestras = n(),
    island = primero_no_na(as.character(island)),
    miteco10 = primero_no_na(miteco10),
    year = primero_no_na(as.character(year)),
    sin_doy = primero_no_na(sin_doy),
    cos_doy = primero_no_na(cos_doy),
    age_mean = media_segura(age),
    prop_male = media_segura(as.numeric(sex == "Male")),
    prop_antibiotic = media_segura(as.numeric(as.character(antibiotic))),
    prop_cavities = media_segura(as.numeric(as.character(cavities))),
    prop_gum = media_segura(as.numeric(as.character(gum))),
    prop_exac6m = media_segura(as.numeric(as.character(exac6m))),
    prop_pool2 = media_segura(as.numeric(as.character(pool) == "2")),
    prop_pool3 = media_segura(as.numeric(as.character(pool) == "3")),
    prop_pool6 = media_segura(as.numeric(as.character(pool) == "6")),
    .groups = "drop"
  ) %>%
  mutate(
    cluster_isla_fecha = as.character(cluster_isla_fecha),
    island = factor(island),
    year = factor(year)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

rownames(meta_cluster) <- meta_cluster$cluster_isla_fecha

cov_cluster <- c(
  "age_mean",
  "prop_male",
  "island",
  "prop_antibiotic",
  "prop_cavities",
  "prop_gum",
  "prop_exac6m",
  "prop_pool2",
  "prop_pool3",
  "prop_pool6",
  "n_muestras",
  "year",
  "sin_doy",
  "cos_doy"
)

ids_c <- dist_labels(d_cluster)
dc <- meta_cluster[ids_c, , drop = FALSE]

vars_c <- c("miteco10", cov_cluster)
cc <- complete.cases(dc[, vars_c, drop = FALSE])
dc <- droplevels(dc[cc, , drop = FALSE])

d_cluster2 <- subdist(d_cluster, rownames(dc))
dc <- dc[match(dist_labels(d_cluster2), rownames(dc)), , drop = FALSE]

cov_c2 <- covariables_activas(dc, cov_cluster)

formula_cluster <- paste("d_cluster2 ~", paste(c("miteco10", cov_c2), collapse = " + "))

set.seed(SEED)

fit_cluster <- adonis2(as.formula(formula_cluster), data = dc, permutations = N_PERM, by = "margin")

tab_c <- as.data.frame(fit_cluster)
tab_c$term <- rownames(tab_c)
fila_c <- tab_c[tab_c$term == "miteco10", , drop = FALSE]

sens_cluster <- tibble(
  n_clusters = nrow(dc),
  permutaciones = N_PERM,
  R2 = fila_c$R2,
  pseudo_F = fila_c$F,
  p = fila_c$`Pr(>F)`,
  formula = formula_cluster
)

write_csv(sens_cluster, file.path(OUT_DIR, "04_SENSIBILIDAD_CLUSTER_ISLA_FECHA_9999.csv"))

# Resumen.

sink(file.path(OUT_DIR, "00_RESUMEN_PERMANOVA_FINAL.txt"))

cat("PERMANOVA FINAL ESTABLE - TFM GEMAS\n")
cat("Semilla:", SEED, "\n")
cat("Permutaciones:", N_PERM, "\n\n")

cat("Tres distancias beta:\n")
print(beta_final)
cat("\n")

cat("PERMDISP Bray-Curtis:\n")
print(permdisp)
cat("\n")

cat("Sensibilidad agregada por isla-fecha:\n")
print(sens_cluster)
cat("\n")

cat("Alcance:\n")
cat("La prueba principal usa 9.999 permutaciones. Los contrastes secundarios se guardan aparte.\n")
cat("Si la entrada es sintética, estas cifras no son resultados del TFM.\n")

sink()

capture.output(sessionInfo(), file = file.path(OUT_DIR, "SESSION_INFO.txt"))

cat("PERMANOVA FINAL TERMINADA\n")
cat("Resultados:\n", OUT_DIR, "\n")

