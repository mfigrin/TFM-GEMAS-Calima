# UNIFRAC_SENSIBILIDAD_GEMAS.R
# Objetivo:
#   Analisis filogenetico de diversidad beta como sensibilidad
#   del analisis principal con Bray-Curtis.
#
# Analisis incluidos:
#   1) PM10/calima 10d -> sdi_exposicion_10d
#   2) Exacerbacion 6M
#   3) Sexo
#
# Distancias:
#   - Weighted UniFrac: objeto filtrado, normalized = TRUE
#   - Unweighted UniFrac: objeto filtrado rarefactado a 10.000
#     lecturas/muestra para reducir el efecto de profundidad
#     sobre presencia/ausencia.
#
# Arbol:
#   Las politomias se resuelven con ape::multi2di(random = FALSE).
#   Las ramas adicionales tienen longitud 0
#
# Salidas:
#   outputs/analisis_microbiota/unifrac_sensibilidad/

# 0. PAQUETES Y OPCIONES

paquetes <- c(
  "phyloseq",
  "ape",
  "vegan",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "writexl"
)

faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]

if (length(faltan) > 0) {
  stop(
    "Faltan paquetes: ",
    paste(faltan, collapse = ", "),
    ". Instálalos antes de ejecutar el script."
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(ape)
  library(vegan)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(writexl)
})

set.seed(20260822)

n_permutaciones <- 999
profundidad_unweighted <- 10000

ruta_ps_ambiental <- file.path(
  "outputs", "analisis_microbiota", "objetos_phyloseq",
  "ps_ambiental_205_filtrado_prev2_total10.rds"
)

ruta_ps_saliva <- file.path(
  "outputs", "analisis_microbiota", "objetos_phyloseq",
  "ps_saliva_206_filtrado_prev2_total10.rds"
)

ruta_salida <- file.path(
  "outputs", "analisis_microbiota", "unifrac_sensibilidad"
)

ruta_tablas <- file.path(ruta_salida, "tablas")
ruta_figuras <- file.path(ruta_salida, "figuras")
ruta_objetos <- file.path(ruta_salida, "objetos")

dir.create(ruta_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_objetos, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ruta_ps_ambiental)) {
  stop("No existe: ", ruta_ps_ambiental)
}

if (!file.exists(ruta_ps_saliva)) {
  stop("No existe: ", ruta_ps_saliva)
}

# 1. FUNCIONES 

convertir_numero <- function(x) {
  suppressWarnings(
    as.numeric(
      gsub(",", ".", trimws(as.character(x)), fixed = TRUE)
    )
  )
}

convertir_logico <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(z))

  out[z %in% c("true", "t", "1", "si", "sí", "yes", "y", "expuesto")] <- TRUE
  out[z %in% c("false", "f", "0", "no", "n", "no_expuesto")] <- FALSE

  out
}

extraer_sample_data <- function(ps_obj) {
  sd_obj <- sample_data(ps_obj)

  sd_df <- data.frame(
    setNames(sd_obj@.Data, names(sd_obj)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(sd_df) <- sample_names(ps_obj)
  sd_df$id_muestra <- rownames(sd_df)

  sd_df %>% relocate(id_muestra)
}

calcular_profundidad <- function(ps_obj) {
  data.frame(
    id_muestra = names(sample_sums(ps_obj)),
    lecturas_totales = as.numeric(sample_sums(ps_obj)),
    stringsAsFactors = FALSE
  )
}

preparar_metadata <- function(ps_obj) {
  meta <- extraer_sample_data(ps_obj) %>%
    left_join(calcular_profundidad(ps_obj), by = "id_muestra")

  if ("Age" %in% names(meta)) {
    meta$Age <- convertir_numero(meta$Age)
  }

  if ("Exacerbations.6M" %in% names(meta)) {
    meta$Exacerbations.6M <- convertir_numero(meta$Exacerbations.6M)

    meta$exacerbacion_6m_factor <- factor(
      ifelse(
        is.na(meta$Exacerbations.6M),
        NA,
        ifelse(
          meta$Exacerbations.6M == 1,
          "Exacerbador",
          "No_exacerbador"
        )
      ),
      levels = c("No_exacerbador", "Exacerbador")
    )
  }

  if ("Sex.letter" %in% names(meta)) {
    meta$Sex.letter <- factor(meta$Sex.letter)
    meta$sexo_factor <- factor(meta$Sex.letter)
  }

  if ("Island" %in% names(meta)) {
    meta$Island <- factor(meta$Island)
  }

  if ("Pool" %in% names(meta)) {
    meta$Pool <- factor(meta$Pool)
  }

  for (v in intersect(
    c("Antibiotic.2mo", "Cavities", "Chewing.gum.30min"),
    names(meta)
  )) {
    meta[[v]] <- convertir_numero(meta[[v]])
    meta[[paste0(v, "_factor")]] <- factor(meta[[v]])
  }

  if ("sdi_exposicion_10d" %in% names(meta)) {
    sdi <- convertir_logico(meta$sdi_exposicion_10d)

    meta$sdi_exposicion_10d_factor <- factor(
      ifelse(
        is.na(sdi),
        NA,
        ifelse(sdi, "Expuesto", "No_expuesto")
      ),
      levels = c("No_expuesto", "Expuesto")
    )
  }

  meta$lecturas_log10 <- log10(meta$lecturas_totales + 1)

  meta
}

resolver_arbol <- function(ps_obj, nombre_dataset) {
  arbol_original <- phy_tree(ps_obj)

  qc_original <- data.frame(
    dataset = nombre_dataset,
    version_arbol = "original",
    n_tips = length(arbol_original$tip.label),
    n_nodos = arbol_original$Nnode,
    binario = ape::is.binary(arbol_original),
    enraizado = ape::is.rooted(arbol_original),
    tips_coinciden_asv = setequal(taxa_names(ps_obj), arbol_original$tip.label),
    tips_duplicados = anyDuplicated(arbol_original$tip.label),
    ramas_na = sum(is.na(arbol_original$edge.length)),
    ramas_cero = sum(arbol_original$edge.length == 0, na.rm = TRUE),
    ramas_negativas = sum(arbol_original$edge.length < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  arbol_binario <- ape::multi2di(
    arbol_original,
    random = FALSE
  )

  qc_binario <- data.frame(
    dataset = nombre_dataset,
    version_arbol = "binario_multi2di",
    n_tips = length(arbol_binario$tip.label),
    n_nodos = arbol_binario$Nnode,
    binario = ape::is.binary(arbol_binario),
    enraizado = ape::is.rooted(arbol_binario),
    tips_coinciden_asv = setequal(taxa_names(ps_obj), arbol_binario$tip.label),
    tips_duplicados = anyDuplicated(arbol_binario$tip.label),
    ramas_na = sum(is.na(arbol_binario$edge.length)),
    ramas_cero = sum(arbol_binario$edge.length == 0, na.rm = TRUE),
    ramas_negativas = sum(arbol_binario$edge.length < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (!ape::is.binary(arbol_binario)) {
    stop("El árbol no quedó binario en ", nombre_dataset)
  }

  if (!ape::is.rooted(arbol_binario)) {
    stop("El árbol no está enraizado en ", nombre_dataset)
  }

  if (!setequal(taxa_names(ps_obj), arbol_binario$tip.label)) {
    stop("ASVs y tips no coinciden en ", nombre_dataset)
  }

  if (any(is.na(arbol_binario$edge.length))) {
    stop("Hay longitudes de rama NA en ", nombre_dataset)
  }

  if (any(arbol_binario$edge.length < 0)) {
    stop("Hay longitudes de rama negativas en ", nombre_dataset)
  }

  ps_out <- ps_obj
  phy_tree(ps_out) <- arbol_binario

  list(
    ps = ps_out,
    qc = bind_rows(qc_original, qc_binario)
  )
}

control_distancia <- function(dist_obj, dataset, distancia) {
  m <- as.matrix(dist_obj)

  data.frame(
    dataset = dataset,
    distancia = distancia,
    n_muestras = nrow(m),
    n_na = sum(is.na(m)),
    n_no_finitos = sum(!is.finite(m)),
    diagonal_max_abs = max(abs(diag(m))),
    asimetria_max_abs = max(abs(m - t(m)), na.rm = TRUE),
    distancia_min = min(m, na.rm = TRUE),
    distancia_max = max(m, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

subdist <- function(dist_obj, ids) {
  m <- as.matrix(dist_obj)
  as.dist(m[ids, ids, drop = FALSE])
}

preparar_datos_modelo <- function(dist_obj, meta, variable, covariables) {
  ids_dist <- labels(dist_obj)

  meta2 <- meta %>%
    filter(id_muestra %in% ids_dist)

  rownames(meta2) <- meta2$id_muestra
  meta2 <- meta2[ids_dist, , drop = FALSE]

  cov_presentes <- covariables[covariables %in% names(meta2)]
  vars_modelo <- unique(c(variable, cov_presentes))

  meta_cc <- meta2 %>%
    select(id_muestra, all_of(vars_modelo)) %>%
    filter(complete.cases(.))

  if (nrow(meta_cc) < 30) {
    return(NULL)
  }

  if (length(unique(meta_cc[[variable]])) < 2) {
    return(NULL)
  }

  rownames(meta_cc) <- meta_cc$id_muestra

  list(
    meta = meta_cc,
    dist = subdist(dist_obj, meta_cc$id_muestra),
    covariables = cov_presentes
  )
}

permanova_ajustada <- function(
  dist_obj,
  meta,
  variable,
  covariables,
  dataset,
  pregunta,
  distancia
) {
  dat <- preparar_datos_modelo(
    dist_obj = dist_obj,
    meta = meta,
    variable = variable,
    covariables = covariables
  )

  if (is.null(dat)) {
    return(NULL)
  }

  formula_txt <- paste(
    "dist_modelo ~",
    paste(c(variable, dat$covariables), collapse = " + ")
  )

  dist_modelo <- dat$dist

  fit <- vegan::adonis2(
    as.formula(formula_txt),
    data = dat$meta,
    permutations = n_permutaciones,
    by = "margin"
  )

  fit_df <- as.data.frame(fit)
  fit_df$termino <- rownames(fit_df)

  fila <- fit_df %>%
    filter(termino == variable)

  if (nrow(fila) != 1) {
    return(NULL)
  }

  data.frame(
    dataset = dataset,
    pregunta = pregunta,
    distancia = distancia,
    variable = variable,
    n_modelo = nrow(dat$meta),
    formula = formula_txt,
    Df = fila$Df,
    SumOfSqs = fila$SumOfSqs,
    R2 = fila$R2,
    F = fila$F,
    p_valor = fila$`Pr(>F)`,
    stringsAsFactors = FALSE
  )
}

dispersion_grupos <- function(
  dist_obj,
  meta,
  variable,
  covariables,
  dataset,
  pregunta,
  distancia
) {
  dat <- preparar_datos_modelo(
    dist_obj = dist_obj,
    meta = meta,
    variable = variable,
    covariables = covariables
  )

  if (is.null(dat)) {
    return(NULL)
  }

  grupo <- factor(dat$meta[[variable]])

  if (nlevels(grupo) < 2) {
    return(NULL)
  }

  bd <- vegan::betadisper(dat$dist, grupo)

  perm <- vegan::permutest(
    bd,
    permutations = n_permutaciones
  )

  data.frame(
    dataset = dataset,
    pregunta = pregunta,
    distancia = distancia,
    variable = variable,
    n = length(grupo),
    n_grupos = nlevels(grupo),
    F = perm$tab[1, "F"],
    p_valor = perm$tab[1, "Pr(>F)"],
    stringsAsFactors = FALSE
  )
}

crear_pcoa <- function(
  dist_obj,
  meta,
  variable,
  covariables,
  dataset,
  pregunta,
  distancia
) {
  dat <- preparar_datos_modelo(
    dist_obj = dist_obj,
    meta = meta,
    variable = variable,
    covariables = covariables
  )

  if (is.null(dat)) {
    return(NULL)
  }

  ord <- ape::pcoa(dat$dist)

  coords <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
  names(coords) <- c("PCoA1", "PCoA2")
  coords$id_muestra <- rownames(coords)

  coords <- coords %>%
    left_join(
      dat$meta %>%
        select(id_muestra, all_of(variable)),
      by = "id_muestra"
    )

  var1 <- 100 * ord$values$Relative_eig[1]
  var2 <- 100 * ord$values$Relative_eig[2]

  p <- ggplot(
    coords,
    aes(
      x = PCoA1,
      y = PCoA2,
      shape = .data[[variable]]
    )
  ) +
    geom_point(
      aes(fill = .data[[variable]]),
      size = 2.5,
      alpha = 0.75
    ) +
    labs(
      title = paste0(
        distancia,
        " - ",
        pregunta
      ),
      x = paste0("PCoA1 (", round(var1, 1), "%)"),
      y = paste0("PCoA2 (", round(var2, 1), "%)"),
      fill = variable,
      shape = variable
    ) +
    theme_minimal(base_size = 12)

  nombre_seguro <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    paste(
      "PCoA",
      dataset,
      distancia,
      variable,
      sep = "_"
    )
  )

  archivo_png <- file.path(
    ruta_figuras,
    paste0(nombre_seguro, ".png")
  )

  ggsave(
    filename = archivo_png,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )

  archivo_csv <- file.path(
    ruta_tablas,
    paste0(nombre_seguro, "_coordenadas.csv")
  )

  readr::write_csv(
    coords,
    archivo_csv,
    na = ""
  )

  data.frame(
    dataset = dataset,
    pregunta = pregunta,
    distancia = distancia,
    variable = variable,
    figura = archivo_png,
    coordenadas = archivo_csv,
    stringsAsFactors = FALSE
  )
}

# 2. CARGAR OBJETOS

ps_ambiental <- readRDS(ruta_ps_ambiental)
ps_saliva <- readRDS(ruta_ps_saliva)

cat("\n============================================\n")
cat("OBJETOS CARGADOS\n")
cat("============================================\n")
cat(
  "Ambiental:",
  nsamples(ps_ambiental),
  "muestras;",
  ntaxa(ps_ambiental),
  "ASVs\n"
)
cat(
  "Saliva:",
  nsamples(ps_saliva),
  "muestras;",
  ntaxa(ps_saliva),
  "ASVs\n"
)

# 3. RESOLVER POLITOMIAS DE FORMA DETERMINISTA

res_arbol_ambiental <- resolver_arbol(
  ps_ambiental,
  "ambiental_205"
)

res_arbol_saliva <- resolver_arbol(
  ps_saliva,
  "saliva_206"
)

ps_ambiental_bin <- res_arbol_ambiental$ps
ps_saliva_bin <- res_arbol_saliva$ps

qc_arbol <- bind_rows(
  res_arbol_ambiental$qc,
  res_arbol_saliva$qc
)

readr::write_csv(
  qc_arbol,
  file.path(ruta_tablas, "01_qc_arboles_unifrac.csv"),
  na = ""
)

saveRDS(
  ps_ambiental_bin,
  file.path(ruta_objetos, "ps_ambiental_205_arbol_binario.rds")
)

saveRDS(
  ps_saliva_bin,
  file.path(ruta_objetos, "ps_saliva_206_arbol_binario.rds")
)


# 4. METADATOS

meta_ambiental <- preparar_metadata(ps_ambiental_bin)
meta_saliva <- preparar_metadata(ps_saliva_bin)


# 5. WEIGHTED UNIFRAC


cat("\nCalculando Weighted UniFrac ambiental...\n")

dist_w_ambiental <- phyloseq::UniFrac(
  ps_ambiental_bin,
  weighted = TRUE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

cat("Weighted UniFrac ambiental: OK\n")

cat("\nCalculando Weighted UniFrac saliva...\n")

dist_w_saliva <- phyloseq::UniFrac(
  ps_saliva_bin,
  weighted = TRUE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

cat("Weighted UniFrac saliva: OK\n")

saveRDS(
  dist_w_ambiental,
  file.path(ruta_objetos, "dist_weighted_unifrac_ambiental_205.rds")
)

saveRDS(
  dist_w_saliva,
  file.path(ruta_objetos, "dist_weighted_unifrac_saliva_206.rds")
)


# 6. UNWEIGHTED UNIFRAC CON RAREFACCION A 10.000

min_ambiental <- min(sample_sums(ps_ambiental_bin))
min_saliva <- min(sample_sums(ps_saliva_bin))

if (min_ambiental < profundidad_unweighted) {
  stop(
    "Hay muestras ambientales con menos de ",
    profundidad_unweighted,
    " lecturas. Mínimo: ",
    min_ambiental
  )
}

if (min_saliva < profundidad_unweighted) {
  stop(
    "Hay muestras de saliva con menos de ",
    profundidad_unweighted,
    " lecturas. Mínimo: ",
    min_saliva
  )
}

cat(
  "\nRarefaccion para Unweighted UniFrac a ",
  profundidad_unweighted,
  " lecturas/muestra...\n",
  sep = ""
)

set.seed(20260822)

ps_ambiental_rare <- phyloseq::rarefy_even_depth(
  ps_ambiental_bin,
  sample.size = profundidad_unweighted,
  rngseed = 20260822,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

set.seed(20260822)

ps_saliva_rare <- phyloseq::rarefy_even_depth(
  ps_saliva_bin,
  sample.size = profundidad_unweighted,
  rngseed = 20260822,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_ambiental_rare <- resolver_arbol(
  ps_ambiental_rare,
  "ambiental_205_rarefactado"
)$ps

ps_saliva_rare <- resolver_arbol(
  ps_saliva_rare,
  "saliva_206_rarefactado"
)$ps

cat("\nCalculando Unweighted UniFrac ambiental...\n")

dist_u_ambiental <- phyloseq::UniFrac(
  ps_ambiental_rare,
  weighted = FALSE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

cat("Unweighted UniFrac ambiental: OK\n")

cat("\nCalculando Unweighted UniFrac saliva...\n")

dist_u_saliva <- phyloseq::UniFrac(
  ps_saliva_rare,
  weighted = FALSE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

cat("Unweighted UniFrac saliva: OK\n")

saveRDS(
  dist_u_ambiental,
  file.path(ruta_objetos, "dist_unweighted_unifrac_ambiental_205_rarefactado_10000.rds")
)

saveRDS(
  dist_u_saliva,
  file.path(ruta_objetos, "dist_unweighted_unifrac_saliva_206_rarefactado_10000.rds")
)


# 7. CONTROL DE CALIDAD DE DISTANCIAS

qc_distancias <- bind_rows(
  control_distancia(
    dist_w_ambiental,
    "ambiental_205",
    "Weighted_UniFrac"
  ),
  control_distancia(
    dist_u_ambiental,
    "ambiental_205",
    "Unweighted_UniFrac_rarefactado_10000"
  ),
  control_distancia(
    dist_w_saliva,
    "saliva_206",
    "Weighted_UniFrac"
  ),
  control_distancia(
    dist_u_saliva,
    "saliva_206",
    "Unweighted_UniFrac_rarefactado_10000"
  )
)

readr::write_csv(
  qc_distancias,
  file.path(ruta_tablas, "02_qc_distancias_unifrac.csv"),
  na = ""
)


# 8. PLAN DE MODELOS

covars_base <- c(
  "Age",
  "Sex.letter",
  "Island",
  "Pool",
  "Antibiotic.2mo_factor",
  "Chewing.gum.30min_factor",
  "Cavities_factor",
  "lecturas_log10"
)

covars_ambiental <- c(
  covars_base,
  "exacerbacion_6m_factor"
)

covars_exacerbacion <- covars_base

covars_sexo <- c(
  "Age",
  "Island",
  "Pool",
  "Antibiotic.2mo_factor",
  "Chewing.gum.30min_factor",
  "Cavities_factor",
  "exacerbacion_6m_factor",
  "lecturas_log10"
)

plan <- tibble::tribble(
  ~dataset, ~pregunta, ~variable,
  "ambiental_205", "PM10_calima_10d", "sdi_exposicion_10d_factor",
  "saliva_206", "Exacerbacion_6M", "exacerbacion_6m_factor",
  "saliva_206", "Sexo", "sexo_factor"
)

readr::write_csv(
  plan,
  file.path(ruta_tablas, "03_plan_unifrac_sensibilidad.csv"),
  na = ""
)


# 9. PERMANOVA, BETADISPER Y PCoA

resultados_permanova <- list()
resultados_dispersion <- list()
figuras_pcoa <- list()

contador <- 1

for (i in seq_len(nrow(plan))) {

  fila <- plan[i, ]

  dataset_i <- fila$dataset
  pregunta_i <- fila$pregunta
  variable_i <- fila$variable

  if (dataset_i == "ambiental_205") {
    meta_i <- meta_ambiental
    distancias_i <- list(
      Weighted_UniFrac = dist_w_ambiental,
      Unweighted_UniFrac_rarefactado_10000 = dist_u_ambiental
    )
    covars_i <- covars_ambiental
  } else {
    meta_i <- meta_saliva
    distancias_i <- list(
      Weighted_UniFrac = dist_w_saliva,
      Unweighted_UniFrac_rarefactado_10000 = dist_u_saliva
    )

    if (pregunta_i == "Exacerbacion_6M") {
      covars_i <- covars_exacerbacion
    } else {
      covars_i <- covars_sexo
    }
  }

  for (nombre_dist in names(distancias_i)) {

    dist_i <- distancias_i[[nombre_dist]]

    cat(
      "\nAnalizando: ",
      pregunta_i,
      " | ",
      nombre_dist,
      "\n",
      sep = ""
    )

    resultados_permanova[[contador]] <- permanova_ajustada(
      dist_obj = dist_i,
      meta = meta_i,
      variable = variable_i,
      covariables = covars_i,
      dataset = dataset_i,
      pregunta = pregunta_i,
      distancia = nombre_dist
    )

    resultados_dispersion[[contador]] <- dispersion_grupos(
      dist_obj = dist_i,
      meta = meta_i,
      variable = variable_i,
      covariables = covars_i,
      dataset = dataset_i,
      pregunta = pregunta_i,
      distancia = nombre_dist
    )

    figuras_pcoa[[contador]] <- crear_pcoa(
      dist_obj = dist_i,
      meta = meta_i,
      variable = variable_i,
      covariables = covars_i,
      dataset = dataset_i,
      pregunta = pregunta_i,
      distancia = nombre_dist
    )

    contador <- contador + 1
  }
}

permanova_final <- bind_rows(resultados_permanova) %>%
  group_by(pregunta) %>%
  mutate(
    p_BH_pregunta = p.adjust(p_valor, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    p_BH_global = p.adjust(p_valor, method = "BH")
  ) %>%
  arrange(pregunta, distancia)

dispersion_final <- bind_rows(resultados_dispersion) %>%
  group_by(pregunta) %>%
  mutate(
    p_BH_pregunta = p.adjust(p_valor, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    p_BH_global = p.adjust(p_valor, method = "BH")
  ) %>%
  arrange(pregunta, distancia)

figuras_final <- bind_rows(figuras_pcoa)

readr::write_csv(
  permanova_final,
  file.path(ruta_tablas, "04_permanova_unifrac_ajustada.csv"),
  na = ""
)

readr::write_csv(
  dispersion_final,
  file.path(ruta_tablas, "05_betadisper_unifrac.csv"),
  na = ""
)

readr::write_csv(
  figuras_final,
  file.path(ruta_tablas, "06_figuras_pcoa_unifrac.csv"),
  na = ""
)


# 10. RESUMEN 

notas <- data.frame(
  aspecto = c(
    "Papel analitico",
    "Weighted UniFrac",
    "Unweighted UniFrac",
    "Arbol",
    "Ramas de longitud cero",
    "PERMANOVA",
    "Dispersion",
    "Interpretacion"
  ),
  detalle = c(
    "Analisis filogenetico de sensibilidad. Bray-Curtis permanece como analisis principal de diversidad beta.",
    "Calculado sobre el objeto filtrado con normalized=TRUE.",
    "Calculado tras rarefaccion reproducible a 10.000 lecturas/muestra para reducir el efecto de profundidad sobre presencia/ausencia.",
    "Politomias resueltas con ape::multi2di(random=FALSE) antes de UniFrac.",
    "Las ramas cero se conservan; no se sustituyen por valores positivos arbitrarios.",
    "adonis2, 999 permutaciones, modelo ajustado y by='margin'.",
    "betadisper/permutest en el mismo conjunto de casos completos del modelo.",
    "Un resultado PERMANOVA significativo debe interpretarse junto con betadisper. Si betadisper es significativo, la separacion puede reflejar diferencias de dispersion ademas de centroides."
  ),
  stringsAsFactors = FALSE
)

writexl::write_xlsx(
  list(
    QC_arbol = qc_arbol,
    QC_distancias = qc_distancias,
    Plan = plan,
    PERMANOVA = permanova_final,
    Betadisper = dispersion_final,
    Figuras_PCoA = figuras_final,
    Notas = notas
  ),
  path = file.path(
    ruta_tablas,
    "07_resumen_UniFrac_sensibilidad.xlsx"
  )
)


# 11. MENSAJE FINAL

cat("\n============================================\n")
cat("ANALISIS UNIFRAC FINALIZADO\n")
cat("============================================\n")
cat("Salida:\n")
cat(ruta_salida, "\n\n")

cat("PERMANOVA:\n")
print(
  permanova_final %>%
    select(
      pregunta,
      distancia,
      n_modelo,
      R2,
      F,
      p_valor,
      p_BH_pregunta,
      p_BH_global
    )
)

cat("\nBetadisper:\n")
print(
  dispersion_final %>%
    select(
      pregunta,
      distancia,
      n,
      F,
      p_valor,
      p_BH_pregunta,
      p_BH_global
    )
)

cat(
  "\nResumen Excel: ",
  file.path(
    ruta_tablas,
    "07_resumen_UniFrac_sensibilidad.xlsx"
  ),
  "\n",
  sep = ""
)
