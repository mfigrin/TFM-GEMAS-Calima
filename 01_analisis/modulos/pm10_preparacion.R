# Modelos complementarios de PM10 y tablas de resultados.
# Los resultados principales se leen de 06; no se recalculan con 999 permutaciones.
# Los coeficientes, ajustes y familias de FDR de PM10 se conservan.

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260830)

PROJECT_DIR <- normalizePath(Sys.getenv("TFM_PROJECT_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)

# Si el proyecto se movió, se intenta usar el directorio de trabajo
# siempre que contenga TFM_GEMAS.Rproj.
if (!dir.exists(PROJECT_DIR)) {
  if (file.exists(file.path(getwd(), "TFM_GEMAS.Rproj"))) {
    PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  } else {
    stop(
      "No se encuentra el proyecto en:\n", PROJECT_DIR,
      "\nAbre TFM_GEMAS.Rproj o modifica PROJECT_DIR al inicio del script."
    )
  }
}

setwd(PROJECT_DIR)

SEED <- 20260830
N_PERM <- 999
ALPHA <- 0.05
FUENTE_TFM <- "Times New Roman"

COLOR_ROJO <- "#E30613"
COLOR_ROJO_OSCURO <- "#B0000B"
COLOR_GRIS <- "#777777"
COLOR_GRIS_CLARO <- "#D9D9D9"
COLOR_NEGRO <- "#222222"

OUT_DIR <- file.path(PROJECT_DIR, "outputs", "22_CIERRE_RESULTADOS_TFM")
TAB_DIR <- file.path(OUT_DIR, "tablas")
FIG_DIR <- file.path(OUT_DIR, "figuras")
OBJ_DIR <- file.path(OUT_DIR, "objetos")
LOG_DIR <- file.path(OUT_DIR, "logs")

invisible(lapply(
  c(OUT_DIR, TAB_DIR, FIG_DIR, OBJ_DIR, LOG_DIR),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# Paquetes.

paquetes <- c(
  "phyloseq",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "ggplot2",
  "vegan",
  "sandwich",
  "lmtest",
  "writexl",
  "scales",
  "patchwork"
)

faltan <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltan) > 0) {
  stop(
    "Faltan paquetes necesarios: ",
    paste(faltan, collapse = ", "),
    "\nInstálalos antes de ejecutar el script."
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(vegan)
  library(sandwich)
  library(lmtest)
  library(writexl)
  library(scales)
  library(patchwork)
})

# Funciones generales.

msg <- function(...) {
  cat(sprintf(...), "\n")
}

registrar <- local({
  x <- tibble(etapa = character(), estado = character(), detalle = character())

  function(etapa = NULL, estado = NULL, detalle = NULL, devolver = FALSE) {
    if (isTRUE(devolver)) {
      return(x)
    }

    x <<- bind_rows(
      x,
      tibble(
        etapa = as.character(etapa),
        estado = as.character(estado),
        detalle = as.character(detalle)
      )
    )

    invisible(NULL)
  }
})

ruta_existente <- function(rutas, etiqueta, obligatorio = TRUE) {
  ok <- rutas[file.exists(rutas)]

  if (length(ok) >= 1) {
    if (length(ok) > 1) {
      msg("AVISO: varias rutas válidas para '%s'. Se usará: %s", etiqueta, ok[1])
    }

    return(normalizePath(ok[1], winslash = "/", mustWork = TRUE))
  }

  if (obligatorio) {
    stop(
      "No se encontró el archivo requerido: ", etiqueta,
      "\nRutas comprobadas:\n- ",
      paste(rutas, collapse = "\n- ")
    )
  }

  NA_character_
}

leer_csv <- function(ruta) {
  readr::read_csv(ruta, show_col_types = FALSE, na = c("", "NA", "NaN", "NULL", "null"))
}

zscore_seguro <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- sd(x, na.rm = TRUE)

  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / s
}

media_segura <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

primero_no_na <- function(x) {
  z <- x[!is.na(x)]

  if (length(z) == 0) {
    return(NA)
  }

  z[1]
}

n_niveles <- function(x) {
  length(unique(x[!is.na(x)]))
}

covariables_activas <- function(data, covariables) {
  covariables[
    vapply(
      covariables,
      function(v) {
        v %in% names(data) && n_niveles(data[[v]]) > 1
      },
      logical(1)
    )
  ]
}

dist_labels <- function(d) {
  labs <- attr(d, "Labels")
  if (is.null(labs)) {
    labs <- labels(d)
  }
  as.character(labs)
}

subdist <- function(d, ids) {
  labs <- dist_labels(d)
  ids <- intersect(as.character(ids), labs)

  if (length(ids) < 2) {
    stop(
      "No hay suficientes IDs comunes para submuestrear la distancia. ",
      "IDs solicitados: ", length(unique(ids)),
      "; etiquetas de la distancia: ", length(labs), "."
    )
  }

  m <- as.matrix(d)
  stats::as.dist(m[ids, ids, drop = FALSE])
}

otu_taxa_en_filas <- function(ps_obj) {
  m <- as(otu_table(ps_obj), "matrix")

  if (!taxa_are_rows(ps_obj)) {
    m <- t(m)
  }

  storage.mode(m) <- "numeric"
  m
}

formatear_mediana_ric <- function(x, digitos = 1) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return("NA")
  }

  q <- quantile(x, probs = c(0.25, 0.50, 0.75), names = FALSE, na.rm = TRUE)

  sprintf(paste0("%.", digitos, "f [%.", digitos, "f–%.", digitos, "f]"), q[2], q[1], q[3])
}

formatear_n_pct <- function(x, condicion = function(z) z == 1) {
  ok <- !is.na(x)

  if (sum(ok) == 0) {
    return("NA")
  }

  n <- sum(condicion(x[ok]), na.rm = TRUE)
  d <- sum(ok)

  sprintf("%d (%.1f %%)", n, 100 * n / d)
}

tema_tfm <- function(base_size = 11) {
  theme_minimal(base_family = FUENTE_TFM, base_size = base_size) +
    theme(
      plot.title = element_text(family = FUENTE_TFM, face = "bold", color = COLOR_NEGRO, hjust = 0),
      plot.subtitle = element_text(family = FUENTE_TFM, color = COLOR_GRIS),
      axis.title = element_text(family = FUENTE_TFM, face = "bold", color = COLOR_NEGRO),
      axis.text = element_text(family = FUENTE_TFM, color = COLOR_NEGRO),
      legend.title = element_text(family = FUENTE_TFM, face = "bold"),
      legend.text = element_text(family = FUENTE_TFM),
      strip.text = element_text(family = FUENTE_TFM, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(12, 12, 12, 12)
    )
}

vcov_cluster <- function(modelo, cluster) {
  sandwich::vcovCL(modelo, cluster = cluster, type = "HC1")
}

extraer_coef_robusto <- function(modelo, cluster, termino) {
  b <- coef(modelo)

  if (!(termino %in% names(b)) || is.na(b[termino])) {
    return(NULL)
  }

  V <- vcov_cluster(modelo, cluster)

  if (!(termino %in% rownames(V))) {
    return(NULL)
  }

  beta <- unname(b[termino])
  se <- sqrt(V[termino, termino])

  if (!is.finite(se) || se <= 0) {
    return(NULL)
  }

  n_clusters <- length(unique(cluster[!is.na(cluster)]))
  gl <- max(n_clusters - 1, 1)
  crit <- qt(0.975, df = gl)
  estadistico <- beta / se
  p <- 2 * pt(abs(estadistico), df = gl, lower.tail = FALSE)

  tibble(
    termino = termino,
    beta = beta,
    se = se,
    IC95_inf = beta - crit * se,
    IC95_sup = beta + crit * se,
    estadistico = estadistico,
    gl_cluster = gl,
    p = p
  )
}

wald_robusto <- function(modelo, cluster, terminos) {
  terminos <- intersect(terminos, names(coef(modelo)))
  terminos <- terminos[!is.na(coef(modelo)[terminos])]

  if (length(terminos) == 0) {
    return(list(chi2 = NA_real_, df = NA_integer_, p = NA_real_))
  }

  V <- vcov_cluster(modelo, cluster)
  V <- V[terminos, terminos, drop = FALSE]
  b <- coef(modelo)[terminos]

  invV <- tryCatch(
    solve(V),
    error = function(e) {
      solve(V + diag(1e-8, nrow(V)))
    }
  )

  chi2 <- as.numeric(t(b) %*% invV %*% b)

  list(chi2 = chi2, df = length(terminos), p = pchisq(chi2, df = length(terminos), lower.tail = FALSE))
}

normalizar_genero <- function(x) {
  y <- tolower(trimws(as.character(x)))
  gsub("[^[:alnum:]]", "", y)
}

# Rutas de entrada.

ruta_master <- ruta_existente(
  c(
    file.path(PROJECT_DIR, "tablas_maestras", "TABLA_MAESTRA_PACIENTES_205_COMPLETA.csv"),
    file.path(PROJECT_DIR, "outputs", "datos_procesados_finales", "TABLA_MAESTRA_PACIENTES_205_COMPLETA.csv")
  ),
  "tabla maestra de 205 pacientes"
)

ruta_ps_amb <- ruta_existente(
  c(
    file.path(
      PROJECT_DIR, "outputs", "analisis_microbiota",
      "reanalisis_completo_MITECO", "objetos",
      "ps_ambiental_205_filtrado_con_MITECO.rds"
    ),
    file.path(
      PROJECT_DIR, "outputs", "analisis_microbiota",
      "objetos_phyloseq",
      "ps_ambiental_205_filtrado_prev2_total10.rds"
    )
  ),
  "objeto phyloseq ambiental filtrado"
)

ruta_ps_saliva206 <- ruta_existente(
  c(
    file.path(
      PROJECT_DIR, "outputs", "analisis_microbiota",
      "objetos_phyloseq",
      "ps_saliva_206_sin_filtrado.rds"
    ),
    file.path(PROJECT_DIR, "objetos_para_rocio", "ps_saliva_206_sin_filtrado.rds")
  ),
  "objeto phyloseq de 206 salivas"
)

ruta_wunifrac <- ruta_existente(
  c(
    file.path(
      PROJECT_DIR, "outputs", "analisis_microbiota",
      "unifrac_sensibilidad", "objetos",
      "dist_weighted_unifrac_ambiental_205.rds"
    )
  ),
  "distancia weighted UniFrac",
  obligatorio = FALSE
)

ruta_uunifrac <- ruta_existente(
  c(
    file.path(
      PROJECT_DIR, "outputs", "analisis_microbiota",
      "unifrac_sensibilidad", "objetos",
      "dist_unweighted_unifrac_ambiental_205_rarefactado_10000.rds"
    )
  ),
  "distancia unweighted UniFrac",
  obligatorio = FALSE
)

ruta_ancom <- ruta_existente(
  c(
    file.path(
      PROJECT_DIR, "outputs", "analisis_microbiota",
      "reanalisis_completo_MITECO", "tablas",
      "23_ANCOMBC2_Genus_todos.csv"
    )
  ),
  "resultados completos ANCOM-BC2"
)

ruta_comp_genus <- ruta_existente(
  c(
    file.path(
      PROJECT_DIR, "outputs", "analisis_microbiota",
      "composicion_taxonomica_v3", "tablas",
      "03_resumen_taxonomico_genus.csv"
    )
  ),
  "resumen taxonómico a género"
)

ruta_prev_genus <- ruta_existente(
  c(
    file.path(
      PROJECT_DIR, "outputs", "analisis_microbiota",
      "composicion_taxonomica_v3", "tablas",
      "13_prevalencia_generos_saliva206.csv"
    )
  ),
  "tabla de prevalencia de géneros"
)

rutas_deseq <- c(
  MITECO10 = file.path(PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION", "DESEQ2_MITECO10.csv"),
  MITECO_dias10 = file.path(
    PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION",
    "DESEQ2_MITECO_DIAS10_Z.csv"
  ),
  NAF10 = file.path(PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION", "DESEQ2_NAF10_Z.csv"),
  MITECO10_x_sexo = file.path(
    PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION",
    "DESEQ2_MITECO10_X_SEXO.csv"
  ),
  MITECO10_x_exac6m = file.path(
    PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION",
    "DESEQ2_MITECO10_X_EXAC6M.csv"
  )
)

if (!all(file.exists(rutas_deseq))) {
  faltantes_deseq <- names(rutas_deseq)[!file.exists(rutas_deseq)]

  stop(
    "Faltan resultados DESeq2: ",
    paste(faltantes_deseq, collapse = ", "),
    "\nEjecuta primero el script 20."
  )
}

# Carga y control de calidad.

master <- leer_csv(ruta_master)

columnas_requeridas <- c(
  "id_muestra",
  "age",
  "sex",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m",
  "fecha",
  "year",
  "doy",
  "sin_doy",
  "cos_doy",
  "cluster_isla_fecha",
  "miteco10",
  "miteco_days_10d",
  "miteco_any_30d",
  "miteco_days_30d",
  "miteco_any_90d",
  "miteco_days_90d",
  "miteco_any_180d",
  "miteco_days_180d",
  "pm10_media_10d",
  "pm10_auc50_10d",
  "n_dias_sdi_10d",
  "pm10_media_180d",
  "pm10_auc50_180d",
  "lecturas_totales",
  "lecturas_log10",
  "Observed",
  "Shannon",
  "Simpson",
  "aitchison_centro_isla"
)

faltan_columnas <- setdiff(columnas_requeridas, names(master))

if (length(faltan_columnas) > 0) {
  stop("Faltan columnas en la tabla maestra:\n- ", paste(faltan_columnas, collapse = "\n- "))
}

if (anyDuplicated(master$id_muestra) > 0) {
  stop("La tabla maestra contiene identificadores duplicados.")
}

if (nrow(master) != 205) {
  stop("Se esperaban 205 pacientes y se encontraron ", nrow(master), ".")
}

meta <- master %>%
  mutate(
    id_muestra = as.character(id_muestra),
    age = as.numeric(age),
    sex = factor(as.character(sex), levels = c("Female", "Male")),
    island = factor(as.character(island), levels = c("La Palma", "Tenerife")),
    pool = factor(pool),
    antibiotic = factor(as.numeric(antibiotic), levels = c(0, 1)),
    cavities = factor(as.numeric(cavities), levels = c(0, 1)),
    gum = factor(as.numeric(gum), levels = c(0, 1)),
    exac6m_num = as.numeric(exac6m),
    exac6m = factor(as.numeric(exac6m), levels = c(0, 1)),
    fecha = as.Date(fecha),
    year = factor(year),
    cluster_isla_fecha = factor(cluster_isla_fecha),
    miteco10 = as.numeric(miteco10),
    intrusion_factor = factor(
      ifelse(miteco10 == 1, "Con intrusión oficial", "Sin intrusión oficial"),
      levels = c("Sin intrusión oficial", "Con intrusión oficial")
    ),
    sexo_es = factor(ifelse(sex == "Male", "Hombres", "Mujeres"), levels = c("Mujeres", "Hombres")),
    menor_18 = age < 18,
    pm10_media_10d_z = zscore_seguro(pm10_media_10d),
    pm10_auc50_10d_z = zscore_seguro(pm10_auc50_10d),
    pm10_media_180d_z = zscore_seguro(pm10_media_180d),
    pm10_auc50_180d_z = zscore_seguro(pm10_auc50_180d)
  )

# IMPORTANTE:
# readr/dplyr devuelven un tibble. Los objetos dist y vegan::adonis2
# se alinean de forma más segura usando row.names de un data.frame base.
# Evitamos así los avisos "Setting row names on a tibble is deprecated"
# y, sobre todo, que se pierdan los IDs al submuestrear las distancias.
meta <- as.data.frame(meta, stringsAsFactors = FALSE)
rownames(meta) <- meta$id_muestra

ps_amb <- readRDS(ruta_ps_amb)
ps_saliva206 <- readRDS(ruta_ps_saliva206)

if (!inherits(ps_amb, "phyloseq")) {
  stop("El objeto ambiental no es de clase phyloseq.")
}

if (!inherits(ps_saliva206, "phyloseq")) {
  stop("El objeto de salivas no es de clase phyloseq.")
}

ids_comunes <- intersect(sample_names(ps_amb), meta$id_muestra)

if (length(ids_comunes) != 205) {
  stop(
    "El objeto phyloseq y la tabla maestra solo comparten ",
    length(ids_comunes), " identificadores; se esperaban 205."
  )
}

ps_amb <- prune_samples(meta$id_muestra, ps_amb)
meta <- meta[
  match(sample_names(ps_amb), meta$id_muestra),
  ,
  drop = FALSE
]

if (!identical(sample_names(ps_amb), meta$id_muestra)) {
  stop("No se pudo alinear el phyloseq con la tabla maestra.")
}

registrar(
  "Carga de datos",
  "OK",
  paste0(
    "205 pacientes; ",
    nsamples(ps_amb), " muestras phyloseq; ",
    ntaxa(ps_amb), " ASV/taxa filtrados."
  )
)

msg("Datos cargados correctamente.")
msg("  Pacientes ambientales: %d", nrow(meta))
msg("  Clusters isla-fecha: %d", n_distinct(meta$cluster_isla_fecha))
msg("  Taxa objeto filtrado: %d", ntaxa(ps_amb))

# Tabla descriptiva de la cohorte.

resumir_columna <- function(data) {
  tibble(
    variable = c(
      "Edad, años",
      "Sexo femenino",
      "Edad menor de 18 años",
      "Tenerife",
      "Exacerbación en los 6 meses previos",
      "Antibióticos en los 2 meses previos",
      "Caries",
      "Consumo de chicle en los 30 minutos previos",
      "Lecturas por muestra",
      "PM10 media D-1 a D-10, µg/m³"
    ),
    valor = c(
      formatear_mediana_ric(data$age, 1),
      formatear_n_pct(data$sex, function(z) z == "Female"),
      formatear_n_pct(data$menor_18, function(z) z),
      formatear_n_pct(data$island, function(z) z == "Tenerife"),
      formatear_n_pct(data$exac6m_num, function(z) z == 1),
      formatear_n_pct(as.numeric(as.character(data$antibiotic)), function(z) z == 1),
      formatear_n_pct(as.numeric(as.character(data$cavities)), function(z) z == 1),
      formatear_n_pct(as.numeric(as.character(data$gum)), function(z) z == 1),
      formatear_mediana_ric(data$lecturas_totales, 0),
      formatear_mediana_ric(data$pm10_media_10d, 1)
    )
  )
}

n_sin <- sum(meta$miteco10 == 0, na.rm = TRUE)
n_con <- sum(meta$miteco10 == 1, na.rm = TRUE)

tabla_cohorte <- resumir_columna(meta) %>%
  rename(`Total (n = 205)` = valor) %>%
  left_join(
    resumir_columna(meta %>% filter(miteco10 == 0)) %>%
      rename(!!paste0("Sin intrusión oficial D-1 a D-10 (n = ", n_sin, ")") := valor),
    by = "variable"
  ) %>%
  left_join(
    resumir_columna(meta %>% filter(miteco10 == 1)) %>%
      rename(!!paste0("Con intrusión oficial D-1 a D-10 (n = ", n_con, ")") := valor),
    by = "variable"
  )

datos_faltantes <- tibble(
  variable = c("Antibióticos", "Caries", "Consumo de chicle"),
  n_ausentes = c(sum(is.na(meta$antibiotic)), sum(is.na(meta$cavities)), sum(is.na(meta$gum)))
)

