# _CIERRE_RESULTADOS_TFM_GEMAS_v3.R
# TFM GEMAS / INTRUSIONES DE POLVO SAHARIANO

# Este script genera una carpeta independiente con los análisis y resultados:
#
#  1. Caracterización de la exposición oficial por ventanas.
#  2. Concordancia entre intrusión oficial y PM10 superficial.
#  3. Reanálisis armonizado de PM10 superficial:
#       - diversidad alfa D-1 a D-10;
#       - interacción PM10 x sexo para Shannon y Simpson;
#       - diversidad beta mediante Bray-Curtis y UniFrac;
#       - sensibilidad agregada por isla-fecha;
#       - carga a 180 días y fenotipo exacerbador.
#  4. Exportación completa de las PERMANOVA marginales.
#  5. Síntesis ANCOM-BC2 frente a DESeq2.
#  6. Tablas y figuras candidatas para el manuscrito.
#
# SALIDA
# outputs/22_CIERRE_RESULTADOS_TFM/
#   ├── tablas/
#   ├── figuras/
#   ├── objetos/
#   └── logs/

# 0. CONFIGURACIÓN


options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260830)

PROJECT_DIR <- "RUTA_LOCAL_NO_DISTRIBUIDA"

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


# 1. PAQUETES

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


# 2. FUNCIONES 


msg <- function(...) {
  cat(sprintf(...), "\n")
}

registrar <- local({
  x <- tibble(
    etapa = character(),
    estado = character(),
    detalle = character()
  )

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
      msg(
        "AVISO: varias rutas válidas para '%s'. Se usará: %s",
        etiqueta,
        ok[1]
      )
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
  readr::read_csv(
    ruta,
    show_col_types = FALSE,
    na = c("", "NA", "NaN", "NULL", "null")
  )
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

  q <- quantile(
    x,
    probs = c(0.25, 0.50, 0.75),
    names = FALSE,
    na.rm = TRUE
  )

  sprintf(
    paste0("%.", digitos, "f [%.", digitos, "f–%.", digitos, "f]"),
    q[2], q[1], q[3]
  )
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
      plot.title = element_text(
        family = FUENTE_TFM,
        face = "bold",
        color = COLOR_NEGRO,
        hjust = 0
      ),
      plot.subtitle = element_text(
        family = FUENTE_TFM,
        color = COLOR_GRIS
      ),
      axis.title = element_text(
        family = FUENTE_TFM,
        face = "bold",
        color = COLOR_NEGRO
      ),
      axis.text = element_text(
        family = FUENTE_TFM,
        color = COLOR_NEGRO
      ),
      legend.title = element_text(
        family = FUENTE_TFM,
        face = "bold"
      ),
      legend.text = element_text(
        family = FUENTE_TFM
      ),
      strip.text = element_text(
        family = FUENTE_TFM,
        face = "bold"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(12, 12, 12, 12)
    )
}

vcov_cluster <- function(modelo, cluster) {
  sandwich::vcovCL(
    modelo,
    cluster = cluster,
    type = "HC1"
  )
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

  list(
    chi2 = chi2,
    df = length(terminos),
    p = pchisq(chi2, df = length(terminos), lower.tail = FALSE)
  )
}

normalizar_genero <- function(x) {
  y <- tolower(trimws(as.character(x)))
  gsub("[^[:alnum:]]", "", y)
}


# 3. RUTAS DE ENTRADA

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
  MITECO10 = file.path(
    PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION",
    "DESEQ2_MITECO10.csv"
  ),
  MITECO_dias10 = file.path(
    PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION",
    "DESEQ2_MITECO_DIAS10_Z.csv"
  ),
  NAF10 = file.path(
    PROJECT_DIR, "outputs", "20_DESEQ2_VALIDACION",
    "DESEQ2_NAF10_Z.csv"
  ),
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


# 4. CARGA Y CONTROL DE CALIDAD

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
  stop(
    "Faltan columnas en la tabla maestra:\n- ",
    paste(faltan_columnas, collapse = "\n- ")
  )
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
    sex = factor(
      as.character(sex),
      levels = c("Female", "Male")
    ),
    island = factor(
      as.character(island),
      levels = c("La Palma", "Tenerife")
    ),
    pool = factor(pool),
    antibiotic = factor(
      as.numeric(antibiotic),
      levels = c(0, 1)
    ),
    cavities = factor(
      as.numeric(cavities),
      levels = c(0, 1)
    ),
    gum = factor(
      as.numeric(gum),
      levels = c(0, 1)
    ),
    exac6m_num = as.numeric(exac6m),
    exac6m = factor(
      as.numeric(exac6m),
      levels = c(0, 1)
    ),
    fecha = as.Date(fecha),
    year = factor(year),
    cluster_isla_fecha = factor(cluster_isla_fecha),
    miteco10 = as.numeric(miteco10),
    intrusion_factor = factor(
      ifelse(
        miteco10 == 1,
        "Con intrusión oficial",
        "Sin intrusión oficial"
      ),
      levels = c(
        "Sin intrusión oficial",
        "Con intrusión oficial"
      )
    ),
    sexo_es = factor(
      ifelse(
        sex == "Male",
        "Hombres",
        "Mujeres"
      ),
      levels = c("Mujeres", "Hombres")
    ),
    menor_18 = age < 18,
    pm10_media_10d_z = zscore_seguro(pm10_media_10d),
    pm10_auc50_10d_z = zscore_seguro(pm10_auc50_10d),
    pm10_media_180d_z = zscore_seguro(pm10_media_180d),
    pm10_auc50_180d_z = zscore_seguro(pm10_auc50_180d)
  )

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


# 5. TABLA DESCRIPTIVA DE LA COHORTE

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
      formatear_n_pct(
        as.numeric(as.character(data$antibiotic)),
        function(z) z == 1
      ),
      formatear_n_pct(
        as.numeric(as.character(data$cavities)),
        function(z) z == 1
      ),
      formatear_n_pct(
        as.numeric(as.character(data$gum)),
        function(z) z == 1
      ),
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
      rename(
        !!paste0("Sin intrusión oficial D-1 a D-10 (n = ", n_sin, ")") := valor
      ),
    by = "variable"
  ) %>%
  left_join(
    resumir_columna(meta %>% filter(miteco10 == 1)) %>%
      rename(
        !!paste0("Con intrusión oficial D-1 a D-10 (n = ", n_con, ")") := valor
      ),
    by = "variable"
  )

datos_faltantes <- tibble(
  variable = c(
    "Antibióticos",
    "Caries",
    "Consumo de chicle"
  ),
  n_ausentes = c(
    sum(is.na(meta$antibiotic)),
    sum(is.na(meta$cavities)),
    sum(is.na(meta$gum))
  )
)


# 6. CARACTERIZACIÓN AMBIENTAL

ventanas <- tibble(
  ventana = factor(
    c(
      "D-1 a D-10",
      "D-1 a D-30",
      "D-1 a D-90",
      "D-1 a D-180"
    ),
    levels = c(
      "D-1 a D-10",
      "D-1 a D-30",
      "D-1 a D-90",
      "D-1 a D-180"
    )
  ),
  variable_any = c(
    "miteco10",
    "miteco_any_30d",
    "miteco_any_90d",
    "miteco_any_180d"
  ),
  variable_days = c(
    "miteco_days_10d",
    "miteco_days_30d",
    "miteco_days_90d",
    "miteco_days_180d"
  )
)

resumen_ventanas <- bind_rows(
  lapply(
    seq_len(nrow(ventanas)),
    function(i) {
      va <- ventanas$variable_any[i]
      vd <- ventanas$variable_days[i]

      exp_bin <- as.numeric(meta[[va]])
      dias <- as.numeric(meta[[vd]])
      expuestos <- exp_bin == 1

      tibble(
        ventana = ventanas$ventana[i],
        n_total = sum(!is.na(exp_bin)),
        n_expuestos = sum(expuestos, na.rm = TRUE),
        porcentaje_expuestos = 100 * mean(expuestos, na.rm = TRUE),
        dias_mediana_total = median(dias, na.rm = TRUE),
        dias_Q1_total = quantile(dias, 0.25, na.rm = TRUE),
        dias_Q3_total = quantile(dias, 0.75, na.rm = TRUE),
        dias_mediana_expuestos = median(dias[expuestos], na.rm = TRUE),
        dias_Q1_expuestos = quantile(dias[expuestos], 0.25, na.rm = TRUE),
        dias_Q3_expuestos = quantile(dias[expuestos], 0.75, na.rm = TRUE)
      )
    }
  )
)

pm10_alto_10d <- meta$n_dias_sdi_10d > 0

tabla_concordancia <- meta %>%
  mutate(
    categoria = case_when(
      miteco10 == 1 & pm10_alto_10d ~
        "Intrusión oficial y PM10 ≥50 µg/m³",
      miteco10 == 1 & !pm10_alto_10d ~
        "Intrusión oficial sin PM10 ≥50 µg/m³",
      miteco10 == 0 & pm10_alto_10d ~
        "PM10 ≥50 µg/m³ sin intrusión oficial",
      miteco10 == 0 & !pm10_alto_10d ~
        "Ausencia de ambas condiciones",
      TRUE ~ NA_character_
    )
  ) %>%
  count(categoria, name = "n") %>%
  mutate(
    porcentaje = 100 * n / sum(n)
  ) %>%
  arrange(
    factor(
      categoria,
      levels = c(
        "Intrusión oficial y PM10 ≥50 µg/m³",
        "Intrusión oficial sin PM10 ≥50 µg/m³",
        "PM10 ≥50 µg/m³ sin intrusión oficial",
        "Ausencia de ambas condiciones"
      )
    )
  )

# La exposición es compartida por los pacientes muestreados en la misma
# isla y fecha. Para describir la relación con PM10 se utiliza una fila
# por combinación isla-fecha.
qc_constancia_cluster <- meta %>%
  group_by(cluster_isla_fecha) %>%
  summarise(
    n_pm10_media = n_distinct(pm10_media_10d),
    n_pm10_auc50 = n_distinct(pm10_auc50_10d),
    n_dias_miteco = n_distinct(miteco_days_10d),
    .groups = "drop"
  )

if (
  any(qc_constancia_cluster$n_pm10_media > 1) ||
    any(qc_constancia_cluster$n_pm10_auc50 > 1) ||
    any(qc_constancia_cluster$n_dias_miteco > 1)
) {
  stop(
    "Alguna exposición ambiental no es constante dentro de un cluster isla-fecha."
  )
}

env_cluster <- meta %>%
  group_by(cluster_isla_fecha) %>%
  summarise(
    island = primero_no_na(as.character(island)),
    fecha = primero_no_na(fecha),
    n_muestras = n(),
    miteco10 = primero_no_na(miteco10),
    miteco_days_10d = primero_no_na(miteco_days_10d),
    pm10_media_10d = primero_no_na(pm10_media_10d),
    pm10_auc50_10d = primero_no_na(pm10_auc50_10d),
    .groups = "drop"
  ) %>%
  mutate(
    island = factor(
      island,
      levels = c("La Palma", "Tenerife")
    )
  )

calcular_correlacion <- function(data, x, y, estrato = "Total") {
  dd <- data %>%
    filter(!is.na(.data[[x]]), !is.na(.data[[y]]))

  if (nrow(dd) < 4) {
    return(
      tibble(
        estrato = estrato,
        x = x,
        y = y,
        n = nrow(dd),
        rho_spearman = NA_real_,
        p = NA_real_
      )
    )
  }

  ct <- suppressWarnings(
    cor.test(
      dd[[x]],
      dd[[y]],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    estrato = estrato,
    x = x,
    y = y,
    n = nrow(dd),
    rho_spearman = unname(ct$estimate),
    p = ct$p.value
  )
}

correlaciones_ambientales <- bind_rows(
  calcular_correlacion(
    env_cluster,
    "miteco_days_10d",
    "pm10_media_10d",
    "Total"
  ),
  calcular_correlacion(
    env_cluster,
    "miteco_days_10d",
    "pm10_auc50_10d",
    "Total"
  ),
  bind_rows(
    lapply(
      levels(env_cluster$island),
      function(isl) {
        bind_rows(
          calcular_correlacion(
            env_cluster %>% filter(island == isl),
            "miteco_days_10d",
            "pm10_media_10d",
            isl
          ),
          calcular_correlacion(
            env_cluster %>% filter(island == isl),
            "miteco_days_10d",
            "pm10_auc50_10d",
            isl
          )
        )
      }
    )
  )
)

registrar(
  "Caracterización ambiental",
  "OK",
  paste0(
    nrow(env_cluster),
    " combinaciones isla-fecha; ",
    sum(meta$miteco10 == 1), "/205 con intrusión D-1 a D-10."
  )
)


# 7. DISTANCIAS MICROBIOLÓGICAS

ps_rel <- transform_sample_counts(
  ps_amb,
  function(x) x / sum(x)
)

dist_bray <- phyloseq::distance(
  ps_rel,
  method = "bray"
)

distancias <- list(
  Bray_Curtis = dist_bray
)

if (!is.na(ruta_wunifrac) && file.exists(ruta_wunifrac)) {
  distancias$Weighted_UniFrac <- readRDS(ruta_wunifrac)
}

if (!is.na(ruta_uunifrac) && file.exists(ruta_uunifrac)) {
  distancias$Unweighted_UniFrac_rarefactado <- readRDS(ruta_uunifrac)
}

for (nm in names(distancias)) {
  n_comunes <- length(
    intersect(labels(distancias[[nm]]), meta$id_muestra)
  )

  if (n_comunes < 200) {
    warning(
      "La distancia ", nm, " solo contiene ", n_comunes,
      " muestras comunes con la tabla maestra."
    )
  }
}

registrar(
  "Distancias microbiológicas",
  "OK",
  paste(names(distancias), collapse = ", ")
)


qc_distancias <- dplyr::bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      labs <- dist_labels(distancias[[nm]])
      tibble(
        distancia = nm,
        n_etiquetas_distancia = length(labs),
        n_ids_metadata = nrow(meta),
        n_ids_comunes = length(intersect(labs, meta$id_muestra)),
        ejemplo_etiqueta_distancia = ifelse(length(labs) > 0, labs[1], NA_character_),
        ejemplo_id_metadata = ifelse(nrow(meta) > 0, meta$id_muestra[1], NA_character_)
      )
    }
  )
)

write_csv(
  qc_distancias,
  file.path(TAB_DIR, "00_QC_alineacion_distancias_metadata.csv"),
  na = ""
)

msg("QC de alineación de distancias:")
print(qc_distancias)


# 8. DIVERSIDAD ALFA: PM10 ARMONIZADO

alpha_outcomes <- c(
  "Observed",
  "Shannon",
  "Simpson",
  "aitchison_centro_isla"
)

covars_alpha <- c(
  "age",
  "sex",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m",
  "lecturas_log10",
  "year",
  "sin_doy",
  "cos_doy"
)

ajustar_alpha <- function(
  outcome,
  exposure,
  data = meta,
  covars = covars_alpha
) {
  vars <- unique(
    c(
      outcome,
      exposure,
      covars,
      "cluster_isla_fecha"
    )
  )

  dd <- data[
    complete.cases(data[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)
  cov2 <- covariables_activas(dd, covars)

  if (n_niveles(dd[[exposure]]) < 2) {
    return(NULL)
  }

  formula_txt <- paste(
    outcome,
    "~",
    paste(c(exposure, cov2), collapse = " + ")
  )

  fit <- lm(
    as.formula(formula_txt),
    data = dd
  )

  z <- extraer_coef_robusto(
    fit,
    dd$cluster_isla_fecha,
    exposure
  )

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    mutate(
      outcome = outcome,
      exposure = exposure,
      n = nrow(dd),
      n_clusters = n_distinct(dd$cluster_isla_fecha),
      formula = formula_txt,
      .before = 1
    )
}

exposiciones_pm10_10d <- c(
  "pm10_media_10d_z",
  "pm10_auc50_10d_z"
)

alpha_pm10_10d <- bind_rows(
  lapply(
    exposiciones_pm10_10d,
    function(e) {
      bind_rows(
        lapply(
          alpha_outcomes,
          function(o) ajustar_alpha(o, e)
        )
      )
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(
    q_BH_outcomes = p.adjust(p, method = "BH")
  ) %>%
  ungroup()

exposiciones_pm10_180d <- c(
  "pm10_media_180d_z",
  "pm10_auc50_180d_z"
)

alpha_pm10_180d <- bind_rows(
  lapply(
    exposiciones_pm10_180d,
    function(e) {
      bind_rows(
        lapply(
          alpha_outcomes,
          function(o) ajustar_alpha(o, e)
        )
      )
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(
    q_BH_outcomes = p.adjust(p, method = "BH")
  ) %>%
  ungroup()


# 9. INTERACCIÓN PM10 x SEXO

ajustar_alpha_interaccion_sexo <- function(
  outcome,
  exposure,
  data = meta
) {
  covars <- setdiff(covars_alpha, "sex")

  vars <- unique(
    c(
      outcome,
      exposure,
      "sex",
      covars,
      "cluster_isla_fecha"
    )
  )

  dd <- data[
    complete.cases(data[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)
  cov2 <- covariables_activas(dd, covars)

  formula_txt <- paste(
    outcome,
    "~",
    paste0(exposure, " * sex"),
    if (length(cov2) > 0) {
      paste0(" + ", paste(cov2, collapse = " + "))
    } else {
      ""
    }
  )

  fit <- lm(
    as.formula(formula_txt),
    data = dd
  )

  terminos_int <- names(coef(fit))[
    grepl(":", names(coef(fit)), fixed = TRUE) &
      grepl(exposure, names(coef(fit)), fixed = TRUE) &
      grepl("sex", names(coef(fit)), fixed = TRUE)
  ]

  w <- wald_robusto(
    fit,
    dd$cluster_isla_fecha,
    terminos_int
  )

  tibble(
    outcome = outcome,
    exposure = exposure,
    modifier = "sex",
    n = nrow(dd),
    n_clusters = n_distinct(dd$cluster_isla_fecha),
    chi2 = w$chi2,
    df = w$df,
    p_interaccion = w$p,
    formula = formula_txt
  )
}

alpha_pm10_x_sexo <- bind_rows(
  lapply(
    exposiciones_pm10_10d,
    function(e) {
      bind_rows(
        lapply(
          c("Shannon", "Simpson"),
          function(o) ajustar_alpha_interaccion_sexo(o, e)
        )
      )
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(
    q_BH_outcomes = p.adjust(
      p_interaccion,
      method = "BH"
    )
  ) %>%
  ungroup()

ajustar_alpha_estratificado_sexo <- function(
  outcome,
  exposure,
  sexo,
  data = meta
) {
  covars <- setdiff(covars_alpha, "sex")

  dd0 <- data %>%
    filter(sex == sexo) %>%
    droplevels()

  z <- ajustar_alpha(
    outcome = outcome,
    exposure = exposure,
    data = dd0,
    covars = covars
  )

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    mutate(
      sex = sexo,
      .before = 1
    )
}

alpha_pm10_sexo_estrat <- bind_rows(
  lapply(
    exposiciones_pm10_10d,
    function(e) {
      bind_rows(
        lapply(
          c("Female", "Male"),
          function(s) {
            bind_rows(
              lapply(
                c("Shannon", "Simpson"),
                function(o) {
                  ajustar_alpha_estratificado_sexo(
                    o,
                    e,
                    s
                  )
                }
              )
            )
          }
        )
      )
    }
  )
)


# 10. PREDICCIONES AJUSTADAS: INTRUSIÓN x SEXO

predicciones_marginales <- function(outcome) {
  covars <- setdiff(covars_alpha, "sex")

  vars <- unique(
    c(
      outcome,
      "miteco10",
      "sex",
      covars,
      "cluster_isla_fecha"
    )
  )

  dd <- meta[
    complete.cases(meta[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)
  cov2 <- covariables_activas(dd, covars)

  formula_txt <- paste(
    outcome,
    "~ miteco10 * sex",
    if (length(cov2) > 0) {
      paste0(" + ", paste(cov2, collapse = " + "))
    } else {
      ""
    }
  )

  fit <- lm(
    as.formula(formula_txt),
    data = dd
  )

  V <- vcov_cluster(
    fit,
    dd$cluster_isla_fecha
  )

  b <- coef(fit)
  validos <- !is.na(b)
  b <- b[validos]
  V <- V[names(b), names(b), drop = FALSE]

  combos <- expand.grid(
    miteco10 = c(0, 1),
    sex = c("Female", "Male"),
    stringsAsFactors = FALSE
  )

  n_clusters <- n_distinct(dd$cluster_isla_fecha)
  crit <- qt(0.975, df = max(n_clusters - 1, 1))

  bind_rows(
    lapply(
      seq_len(nrow(combos)),
      function(i) {
        nd <- dd
        nd$miteco10 <- combos$miteco10[i]
        nd$sex <- factor(
          combos$sex[i],
          levels = levels(dd$sex)
        )

        X <- model.matrix(
          delete.response(terms(fit)),
          data = nd
        )

        faltan_X <- setdiff(names(b), colnames(X))

        if (length(faltan_X) > 0) {
          extra <- matrix(
            0,
            nrow = nrow(X),
            ncol = length(faltan_X),
            dimnames = list(NULL, faltan_X)
          )
          X <- cbind(X, extra)
        }

        X <- X[, names(b), drop = FALSE]
        L <- colMeans(X)

        estimacion <- as.numeric(sum(L * b))
        se <- sqrt(as.numeric(t(L) %*% V %*% L))

        tibble(
          outcome = outcome,
          miteco10 = combos$miteco10[i],
          exposicion = ifelse(
            combos$miteco10[i] == 1,
            "Con intrusión oficial",
            "Sin intrusión oficial"
          ),
          sex = combos$sex[i],
          sexo = ifelse(
            combos$sex[i] == "Male",
            "Hombres",
            "Mujeres"
          ),
          n_modelo = nrow(dd),
          n_clusters = n_clusters,
          estimacion_ajustada = estimacion,
          se = se,
          IC95_inf = estimacion - crit * se,
          IC95_sup = estimacion + crit * se
        )
      }
    )
  )
}

predicciones_sexo <- bind_rows(
  predicciones_marginales("Shannon"),
  predicciones_marginales("Simpson")
)


# 11. PERMANOVA COMPLETA

covars_beta <- c(
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

permanova_completa <- function(
  dist_obj,
  distancia,
  exposure,
  data = meta,
  covars = covars_beta,
  nperm = N_PERM
) {
  labs_dist <- dist_labels(dist_obj)

  ids <- intersect(
    labs_dist,
    as.character(data$id_muestra)
  )

  if (length(ids) < 4) {
    stop(
      "PERMANOVA: no se pudieron alinear las muestras para ",
      distancia, " / ", exposure, ". ",
      "IDs comunes encontrados: ", length(ids), "."
    )
  }

  idx <- match(ids, data$id_muestra)
  dd <- data[idx, , drop = FALSE]
  dd <- as.data.frame(dd, stringsAsFactors = FALSE)

  vars <- unique(c(exposure, covars))
  faltan_vars <- setdiff(vars, names(dd))

  if (length(faltan_vars) > 0) {
    stop(
      "PERMANOVA ", distancia, " / ", exposure,
      ": faltan variables: ",
      paste(faltan_vars, collapse = ", ")
    )
  }

  cc <- complete.cases(dd[, vars, drop = FALSE])
  dd <- dd[cc, , drop = FALSE]
  dd <- droplevels(dd)

  if (nrow(dd) < 4) {
    stop(
      "PERMANOVA ", distancia, " / ", exposure,
      ": quedan menos de 4 casos completos (n = ", nrow(dd), ")."
    )
  }

  # Submuestreamos la matriz de distancia con los IDs reales, no con
  # rownames de un tibble.
  d2 <- subdist(
    dist_obj,
    dd$id_muestra
  )

  # Reordenamos metadata exactamente como la matriz de distancia.
  labs_d2 <- dist_labels(d2)
  dd <- dd[
    match(labs_d2, dd$id_muestra),
    ,
    drop = FALSE
  ]
  rownames(dd) <- dd$id_muestra

  if (!identical(rownames(dd), labs_d2)) {
    stop(
      "PERMANOVA ", distancia, " / ", exposure,
      ": fallo de alineación entre metadata y matriz de distancia."
    )
  }

  cov2 <- covariables_activas(dd, covars)

  if (n_niveles(dd[[exposure]]) < 2) {
    warning(
      "PERMANOVA omitida para ", distancia, " / ", exposure,
      ": la exposición tiene menos de dos valores distintos tras casos completos."
    )
    return(NULL)
  }

  formula_txt <- paste(
    "d2 ~",
    paste(c(exposure, cov2), collapse = " + ")
  )

  set.seed(SEED)

  fit <- vegan::adonis2(
    as.formula(formula_txt),
    data = dd,
    permutations = nperm,
    by = "margin"
  )

  tab <- as.data.frame(fit)
  tab$term <- rownames(tab)

  # vegan puede cambiar ligeramente los nombres de columnas entre versiones.
  col_p <- intersect(c("Pr(>F)", "Pr..F."), names(tab))
  if (length(col_p) == 0) {
    stop(
      "PERMANOVA ", distancia, " / ", exposure,
      ": no se encontró la columna de p en adonis2. Columnas disponibles: ",
      paste(names(tab), collapse = ", ")
    )
  }
  col_p <- col_p[1]

  out <- tibble(
    distancia = rep(distancia, nrow(tab)),
    exposure = rep(exposure, nrow(tab)),
    n = rep(nrow(dd), nrow(tab)),
    term = as.character(tab$term),
    Df = as.numeric(tab$Df),
    SumOfSqs = as.numeric(tab$SumOfSqs),
    R2 = as.numeric(tab$R2),
    F = as.numeric(tab$F),
    p = as.numeric(tab[[col_p]]),
    formula = rep(formula_txt, nrow(tab))
  )

  idx_pruebas <- !out$term %in% c("Residual", "Total") &
    !is.na(out$p)

  out$q_BH_modelo <- NA_real_
  out$q_BH_modelo[idx_pruebas] <- p.adjust(
    out$p[idx_pruebas],
    method = "BH"
  )

  out
}

exposiciones_beta <- c(
  "miteco10",
  "pm10_media_10d_z",
  "pm10_auc50_10d_z"
)

resultados_permanova <- list()

for (nm in names(distancias)) {
  for (e in exposiciones_beta) {
    msg("PERMANOVA: %s ~ %s", nm, e)

    tmp <- tryCatch(
      permanova_completa(
        distancias[[nm]],
        nm,
        e
      ),
      error = function(err) {
        registrar(
          paste0("PERMANOVA ", nm, " / ", e),
          "ERROR",
          conditionMessage(err)
        )
        stop(err)
      }
    )

    if (!is.null(tmp) && nrow(tmp) > 0) {
      resultados_permanova[[paste(nm, e, sep = "__")]] <- tmp
    }
  }
}

permanova_todas <- dplyr::bind_rows(resultados_permanova)

if (
  nrow(permanova_todas) == 0 ||
  !all(c("term", "exposure", "distancia", "p") %in% names(permanova_todas))
) {
  stop(
    "No se generó una tabla PERMANOVA válida. ",
    "Columnas disponibles: ",
    paste(names(permanova_todas), collapse = ", "),
    ". Revisa los mensajes inmediatamente anteriores."
  )
}

beta_exposiciones <- permanova_todas %>%
  dplyr::filter(.data$term == .data$exposure) %>%
  dplyr::group_by(.data$exposure) %>%
  dplyr::mutate(
    q_BH_distancias = p.adjust(
      .data$p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

if (nrow(beta_exposiciones) == 0) {
  stop(
    "Las PERMANOVA se ejecutaron, pero no se localizaron las filas ",
    "correspondientes a las exposiciones. Términos disponibles: ",
    paste(unique(permanova_todas$term), collapse = ", ")
  )
}

permanova_bray_miteco <- permanova_todas %>%
  filter(
    distancia == "Bray_Curtis",
    exposure == "miteco10"
  )

permanova_bray_pm10 <- permanova_todas %>%
  filter(
    distancia == "Bray_Curtis",
    exposure == "pm10_media_10d_z"
  )


# 12. PERMDISP PARA LA EXPOSICIÓN BINARIA


ids_bray <- intersect(
  dist_labels(dist_bray),
  meta$id_muestra
)

grupo_bray <- meta[
  match(ids_bray, meta$id_muestra),
  "intrusion_factor",
  drop = TRUE
]

d_bray_disp <- subdist(dist_bray, ids_bray)
grupo_bray <- grupo_bray[
  match(dist_labels(d_bray_disp), ids_bray)
]

bd <- vegan::betadisper(
  d_bray_disp,
  grupo_bray
)

set.seed(SEED)

perm_bd <- vegan::permutest(
  bd,
  permutations = N_PERM
)

betadisper_miteco <- tibble(
  distancia = "Bray_Curtis",
  grupo = "Intrusión oficial D-1 a D-10",
  n = length(grupo_bray),
  F = perm_bd$tab[1, "F"],
  p_betadisper = perm_bd$tab[1, "Pr(>F)"]
)


# 13. SENSIBILIDAD BRAY AGREGADA POR ISLA-FECHA

otu_rel <- otu_taxa_en_filas(ps_rel)

clusters <- as.character(meta$cluster_isla_fecha)
names(clusters) <- meta$id_muestra

ids_otu <- colnames(otu_rel)
clusters_otu <- clusters[ids_otu]
niveles_cluster <- unique(clusters_otu)

perfil_cluster <- sapply(
  niveles_cluster,
  function(cl) {
    ids <- ids_otu[clusters_otu == cl]

    if (length(ids) == 1) {
      otu_rel[, ids, drop = TRUE]
    } else {
      rowMeans(
        otu_rel[, ids, drop = FALSE]
      )
    }
  }
)

colnames(perfil_cluster) <- niveles_cluster

dist_cluster_bray <- vegan::vegdist(
  t(perfil_cluster),
  method = "bray"
)

meta_cluster <- meta %>%
  group_by(cluster_isla_fecha) %>%
  summarise(
    n_muestras = n(),
    island = primero_no_na(as.character(island)),
    fecha = primero_no_na(fecha),
    miteco10 = primero_no_na(miteco10),
    pm10_media_10d_z = primero_no_na(pm10_media_10d_z),
    pm10_auc50_10d_z = primero_no_na(pm10_auc50_10d_z),
    year = primero_no_na(as.character(year)),
    sin_doy = primero_no_na(sin_doy),
    cos_doy = primero_no_na(cos_doy),
    age_mean = media_segura(age),
    prop_male = media_segura(as.numeric(sex == "Male")),
    prop_antibiotic = media_segura(
      as.numeric(as.character(antibiotic))
    ),
    prop_cavities = media_segura(
      as.numeric(as.character(cavities))
    ),
    prop_gum = media_segura(
      as.numeric(as.character(gum))
    ),
    prop_exac6m = media_segura(exac6m_num),
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

covars_cluster <- c(
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

permanova_cluster <- function(exposure) {
  ids <- dist_labels(dist_cluster_bray)
  dd <- meta_cluster[ids, , drop = FALSE]

  vars <- unique(c(exposure, covars_cluster))
  cc <- complete.cases(dd[, vars, drop = FALSE])
  dd <- droplevels(dd[cc, , drop = FALSE])

  d2 <- subdist(
    dist_cluster_bray,
    rownames(dd)
  )

  dd <- dd[
    dist_labels(d2),
    ,
    drop = FALSE
  ]

  cov2 <- covariables_activas(dd, covars_cluster)

  formula_txt <- paste(
    "d2 ~",
    paste(c(exposure, cov2), collapse = " + ")
  )

  set.seed(SEED)

  fit <- vegan::adonis2(
    as.formula(formula_txt),
    data = dd,
    permutations = N_PERM,
    by = "margin"
  )

  tab <- as.data.frame(fit)
  tab$term <- rownames(tab)
  fila <- tab[tab$term == exposure, , drop = FALSE]

  tibble(
    distancia = "Bray_cluster_isla_fecha",
    exposure = exposure,
    n_clusters = nrow(dd),
    R2 = fila$R2,
    F = fila$F,
    p = fila$`Pr(>F)`,
    formula = formula_txt
  )
}

beta_cluster <- bind_rows(
  lapply(
    c(
      "miteco10",
      "pm10_media_10d_z",
      "pm10_auc50_10d_z"
    ),
    permanova_cluster
  )
) %>%
  mutate(
    q_BH = p.adjust(
      p,
      method = "BH"
    )
  )


# 14. PM10 A 180 DÍAS Y FENOTIPO EXACERBADOR

ajustar_exacerbacion_outcome <- function(
  exposure,
  incluir_antibiotico = TRUE
) {
  covars <- c(
    "age",
    "sex",
    "island",
    "pool",
    "cavities",
    "gum",
    "year",
    "sin_doy",
    "cos_doy"
  )

  if (incluir_antibiotico) {
    covars <- c(covars, "antibiotic")
  }

  vars <- unique(
    c(
      "exac6m_num",
      exposure,
      covars,
      "cluster_isla_fecha"
    )
  )

  dd <- meta[
    complete.cases(meta[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)
  cov2 <- covariables_activas(dd, covars)

  formula_txt <- paste(
    "exac6m_num ~",
    paste(c(exposure, cov2), collapse = " + ")
  )

  fit <- glm(
    as.formula(formula_txt),
    family = binomial(),
    data = dd
  )

  z <- extraer_coef_robusto(
    fit,
    dd$cluster_isla_fecha,
    exposure
  )

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    transmute(
      exposure = exposure,
      modelo = ifelse(
        incluir_antibiotico,
        "Con ajuste por antibióticos",
        "Sin ajuste por antibióticos"
      ),
      n = nrow(dd),
      n_clusters = n_distinct(dd$cluster_isla_fecha),
      beta = beta,
      se = se,
      OR_por_1_DE = exp(beta),
      OR_IC95_inf = exp(IC95_inf),
      OR_IC95_sup = exp(IC95_sup),
      p = p,
      formula = formula_txt
    )
}

exac_pm10_180d <- bind_rows(
  lapply(
    exposiciones_pm10_180d,
    function(e) {
      bind_rows(
        ajustar_exacerbacion_outcome(e, TRUE),
        ajustar_exacerbacion_outcome(e, FALSE)
      )
    }
  )
) %>%
  group_by(modelo) %>%
  mutate(
    q_BH_exposiciones = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()


# 15. SÍNTESIS ANCOM-BC2 FRENTE A DESEQ2

ancom <- leer_csv(ruta_ancom)

deseq_lista <- lapply(
  names(rutas_deseq),
  function(nm) {
    leer_csv(rutas_deseq[[nm]]) %>%
      mutate(
        analisis = nm
      )
  }
)

names(deseq_lista) <- names(rutas_deseq)

taxones_apriori <- c(
  "Corynebacterium",
  "Dolosigranulum",
  "Moraxella",
  "Streptococcus",
  "Veillonella",
  "Neisseria",
  "Bifidobacterium",
  "Capnocytophaga",
  "Absconditabacteriales_(SR1)"
)

keys_apriori <- normalizar_genero(taxones_apriori)
key_tannerella <- normalizar_genero("Tannerella")

clasificar_senal <- function(
  ancom_q,
  ancom_robusto,
  deseq_fdr,
  concordante
) {
  ancom_sig <- !is.na(ancom_q) && ancom_q < ALPHA
  deseq_sig <- !is.na(deseq_fdr) && deseq_fdr < ALPHA
  robusto <- isTRUE(ancom_robusto)

  if (robusto) {
    return("Robusta según ANCOM-BC2")
  }

  if (ancom_sig && deseq_sig && isTRUE(concordante)) {
    return("Concordante en ambos métodos, no robusta a pseudoconteos")
  }

  if (ancom_sig && deseq_sig && !isTRUE(concordante)) {
    return("Significativa en ambos métodos con dirección discordante")
  }

  if (ancom_sig && !deseq_sig) {
    return("Señal ANCOM-BC2 no robusta y no confirmada por DESeq2")
  }

  if (!ancom_sig && deseq_sig) {
    return("Señal dependiente de DESeq2")
  }

  if (isTRUE(concordante)) {
    return("Dirección concordante sin evidencia tras FDR")
  }

  "Sin evidencia consistente"
}

integrar_metodos <- function(analisis_id) {
  a <- ancom %>%
    filter(analisis == analisis_id) %>%
    transmute(
      taxon_id,
      Genus_ANCOM = Genus,
      beta_ANCOM = beta,
      p_ANCOM = p,
      q_ANCOM = q,
      passed_sensitivity = passed_sensitivity,
      robusto_ANCOM = robusto
    )

  d <- deseq_lista[[analisis_id]] %>%
    transmute(
      taxon_id,
      Genus_DESeq2 = coalesce(Genus_label, Genus),
      baseMean_DESeq2 = baseMean,
      log2FC_DESeq2 = log2FC,
      p_DESeq2 = P.Value,
      FDR_DESeq2 = FDR
    )

  full_join(a, d, by = "taxon_id") %>%
    mutate(
      analisis = analisis_id,
      Genus = coalesce(Genus_ANCOM, Genus_DESeq2),
      genus_key = normalizar_genero(Genus),
      direccion_concordante = case_when(
        is.na(beta_ANCOM) | is.na(log2FC_DESeq2) ~ NA,
        beta_ANCOM == 0 | log2FC_DESeq2 == 0 ~ NA,
        TRUE ~ sign(beta_ANCOM) == sign(log2FC_DESeq2)
      ),
      taxon_apriori = genus_key %in% keys_apriori,
      tannerella_posthoc = genus_key == key_tannerella,
      significativo_ANCOM = !is.na(q_ANCOM) & q_ANCOM < ALPHA,
      significativo_DESeq2 = !is.na(FDR_DESeq2) & FDR_DESeq2 < ALPHA,
      seleccionado = taxon_apriori |
        tannerella_posthoc |
        significativo_ANCOM |
        significativo_DESeq2
    ) %>%
    rowwise() %>%
    mutate(
      clasificacion = clasificar_senal(
        q_ANCOM,
        robusto_ANCOM,
        FDR_DESeq2,
        direccion_concordante
      ),
      motivo_seleccion = paste(
        c(
          if (taxon_apriori) "A priori" else NULL,
          if (tannerella_posthoc) "Post hoc: Tannerella" else NULL,
          if (significativo_ANCOM) "q ANCOM-BC2 <0,05" else NULL,
          if (significativo_DESeq2) "FDR DESeq2 <0,05" else NULL
        ),
        collapse = "; "
      )
    ) %>%
    ungroup()
}

comparacion_metodos_completa <- bind_rows(
  lapply(
    names(rutas_deseq),
    integrar_metodos
  )
)

comparacion_metodos_seleccion <- comparacion_metodos_completa %>%
  filter(seleccionado) %>%
  arrange(
    factor(
      analisis,
      levels = names(rutas_deseq)
    ),
    desc(significativo_DESeq2),
    desc(significativo_ANCOM),
    Genus
  )

resumen_metodos <- comparacion_metodos_completa %>%
  group_by(analisis) %>%
  summarise(
    n_taxa_ANCOM = sum(!is.na(beta_ANCOM)),
    n_taxa_DESeq2 = sum(!is.na(log2FC_DESeq2)),
    n_q_ANCOM_005 = sum(significativo_ANCOM, na.rm = TRUE),
    n_robustos_ANCOM = sum(robusto_ANCOM, na.rm = TRUE),
    n_FDR_DESeq2_005 = sum(significativo_DESeq2, na.rm = TRUE),
    n_concordantes_entre_comunes = sum(
      direccion_concordante %in% TRUE,
      na.rm = TRUE
    ),
    n_comunes_con_direccion = sum(
      !is.na(direccion_concordante)
    ),
    porcentaje_concordancia_direccion = ifelse(
      n_comunes_con_direccion > 0,
      100 * n_concordantes_entre_comunes / n_comunes_con_direccion,
      NA_real_
    ),
    .groups = "drop"
  )

registrar(
  "Comparación ANCOM-BC2/DESeq2",
  "OK",
  paste0(
    nrow(comparacion_metodos_seleccion),
    " filas seleccionadas mediante reglas predefinidas."
  )
)


# 16. COMPOSICIÓN TAXONÓMICA Y MICROBIOMA NÚCLEO

composicion_genus <- leer_csv(ruta_comp_genus)
prevalencia_genus <- leer_csv(ruta_prev_genus)

top15_genus <- composicion_genus %>%
  arrange(desc(abundancia_relativa_media)) %>%
  slice_head(n = 15) %>%
  mutate(
    taxon = factor(
      taxon,
      levels = rev(taxon)
    )
  )

prevalencia_plot <- prevalencia_genus %>%
  mutate(
    core = case_when(
      prevalencia_prop >= 0.95 ~ "Núcleo ≥95 %",
      prevalencia_prop >= 0.90 ~ "Núcleo 90–<95 %",
      TRUE ~ "Fuera del núcleo"
    ),
    core = factor(
      core,
      levels = c(
        "Fuera del núcleo",
        "Núcleo 90–<95 %",
        "Núcleo ≥95 %"
      )
    ),
    etiquetar = taxon %in% c(
      head(
        composicion_genus %>%
          arrange(desc(abundancia_relativa_media)) %>%
          pull(taxon),
        8
      ),
      taxones_apriori
    )
  )

resumen_core <- tibble(
  umbral_prevalencia = c("≥90 %", "≥95 %"),
  n_generos = c(
    sum(prevalencia_genus$prevalencia_prop >= 0.90),
    sum(prevalencia_genus$prevalencia_prop >= 0.95)
  )
)


# 17. FIGURAS


# 17A. Frecuencia de exposición por ventanas
fig_ambiental_A <- ggplot(
  resumen_ventanas,
  aes(
    x = ventana,
    y = porcentaje_expuestos
  )
) +
  geom_col(
    fill = COLOR_ROJO,
    width = 0.66
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f %%\n(%d/%d)",
        porcentaje_expuestos,
        n_expuestos,
        n_total
      )
    ),
    vjust = -0.25,
    family = FUENTE_TFM,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 108),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, " %")
  ) +
  labs(
    title = "A. Frecuencia de exposición según la ventana temporal",
    x = NULL,
    y = "Participantes con al menos una intrusión oficial"
  ) +
  tema_tfm()

rho_media <- correlaciones_ambientales %>%
  filter(
    estrato == "Total",
    y == "pm10_media_10d"
  ) %>%
  slice(1)

fig_ambiental_B <- ggplot(
  env_cluster,
  aes(
    x = miteco_days_10d,
    y = pm10_media_10d,
    color = island
  )
) +
  geom_point(
    alpha = 0.75,
    size = 2.2
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.7,
    color = COLOR_NEGRO,
    fill = COLOR_GRIS_CLARO
  ) +
  scale_color_manual(
    values = c(
      "La Palma" = COLOR_GRIS,
      "Tenerife" = COLOR_ROJO
    )
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "Spearman ρ = %.2f; p = %s\nn = %d combinaciones isla-fecha",
      rho_media$rho_spearman,
      format.pval(rho_media$p, digits = 2, eps = 0.001),
      rho_media$n
    ),
    hjust = 1.05,
    vjust = 1.3,
    family = FUENTE_TFM,
    size = 3.5
  ) +
  labs(
    title = "B. Intrusión oficial y concentración superficial de PM10",
    x = "Número de días con intrusión oficial en D-1 a D-10",
    y = "PM10 media D-1 a D-10 (µg/m³)",
    color = "Isla"
  ) +
  tema_tfm()

ggsave(
  file.path(FIG_DIR, "Figura_ambiental_A_frecuencia_ventanas.png"),
  fig_ambiental_A,
  width = 7.2,
  height = 5.2,
  dpi = 320,
  bg = "white"
)

ggsave(
  file.path(FIG_DIR, "Figura_ambiental_B_MITECO_vs_PM10.png"),
  fig_ambiental_B,
  width = 7.2,
  height = 5.2,
  dpi = 320,
  bg = "white"
)

fig_ambiental <- patchwork::wrap_plots(
  list(fig_ambiental_A, fig_ambiental_B),
  ncol = 1,
  heights = c(0.85, 1.15)
)

ggsave(
  file.path(FIG_DIR, "Figura_2_exposicion_oficial_y_PM10_superficial.png"),
  fig_ambiental,
  width = 8.2,
  height = 10.3,
  dpi = 320,
  bg = "white"
)


# 17B. Composición taxonómica y núcleo
fig_comp_A <- ggplot(
  top15_genus,
  aes(
    x = abundancia_relativa_media * 100,
    y = taxon
  )
) +
  geom_col(
    fill = COLOR_ROJO,
    width = 0.72
  ) +
  labs(
    title = "A. Géneros con mayor abundancia relativa media",
    x = "Abundancia relativa media (%)",
    y = NULL
  ) +
  tema_tfm()

fig_comp_B <- ggplot(
  prevalencia_plot,
  aes(
    x = prevalencia_prop * 100,
    y = abundancia_relativa_media * 100,
    color = core
  )
) +
  geom_vline(
    xintercept = c(90, 95),
    linetype = c("dashed", "solid"),
    color = c(COLOR_GRIS, COLOR_ROJO),
    linewidth = 0.6
  ) +
  geom_point(
    alpha = 0.75,
    size = 2
  ) +
  scale_y_log10(
    labels = scales::label_number(
      accuracy = 0.001,
      suffix = " %"
    )
  ) +
  scale_color_manual(
    values = c(
      "Fuera del núcleo" = COLOR_GRIS_CLARO,
      "Núcleo 90–<95 %" = COLOR_GRIS,
      "Núcleo ≥95 %" = COLOR_ROJO
    )
  ) +
  labs(
    title = "B. Prevalencia y abundancia de los géneros bacterianos",
    x = "Prevalencia en las 206 muestras (%)",
    y = "Abundancia relativa media, escala logarítmica",
    color = "Clasificación"
  ) +
  tema_tfm()

datos_etiquetas <- prevalencia_plot %>%
  filter(etiquetar) %>%
  arrange(desc(abundancia_relativa_media)) %>%
  slice_head(n = 12)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  fig_comp_B <- fig_comp_B +
    ggrepel::geom_text_repel(
      data = datos_etiquetas,
      aes(label = taxon),
      family = FUENTE_TFM,
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      show.legend = FALSE
    )
} else {
  fig_comp_B <- fig_comp_B +
    geom_text(
      data = datos_etiquetas %>% slice_head(n = 8),
      aes(label = taxon),
      family = FUENTE_TFM,
      size = 2.8,
      nudge_x = 1.2,
      check_overlap = TRUE,
      show.legend = FALSE
    )
}

ggsave(
  file.path(FIG_DIR, "Figura_composicion_A_top15_generos.png"),
  fig_comp_A,
  width = 7.2,
  height = 5.8,
  dpi = 320,
  bg = "white"
)

ggsave(
  file.path(FIG_DIR, "Figura_composicion_B_prevalencia_core.png"),
  fig_comp_B,
  width = 7.2,
  height = 5.8,
  dpi = 320,
  bg = "white"
)

fig_comp <- patchwork::wrap_plots(
  list(fig_comp_A, fig_comp_B),
  ncol = 1,
  heights = c(0.9, 1.1)
)

ggsave(
  file.path(FIG_DIR, "Figura_3_composicion_taxonomica_y_microbioma_nucleo.png"),
  fig_comp,
  width = 8.2,
  height = 10.5,
  dpi = 320,
  bg = "white"
)


# 17C. Diversidad ajustada por sexo
predicciones_plot <- predicciones_sexo %>%
  mutate(
    sexo = factor(
      sexo,
      levels = c("Mujeres", "Hombres")
    ),
    exposicion = factor(
      exposicion,
      levels = c(
        "Sin intrusión oficial",
        "Con intrusión oficial"
      )
    )
  )

crear_fig_sexo <- function(outcome_name, ylab) {
  ggplot(
    predicciones_plot %>% filter(.data$outcome == outcome_name),
    aes(
      x = exposicion,
      y = estimacion_ajustada,
      color = sexo,
      group = sexo
    )
  ) +
    geom_line(
      linewidth = 0.75,
      position = position_dodge(width = 0.12)
    ) +
    geom_point(
      size = 2.8,
      position = position_dodge(width = 0.12)
    ) +
    geom_errorbar(
      aes(
        ymin = IC95_inf,
        ymax = IC95_sup
      ),
      width = 0.08,
      linewidth = 0.65,
      position = position_dodge(width = 0.12)
    ) +
    scale_color_manual(
      values = c(
        "Mujeres" = COLOR_GRIS,
        "Hombres" = COLOR_ROJO
      )
    ) +
    labs(
      title = outcome_name,
      x = NULL,
      y = ylab,
      color = "Sexo"
    ) +
    tema_tfm() +
    theme(
      axis.text.x = element_text(
        angle = 12,
        hjust = 1
      )
    )
}

fig_sexo_shannon <- crear_fig_sexo(
  "Shannon",
  "Índice de Shannon ajustado"
)

fig_sexo_simpson <- crear_fig_sexo(
  "Simpson",
  "Índice de Simpson ajustado"
)

ggsave(
  file.path(FIG_DIR, "Figura_sexo_A_Shannon_ajustado.png"),
  fig_sexo_shannon,
  width = 7.2,
  height = 5.2,
  dpi = 320,
  bg = "white"
)

ggsave(
  file.path(FIG_DIR, "Figura_sexo_B_Simpson_ajustado.png"),
  fig_sexo_simpson,
  width = 7.2,
  height = 5.2,
  dpi = 320,
  bg = "white"
)

fig_sexo <- patchwork::wrap_plots(
  list(fig_sexo_shannon, fig_sexo_simpson),
  ncol = 2,
  guides = "collect"
)

ggsave(
  file.path(FIG_DIR, "Figura_4_diversidad_ajustada_intrusion_por_sexo.png"),
  fig_sexo,
  width = 10.5,
  height = 5.7,
  dpi = 320,
  bg = "white"
)


# 17D. PCoA Bray-Curtis candidata
pcoa <- cmdscale(
  dist_bray,
  eig = TRUE,
  k = 2,
  add = TRUE
)

eig_positivos <- pcoa$eig[pcoa$eig > 0]
var1 <- ifelse(
  length(eig_positivos) > 0,
  100 * pcoa$eig[1] / sum(eig_positivos),
  NA_real_
)
var2 <- ifelse(
  length(eig_positivos) > 1,
  100 * pcoa$eig[2] / sum(eig_positivos),
  NA_real_
)

pcoa_df <- as.data.frame(pcoa$points)
names(pcoa_df)[1:2] <- c("PCoA1", "PCoA2")
pcoa_df$id_muestra <- rownames(pcoa_df)

pcoa_df <- pcoa_df %>%
  left_join(
    meta %>%
      select(
        id_muestra,
        intrusion_factor
      ),
    by = "id_muestra"
  )

fila_miteco_bray <- beta_exposiciones %>%
  filter(
    distancia == "Bray_Curtis",
    exposure == "miteco10"
  ) %>%
  slice(1)

fig_pcoa <- ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = intrusion_factor
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.9
  ) +
  stat_ellipse(
    type = "norm",
    linewidth = 0.7,
    level = 0.80,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Sin intrusión oficial" = COLOR_GRIS,
      "Con intrusión oficial" = COLOR_ROJO
    )
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.04,
    vjust = 1.25,
    family = FUENTE_TFM,
    size = 3.3,
    label = sprintf(
      "PERMANOVA ajustada: R² = %.3f; p = %s\nPERMDISP: p = %s",
      fila_miteco_bray$R2,
      format.pval(fila_miteco_bray$p, digits = 2, eps = 0.001),
      format.pval(
        betadisper_miteco$p_betadisper,
        digits = 2,
        eps = 0.001
      )
    )
  ) +
  labs(
    title = "PCoA de la diversidad beta de Bray-Curtis",
    x = sprintf("PCoA1 (%.1f %% de la variación)", var1),
    y = sprintf("PCoA2 (%.1f %% de la variación)", var2),
    color = "Exposición reciente"
  ) +
  tema_tfm()

ggsave(
  file.path(FIG_DIR, "Figura_candidata_PCoA_Bray_intrusion_oficial.png"),
  fig_pcoa,
  width = 7.5,
  height = 6.0,
  dpi = 320,
  bg = "white"
)


# 18. TABLA INTEGRADA DE INDICADORES AMBIENTALES

ruta_alpha_miteco <- file.path(
  PROJECT_DIR, "outputs", "analisis_microbiota",
  "reanalisis_completo_MITECO", "tablas",
  "02_alpha_global_MITECO.csv"
)

ruta_beta_miteco <- file.path(
  PROJECT_DIR, "outputs", "analisis_microbiota",
  "reanalisis_completo_MITECO", "tablas",
  "12_beta_global_MITECO.csv"
)

alpha_miteco_existente <- leer_csv(ruta_alpha_miteco)
beta_miteco_existente <- leer_csv(ruta_beta_miteco)

etiquetas_exposicion <- c(
  miteco10 = "Intrusión oficial D-1 a D-10",
  miteco_days10_z = "Número de días de intrusión D-1 a D-10, por 1 DE",
  naf10_z = "Aporte natural africano D-1 a D-10, por 1 DE",
  pm10_media_10d_z = "PM10 superficial media D-1 a D-10, por 1 DE",
  pm10_auc50_10d_z = "Carga superficial sobre 50 µg/m³ D-1 a D-10, por 1 DE"
)

alpha_integrada <- bind_rows(
  alpha_miteco_existente %>%
    select(
      outcome,
      exposure,
      n,
      beta,
      IC95_inf,
      IC95_sup,
      p,
      q = p_BH_outcomes
    ),
  alpha_pm10_10d %>%
    select(
      outcome,
      exposure,
      n,
      beta,
      IC95_inf,
      IC95_sup,
      p,
      q = q_BH_outcomes
    )
) %>%
  mutate(
    indicador = unname(etiquetas_exposicion[exposure]),
    tipo_resultado = "Diversidad alfa",
    metrica = outcome
  )

beta_integrada <- bind_rows(
  beta_miteco_existente %>%
    select(
      distancia,
      exposure,
      n,
      R2,
      F,
      p,
      q = p_BH_distancias
    ),
  beta_exposiciones %>%
    filter(
      exposure %in% c(
        "pm10_media_10d_z",
        "pm10_auc50_10d_z"
      )
    ) %>%
    select(
      distancia,
      exposure,
      n,
      R2,
      F,
      p,
      q = q_BH_distancias
    )
) %>%
  mutate(
    indicador = unname(etiquetas_exposicion[exposure]),
    tipo_resultado = "Diversidad beta",
    metrica = distancia
  )

tabla_indicadores_ambientales <- bind_rows(
  alpha_integrada %>%
    transmute(
      indicador,
      tipo_resultado,
      metrica,
      n,
      estimador = beta,
      IC95_inf,
      IC95_sup,
      F = NA_real_,
      p,
      q
    ),
  beta_integrada %>%
    transmute(
      indicador,
      tipo_resultado,
      metrica,
      n,
      estimador = R2,
      IC95_inf = NA_real_,
      IC95_sup = NA_real_,
      F,
      p,
      q
    )
)


# 19. EXPORTACIÓN

write_csv(
  tabla_cohorte,
  file.path(TAB_DIR, "01_caracteristicas_cohorte_segun_intrusion.csv"),
  na = ""
)

write_csv(
  datos_faltantes,
  file.path(TAB_DIR, "01b_datos_faltantes_tabla_cohorte.csv"),
  na = ""
)

write_csv(
  resumen_ventanas,
  file.path(TAB_DIR, "02_frecuencia_intrusion_por_ventanas.csv"),
  na = ""
)

write_csv(
  tabla_concordancia,
  file.path(TAB_DIR, "03_concordancia_intrusion_oficial_PM10.csv"),
  na = ""
)

write_csv(
  env_cluster,
  file.path(TAB_DIR, "04_exposicion_por_cluster_isla_fecha.csv"),
  na = ""
)

write_csv(
  correlaciones_ambientales,
  file.path(TAB_DIR, "05_correlaciones_intrusion_PM10_cluster.csv"),
  na = ""
)

write_csv(
  alpha_pm10_10d,
  file.path(TAB_DIR, "06_alpha_PM10_D1_D10.csv"),
  na = ""
)

write_csv(
  alpha_pm10_x_sexo,
  file.path(TAB_DIR, "07_alpha_PM10_interaccion_sexo.csv"),
  na = ""
)

write_csv(
  alpha_pm10_sexo_estrat,
  file.path(TAB_DIR, "08_alpha_PM10_estratificado_sexo.csv"),
  na = ""
)

write_csv(
  predicciones_sexo,
  file.path(TAB_DIR, "09_estimaciones_ajustadas_intrusion_sexo.csv"),
  na = ""
)

write_csv(
  alpha_pm10_180d,
  file.path(TAB_DIR, "10_alpha_PM10_180d.csv"),
  na = ""
)

write_csv(
  beta_exposiciones,
  file.path(TAB_DIR, "11_beta_exposiciones_principales.csv"),
  na = ""
)

write_csv(
  permanova_todas,
  file.path(TAB_DIR, "12_PERMANOVA_todos_terminos.csv"),
  na = ""
)

write_csv(
  permanova_bray_miteco,
  file.path(TAB_DIR, "13_PERMANOVA_Bray_covariables_intrusion.csv"),
  na = ""
)

write_csv(
  permanova_bray_pm10,
  file.path(TAB_DIR, "14_PERMANOVA_Bray_covariables_PM10_media.csv"),
  na = ""
)

write_csv(
  betadisper_miteco,
  file.path(TAB_DIR, "15_PERMDISP_Bray_intrusion.csv"),
  na = ""
)

write_csv(
  beta_cluster,
  file.path(TAB_DIR, "16_sensibilidad_Bray_cluster_isla_fecha.csv"),
  na = ""
)

write_csv(
  exac_pm10_180d,
  file.path(TAB_DIR, "17_PM10_180d_y_fenotipo_exacerbador.csv"),
  na = ""
)

write_csv(
  comparacion_metodos_completa,
  file.path(TAB_DIR, "18_ANCOMBC2_DESeq2_comparacion_completa.csv"),
  na = ""
)

write_csv(
  comparacion_metodos_seleccion,
  file.path(TAB_DIR, "19_ANCOMBC2_DESeq2_senales_seleccionadas.csv"),
  na = ""
)

write_csv(
  resumen_metodos,
  file.path(TAB_DIR, "20_ANCOMBC2_DESeq2_resumen.csv"),
  na = ""
)

write_csv(
  top15_genus,
  file.path(TAB_DIR, "21_top15_generos.csv"),
  na = ""
)

write_csv(
  prevalencia_plot,
  file.path(TAB_DIR, "22_prevalencia_y_core_generos.csv"),
  na = ""
)

write_csv(
  resumen_core,
  file.path(TAB_DIR, "23_resumen_microbioma_nucleo.csv"),
  na = ""
)

write_csv(
  pcoa_df,
  file.path(TAB_DIR, "24_coordenadas_PCoA_Bray.csv"),
  na = ""
)

write_csv(
  tabla_indicadores_ambientales,
  file.path(TAB_DIR, "25_tabla_integrada_indicadores_ambientales.csv"),
  na = ""
)

write_csv(
  qc_constancia_cluster,
  file.path(TAB_DIR, "26_QC_constancia_exposicion_cluster.csv"),
  na = ""
)

manifest <- registrar(devolver = TRUE)

write_csv(
  manifest,
  file.path(TAB_DIR, "27_manifest_ejecucion.csv"),
  na = ""
)

# Libro Excel de apoyo para la redacción.
writexl::write_xlsx(
  list(
    Cohorte = tabla_cohorte,
    Faltantes = datos_faltantes,
    Ventanas = resumen_ventanas,
    Concordancia = tabla_concordancia,
    Correlaciones = correlaciones_ambientales,
    Alpha_PM10_10d = alpha_pm10_10d,
    PM10_x_sexo = alpha_pm10_x_sexo,
    PM10_sexo_estrat = alpha_pm10_sexo_estrat,
    Ajustadas_sexo = predicciones_sexo,
    Alpha_PM10_180d = alpha_pm10_180d,
    Beta_exposiciones = beta_exposiciones,
    PERMANOVA_Bray_intr = permanova_bray_miteco,
    PERMANOVA_Bray_PM10 = permanova_bray_pm10,
    PERMDISP = betadisper_miteco,
    Beta_cluster = beta_cluster,
    PM10_180_exac = exac_pm10_180d,
    Taxa_seleccionados = comparacion_metodos_seleccion,
    Metodos_resumen = resumen_metodos,
    Top15_generos = top15_genus,
    Core_resumen = resumen_core,
    Indicadores = tabla_indicadores_ambientales,
    Manifest = manifest
  ),
  file.path(TAB_DIR, "28_RESUMEN_CIERRE_RESULTADOS_TFM.xlsx")
)

saveRDS(
  list(
    meta = meta,
    env_cluster = env_cluster,
    distancias = distancias,
    dist_cluster_bray = dist_cluster_bray,
    alpha_pm10_10d = alpha_pm10_10d,
    alpha_pm10_x_sexo = alpha_pm10_x_sexo,
    alpha_pm10_180d = alpha_pm10_180d,
    permanova_todas = permanova_todas,
    comparacion_metodos = comparacion_metodos_completa,
    predicciones_sexo = predicciones_sexo
  ),
  file.path(OBJ_DIR, "objetos_cierre_resultados_TFM.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG_DIR, "SESSION_INFO.txt")
)


# 20. INFORME DE CONTROL Y FIN


rho_txt <- correlaciones_ambientales %>%
  filter(
    estrato == "Total",
    y == "pm10_media_10d"
  ) %>%
  slice(1)

pm10_beta_bray <- beta_exposiciones %>%
  filter(
    distancia == "Bray_Curtis",
    exposure == "pm10_media_10d_z"
  ) %>%
  slice(1)

sink(file.path(LOG_DIR, "00_RESUMEN_QC_CIERRE_RESULTADOS.txt"))

cat("CIERRE DE RESULTADOS TFM GEMAS\n")
cat("================================\n\n")
cat("Fecha de ejecución:", format(Sys.time()), "\n")
cat("Proyecto:", PROJECT_DIR, "\n\n")

cat("COHORTE\n")
cat("-------\n")
cat("Pacientes:", nrow(meta), "\n")
cat("La Palma:", sum(meta$island == "La Palma"), "\n")
cat("Tenerife:", sum(meta$island == "Tenerife"), "\n")
cat("Mujeres:", sum(meta$sex == "Female"), "\n")
cat("Hombres:", sum(meta$sex == "Male"), "\n")
cat("Menores de 18 años:", sum(meta$menor_18), "\n")
cat("Exacerbadores:", sum(meta$exac6m_num == 1), "\n")
cat("No exacerbadores:", sum(meta$exac6m_num == 0), "\n\n")

cat("EXPOSICIÓN OFICIAL\n")
cat("------------------\n")
for (i in seq_len(nrow(resumen_ventanas))) {
  cat(
    as.character(resumen_ventanas$ventana[i]), ": ",
    resumen_ventanas$n_expuestos[i], "/",
    resumen_ventanas$n_total[i], " (",
    round(resumen_ventanas$porcentaje_expuestos[i], 1), "%)\n",
    sep = ""
  )
}
cat("\n")

cat("RELACIÓN INTRUSIÓN-PM10\n")
cat("-----------------------\n")
cat("Clusters isla-fecha:", nrow(env_cluster), "\n")
cat(
  "Spearman días de intrusión vs PM10 media: rho = ",
  round(rho_txt$rho_spearman, 3),
  "; p = ",
  format.pval(rho_txt$p, digits = 3, eps = 0.001),
  "\n\n",
  sep = ""
)

cat("REANÁLISIS PM10\n")
cat("----------------\n")
cat("Modelos alfa D-1 a D-10:", nrow(alpha_pm10_10d), "\n")
cat("Interacciones PM10 x sexo:", nrow(alpha_pm10_x_sexo), "\n")
cat("Distancias beta incluidas:", paste(names(distancias), collapse = ", "), "\n")
if (nrow(pm10_beta_bray) == 1) {
  cat(
    "PM10 media -> Bray-Curtis: R2 = ",
    round(pm10_beta_bray$R2, 4),
    "; p = ",
    format.pval(pm10_beta_bray$p, digits = 3, eps = 0.001),
    "; q = ",
    format.pval(pm10_beta_bray$q_BH_distancias, digits = 3, eps = 0.001),
    "\n",
    sep = ""
  )
}
cat("\n")

cat("COMPARACIÓN TAXONÓMICA\n")
cat("----------------------\n")
print(resumen_metodos)
cat("\n")

cat("MICROBIOMA NÚCLEO\n")
cat("-----------------\n")
cat(
  "Géneros con prevalencia >=90%: ",
  resumen_core$n_generos[resumen_core$umbral_prevalencia == "≥90 %"],
  "\n",
  sep = ""
)
cat(
  "Géneros con prevalencia >=95%: ",
  resumen_core$n_generos[resumen_core$umbral_prevalencia == "≥95 %"],
  "\n\n",
  sep = ""
)

cat("ARCHIVOS PRINCIPALES PARA REDACTAR\n")
cat("----------------------------------\n")
cat("tablas/01_caracteristicas_cohorte_segun_intrusion.csv\n")
cat("tablas/02_frecuencia_intrusion_por_ventanas.csv\n")
cat("tablas/03_concordancia_intrusion_oficial_PM10.csv\n")
cat("tablas/06_alpha_PM10_D1_D10.csv\n")
cat("tablas/11_beta_exposiciones_principales.csv\n")
cat("tablas/13_PERMANOVA_Bray_covariables_intrusion.csv\n")
cat("tablas/19_ANCOMBC2_DESeq2_senales_seleccionadas.csv\n")
cat("tablas/28_RESUMEN_CIERRE_RESULTADOS_TFM.xlsx\n")
cat("figuras/Figura_2_exposicion_oficial_y_PM10_superficial.png\n")
cat("figuras/Figura_3_composicion_taxonomica_y_microbioma_nucleo.png\n")
cat("figuras/Figura_4_diversidad_ajustada_intrusion_por_sexo.png\n")
cat("figuras/Figura_candidata_PCoA_Bray_intrusion_oficial.png\n\n")

cat("NOTA\n")
cat("----\n")
cat(
  "Los análisis de este script son de cierre y sensibilidad. ",
  "No deben interpretarse como evidencia causal.\n"
)

sink()

msg("============================================================")
msg("SCRIPT DE CIERRE FINALIZADO")
msg("Resultados guardados en:")
msg("%s", OUT_DIR)
msg("")
msg("Pásame primero estos archivos:")
msg("  logs/00_RESUMEN_QC_CIERRE_RESULTADOS.txt")
msg("  tablas/28_RESUMEN_CIERRE_RESULTADOS_TFM.xlsx")
msg("  tablas/19_ANCOMBC2_DESeq2_senales_seleccionadas.csv")
msg("  figuras/Figura_2_exposicion_oficial_y_PM10_superficial.png")
msg("============================================================")
