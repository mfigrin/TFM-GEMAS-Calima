# REANALISIS_COMPLETO_MITECO_GEMAS.R
# TFM GEMAS

# Este script reemplaza, para la inferencia ambiental principal,
# la definicion operativa previa PM10 >=50 ug/m3 por la presencia
# oficial de episodios africanos MITECO.

# EXPOSICIONES AMBIENTALES
# ------------------------
# PRINCIPAL:
#   miteco10 = >=1 dia oficial MITECO en D-1..D-10.
#
# SECUNDARIA DE CARGA:
#   dias MITECO D-1..D-10, por 1 DE.
#
# SENSIBILIDAD DE INTENSIDAD:
#   aporte africano MITECO ("Descuento") medio D-1..D-10, por 1 DE.
#
# ANALISIS 
#  1. QC de cohorte y exposicion.
#  2. Diversidad alfa: Observed, Shannon, Simpson.
#  3. Distancia Aitchison al centro de isla.
#  4. Diversidad beta: Bray-Curtis.
#  5. Weighted y unweighted UniFrac si existen los objetos 08B.
#  6. Betadisper.
#  7. Sensibilidad beta agregada por cluster isla-fecha.
#  8. Modificacion por sexo, isla, exacerbacion y regimen estacional.
#  9. Ventanas temporales D1-10 / D11-30 / D31-90.
# 10. Dosis-respuesta no lineal para dias MITECO y aporte NAF.
# 11. Core microbiome 90% y 95%.
# 12. Composicion taxonomica descriptiva a genero.
# 13. ANCOM-BC2 a genero:
#       - MITECO10
#       - dias MITECO10
#       - aporte NAF10
#       - MITECO10 x sexo
#       - MITECO10 x exacerbacion
# 14. Sensibilidad ANCOM-BC2 Family / Genus / ASV para:
#       - MITECO10
#       - MITECO10 x sexo
# 15. Fenotipo exacerbador 6M:
#       - asociacion microbiota-exacerbacion actualizada con carga
#         MITECO180 y calendario;
#       - carga MITECO180 -> exacerbacion (con/sin antibiotico).
# 16. Asociacion de los estados comunitarios ya definidos (si existen)
#     con MITECO. El clustering NO se recalcula porque es independiente
#     de la exposicion y su silhouette previo fue bajo.
# Deben haberse ejecutado:
#   07_preparacion_phyloseq_analitico.R
#   30_INTEGRAR_MITECO_PM10_GEMAS_v3.R

# 0. CONFIGURACION

set.seed(20260826)

N_PERM <- 999
ALPHA_FDR <- 0.05
PREV_ANCOM <- 0.10
PSEUDOCOUNT_CLR <- 0.5

# Si ANCOM-BC2 a ASV resulta demasiado lento en el primer ensayo,
# se puede poner FALSE para obtener el borrador preliminar y
# reactivarlo antes de la entrega final.
EJECUTAR_ANCOM_ASV <- TRUE

# Si no existen los objetos UniFrac creados previamente, el script
# continuara sin detenerse y lo registrara en el manifest.
EJECUTAR_UNIFRAC_SI_EXISTE <- TRUE


# 1. PAQUETES

paquetes <- c(
  "phyloseq",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "ggplot2",
  "vegan",
  "ape",
  "writexl",
  "sandwich",
  "lmtest",
  "ANCOMBC",
  "splines"
)

faltan <- paquetes[
  !vapply(
    paquetes,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(faltan) > 0) {
  stop(
    "Faltan paquetes: ",
    paste(faltan, collapse = ", "),
    "\nInstalalos antes de ejecutar el reanalisis."
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
  library(ape)
  library(writexl)
  library(sandwich)
  library(lmtest)
  library(ANCOMBC)
})


# 2. RUTAS

localizar_archivo_unico <- function(
  nombre,
  raiz = "outputs",
  obligatorio = TRUE
) {
  encontrados <- list.files(
    raiz,
    pattern = paste0(
      "^",
      gsub(
        "([.(){}+*?^$|\\\\])",
        "\\\\\\1",
        nombre
      ),
      "$"
    ),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )

  if (length(encontrados) == 1) {
    cat(
      "Encontrado: ",
      nombre,
      " -> ",
      encontrados,
      "\n",
      sep = ""
    )

    return(encontrados)
  }

  if (length(encontrados) > 1) {
    stop(
      "Se encontraron varias copias de '",
      nombre,
      "'. No puedo elegir de forma segura:\n- ",
      paste(encontrados, collapse = "\n- ")
    )
  }

  if (obligatorio) {
    stop(
      "No se encontro '",
      nombre,
      "' dentro de ",
      raiz,
      "/.\n",
      "Ejecuta primero el script que genera ese objeto o comprueba su nombre real."
    )
  }

  NA_character_
}


ruta_ps_amb_f <- localizar_archivo_unico(
  "ps_ambiental_205_filtrado_prev2_total10.rds"
)

ruta_ps_amb_u <- localizar_archivo_unico(
  "ps_ambiental_205_sin_filtrado.rds"
)

ruta_ps_saliva_u <- localizar_archivo_unico(
  "ps_saliva_206_sin_filtrado.rds"
)

ruta_miteco_muestra <- localizar_archivo_unico(
  "04_exposicion_MITECO_PM10_por_muestra.csv"
)

ruta_miteco_traza <- localizar_archivo_unico(
  "02_trazabilidad_MITECO_PM10_D1_D180.csv"
)

ruta_wunifrac <- localizar_archivo_unico(
  "dist_weighted_unifrac_ambiental_205.rds",
  obligatorio = FALSE
)

ruta_uunifrac <- localizar_archivo_unico(
  "dist_unweighted_unifrac_ambiental_205_rarefactado_10000.rds",
  obligatorio = FALSE
)

ruta_states <- localizar_archivo_unico(
  "02_asignacion_estados_comunitarios.csv",
  obligatorio = FALSE
)

ruta_salida <- file.path(
  "outputs",
  "analisis_microbiota",
  "reanalisis_completo_MITECO"
)

ruta_tablas <- file.path(
  ruta_salida,
  "tablas"
)

ruta_figuras <- file.path(
  ruta_salida,
  "figuras"
)

ruta_objetos <- file.path(
  ruta_salida,
  "objetos"
)

for (r in c(
  ruta_salida,
  ruta_tablas,
  ruta_figuras,
  ruta_objetos
)) {
  dir.create(
    r,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# 3. FUNCIONES

a_num <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  suppressWarnings(
    as.numeric(
      gsub(
        ",",
        ".",
        trimws(as.character(x)),
        fixed = TRUE
      )
    )
  )
}


a_binario <- function(x) {
  if (is.logical(x)) {
    return(as.numeric(x))
  }

  y <- tolower(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA_real_,
    length(y)
  )

  out[
    y %in% c(
      "1",
      "true",
      "t",
      "si",
      "sí",
      "yes",
      "y",
      "expuesto",
      "exposed",
      "exacerbador",
      "male",
      "m",
      "masculino",
      "hombre"
    )
  ] <- 1

  out[
    y %in% c(
      "0",
      "false",
      "f",
      "no",
      "n",
      "no_expuesto",
      "non-exposed",
      "nonexposed",
      "no_exacerbador",
      "female",
      "f",
      "femenino",
      "mujer"
    )
  ] <- 0

  suppressWarnings({
    z <- as.numeric(y)
  })

  idx <- is.na(out) &
    !is.na(z) &
    z %in% c(0, 1)

  out[idx] <- z[idx]

  out
}


limpiar_texto <- function(x) {
  y <- trimws(
    as.character(x)
  )

  y[
    y %in% c(
      "",
      "NA",
      "NaN",
      "NULL",
      "null"
    )
  ] <- NA_character_

  y
}


normalizar_isla <- function(x) {
  y <- limpiar_texto(x)

  y[
    y %in% c(
      "La_Palma",
      "La palma",
      "LA PALMA"
    )
  ] <- "La Palma"

  y[
    y %in% c(
      "TENERIFE",
      "tenerife"
    )
  ] <- "Tenerife"

  y
}


zscore <- function(x) {
  s <- sd(
    x,
    na.rm = TRUE
  )

  if (
    is.na(s) ||
      s == 0
  ) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) /
    s
}


sample_data_df <- function(ps_obj) {
  sd_obj <- sample_data(ps_obj)

  out <- data.frame(
    setNames(
      sd_obj@.Data,
      names(sd_obj)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(out) <- sample_names(ps_obj)
  out$id_muestra <- rownames(out)

  out
}


otu_taxa_rows <- function(ps_obj) {
  m <- as(
    otu_table(ps_obj),
    "matrix"
  )

  if (!taxa_are_rows(ps_obj)) {
    m <- t(m)
  }

  m
}


tax_df <- function(ps_obj) {
  tt <- as.data.frame(
    tax_table(ps_obj),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (v in names(tt)) {
    tt[[v]] <- limpiar_texto(
      tt[[v]]
    )
  }

  tt$taxon_id <- rownames(tt)
  tt
}


n_niveles <- function(x) {
  if (is.factor(x)) {
    return(
      nlevels(
        droplevels(x)
      )
    )
  }

  length(
    unique(
      x[
        !is.na(x)
      ]
    )
  )
}


covars_activas <- function(
  d,
  vars
) {
  vars[
    vapply(
      vars,
      function(v) {
        v %in% names(d) &&
          n_niveles(
            d[[v]]
          ) >= 2
      },
      logical(1)
    )
  ]
}


subdist <- function(
  dist_obj,
  ids
) {
  m <- as.matrix(
    dist_obj
  )

  ids <- intersect(
    ids,
    rownames(m)
  )

  as.dist(
    m[
      ids,
      ids,
      drop = FALSE
    ]
  )
}


vcov_cluster <- function(
  fit,
  cluster
) {
  sandwich::vcovCL(
    fit,
    cluster = cluster,
    type = "HC1"
  )
}


extraer_coef_robusto <- function(
  fit,
  cluster,
  termino
) {
  ct <- lmtest::coeftest(
    fit,
    vcov. = vcov_cluster(
      fit,
      cluster
    )
  )

  if (!termino %in% rownames(ct)) {
    return(NULL)
  }

  beta <- ct[
    termino,
    "Estimate"
  ]

  se <- ct[
    termino,
    "Std. Error"
  ]

  tibble(
    termino = termino,
    beta = beta,
    se = se,
    IC95_inf = beta -
      1.96 *
      se,
    IC95_sup = beta +
      1.96 *
      se,
    p = ct[
      termino,
      ncol(ct)
    ]
  )
}


wald_terms <- function(
  fit,
  cluster,
  terms
) {
  b <- coef(
    fit
  )

  V <- vcov_cluster(
    fit,
    cluster
  )

  terms <- intersect(
    terms,
    names(b)
  )

  terms <- terms[
    !is.na(
      b[
        terms
      ]
    )
  ]

  if (length(terms) == 0) {
    return(
      tibble(
        chi2 = NA_real_,
        df = 0L,
        p = NA_real_
      )
    )
  }

  bb <- b[
    terms
  ]

  VV <- V[
    terms,
    terms,
    drop = FALSE
  ]

  inv <- tryCatch(
    solve(
      VV
    ),
    error = function(e) {
      qr.solve(
        VV
      )
    }
  )

  W <- as.numeric(
    t(bb) %*%
      inv %*%
      bb
  )

  tibble(
    chi2 = W,
    df = length(bb),
    p = pchisq(
      W,
      df = length(bb),
      lower.tail = FALSE
    )
  )
}


buscar_terminos_interaccion <- function(
  fit,
  token1,
  token2
) {
  nm <- names(
    coef(fit)
  )

  nm[
    grepl(
      token1,
      nm,
      fixed = TRUE
    ) &
      grepl(
        token2,
        nm,
        fixed = TRUE
      ) &
      grepl(
        ":",
        nm,
        fixed = TRUE
      )
  ]
}


manifest <- tibble(
  bloque = character(),
  estado = character(),
  detalle = character()
)


registrar <- function(
  bloque,
  estado,
  detalle = ""
) {
  manifest <<- bind_rows(
    manifest,
    tibble(
      bloque = bloque,
      estado = estado,
      detalle = detalle
    )
  )
}


# 4. CARGA E INTEGRACION DE DATOS

ps_amb_f <- readRDS(
  ruta_ps_amb_f
)

ps_amb_u <- readRDS(
  ruta_ps_amb_u
)

ps_saliva_u <- readRDS(
  ruta_ps_saliva_u
)

miteco <- readr::read_csv(
  ruta_miteco_muestra,
  show_col_types = FALSE
) %>%
  mutate(
    fecha_muestra = as.Date(
      fecha_muestra
    )
  )

traza <- readr::read_csv(
  ruta_miteco_traza,
  show_col_types = FALSE
) %>%
  mutate(
    fecha = as.Date(
      fecha
    ),
    fecha_muestra = as.Date(
      fecha_muestra
    ),
    miteco_africano = as.integer(
      tidyr::replace_na(
        miteco_africano,
        0
      )
    )
  )


# 5. VENTANAS TEMPORALES OFICIALES MITECO


ventanas <- traza %>%
  group_by(
    id_muestra
  ) %>%
  summarise(
    miteco_days_D1_3 = sum(
      miteco_africano[
        between(
          lag,
          1,
          3
        )
      ],
      na.rm = TRUE
    ),

    miteco_days_D4_10 = sum(
      miteco_africano[
        between(
          lag,
          4,
          10
        )
      ],
      na.rm = TRUE
    ),

    miteco_days_D1_10_check = sum(
      miteco_africano[
        between(
          lag,
          1,
          10
        )
      ],
      na.rm = TRUE
    ),

    miteco_days_D11_30 = sum(
      miteco_africano[
        between(
          lag,
          11,
          30
        )
      ],
      na.rm = TRUE
    ),

    miteco_days_D31_90 = sum(
      miteco_africano[
        between(
          lag,
          31,
          90
        )
      ],
      na.rm = TRUE
    ),

    miteco_days_D1_30 = sum(
      miteco_africano[
        between(
          lag,
          1,
          30
        )
      ],
      na.rm = TRUE
    ),

    miteco_days_D1_90 = sum(
      miteco_africano[
        between(
          lag,
          1,
          90
        )
      ],
      na.rm = TRUE
    ),

    miteco_recency_D1_10 = ifelse(
      any(
        miteco_africano[
          between(
            lag,
            1,
            10
          )
        ] == 1,
        na.rm = TRUE
      ),
      min(
        lag[
          between(
            lag,
            1,
            10
          ) &
            miteco_africano == 1
        ],
        na.rm = TRUE
      ),
      NA_real_
    ),

    .groups = "drop"
  )


miteco <- miteco %>%
  left_join(
    ventanas,
    by = "id_muestra"
  ) %>%
  mutate(
    # Auditoria: la suma diaria debe reproducir script 30.
    qc_dias10_coincide = miteco_days_D1_10_check ==
      miteco_days_10d,

    miteco10 = as.numeric(
      miteco_any_10d
    ),

    miteco10_factor = factor(
      ifelse(
        miteco10 == 1,
        "MITECO_si",
        "MITECO_no"
      ),
      levels = c(
        "MITECO_no",
        "MITECO_si"
      )
    ),

    miteco_days10_z = zscore(
      miteco_days_10d
    ),

    naf10_z = zscore(
      naf_mean_calendar_10d
    ),

    miteco_days180_z = zscore(
      miteco_days_180d
    ),

    recent_D1_10_z = zscore(
      miteco_days_10d
    ),

    residual_D11_30_z = zscore(
      miteco_days_D11_30
    ),

    prior_D31_90_z = zscore(
      miteco_days_D31_90
    ),

    grado_miteco10 = factor(
      case_when(
        miteco_days_10d == 0 ~
          "0_sin_MITECO",

        miteco_days_10d >= 1 &
          miteco_days_10d <= 3 ~
          "1_3_dias",

        miteco_days_10d >= 4 ~
          "4plus_dias",

        TRUE ~ NA_character_
      ),
      levels = c(
        "0_sin_MITECO",
        "1_3_dias",
        "4plus_dias"
      )
    ),

    patron_miteco10 = factor(
      case_when(
        miteco_days_D1_3 == 0 &
          miteco_days_D4_10 == 0 ~
          "Sin_MITECO_10d",

        miteco_days_D1_3 > 0 &
          miteco_days_D4_10 == 0 ~
          "Solo_D1_3",

        miteco_days_D1_3 == 0 &
          miteco_days_D4_10 > 0 ~
          "Solo_D4_10",

        miteco_days_D1_3 > 0 &
          miteco_days_D4_10 > 0 ~
          "D1_3_y_D4_10",

        TRUE ~ NA_character_
      ),
      levels = c(
        "Sin_MITECO_10d",
        "Solo_D1_3",
        "Solo_D4_10",
        "D1_3_y_D4_10"
      )
    )
  )


if (!all(
  miteco$qc_dias10_coincide,
  na.rm = TRUE
)) {
  stop(
    "La reconstruccion MITECO D1-D10 no reproduce miteco_days_10d del script 30."
  )
}

registrar(
  "Integracion MITECO",
  "OK",
  paste0(
    "n=",
    nrow(miteco),
    "; MITECO10+=",
    sum(
      miteco$miteco10 == 1,
      na.rm = TRUE
    )
  )
)


# 6. METADATA ANALITICA

meta_base <- sample_data_df(
  ps_amb_f
)

if (
  "fecha_muestra" %in% names(
    meta_base
  )
) {
  names(
    meta_base
  )[
    names(
      meta_base
    ) == "fecha_muestra"
  ] <- "fecha_muestra_phyloseq"
}

meta_base <- meta_base %>%
  mutate(
    age = a_num(
      Age
    ),

    sex_num = a_binario(
      Sex.letter
    ),

    sex = factor(
      ifelse(
        sex_num == 1,
        "Male",
        ifelse(
          sex_num == 0,
          "Female",
          NA_character_
        )
      ),
      levels = c(
        "Female",
        "Male"
      )
    ),

    island = factor(
      normalizar_isla(
        Island
      ),
      levels = c(
        "La Palma",
        "Tenerife"
      )
    ),

    pool = factor(
      Pool
    ),

    antibiotic = factor(
      a_binario(
        Antibiotic.2mo
      ),
      levels = c(
        0,
        1
      )
    ),

    cavities = factor(
      a_binario(
        Cavities
      ),
      levels = c(
        0,
        1
      )
    ),

    gum = factor(
      a_binario(
        Chewing.gum.30min
      ),
      levels = c(
        0,
        1
      )
    ),

    exac6m_num = a_num(
      Exacerbations.6M
    ),

    exac6m_num = ifelse(
      exac6m_num %in% c(
        0,
        1
      ),
      exac6m_num,
      NA_real_
    ),

    exac6m = factor(
      exac6m_num,
      levels = c(
        0,
        1
      )
    )
  ) %>%
  left_join(
    miteco,
    by = "id_muestra"
  )

# QC DE LA UNION DE FECHAS

if (
  !"fecha_muestra" %in% names(
    meta_base
  )
) {
  stop(
    "Tras unir MITECO no existe fecha_muestra. ",
    "Columnas disponibles: ",
    paste(
      names(
        meta_base
      ),
      collapse = ", "
    )
  )
}

meta_base$fecha_muestra <- as.Date(
  meta_base$fecha_muestra
)

if (
  "fecha_muestra_phyloseq" %in% names(
    meta_base
  )
) {
  meta_base$fecha_muestra_phyloseq <- as.Date(
    meta_base$fecha_muestra_phyloseq
  )

  n_mismatch_fecha_phyloseq <- sum(
    !is.na(
      meta_base$fecha_muestra
    ) &
      !is.na(
        meta_base$fecha_muestra_phyloseq
      ) &
      meta_base$fecha_muestra !=
        meta_base$fecha_muestra_phyloseq
  )
} else {
  n_mismatch_fecha_phyloseq <- NA_integer_
}

fecha_original_phyloseq <- as.Date(
  meta_base$Date
)

n_mismatch_fecha_original <- sum(
  !is.na(
    meta_base$fecha_muestra
  ) &
    !is.na(
      fecha_original_phyloseq
    ) &
    meta_base$fecha_muestra !=
      fecha_original_phyloseq
)

n_sin_fecha_miteco <- sum(
  is.na(
    meta_base$fecha_muestra
  )
)

qc_integracion_fechas <- tibble(
  n_muestras = nrow(
    meta_base
  ),
  n_sin_fecha_MITECO = n_sin_fecha_miteco,
  discrepancias_fecha_vs_phyloseq_ambiental =
    n_mismatch_fecha_phyloseq,
  discrepancias_fecha_vs_Date_original =
    n_mismatch_fecha_original
)

print(
  qc_integracion_fechas
)

if (
  n_sin_fecha_miteco > 0
) {
  stop(
    "Hay ",
    n_sin_fecha_miteco,
    " muestras sin fecha tras la union MITECO."
  )
}

if (
  n_mismatch_fecha_original > 0
) {
  stop(
    "La fecha MITECO discrepa de Date original del phyloseq en ",
    n_mismatch_fecha_original,
    " muestras. No se continua hasta auditarlo."
  )
}

if (
  !is.na(
    n_mismatch_fecha_phyloseq
  ) &&
    n_mismatch_fecha_phyloseq > 0
) {
  warning(
    "La fecha MITECO discrepa de la antigua fecha_muestra ambiental en ",
    n_mismatch_fecha_phyloseq,
    " muestras, aunque coincide con Date original."
  )
}

meta_base$lecturas_totales <- as.numeric(
  sample_sums(
    ps_amb_f
  )[
    meta_base$id_muestra
  ]
)

meta_base$lecturas_log10 <- log10(
  meta_base$lecturas_totales +
    1
)

meta_base <- meta_base %>%
  mutate(
    fecha = as.Date(
      fecha_muestra
    ),

    year = factor(
      format(
        fecha,
        "%Y"
      )
    ),

    doy = as.integer(
      format(
        fecha,
        "%j"
      )
    ),

    sin_doy = sin(
      2 *
        pi *
        doy /
        365.25
    ),

    cos_doy = cos(
      2 *
        pi *
        doy /
        365.25
    ),

    cluster_isla_fecha = factor(
      paste(
        as.character(
          island
        ),
        fecha,
        sep = "__"
      )
    ),

    regimen_strict = factor(
      case_when(
        format(
          fecha,
          "%m"
        ) %in% c(
          "11",
          "12",
          "01"
        ) ~
          "Winter_proxy_NDJ",

        format(
          fecha,
          "%m"
        ) %in% c(
          "07",
          "08"
        ) ~
          "Summer_proxy_JA",

        TRUE ~ NA_character_
      ),
      levels = c(
        "Winter_proxy_NDJ",
        "Summer_proxy_JA"
      )
    ),

    miteco_sex_int = miteco10 *
      sex_num,

    miteco_exac_int = miteco10 *
      exac6m_num
  )


# 7. ACTUALIZAR PHYLOSEQ CON METADATA MITECO

actualizar_ps <- function(
  ps_obj,
  meta
) {
  ids <- intersect(
    sample_names(
      ps_obj
    ),
    meta$id_muestra
  )

  p <- prune_samples(
    ids,
    ps_obj
  )

  m <- meta %>%
    filter(
      id_muestra %in% ids
    ) %>%
    as.data.frame(
      stringsAsFactors = FALSE
    )

  rownames(m) <- m$id_muestra

  m <- m[
    sample_names(
      p
    ),
    ,
    drop = FALSE
  ]

  sample_data(
    p
  ) <- sample_data(
    m
  )

  p
}


ps_amb_f_m <- actualizar_ps(
  ps_amb_f,
  meta_base
)

ps_amb_u_m <- actualizar_ps(
  ps_amb_u,
  meta_base
)

saveRDS(
  ps_amb_f_m,
  file.path(
    ruta_objetos,
    "ps_ambiental_205_filtrado_con_MITECO.rds"
  )
)

saveRDS(
  ps_amb_u_m,
  file.path(
    ruta_objetos,
    "ps_ambiental_205_sin_filtrar_con_MITECO.rds"
  )
)

# 8. QC DESCRIPTIVO EXPOSICION

qc_exposicion <- meta_base %>%
  summarise(
    n = n(),

    miteco10_si = sum(
      miteco10 == 1,
      na.rm = TRUE
    ),

    miteco10_no = sum(
      miteco10 == 0,
      na.rm = TRUE
    ),

    miteco30_si = sum(
      miteco_days_D1_30 > 0,
      na.rm = TRUE
    ),

    miteco30_no = sum(
      miteco_days_D1_30 == 0,
      na.rm = TRUE
    ),

    miteco90_si = sum(
      miteco_days_D1_90 > 0,
      na.rm = TRUE
    ),

    miteco90_no = sum(
      miteco_days_D1_90 == 0,
      na.rm = TRUE
    ),

    mediana_dias_MITECO_D1_10 = median(
      miteco_days_10d,
      na.rm = TRUE
    ),

    mediana_dias_MITECO_D11_30 = median(
      miteco_days_D11_30,
      na.rm = TRUE
    ),

    mediana_dias_MITECO_D31_90 = median(
      miteco_days_D31_90,
      na.rm = TRUE
    )
  )


qc_grado <- meta_base %>%
  count(
    grado_miteco10,
    name = "n"
  )


qc_patron <- meta_base %>%
  count(
    patron_miteco10,
    name = "n"
  )


qc_sexo <- meta_base %>%
  count(
    sex,
    miteco10_factor,
    name = "n"
  )


qc_isla <- meta_base %>%
  count(
    island,
    miteco10_factor,
    name = "n"
  )


qc_regimen <- meta_base %>%
  filter(
    !is.na(
      regimen_strict
    )
  ) %>%
  count(
    regimen_strict,
    miteco10_factor,
    name = "n"
  )


write_csv(
  qc_integracion_fechas,
  file.path(
    ruta_tablas,
    "00_QC_integracion_fechas_MITECO.csv"
  ),
  na = ""
)

write_csv(
  qc_exposicion,
  file.path(
    ruta_tablas,
    "00_QC_exposicion_MITECO.csv"
  ),
  na = ""
)

write_csv(
  qc_grado,
  file.path(
    ruta_tablas,
    "00b_QC_grado_dias_MITECO10.csv"
  ),
  na = ""
)

write_csv(
  qc_patron,
  file.path(
    ruta_tablas,
    "00c_QC_patron_temporal_MITECO10.csv"
  ),
  na = ""
)

write_csv(
  qc_regimen,
  file.path(
    ruta_tablas,
    "00d_QC_winter_summer_MITECO.csv"
  ),
  na = ""
)


# 9. DIVERSIDAD ALFA + AITCHISON

alpha <- estimate_richness(
  ps_amb_f_m,
  measures = c(
    "Observed",
    "Shannon",
    "Simpson"
  )
)

alpha$id_muestra <- rownames(
  alpha
)

# Aitchison a centro de isla:
# misma definicion usada en el analisis exploratorio previo.
psg_ait <- tax_glom(
  ps_amb_u_m,
  taxrank = "Genus",
  NArm = TRUE
)

otu_g <- otu_taxa_rows(
  psg_ait
)

prev_g <- rowMeans(
  otu_g > 0
)

keep_g <- names(
  prev_g
)[
  prev_g >= 0.10
]

psg_ait <- prune_taxa(
  keep_g,
  psg_ait
)

otu_g <- otu_taxa_rows(
  psg_ait
)

x_clr <- t(
  otu_g +
    PSEUDOCOUNT_CLR
)

clr <- log(
  x_clr
) -
  rowMeans(
    log(
      x_clr
    )
  )

meta_ait <- meta_base[
  match(
    rownames(clr),
    meta_base$id_muestra
  ),
  ,
  drop = FALSE
]

centros_isla <- lapply(
  levels(
    droplevels(
      meta_ait$island
    )
  ),
  function(isl) {
    ids <- meta_ait$id_muestra[
      meta_ait$island == isl
    ]

    colMeans(
      clr[
        ids,
        ,
        drop = FALSE
      ]
    )
  }
)

names(
  centros_isla
) <- levels(
  droplevels(
    meta_ait$island
  )
)

aitchison_isla <- vapply(
  meta_ait$id_muestra,
  function(id) {
    isl <- as.character(
      meta_ait$island[
        meta_ait$id_muestra == id
      ][1]
    )

    sqrt(
      sum(
        (
          clr[
            id,
            ,
            drop = TRUE
          ] -
            centros_isla[[isl]]
        )^2
      )
    )
  },
  numeric(1)
)


metricas <- meta_base %>%
  left_join(
    alpha,
    by = "id_muestra"
  ) %>%
  mutate(
    aitchison_centro_isla = aitchison_isla[
      id_muestra
    ]
  )


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


ajustar_alpha_cont <- function(
  outcome,
  exposure,
  data = metricas,
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
    complete.cases(
      data[
        ,
        vars,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  cov2 <- covars_activas(
    dd,
    covars
  )

  if (
    n_niveles(
      dd[[exposure]]
    ) < 2
  ) {
    return(NULL)
  }

  fit <- lm(
    as.formula(
      paste(
        outcome,
        "~",
        paste(
          c(
            exposure,
            cov2
          ),
          collapse = " + "
        )
      )
    ),
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
      n_clusters = n_distinct(
        dd$cluster_isla_fecha
      ),
      .before = 1
    )
}


exposiciones_alpha <- c(
  "miteco10",
  "miteco_days10_z",
  "naf10_z"
)


alpha_global <- bind_rows(
  lapply(
    exposiciones_alpha,
    function(e) {
      bind_rows(
        lapply(
          alpha_outcomes,
          function(o) {
            ajustar_alpha_cont(
              o,
              e
            )
          }
        )
      )
    }
  )
) %>%
  group_by(
    exposure
  ) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()


# 10. ALFA: GRADO / PATRON TEMPORAL


ajustar_alpha_factor <- function(
  outcome,
  exposure,
  data = metricas,
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
    complete.cases(
      data[
        ,
        vars,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  if (
    n_niveles(
      dd[[exposure]]
    ) < 2
  ) {
    return(NULL)
  }

  cov2 <- covars_activas(
    dd,
    covars
  )

  fit <- lm(
    as.formula(
      paste(
        outcome,
        "~",
        paste(
          c(
            exposure,
            cov2
          ),
          collapse = " + "
        )
      )
    ),
    data = dd
  )

  terms <- names(
    coef(
      fit
    )
  )

  terms <- terms[
    grepl(
      paste0(
        "^",
        exposure
      ),
      terms
    )
  ]

  w <- wald_terms(
    fit,
    dd$cluster_isla_fecha,
    terms
  )

  tibble(
    outcome = outcome,
    exposure = exposure,
    n = nrow(dd),
    n_clusters = n_distinct(
      dd$cluster_isla_fecha
    ),
    chi2 = w$chi2,
    df = w$df,
    p = w$p
  )
}


alpha_factor <- bind_rows(
  lapply(
    c(
      "grado_miteco10",
      "patron_miteco10"
    ),
    function(e) {
      bind_rows(
        lapply(
          alpha_outcomes,
          function(o) {
            ajustar_alpha_factor(
              o,
              e
            )
          }
        )
      )
    }
  )
) %>%
  group_by(
    exposure
  ) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()


# 11. ALFA: INTERACCIONES

ajustar_alpha_interaccion <- function(
  outcome,
  exposure,
  modifier,
  data = metricas,
  covars = covars_alpha
) {
  cov2 <- setdiff(
    covars,
    modifier
  )

  vars <- unique(
    c(
      outcome,
      exposure,
      modifier,
      cov2,
      "cluster_isla_fecha"
    )
  )

  dd <- data[
    complete.cases(
      data[
        ,
        vars,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  if (
    n_niveles(
      dd[[exposure]]
    ) < 2 ||
      n_niveles(
        dd[[modifier]]
      ) < 2
  ) {
    return(NULL)
  }

  cov2 <- covars_activas(
    dd,
    cov2
  )

  fit <- lm(
    as.formula(
      paste(
        outcome,
        "~",
        paste0(
          exposure,
          " * ",
          modifier
        ),
        if (
          length(cov2) > 0
        ) {
          paste0(
            " + ",
            paste(
              cov2,
              collapse = " + "
            )
          )
        } else {
          ""
        }
      )
    ),
    data = dd
  )

  terms <- buscar_terminos_interaccion(
    fit,
    exposure,
    modifier
  )

  w <- wald_terms(
    fit,
    dd$cluster_isla_fecha,
    terms
  )

  tibble(
    outcome = outcome,
    exposure = exposure,
    modifier = modifier,
    n = nrow(dd),
    n_clusters = n_distinct(
      dd$cluster_isla_fecha
    ),
    chi2 = w$chi2,
    df = w$df,
    p_interaccion = w$p
  )
}


alpha_int_sexo <- bind_rows(
  lapply(
    c(
      "miteco10",
      "miteco_days10_z",
      "naf10_z"
    ),
    function(e) {
      bind_rows(
        lapply(
          alpha_outcomes,
          function(o) {
            ajustar_alpha_interaccion(
              o,
              e,
              "sex"
            )
          }
        )
      )
    }
  )
) %>%
  group_by(
    exposure
  ) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p_interaccion,
      method = "BH"
    )
  ) %>%
  ungroup()


alpha_int_isla <- bind_rows(
  lapply(
    alpha_outcomes,
    function(o) {
      ajustar_alpha_interaccion(
        o,
        "miteco10",
        "island"
      )
    }
  )
) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p_interaccion,
      method = "BH"
    )
  )


alpha_int_exac <- bind_rows(
  lapply(
    alpha_outcomes,
    function(o) {
      ajustar_alpha_interaccion(
        o,
        "miteco10",
        "exac6m"
      )
    }
  )
) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p_interaccion,
      method = "BH"
    )
  )


metricas_regimen <- metricas %>%
  filter(
    !is.na(
      regimen_strict
    )
  ) %>%
  droplevels()


alpha_int_regimen <- bind_rows(
  lapply(
    alpha_outcomes,
    function(o) {
      ajustar_alpha_interaccion(
        o,
        "miteco10",
        "regimen_strict",
        data = metricas_regimen
      )
    }
  )
) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p_interaccion,
      method = "BH"
    )
  )


# 12. ALFA: ESTRATIFICADO POR SEXO

alpha_sexo_estrat <- bind_rows(
  lapply(
    c(
      "Female",
      "Male"
    ),
    function(s) {
      d <- metricas %>%
        filter(
          sex == s
        ) %>%
        droplevels()

      cov_s <- setdiff(
        covars_alpha,
        "sex"
      )

      bind_rows(
        lapply(
          c(
            "miteco10",
            "miteco_days10_z"
          ),
          function(e) {
            bind_rows(
              lapply(
                alpha_outcomes,
                function(o) {
                  z <- ajustar_alpha_cont(
                    o,
                    e,
                    data = d,
                    covars = cov_s
                  )

                  if (is.null(z)) {
                    return(NULL)
                  }

                  z %>%
                    mutate(
                      sex = s,
                      .before = 1
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


# 13. VENTANAS TEMPORALES SIMULTANEAS

windows <- c(
  "recent_D1_10_z",
  "residual_D11_30_z",
  "prior_D31_90_z"
)


ajustar_alpha_tres_ventanas <- function(
  outcome
) {
  vars <- unique(
    c(
      outcome,
      windows,
      covars_alpha,
      "cluster_isla_fecha"
    )
  )

  dd <- metricas[
    complete.cases(
      metricas[
        ,
        vars,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  cov2 <- covars_activas(
    dd,
    covars_alpha
  )

  fit <- lm(
    as.formula(
      paste(
        outcome,
        "~",
        paste(
          c(
            windows,
            cov2
          ),
          collapse = " + "
        )
      )
    ),
    data = dd
  )

  bind_rows(
    lapply(
      windows,
      function(w) {
        z <- extraer_coef_robusto(
          fit,
          dd$cluster_isla_fecha,
          w
        )

        if (is.null(z)) {
          return(NULL)
        }

        z %>%
          mutate(
            outcome = outcome,
            window = w,
            n = nrow(dd),
            n_clusters = n_distinct(
              dd$cluster_isla_fecha
            ),
            .before = 1
          )
      }
    )
  )
}


alpha_ventanas <- bind_rows(
  lapply(
    alpha_outcomes,
    ajustar_alpha_tres_ventanas
  )
) %>%
  group_by(
    window
  ) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()


ajustar_alpha_ventanas_sexo <- function(
  outcome
) {
  cov2 <- setdiff(
    covars_alpha,
    "sex"
  )

  vars <- unique(
    c(
      outcome,
      windows,
      "sex",
      cov2,
      "cluster_isla_fecha"
    )
  )

  dd <- metricas[
    complete.cases(
      metricas[
        ,
        vars,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  cov2 <- covars_activas(
    dd,
    cov2
  )

  fit <- lm(
    as.formula(
      paste(
        outcome,
        "~ (",
        paste(
          windows,
          collapse = " + "
        ),
        ") * sex",
        if (
          length(cov2) > 0
        ) {
          paste0(
            " + ",
            paste(
              cov2,
              collapse = " + "
            )
          )
        } else {
          ""
        }
      )
    ),
    data = dd
  )

  bind_rows(
    lapply(
      windows,
      function(w) {
        terms <- buscar_terminos_interaccion(
          fit,
          w,
          "sex"
        )

        wt <- wald_terms(
          fit,
          dd$cluster_isla_fecha,
          terms
        )

        tibble(
          outcome = outcome,
          window = w,
          n = nrow(dd),
          n_clusters = n_distinct(
            dd$cluster_isla_fecha
          ),
          chi2 = wt$chi2,
          df = wt$df,
          p_interaccion = wt$p
        )
      }
    )
  )
}


alpha_ventanas_sexo <- bind_rows(
  lapply(
    alpha_outcomes,
    ajustar_alpha_ventanas_sexo
  )
) %>%
  group_by(
    window
  ) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p_interaccion,
      method = "BH"
    )
  ) %>%
  ungroup()


# 14. DOSIS-RESPUESTA NO LINEAL

analizar_spline <- function(
  outcome,
  exposure_raw,
  data = metricas
) {
  vars <- unique(
    c(
      outcome,
      exposure_raw,
      covars_alpha
    )
  )

  dd <- data[
    complete.cases(
      data[
        ,
        vars,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  cov2 <- covars_activas(
    dd,
    covars_alpha
  )

  rhs_cov <- if (
    length(cov2) > 0
  ) {
    paste(
      cov2,
      collapse = " + "
    )
  } else {
    "1"
  }

  m0 <- lm(
    as.formula(
      paste(
        outcome,
        "~",
        rhs_cov
      )
    ),
    data = dd
  )

  m1 <- lm(
    as.formula(
      paste(
        outcome,
        "~",
        exposure_raw,
        "+",
        rhs_cov
      )
    ),
    data = dd
  )

  m2 <- lm(
    as.formula(
      paste(
        outcome,
        "~ splines::ns(",
        exposure_raw,
        ", df=3) +",
        rhs_cov
      )
    ),
    data = dd
  )

  a01 <- anova(
    m0,
    m1
  )

  a12 <- anova(
    m1,
    m2
  )

  tibble(
    outcome = outcome,
    exposure = exposure_raw,
    n = nrow(dd),
    p_global_lineal = a01$`Pr(>F)`[2],
    p_no_lineal_sobre_lineal = a12$`Pr(>F)`[2],
    AIC_lineal = AIC(m1),
    AIC_spline = AIC(m2)
  )
}


spline_res <- bind_rows(
  lapply(
    c(
      "miteco_days_10d",
      "naf_mean_calendar_10d"
    ),
    function(e) {
      bind_rows(
        lapply(
          c(
            "Observed",
            "Shannon",
            "Simpson"
          ),
          function(o) {
            analizar_spline(
              o,
              e
            )
          }
        )
      )
    }
  )
) %>%
  group_by(
    exposure
  ) %>%
  mutate(
    q_no_lineal = p.adjust(
      p_no_lineal_sobre_lineal,
      method = "BH"
    )
  ) %>%
  ungroup()



# 15. BETA-DIVERSIDAD: DISTANCIAS

ps_rel <- transform_sample_counts(
  ps_amb_f_m,
  function(x) {
    x /
      sum(x)
  }
)

dist_bray <- phyloseq::distance(
  ps_rel,
  method = "bray"
)

distancias <- list(
  Bray_Curtis = dist_bray
)

if (
  EJECUTAR_UNIFRAC_SI_EXISTE &&
    !is.na(ruta_wunifrac) &&
    file.exists(
      ruta_wunifrac
    )
) {
  distancias$Weighted_UniFrac <- readRDS(
    ruta_wunifrac
  )

  registrar(
    "Weighted UniFrac",
    "OK",
    ruta_wunifrac
  )
} else {
  registrar(
    "Weighted UniFrac",
    "OMITIDO",
    "No existe el objeto 08B."
  )
}

if (
  EJECUTAR_UNIFRAC_SI_EXISTE &&
    !is.na(ruta_uunifrac) &&
    file.exists(
      ruta_uunifrac
    )
) {
  distancias$Unweighted_UniFrac_rarefactado <- readRDS(
    ruta_uunifrac
  )

  registrar(
    "Unweighted UniFrac",
    "OK",
    ruta_uunifrac
  )
} else {
  registrar(
    "Unweighted UniFrac",
    "OMITIDO",
    "No existe el objeto 08B."
  )
}


# 16. PERMANOVA GENERAL

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


permanova_cont <- function(
  dist_obj,
  distancia,
  exposure,
  data = meta_base,
  covars = covars_beta
) {
  ids <- intersect(
    labels(
      dist_obj
    ),
    data$id_muestra
  )

  dd <- data[
    match(
      ids,
      data$id_muestra
    ),
    ,
    drop = FALSE
  ]

  rownames(dd) <- dd$id_muestra

  vars <- unique(
    c(
      exposure,
      covars
    )
  )

  cc <- complete.cases(
    dd[
      ,
      vars,
      drop = FALSE
    ]
  )

  dd <- dd[
    cc,
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  d2 <- subdist(
    dist_obj,
    rownames(dd)
  )

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  cov2 <- covars_activas(
    dd,
    covars
  )

  if (
    n_niveles(
      dd[[exposure]]
    ) < 2
  ) {
    return(NULL)
  }

  ftxt <- paste(
    "d2 ~",
    paste(
      c(
        exposure,
        cov2
      ),
      collapse = " + "
    )
  )

  fit <- vegan::adonis2(
    as.formula(
      ftxt
    ),
    data = dd,
    permutations = N_PERM,
    by = "margin"
  )

  tab <- as.data.frame(
    fit
  )

  tab$term <- rownames(
    tab
  )

  fila <- tab[
    tab$term == exposure,
    ,
    drop = FALSE
  ]

  if (
    nrow(
      fila
    ) != 1
  ) {
    return(NULL)
  }

  tibble(
    distancia = distancia,
    exposure = exposure,
    n = nrow(dd),
    R2 = fila$R2,
    F = fila$F,
    p = fila$`Pr(>F)`,
    formula = ftxt
  )
}


beta_global <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      d <- distancias[[nm]]

      bind_rows(
        lapply(
          c(
            "miteco10",
            "miteco_days10_z",
            "naf10_z"
          ),
          function(e) {
            permanova_cont(
              d,
              nm,
              e
            )
          }
        )
      )
    }
  )
) %>%
  group_by(
    exposure
  ) %>%
  mutate(
    p_BH_distancias = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()



# 17. PERMANOVA INTERACCIONES

permanova_interaccion <- function(
  dist_obj,
  distancia,
  exposure,
  modifier,
  data = meta_base,
  covars = covars_beta
) {
  ids <- intersect(
    labels(
      dist_obj
    ),
    data$id_muestra
  )

  dd <- data[
    match(
      ids,
      data$id_muestra
    ),
    ,
    drop = FALSE
  ]

  rownames(dd) <- dd$id_muestra

  cov2 <- setdiff(
    covars,
    modifier
  )

  vars <- unique(
    c(
      exposure,
      modifier,
      cov2
    )
  )

  cc <- complete.cases(
    dd[
      ,
      vars,
      drop = FALSE
    ]
  )

  dd <- droplevels(
    dd[
      cc,
      ,
      drop = FALSE
    ]
  )

  d2 <- subdist(
    dist_obj,
    rownames(dd)
  )

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  cov2 <- covars_activas(
    dd,
    cov2
  )

  if (
    n_niveles(
      dd[[exposure]]
    ) < 2 ||
      n_niveles(
        dd[[modifier]]
      ) < 2
  ) {
    return(NULL)
  }

  ftxt <- paste(
    "d2 ~",
    paste0(
      exposure,
      " * ",
      modifier
    ),
    if (
      length(cov2) > 0
    ) {
      paste0(
        " + ",
        paste(
          cov2,
          collapse = " + "
        )
      )
    } else {
      ""
    }
  )

  fit <- vegan::adonis2(
    as.formula(
      ftxt
    ),
    data = dd,
    permutations = N_PERM,
    by = "margin"
  )

  tab <- as.data.frame(
    fit
  )

  tab$term <- rownames(
    tab
  )

  idx <- grepl(
    exposure,
    tab$term,
    fixed = TRUE
  ) &
    grepl(
      modifier,
      tab$term,
      fixed = TRUE
    ) &
    grepl(
      ":",
      tab$term,
      fixed = TRUE
    )

  fila <- tab[
    idx,
    ,
    drop = FALSE
  ]

  if (
    nrow(
      fila
    ) != 1
  ) {
    return(
      tibble(
        distancia = distancia,
        exposure = exposure,
        modifier = modifier,
        n = nrow(dd),
        R2 = NA_real_,
        F = NA_real_,
        p = NA_real_,
        formula = ftxt
      )
    )
  }

  tibble(
    distancia = distancia,
    exposure = exposure,
    modifier = modifier,
    n = nrow(dd),
    R2 = fila$R2,
    F = fila$F,
    p = fila$`Pr(>F)`,
    formula = ftxt
  )
}


beta_int_sexo <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_interaccion(
        distancias[[nm]],
        nm,
        "miteco10",
        "sex"
      )
    }
  )
) %>%
  mutate(
    p_BH_distancias = p.adjust(
      p,
      method = "BH"
    )
  )


beta_int_isla <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_interaccion(
        distancias[[nm]],
        nm,
        "miteco10",
        "island"
      )
    }
  )
) %>%
  mutate(
    p_BH_distancias = p.adjust(
      p,
      method = "BH"
    )
  )


beta_int_exac <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_interaccion(
        distancias[[nm]],
        nm,
        "miteco10",
        "exac6m"
      )
    }
  )
) %>%
  mutate(
    p_BH_distancias = p.adjust(
      p,
      method = "BH"
    )
  )


beta_int_regimen <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_interaccion(
        distancias[[nm]],
        nm,
        "miteco10",
        "regimen_strict",
        data = meta_base %>%
          filter(
            !is.na(
              regimen_strict
            )
          ) %>%
          droplevels()
      )
    }
  )
) %>%
  mutate(
    p_BH_distancias = p.adjust(
      p,
      method = "BH"
    )
  )



# 18. BETA: VENTANAS TEMPORALES SIMULTANEAS


permanova_ventanas <- function(
  dist_obj,
  distancia
) {
  ids <- labels(
    dist_obj
  )

  dd <- meta_base[
    match(
      ids,
      meta_base$id_muestra
    ),
    ,
    drop = FALSE
  ]

  rownames(dd) <- dd$id_muestra

  vars <- unique(
    c(
      windows,
      covars_beta
    )
  )

  cc <- complete.cases(
    dd[
      ,
      vars,
      drop = FALSE
    ]
  )

  dd <- droplevels(
    dd[
      cc,
      ,
      drop = FALSE
    ]
  )

  d2 <- subdist(
    dist_obj,
    rownames(dd)
  )

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  cov2 <- covars_activas(
    dd,
    covars_beta
  )

  ftxt <- paste(
    "d2 ~",
    paste(
      c(
        windows,
        cov2
      ),
      collapse = " + "
    )
  )

  fit <- vegan::adonis2(
    as.formula(
      ftxt
    ),
    data = dd,
    permutations = N_PERM,
    by = "margin"
  )

  tab <- as.data.frame(
    fit
  )

  tab$term <- rownames(
    tab
  )

  bind_rows(
    lapply(
      windows,
      function(w) {
        fila <- tab[
          tab$term == w,
          ,
          drop = FALSE
        ]

        if (
          nrow(
            fila
          ) != 1
        ) {
          return(NULL)
        }

        tibble(
          distancia = distancia,
          window = w,
          n = nrow(dd),
          R2 = fila$R2,
          F = fila$F,
          p = fila$`Pr(>F)`,
          formula = ftxt
        )
      }
    )
  )
}


beta_ventanas <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_ventanas(
        distancias[[nm]],
        nm
      )
    }
  )
) %>%
  group_by(
    window
  ) %>%
  mutate(
    p_BH_distancias = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()



# 19. BETA: DISPERSIÓN


betadisper_test <- function(
  dist_obj,
  distancia,
  grupo,
  data = meta_base
) {
  ids <- labels(
    dist_obj
  )

  dd <- data[
    match(
      ids,
      data$id_muestra
    ),
    ,
    drop = FALSE
  ]

  rownames(dd) <- dd$id_muestra

  cc <- !is.na(
    dd[[grupo]]
  )

  dd <- droplevels(
    dd[
      cc,
      ,
      drop = FALSE
    ]
  )

  d2 <- subdist(
    dist_obj,
    rownames(dd)
  )

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  g <- factor(
    dd[[grupo]]
  )

  if (
    nlevels(
      droplevels(g)
    ) < 2
  ) {
    return(NULL)
  }

  bd <- vegan::betadisper(
    d2,
    group = g
  )

  perm <- vegan::permutest(
    bd,
    permutations = N_PERM
  )

  tibble(
    distancia = distancia,
    grupo = grupo,
    n = nrow(dd),
    p_betadisper = perm$tab[1, "Pr(>F)"]
  )
}


meta_base$miteco_sex_group <- factor(
  ifelse(
    is.na(
      meta_base$miteco10
    ) |
      is.na(
        meta_base$sex
      ),
    NA_character_,
    paste0(
      as.character(
        meta_base$sex
      ),
      "_",
      ifelse(
        meta_base$miteco10 == 1,
        "MITECO+",
        "MITECO-"
      )
    )
  )
)


betadisper_res <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      bind_rows(
        betadisper_test(
          distancias[[nm]],
          nm,
          "miteco10_factor"
        ),
        betadisper_test(
          distancias[[nm]],
          nm,
          "miteco_sex_group"
        )
      )
    }
  )
) %>%
  group_by(
    grupo
  ) %>%
  mutate(
    p_BH_distancias = p.adjust(
      p_betadisper,
      method = "BH"
    )
  ) %>%
  ungroup()


# 20. BETA: SENSIBILIDAD CLUSTER ISLA-FECHA, BRAY

otu_rel <- otu_taxa_rows(
  ps_rel
)

clusters <- as.character(
  meta_base$cluster_isla_fecha
)

names(
  clusters
) <- meta_base$id_muestra

ids_otu <- colnames(
  otu_rel
)

clusters_otu <- clusters[
  ids_otu
]

cluster_levels <- unique(
  clusters_otu
)

perfil_cluster <- sapply(
  cluster_levels,
  function(cl) {
    ids <- ids_otu[
      clusters_otu == cl
    ]

    if (
      length(ids) == 1
    ) {
      otu_rel[
        ,
        ids,
        drop = TRUE
      ]
    } else {
      rowMeans(
        otu_rel[
          ,
          ids,
          drop = FALSE
        ]
      )
    }
  }
)

colnames(
  perfil_cluster
) <- cluster_levels

dist_cluster_bray <- vegan::vegdist(
  t(
    perfil_cluster
  ),
  method = "bray"
)


qc_cluster_constancia <- meta_base %>%
  group_by(
    cluster_isla_fecha
  ) %>%
  summarise(
    n_miteco10 = n_distinct(
      miteco10[
        !is.na(
          miteco10
        )
      ]
    ),
    n_days10 = n_distinct(
      miteco_days_10d[
        !is.na(
          miteco_days_10d
        )
      ]
    ),
    .groups = "drop"
  )

if (
  any(
    qc_cluster_constancia$n_miteco10 > 1,
    na.rm = TRUE
  ) ||
    any(
      qc_cluster_constancia$n_days10 > 1,
      na.rm = TRUE
    )
) {
  stop(
    "La exposicion MITECO no es constante dentro de algun cluster isla-fecha."
  )
}

meta_cluster <- meta_base %>%
  group_by(
    cluster_isla_fecha
  ) %>%
  summarise(
    n_muestras = n(),

    island = first(
      as.character(
        island
      )
    ),

    fecha = first(
      fecha
    ),

    miteco10 = first(
      miteco10
    ),

    miteco_days10_z = first(
      miteco_days10_z
    ),

    naf10_z = first(
      naf10_z
    ),

    recent_D1_10_z = first(
      recent_D1_10_z
    ),

    residual_D11_30_z = first(
      residual_D11_30_z
    ),

    prior_D31_90_z = first(
      prior_D31_90_z
    ),

    year = first(
      as.character(
        year
      )
    ),

    sin_doy = first(
      sin_doy
    ),

    cos_doy = first(
      cos_doy
    ),

    age_mean = mean(
      age,
      na.rm = TRUE
    ),

    prop_male = mean(
      sex_num,
      na.rm = TRUE
    ),

    prop_antibiotic = mean(
      as.numeric(
        as.character(
          antibiotic
        )
      ),
      na.rm = TRUE
    ),

    prop_cavities = mean(
      as.numeric(
        as.character(
          cavities
        )
      ),
      na.rm = TRUE
    ),

    prop_gum = mean(
      as.numeric(
        as.character(
          gum
        )
      ),
      na.rm = TRUE
    ),

    prop_exac6m = mean(
      exac6m_num,
      na.rm = TRUE
    ),

    prop_pool2 = mean(
      as.character(
        pool
      ) == "2",
      na.rm = TRUE
    ),

    prop_pool3 = mean(
      as.character(
        pool
      ) == "3",
      na.rm = TRUE
    ),

    prop_pool6 = mean(
      as.character(
        pool
      ) == "6",
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    cluster_isla_fecha = as.character(
      cluster_isla_fecha
    ),

    island = factor(
      island
    ),

    year = factor(
      year
    )
  ) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )

rownames(
  meta_cluster
) <- meta_cluster$cluster_isla_fecha


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


permanova_cluster <- function(
  exposure
) {
  ids <- labels(
    dist_cluster_bray
  )

  dd <- meta_cluster[
    ids,
    ,
    drop = FALSE
  ]

  vars <- unique(
    c(
      exposure,
      covars_cluster
    )
  )

  cc <- complete.cases(
    dd[
      ,
      vars,
      drop = FALSE
    ]
  )

  dd <- droplevels(
    dd[
      cc,
      ,
      drop = FALSE
    ]
  )

  d2 <- subdist(
    dist_cluster_bray,
    rownames(dd)
  )

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  cov2 <- covars_activas(
    dd,
    covars_cluster
  )

  ftxt <- paste(
    "d2 ~",
    paste(
      c(
        exposure,
        cov2
      ),
      collapse = " + "
    )
  )

  fit <- vegan::adonis2(
    as.formula(
      ftxt
    ),
    data = dd,
    permutations = N_PERM,
    by = "margin"
  )

  tab <- as.data.frame(
    fit
  )

  tab$term <- rownames(
    tab
  )

  fila <- tab[
    tab$term == exposure,
    ,
    drop = FALSE
  ]

  tibble(
    distancia = "Bray_cluster_isla_fecha",
    exposure = exposure,
    n_clusters = nrow(dd),
    R2 = fila$R2,
    F = fila$F,
    p = fila$`Pr(>F)`,
    formula = ftxt
  )
}


beta_cluster <- bind_rows(
  lapply(
    c(
      "miteco10",
      "miteco_days10_z",
      "naf10_z"
    ),
    permanova_cluster
  )
) %>%
  mutate(
    p_BH = p.adjust(
      p,
      method = "BH"
    )
  )


# 21. CORE MICROBIOME

ps_saliva_g <- tax_glom(
  ps_saliva_u,
  taxrank = "Genus",
  NArm = TRUE
)

otu_sal_g <- otu_taxa_rows(
  ps_saliva_g
)

prev_core <- rowMeans(
  otu_sal_g > 0
)

core90_ids <- names(
  prev_core
)[
  prev_core >= 0.90
]

core95_ids <- names(
  prev_core
)[
  prev_core >= 0.95
]


metricas_core <- function(
  ps_obj,
  core_ids,
  nombre_core
) {
  pg <- tax_glom(
    ps_obj,
    taxrank = "Genus",
    NArm = TRUE
  )

  ids <- intersect(
    core_ids,
    taxa_names(
      pg
    )
  )

  o <- otu_taxa_rows(
    pg
  )

  if (
    length(ids) == 0
  ) {
    stop(
      "No hay taxa core tras intersectar para ",
      nombre_core
    )
  }

  present <- colSums(
    o[
      ids,
      ,
      drop = FALSE
    ] >
      0
  )

  core_reads <- colSums(
    o[
      ids,
      ,
      drop = FALSE
    ]
  )

  total <- colSums(
    o
  )

  tibble(
    id_muestra = colnames(o),
    core = nombre_core,
    n_core_total = length(ids),
    n_core_present = as.numeric(
      present
    ),
    frac_core_present = as.numeric(
      present
    ) /
      length(ids),
    prop_reads_core = as.numeric(
      core_reads
    ) /
      as.numeric(total)
  )
}


core_metrics <- bind_rows(
  metricas_core(
    ps_amb_u_m,
    core90_ids,
    "core90"
  ),
  metricas_core(
    ps_amb_u_m,
    core95_ids,
    "core95"
  )
) %>%
  left_join(
    meta_base,
    by = "id_muestra"
  )


core_outcomes <- c(
  "n_core_present",
  "prop_reads_core"
)


ajustar_core <- function(
  d,
  outcome,
  exposure
) {
  z <- ajustar_alpha_cont(
    outcome,
    exposure,
    data = d,
    covars = covars_alpha
  )

  if (is.null(z)) {
    return(NULL)
  }

  z
}


core_global <- bind_rows(
  lapply(
    c(
      "core90",
      "core95"
    ),
    function(cx) {
      d <- core_metrics %>%
        filter(
          core == cx
        )

      bind_rows(
        lapply(
          c(
            "miteco10",
            "miteco_days10_z"
          ),
          function(e) {
            bind_rows(
              lapply(
                core_outcomes,
                function(o) {
                  z <- ajustar_core(
                    d,
                    o,
                    e
                  )

                  if (is.null(z)) {
                    return(NULL)
                  }

                  z %>%
                    mutate(
                      core = cx,
                      .before = 1
                    )
                }
              )
            )
          }
        )
      )
    }
  )
) %>%
  group_by(
    core,
    exposure
  ) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()


core_int_sexo <- bind_rows(
  lapply(
    c(
      "core90",
      "core95"
    ),
    function(cx) {
      d <- core_metrics %>%
        filter(
          core == cx
        )

      bind_rows(
        lapply(
          core_outcomes,
          function(o) {
            z <- ajustar_alpha_interaccion(
              o,
              "miteco10",
              "sex",
              data = d
            )

            if (is.null(z)) {
              return(NULL)
            }

            z %>%
              mutate(
                core = cx,
                .before = 1
              )
          }
        )
      )
    }
  )
) %>%
  group_by(
    core
  ) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p_interaccion,
      method = "BH"
    )
  ) %>%
  ungroup()


# 22. COMPOSICION TAXONOMICA DESCRIPTIVA

ps_gen_desc <- tax_glom(
  ps_amb_u_m,
  taxrank = "Genus",
  NArm = TRUE
)

otu_desc <- otu_taxa_rows(
  ps_gen_desc
)

tt_desc <- tax_df(
  ps_gen_desc
)

rel_desc <- sweep(
  otu_desc,
  2,
  colSums(
    otu_desc
  ),
  "/"
)

genus_names <- tt_desc$Genus
genus_names[
  is.na(
    genus_names
  )
] <- "Genus_no_asignado"

genus_names <- make.unique(
  genus_names,
  sep = "__"
)

rownames(
  rel_desc
) <- genus_names

ra_mean <- rowMeans(
  rel_desc
)

top15 <- names(
  sort(
    ra_mean,
    decreasing = TRUE
  )
)[
  seq_len(
    min(
      15,
      length(
        ra_mean
      )
    )
  )
]


comp_long <- as.data.frame(
  t(
    rel_desc[
      top15,
      ,
      drop = FALSE
    ]
  )
)

comp_long$id_muestra <- rownames(
  comp_long
)

comp_long <- comp_long %>%
  left_join(
    meta_base %>%
      select(
        id_muestra,
        miteco10_factor,
        sex
      ),
    by = "id_muestra"
  ) %>%
  pivot_longer(
    cols = all_of(
      top15
    ),
    names_to = "Genus",
    values_to = "abundancia_relativa"
  )


comp_resumen <- comp_long %>%
  group_by(
    miteco10_factor,
    sex,
    Genus
  ) %>%
  summarise(
    n = n(),
    media = mean(
      abundancia_relativa,
      na.rm = TRUE
    ),
    mediana = median(
      abundancia_relativa,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 23. ANCOM-BC2: FUNCIONES
# ------------------------------------------------------------

prep_meta_ancom <- function(
  ps_obj
) {
  m <- sample_data_df(
    ps_obj
  )

  # Las variables ya estan integradas porque ps_amb_u_m se actualizo.
  m$age <- a_num(
    m$age
  )

  m$sex_num <- a_num(
    m$sex_num
  )

  m$miteco10 <- a_num(
    m$miteco10
  )

  m$miteco_days10_z <- a_num(
    m$miteco_days10_z
  )

  m$naf10_z <- a_num(
    m$naf10_z
  )

  m$exac6m_num <- a_num(
    m$exac6m_num
  )

  m$miteco_sex_int <- m$miteco10 *
    m$sex_num

  m$miteco_exac_int <- m$miteco10 *
    m$exac6m_num

  m$island <- factor(
    m$island
  )

  m$pool <- factor(
    m$pool
  )

  m$antibiotic <- factor(
    m$antibiotic
  )

  m$cavities <- factor(
    m$cavities
  )

  m$gum <- factor(
    m$gum
  )

  m$year <- factor(
    m$year
  )

  m
}


preparar_nivel_ancom <- function(
  ps_obj,
  nivel
) {
  if (
    nivel == "ASV"
  ) {
    p <- ps_obj
  } else {
    if (
      !nivel %in%
        rank_names(
          ps_obj
        )
    ) {
      stop(
        "No existe el rango ",
        nivel
      )
    }

    p <- tax_glom(
      ps_obj,
      taxrank = nivel,
      NArm = TRUE
    )
  }

  o <- otu_taxa_rows(
    p
  )

  nmin <- ceiling(
    PREV_ANCOM *
      ncol(o)
  )

  keep <- rownames(
    o
  )[
    rowSums(
      o > 0
    ) >=
      nmin
  ]

  p <- prune_taxa(
    keep,
    p
  )

  p <- prune_taxa(
    taxa_sums(
      p
    ) >
      0,
    p
  )

  p
}


ejecutar_ancom <- function(
  ps_obj,
  nivel,
  analisis,
  main_var,
  covars
) {
  meta <- prep_meta_ancom(
    ps_obj
  )

  vars <- c(
    main_var,
    covars
  )

  vars <- vars[
    vars %in% names(
      meta
    )
  ]

  cc <- complete.cases(
    meta[
      ,
      vars,
      drop = FALSE
    ]
  )

  meta2 <- droplevels(
    meta[
      cc,
      ,
      drop = FALSE
    ]
  )

  cov2 <- covars_activas(
    meta2,
    covars
  )

  vars_final <- c(
    main_var,
    cov2
  )

  cc <- complete.cases(
    meta[
      ,
      vars_final,
      drop = FALSE
    ]
  )

  ids <- meta$id_muestra[
    cc
  ]

  p <- prune_samples(
    ids,
    ps_obj
  )

  p <- prune_taxa(
    taxa_sums(
      p
    ) >
      0,
    p
  )

  p <- preparar_nivel_ancom(
    p,
    nivel
  )

  meta3 <- prep_meta_ancom(
    p
  )

  rownames(
    meta3
  ) <- meta3$id_muestra

  sample_data(
    p
  ) <- sample_data(
    meta3
  )

  fix_formula <- paste(
    c(
      main_var,
      cov2
    ),
    collapse = " + "
  )

  fit <- ANCOMBC::ancombc2(
    data = p,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = fix_formula,
    rand_formula = NULL,
    p_adj_method = "BH",
    pseudo_sens = TRUE,
    prv_cut = 0,
    lib_cut = 0,
    s0_perc = 0.05,
    group = NULL,
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = ALPHA_FDR,
    n_cl = 1,
    verbose = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE
  )

  saveRDS(
    fit,
    file.path(
      ruta_objetos,
      paste0(
        "ANCOMBC2_",
        analisis,
        "_",
        nivel,
        ".rds"
      )
    )
  )

  rr <- as.data.frame(
    fit$res,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (
    "taxon" %in%
      names(rr)
  ) {
    rr$taxon_id <- as.character(
      rr$taxon
    )
  } else {
    rr$taxon_id <- rownames(
      rr
    )
  }

  lfc_cols <- grep(
    "^lfc_",
    names(rr),
    value = TRUE
  )

  lc <- lfc_cols[
    grepl(
      main_var,
      lfc_cols,
      fixed = TRUE
    )
  ]

  if (
    length(lc) != 1
  ) {
    stop(
      "No se identifica de forma unica el coeficiente ",
      main_var,
      " en ",
      analisis,
      " ",
      nivel,
      ". Candidatos: ",
      paste(
        lc,
        collapse = ", "
      )
    )
  }

  suf <- sub(
    "^lfc_",
    "",
    lc
  )

  getcol <- function(
    pref
  ) {
    z <- paste0(
      pref,
      "_",
      suf
    )

    if (
      z %in%
        names(rr)
    ) {
      z
    } else {
      NA_character_
    }
  }

  cq <- getcol(
    "q"
  )

  cp <- getcol(
    "p"
  )

  cd <- getcol(
    "diff"
  )

  cpass <- getcol(
    "passed_ss"
  )

  crob <- getcol(
    "diff_robust"
  )

  if (
    is.na(cq) ||
      is.na(cp)
  ) {
    stop(
      "No se encuentran columnas p/q en ",
      analisis,
      " ",
      nivel
    )
  }

  out <- data.frame(
    analisis = analisis,
    nivel = nivel,
    n_muestras = nsamples(
      p
    ),
    n_taxa = ntaxa(
      p
    ),
    taxon_id = rr$taxon_id,
    beta = rr[[lc]],
    p = rr[[cp]],
    q = rr[[cq]],
    diff_q = if (
      !is.na(cd)
    ) {
      rr[[cd]]
    } else {
      NA
    },
    passed_sensitivity = if (
      !is.na(cpass)
    ) {
      rr[[cpass]]
    } else {
      NA
    },
    diff_robust = if (
      !is.na(crob)
    ) {
      rr[[crob]]
    } else {
      NA
    },
    stringsAsFactors = FALSE
  )

  tt <- as.data.frame(
    tax_table(
      p
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  tt$taxon_id <- rownames(
    tt
  )

  out <- left_join(
    out,
    tt,
    by = "taxon_id"
  )

  out$significativo_q <- !is.na(
    out$q
  ) &
    out$q <
      ALPHA_FDR

  out$robusto <- ifelse(
    !is.na(
      out$diff_robust
    ),
    as.logical(
      out$diff_robust
    ),
    out$significativo_q &
      (
        is.na(
          out$passed_sensitivity
        ) |
          as.logical(
            out$passed_sensitivity
          )
      )
  )

  out
}


cov_ancom_global <- c(
  "age",
  "sex_num",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m_num",
  "year",
  "sin_doy",
  "cos_doy"
)


cov_ancom_sexint <- c(
  "miteco10",
  "sex_num",
  "age",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m_num",
  "year",
  "sin_doy",
  "cos_doy"
)


cov_ancom_exacint <- c(
  "miteco10",
  "exac6m_num",
  "age",
  "sex_num",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "year",
  "sin_doy",
  "cos_doy"
)


# 24. ANCOM-BC2 GENERO

ancom_genus <- list()

ancom_genus[[1]] <- ejecutar_ancom(
  ps_amb_u_m,
  "Genus",
  "MITECO10",
  "miteco10",
  cov_ancom_global
)

ancom_genus[[2]] <- ejecutar_ancom(
  ps_amb_u_m,
  "Genus",
  "MITECO_dias10",
  "miteco_days10_z",
  cov_ancom_global
)

ancom_genus[[3]] <- ejecutar_ancom(
  ps_amb_u_m,
  "Genus",
  "NAF10",
  "naf10_z",
  cov_ancom_global
)

ancom_genus[[4]] <- ejecutar_ancom(
  ps_amb_u_m,
  "Genus",
  "MITECO10_x_sexo",
  "miteco_sex_int",
  cov_ancom_sexint
)

ancom_genus[[5]] <- ejecutar_ancom(
  ps_amb_u_m,
  "Genus",
  "MITECO10_x_exac6m",
  "miteco_exac_int",
  cov_ancom_exacint
)

ancom_genus_all <- bind_rows(
  ancom_genus
)

ancom_genus_resumen <- ancom_genus_all %>%
  group_by(
    analisis,
    nivel
  ) %>%
  summarise(
    n_muestras = first(
      n_muestras
    ),

    n_taxa = first(
      n_taxa
    ),

    n_q005 = sum(
      significativo_q,
      na.rm = TRUE
    ),

    n_robustos = sum(
      robusto,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

registrar(
  "ANCOM-BC2 Genus",
  "OK",
  paste0(
    "Modelos=",
    n_distinct(
      ancom_genus_all$analisis
    )
  )
)


# 25. ANCOM-BC2 SENSIBILIDAD FAMILY / GENUS / ASV


niveles_ancom <- c(
  "Family",
  "Genus"
)

if (
  EJECUTAR_ANCOM_ASV
) {
  niveles_ancom <- c(
    niveles_ancom,
    "ASV"
  )
}

ancom_sens <- list()
kk <- 1

for (niv in niveles_ancom) {
  ancom_sens[[kk]] <- ejecutar_ancom(
    ps_amb_u_m,
    niv,
    "MITECO10_sensibilidad",
    "miteco10",
    cov_ancom_global
  )

  kk <- kk +
    1

  ancom_sens[[kk]] <- ejecutar_ancom(
    ps_amb_u_m,
    niv,
    "MITECO10_x_sexo_sensibilidad",
    "miteco_sex_int",
    cov_ancom_sexint
  )

  kk <- kk +
    1
}

ancom_sens_all <- bind_rows(
  ancom_sens
)

ancom_sens_resumen <- ancom_sens_all %>%
  group_by(
    analisis,
    nivel
  ) %>%
  summarise(
    n_muestras = first(
      n_muestras
    ),

    n_taxa = first(
      n_taxa
    ),

    n_q005 = sum(
      significativo_q,
      na.rm = TRUE
    ),

    n_robustos = sum(
      robusto,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

registrar(
  "ANCOM-BC2 sensibilidad taxonomica",
  "OK",
  paste(
    niveles_ancom,
    collapse = "/"
  )
)


# 26. EXACERBACION 6M

# 26A. Microbiota alfa ~ exacerbacion, añadiendo carga MITECO 180d como covariable ambiental y manteniendo sensibilidad con/sin ATB.

ajustar_exac_alpha <- function(
  outcome,
  incluir_antibiotico = TRUE
) {
  cov_exac <- c(
    "age",
    "sex",
    "island",
    "pool",
    "cavities",
    "gum",
    "lecturas_log10",
    "miteco_days180_z",
    "year",
    "sin_doy",
    "cos_doy"
  )

  if (
    incluir_antibiotico
  ) {
    cov_exac <- c(
      cov_exac,
      "antibiotic"
    )
  }

  vars <- unique(
    c(
      outcome,
      "exac6m",
      cov_exac,
      "cluster_isla_fecha"
    )
  )

  dd <- metricas[
    complete.cases(
      metricas[
        ,
        vars,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  cov2 <- covars_activas(
    dd,
    cov_exac
  )

  fit <- lm(
    as.formula(
      paste(
        outcome,
        "~ exac6m +",
        paste(
          cov2,
          collapse = " + "
        )
      )
    ),
    data = dd
  )

  term <- names(
    coef(
      fit
    )
  )

  term <- term[
    grepl(
      "^exac6m",
      term
    )
  ]

  if (
    length(term) != 1
  ) {
    return(NULL)
  }

  z <- extraer_coef_robusto(
    fit,
    dd$cluster_isla_fecha,
    term
  )

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    mutate(
      outcome = outcome,
      modelo = ifelse(
        incluir_antibiotico,
        "Con_antibiotico",
        "Sin_antibiotico"
      ),
      n = nrow(dd),
      n_clusters = n_distinct(
        dd$cluster_isla_fecha
      ),
      .before = 1
    )
}


exac_alpha <- bind_rows(
  lapply(
    alpha_outcomes,
    function(o) {
      bind_rows(
        ajustar_exac_alpha(
          o,
          TRUE
        ),
        ajustar_exac_alpha(
          o,
          FALSE
        )
      )
    }
  )
) %>%
  group_by(
    modelo
  ) %>%
  mutate(
    p_BH_outcomes = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()


# 26B. Carga ambiental MITECO180 -> exacerbacion 6M.
ajustar_exac_outcome <- function(
  incluir_antibiotico = TRUE
) {
  cov_exac <- c(
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

  if (
    incluir_antibiotico
  ) {
    cov_exac <- c(
      cov_exac,
      "antibiotic"
    )
  }

  vars <- unique(
    c(
      "exac6m",
      "miteco_days180_z",
      cov_exac,
      "cluster_isla_fecha"
    )
  )

  dd <- metricas[
    complete.cases(
      metricas[
        ,
        vars,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  dd <- droplevels(
    dd
  )

  cov2 <- covars_activas(
    dd,
    cov_exac
  )

  fit <- glm(
    as.formula(
      paste(
        "exac6m ~ miteco_days180_z +",
        paste(
          cov2,
          collapse = " + "
        )
      )
    ),
    family = binomial(),
    data = dd
  )

  z <- extraer_coef_robusto(
    fit,
    dd$cluster_isla_fecha,
    "miteco_days180_z"
  )

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    transmute(
      modelo = ifelse(
        incluir_antibiotico,
        "Con_antibiotico",
        "Sin_antibiotico"
      ),
      n = nrow(dd),
      n_clusters = n_distinct(
        dd$cluster_isla_fecha
      ),
      beta = beta,
      se = se,
      OR_por_1DE = exp(
        beta
      ),
      OR_IC95_inf = exp(
        IC95_inf
      ),
      OR_IC95_sup = exp(
        IC95_sup
      ),
      p = p
    )
}


exac_miteco180 <- bind_rows(
  ajustar_exac_outcome(
    TRUE
  ),
  ajustar_exac_outcome(
    FALSE
  )
)


# 27. ESTADOS COMUNITARIOS

states_res <- tibble(
  disponible = FALSE,
  prueba = NA_character_,
  p = NA_real_,
  nota = "No se encontro el archivo previo de estados comunitarios."
)

if (
  !is.na(ruta_states) &&
    file.exists(
      ruta_states
    )
) {
  st <- readr::read_csv(
    ruta_states,
    show_col_types = FALSE
  )

  # localizar columna ID de forma tolerante
  id_col <- intersect(
    c(
      "id_muestra",
      "Sample",
      "sample"
    ),
    names(st)
  )

  state_col <- intersect(
    c(
      "community_state",
      "estado",
      "cluster"
    ),
    names(st)
  )

  if (
    length(id_col) >= 1 &&
      length(state_col) >= 1
  ) {
    st2 <- st %>%
      transmute(
        id_muestra = .data[[id_col[1]]],
        community_state = factor(
          .data[[state_col[1]]]
        )
      ) %>%
      left_join(
        meta_base %>%
          select(
            id_muestra,
            miteco10_factor
          ),
        by = "id_muestra"
      ) %>%
      filter(
        !is.na(
          community_state
        ),
        !is.na(
          miteco10_factor
        )
      )

    tab_state <- table(
      st2$community_state,
      st2$miteco10_factor
    )

    test_state <- if (
      any(
        chisq.test(
          tab_state,
          simulate.p.value = FALSE
        )$expected <
          5
      )
    ) {
      fisher.test(
        tab_state
      )
    } else {
      chisq.test(
        tab_state
      )
    }

    states_res <- tibble(
      disponible = TRUE,
      prueba = class(
        test_state
      )[1],
      p = test_state$p.value,
      nota = "Asociacion simple; los estados previos tienen baja separacion (silhouette) y se consideran exploratorios."
    )

    registrar(
      "Estados comunitarios",
      "OK",
      paste0(
        "p simple=",
        signif(
          test_state$p.value,
          4
        )
      )
    )
  } else {
    registrar(
      "Estados comunitarios",
      "OMITIDO",
      "No se localizaron columnas ID/estado en el output previo."
    )
  }
} else {
  registrar(
    "Estados comunitarios",
    "OMITIDO",
    "No existe output previo."
  )
}


# 28. FIGURAS DESCRIPTIVAS PRINCIPALES


fig_alpha <- metricas %>%
  filter(
    !is.na(
      miteco10_factor
    ),
    !is.na(
      sex
    )
  ) %>%
  ggplot(
    aes(
      x = miteco10_factor,
      y = Shannon,
      fill = sex
    )
  ) +
  geom_boxplot(
    outlier.alpha = 0.35
  ) +
  facet_wrap(
    ~ sex
  ) +
  labs(
    title = "Diversidad Shannon segun episodio MITECO D-1 a D-10 y sexo",
    x = NULL,
    y = "Shannon"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(
    ruta_figuras,
    "01_Shannon_MITECO10_por_sexo.png"
  ),
  fig_alpha,
  width = 9,
  height = 5.5,
  dpi = 300
)


fig_days <- metricas %>%
  filter(
    !is.na(
      sex
    )
  ) %>%
  ggplot(
    aes(
      x = miteco_days_10d,
      y = Shannon
    )
  ) +
  geom_point(
    alpha = 0.45
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  facet_wrap(
    ~ sex
  ) +
  labs(
    title = "Shannon y numero de dias MITECO en D-1 a D-10",
    x = "N dias oficiales MITECO",
    y = "Shannon"
  ) +
  theme_minimal()

ggsave(
  file.path(
    ruta_figuras,
    "02_Shannon_dias_MITECO10_por_sexo.png"
  ),
  fig_days,
  width = 9,
  height = 5.5,
  dpi = 300
)


# 29. EXPORTAR TABLAS

write_csv(
  metricas,
  file.path(
    ruta_tablas,
    "01_dataset_metricas_MITECO_analitico.csv"
  ),
  na = ""
)

write_csv(
  alpha_global,
  file.path(
    ruta_tablas,
    "02_alpha_global_MITECO.csv"
  ),
  na = ""
)

write_csv(
  alpha_factor,
  file.path(
    ruta_tablas,
    "03_alpha_grado_patron_MITECO.csv"
  ),
  na = ""
)

write_csv(
  alpha_int_sexo,
  file.path(
    ruta_tablas,
    "04_alpha_interaccion_MITECO_sexo.csv"
  ),
  na = ""
)

write_csv(
  alpha_sexo_estrat,
  file.path(
    ruta_tablas,
    "05_alpha_estratificado_sexo.csv"
  ),
  na = ""
)

write_csv(
  alpha_int_isla,
  file.path(
    ruta_tablas,
    "06_alpha_interaccion_MITECO_isla.csv"
  ),
  na = ""
)

write_csv(
  alpha_int_exac,
  file.path(
    ruta_tablas,
    "07_alpha_interaccion_MITECO_exac6m.csv"
  ),
  na = ""
)

write_csv(
  alpha_int_regimen,
  file.path(
    ruta_tablas,
    "08_alpha_interaccion_MITECO_winter_summer.csv"
  ),
  na = ""
)

write_csv(
  alpha_ventanas,
  file.path(
    ruta_tablas,
    "09_alpha_ventanas_D1_10_D11_30_D31_90.csv"
  ),
  na = ""
)

write_csv(
  alpha_ventanas_sexo,
  file.path(
    ruta_tablas,
    "10_alpha_ventanas_interaccion_sexo.csv"
  ),
  na = ""
)

write_csv(
  spline_res,
  file.path(
    ruta_tablas,
    "11_alpha_dosis_respuesta_splines_MITECO.csv"
  ),
  na = ""
)

write_csv(
  beta_global,
  file.path(
    ruta_tablas,
    "12_beta_global_MITECO.csv"
  ),
  na = ""
)

write_csv(
  beta_int_sexo,
  file.path(
    ruta_tablas,
    "13_beta_interaccion_MITECO_sexo.csv"
  ),
  na = ""
)

write_csv(
  beta_int_isla,
  file.path(
    ruta_tablas,
    "14_beta_interaccion_MITECO_isla.csv"
  ),
  na = ""
)

write_csv(
  beta_int_exac,
  file.path(
    ruta_tablas,
    "15_beta_interaccion_MITECO_exac6m.csv"
  ),
  na = ""
)

write_csv(
  beta_int_regimen,
  file.path(
    ruta_tablas,
    "16_beta_interaccion_MITECO_winter_summer.csv"
  ),
  na = ""
)

write_csv(
  beta_ventanas,
  file.path(
    ruta_tablas,
    "17_beta_ventanas_temporales_MITECO.csv"
  ),
  na = ""
)

write_csv(
  betadisper_res,
  file.path(
    ruta_tablas,
    "18_betadisper_MITECO.csv"
  ),
  na = ""
)

write_csv(
  beta_cluster,
  file.path(
    ruta_tablas,
    "19_beta_sensibilidad_cluster_isla_fecha.csv"
  ),
  na = ""
)

write_csv(
  qc_cluster_constancia,
  file.path(
    ruta_tablas,
    "19b_QC_constancia_MITECO_cluster_isla_fecha.csv"
  ),
  na = ""
)

write_csv(
  core_global,
  file.path(
    ruta_tablas,
    "20_core_microbiome_global_MITECO.csv"
  ),
  na = ""
)

write_csv(
  core_int_sexo,
  file.path(
    ruta_tablas,
    "21_core_microbiome_interaccion_sexo.csv"
  ),
  na = ""
)

write_csv(
  comp_resumen,
  file.path(
    ruta_tablas,
    "22_composicion_top15_generos_MITECO_sexo.csv"
  ),
  na = ""
)

write_csv(
  ancom_genus_all,
  file.path(
    ruta_tablas,
    "23_ANCOMBC2_Genus_todos.csv"
  ),
  na = ""
)

write_csv(
  ancom_genus_resumen,
  file.path(
    ruta_tablas,
    "24_ANCOMBC2_Genus_resumen.csv"
  ),
  na = ""
)

write_csv(
  ancom_genus_all %>%
    filter(
      robusto
    ),
  file.path(
    ruta_tablas,
    "25_ANCOMBC2_Genus_robustos.csv"
  ),
  na = ""
)

write_csv(
  ancom_sens_all,
  file.path(
    ruta_tablas,
    "26_ANCOMBC2_sensibilidad_taxonomica_todos.csv"
  ),
  na = ""
)

write_csv(
  ancom_sens_resumen,
  file.path(
    ruta_tablas,
    "27_ANCOMBC2_sensibilidad_taxonomica_resumen.csv"
  ),
  na = ""
)

write_csv(
  exac_alpha,
  file.path(
    ruta_tablas,
    "28_exacerbacion6M_alpha_ajustada_MITECO.csv"
  ),
  na = ""
)

write_csv(
  exac_miteco180,
  file.path(
    ruta_tablas,
    "29_exacerbacion6M_carga_MITECO180.csv"
  ),
  na = ""
)

write_csv(
  states_res,
  file.path(
    ruta_tablas,
    "30_estados_comunitarios_vs_MITECO.csv"
  ),
  na = ""
)


# 30. RESUMEN DE HALLAZGOS PARA REDACCION

resumen_redaccion <- bind_rows(
  alpha_global %>%
    transmute(
      bloque = "Alpha_global",
      analisis = paste(
        outcome,
        exposure,
        sep = " ~ "
      ),
      estimador = beta,
      p = p,
      q = p_BH_outcomes
    ),

  alpha_int_sexo %>%
    transmute(
      bloque = "Alpha_interaccion_sexo",
      analisis = paste(
        outcome,
        exposure,
        "x sexo"
      ),
      estimador = chi2,
      p = p_interaccion,
      q = p_BH_outcomes
    ),

  alpha_ventanas %>%
    transmute(
      bloque = "Alpha_ventanas",
      analisis = paste(
        outcome,
        window,
        sep = " ~ "
      ),
      estimador = beta,
      p = p,
      q = p_BH_outcomes
    ),

  beta_global %>%
    transmute(
      bloque = "Beta_global",
      analisis = paste(
        distancia,
        exposure
      ),
      estimador = R2,
      p = p,
      q = p_BH_distancias
    ),

  beta_int_sexo %>%
    transmute(
      bloque = "Beta_interaccion_sexo",
      analisis = paste(
        distancia,
        "MITECO10 x sexo"
      ),
      estimador = R2,
      p = p,
      q = p_BH_distancias
    ),

  core_global %>%
    transmute(
      bloque = "Core",
      analisis = paste(
        core,
        outcome,
        exposure
      ),
      estimador = beta,
      p = p,
      q = p_BH_outcomes
    )
)


write_csv(
  resumen_redaccion,
  file.path(
    ruta_tablas,
    "31_resumen_resultados_para_redaccion.csv"
  ),
  na = ""
)


# 33.  FINAL

cat(
  "\n============================================================\n"
)

cat(
  "REANALISIS COMPLETO MITECO FINALIZADO\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "QC exposicion:\n"
)

print(
  qc_exposicion
)

cat(
  "\nAlpha global:\n"
)

print(
  alpha_global
)

cat(
  "\nAlpha MITECO x sexo:\n"
)

print(
  alpha_int_sexo
)

cat(
  "\nBeta global:\n"
)

print(
  beta_global
)

cat(
  "\nBeta MITECO x sexo:\n"
)

print(
  beta_int_sexo
)

cat(
  "\nANCOM-BC2 resumen Genero:\n"
)

print(
  ancom_genus_resumen
)

cat(
  "\nANCOM-BC2 sensibilidad taxonomica:\n"
)

print(
  ancom_sens_resumen
)

cat(
  "\nManifest:\n"
)

print(
  manifest
)

cat(
  "\nArchivo maestro:\n",
  file.path(
    ruta_tablas,
    "33_RESUMEN_REANALISIS_COMPLETO_MITECO_GEMAS.xlsx"
  ),
  "\n"
)
