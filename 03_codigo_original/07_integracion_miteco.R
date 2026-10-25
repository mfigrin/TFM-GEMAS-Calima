# INTEGRAR_MITECO_PM10_GEMAS.R
# TFM GEMAS
#
# OBJETIVO
# Integrar dos dimensiones ambientales complementarias:
#
# 1) PRESENCIA OFICIAL DE EPISODIO AFRICANO (MITECO)
#    - calendario oficial diario para Canarias 2016-2020
#    - independiente del umbral PM10 >=50
#
# 2) IMPACTO / INTENSIDAD EN SUPERFICIE
#    a) PM10 local ya usado en el TFM
#    b) aporte africano diario oficial ("Descuento") de MITECO:
#       - Tenerife: El Rio
#       - La Palma: Echedo-Valverde
#         (proxy MITECO para El Hierro, La Palma y La Gomera)
#
# La salida permite distinguir:
# - Sin episodio MITECO
# - MITECO sin PM10 local >=50
# - MITECO + PM10 local >=50
# - PM10 >=50 sin MITECO (discordancia/QC)

# 0. PAQUETES

paquetes <- c(
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "stringr",
  "readxl",
  "writexl"
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
    "\nInstalalos antes de ejecutar."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(writexl)
})

# 1. RUTAS

ruta_miteco <- file.path(
  "data_raw",
  "CALIMA",
  "CALIMA MITECO 2016-2020"
)

ruta_calendario <- file.path(
  ruta_miteco,
  "00_MITECO_Canarias_calendario_africanos_2016_2020.csv"
)

ruta_traza_pm10 <- file.path(
  "outputs",
  "analisis_microbiota",
  "patrones_exposicion_pm10",
  "tablas",
  "01_trazabilidad_PM10_D1_D180.csv"
)

ruta_caract_pm10 <- file.path(
  "outputs",
  "analisis_microbiota",
  "patrones_exposicion_pm10",
  "tablas",
  "03_caracteristicas_episodios_PM10_por_muestra.csv"
)

for (f in c(
  ruta_calendario,
  ruta_traza_pm10,
  ruta_caract_pm10
)) {
  if (!file.exists(f)) {
    stop(
      "No existe: ",
      f
    )
  }
}

ruta_salida <- file.path(
  "outputs",
  "analisis_microbiota",
  "miteco_integrado"
)

ruta_tablas <- file.path(
  ruta_salida,
  "tablas"
)

ruta_objetos <- file.path(
  ruta_salida,
  "objetos"
)

dir.create(
  ruta_tablas,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ruta_objetos,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. FUNCIONES

normalizar_texto <- function(x) {

  y <- as.character(
    x
  )

  y <- iconv(
    y,
    from = "",
    to = "ASCII//TRANSLIT"
  )

  y <- tolower(
    trimws(
      y
    )
  )

  y <- gsub(
    "[^a-z0-9]+",
    " ",
    y
  )

  trimws(
    y
  )
}


a_num <- function(x) {

  if (is.numeric(x)) {
    return(
      as.numeric(
        x
      )
    )
  }

  y <- trimws(
    as.character(
      x
    )
  )

  y[
    y == ""
  ] <- NA_character_

  y <- gsub(
    ",",
    ".",
    y,
    fixed = TRUE
  )

  suppressWarnings(
    as.numeric(
      y
    )
  )
}


a_logico <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(
    trimws(
      as.character(
        x
      )
    )
  )

  out <- rep(
    NA,
    length(
      y
    )
  )

  out[
    y %in% c(
      "TRUE",
      "T",
      "1",
      "SI",
      "SÍ",
      "YES"
    )
  ] <- TRUE

  out[
    y %in% c(
      "FALSE",
      "F",
      "0",
      "NO"
    )
  ] <- FALSE

  as.logical(
    out
  )
}


parse_fecha_excel <- function(x) {

  if (inherits(x, "Date")) {
    return(
      as.Date(
        x
      )
    )
  }

  if (
    inherits(x, "POSIXct") ||
      inherits(x, "POSIXt")
  ) {
    return(
      as.Date(
        x
      )
    )
  }

  # Intentar fecha textual ISO / habitual.
  y <- trimws(
    as.character(
      x
    )
  )

  out <- suppressWarnings(
    as.Date(
      y
    )
  )

  # Excel serial si no fue posible.
  idx <- is.na(out)

  if (any(idx)) {

    num <- suppressWarnings(
      as.numeric(
        y[
          idx
        ]
      )
    )

    out[
      idx
    ] <- as.Date(
      num,
      origin = "1899-12-30"
    )
  }

  out
}


encontrar_bloque_metrica <- function(
  metricas,
  tipo = c(
    "pm10",
    "perc40",
    "descuento"
  )
) {

  tipo <- match.arg(
    tipo
  )

  if (
    tipo == "pm10"
  ) {
    idx <- which(
      grepl(
        "^pm10$",
        metricas
      )
    )
  }

  if (
    tipo == "perc40"
  ) {
    idx <- which(
      grepl(
        "^perc ?40$",
        metricas
      ) |
        grepl(
          "percentil.*40",
          metricas
        )
    )
  }

  if (
    tipo == "descuento"
  ) {
    idx <- which(
      grepl(
        "descuento",
        metricas
      )
    )
  }

  idx
}


localizar_estacion_en_bloque_pm10 <- function(
  estaciones,
  idx_pm10,
  patron_estacion
) {

  estaciones_pm10 <- estaciones[
    idx_pm10
  ]

  pos <- which(
    grepl(
      patron_estacion,
      estaciones_pm10
    )
  )

  if (
    length(
      pos
    ) != 1
  ) {
    stop(
      "No se encontro una posicion unica para la estacion '",
      patron_estacion,
      "' dentro del bloque PM10. Posiciones: ",
      paste(
        pos,
        collapse = ", "
      ),
      "\nEstaciones PM10 detectadas: ",
      paste(
        estaciones_pm10,
        collapse = " | "
      )
    )
  }

  pos
}


leer_canarias_miteco <- function(
  archivo
) {

  cat(
    "Leyendo: ",
    basename(
      archivo
    ),
    "\n",
    sep = ""
  )

  raw <- readxl::read_excel(
    archivo,
    sheet = "CANARIAS",
    col_names = FALSE
  )

  if (nrow(raw) < 10) {
    stop(
      "Hoja CANARIAS inesperadamente corta en ",
      basename(
        archivo
      )
    )
  }

  estaciones <- normalizar_texto(
    unlist(
      raw[
        1,
        ,
        drop = TRUE
      ]
    )
  )

  metricas <- normalizar_texto(
    unlist(
      raw[
        2,
        ,
        drop = TRUE
      ]
    )
  )



  idx_pm10 <- encontrar_bloque_metrica(
    metricas,
    "pm10"
  )

  idx_perc <- encontrar_bloque_metrica(
    metricas,
    "perc40"
  )

  idx_desc <- encontrar_bloque_metrica(
    metricas,
    "descuento"
  )

  if (
    length(
      idx_pm10
    ) == 0
  ) {
    stop(
      "No se detecto el bloque PM10 en ",
      basename(
        archivo
      ),
      ". Metricas leidas: ",
      paste(
        metricas,
        collapse = " | "
      )
    )
  }

  tiene_perc <- length(
    idx_perc
  ) > 0

  if (
    length(
      idx_desc
    ) == 0
  ) {
    stop(
      "No se detecto el bloque Descuento en ",
      basename(
        archivo
      ),
      ". Metricas leidas: ",
      paste(
        metricas,
        collapse = " | "
      )
    )
  }

  if (
    length(
      idx_desc
    ) !=
      length(
        idx_pm10
      )
  ) {
    stop(
      "Los bloques PM10 y Descuento no tienen la misma longitud en ",
      basename(
        archivo
      ),
      ": PM10=",
      length(
        idx_pm10
      ),
      "; Descuento=",
      length(
        idx_desc
      )
    )
  }

  if (
    tiene_perc &&
      length(
        idx_perc
      ) !=
        length(
          idx_pm10
        )
  ) {
    warning(
      "Bloque PERC40 con longitud distinta en ",
      basename(
        archivo
      ),
      ". Se ignorara PERC40 para ese año."
    )

    tiene_perc <- FALSE
  }

  pos_elrio <- localizar_estacion_en_bloque_pm10(
    estaciones,
    idx_pm10,
    "el rio"
  )

  pos_echedo <- localizar_estacion_en_bloque_pm10(
    estaciones,
    idx_pm10,
    "echedo"
  )

  i_elrio_pm10 <- idx_pm10[
    pos_elrio
  ]

  i_elrio_desc <- idx_desc[
    pos_elrio
  ]

  i_echedo_pm10 <- idx_pm10[
    pos_echedo
  ]

  i_echedo_desc <- idx_desc[
    pos_echedo
  ]

  i_elrio_perc <- if (
    tiene_perc
  ) {
    idx_perc[
      pos_elrio
    ]
  } else {
    NA_integer_
  }

  i_echedo_perc <- if (
    tiene_perc
  ) {
    idx_perc[
      pos_echedo
    ]
  } else {
    NA_integer_
  }

  cat(
    "  Bloques detectados: PM10=",
    length(
      idx_pm10
    ),
    "; PERC40=",
    length(
      idx_perc
    ),
    "; Descuento=",
    length(
      idx_desc
    ),
    "
",
    sep = ""
  )

  cat(
    "  El Rio columnas: PM10=",
    i_elrio_pm10,
    "; PERC40=",
    ifelse(
      is.na(
        i_elrio_perc
      ),
      "NA",
      i_elrio_perc
    ),
    "; Descuento=",
    i_elrio_desc,
    "
",
    sep = ""
  )

  cat(
    "  Echedo columnas: PM10=",
    i_echedo_pm10,
    "; PERC40=",
    ifelse(
      is.na(
        i_echedo_perc
      ),
      "NA",
      i_echedo_perc
    ),
    "; Descuento=",
    i_echedo_desc,
    "
",
    sep = ""
  )

  fecha <- parse_fecha_excel(
    raw[[1]][
      -c(
        1,
        2
      )
    ]
  )

  extraer_estacion <- function(
    nombre,
    idx_pm10,
    idx_perc,
    idx_desc
  ) {

    pm10 <- a_num(
      raw[[idx_pm10]][
        -c(
          1,
          2
        )
      ]
    )

    perc40 <- if (
      is.na(
        idx_perc
      )
    ) {
      rep(
        NA_real_,
        length(
          fecha
        )
      )
    } else {
      a_num(
        raw[[idx_perc]][
          -c(
            1,
            2
          )
        ]
      )
    }

    aporte <- a_num(
      raw[[idx_desc]][
        -c(
          1,
          2
        )
      ]
    )

    # QC físico:
    # el aporte natural no puede exceder el PM10 observado.
    aporte[
      !is.na(
        aporte
      ) &
        (
          is.na(
            pm10
          ) |
            aporte <
              -0.001 |
            aporte >
              pm10 +
                0.5
        )
    ] <- NA_real_

    tibble(
      fecha = fecha,
      estacion_miteco = nombre,
      pm10_miteco = pm10,
      perc40_miteco = perc40,
      aporte_naf_miteco = aporte
    ) %>%
      filter(
        !is.na(
          fecha
        )
      )
  }

  bind_rows(
    extraer_estacion(
      "El_Rio",
      i_elrio_pm10,
      i_elrio_perc,
      i_elrio_desc
    ),
    extraer_estacion(
      "Echedo_Valverde",
      i_echedo_pm10,
      i_echedo_perc,
      i_echedo_desc
    )
  )
}


contar_episodios <- function(
  ids
) {

  ids <- ids[
    !is.na(
      ids
    )
  ]

  if (length(ids) == 0) {
    return(
      0L
    )
  }

  length(
    unique(
      ids
    )
  )
}


# 3. CALENDARIO OFICIAL MITECO

cal <- readr::read_csv(
  ruta_calendario,
  show_col_types = FALSE
) %>%
  mutate(
    fecha = as.Date(
      fecha
    ),
    miteco_africano = as.integer(
      miteco_africano
    )
  )

if (
  min(
    cal$fecha
  ) != as.Date(
    "2016-01-01"
  ) ||
    max(
      cal$fecha
    ) != as.Date(
      "2020-12-31"
    )
) {
  stop(
    "El calendario MITECO no cubre exactamente 2016-2020."
  )
}

qc_anual <- cal %>%
  group_by(
    anio
  ) %>%
  summarise(
    dias_miteco = sum(
      miteco_africano,
      na.rm = TRUE
    ),
    dias_anio = n(),
    porcentaje = 100 *
      dias_miteco /
      dias_anio,
    .groups = "drop"
  )

esperado <- tibble(
  anio = 2016:2020,
  dias_esperados = c(
    154,
    201,
    129,
    149,
    171
  ),
  pct_reportado = c(
    42,
    55,
    35,
    41,
    47
  )
)

qc_anual <- qc_anual %>%
  left_join(
    esperado,
    by = "anio"
  ) %>%
  mutate(
    coincide_dias = dias_miteco ==
      dias_esperados
  )

if (!all(qc_anual$coincide_dias)) {
  stop(
    "El calendario MITECO no reproduce los dias anuales esperados."
  )
}


# 4. APORTES NAF DE LAS HOJAS MITECO

archivos_descuentos <- list.files(
  ruta_miteco,
  pattern = "^pm10-descuentos_todas_estaciones_20(16|17|18|19|20)_.*\\.xls[x]?$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(archivos_descuentos) != 5) {
  stop(
    "Se esperaban 5 ficheros pm10-descuentos 2016-2020 y se encontraron ",
    length(
      archivos_descuentos
    )
  )
}

miteco_estaciones <- bind_rows(
  lapply(
    sort(
      archivos_descuentos
    ),
    leer_canarias_miteco
  )
) %>%
  left_join(
    cal %>%
      select(
        fecha,
        miteco_africano,
        episodio_miteco_id
      ),
    by = "fecha"
  ) %>%
  mutate(
    # Fuera de episodio oficial el aporte NAF se conoce como 0.
    # Dentro de episodio, si la estacion carece de dato se mantiene NA.
    aporte_naf_usable = case_when(
      miteco_africano == 0 ~ 0,
      miteco_africano == 1 &
        !is.na(
          aporte_naf_miteco
        ) ~ aporte_naf_miteco,
      TRUE ~ NA_real_
    )
  )

readr::write_csv(
  miteco_estaciones,
  file.path(
    ruta_tablas,
    "01_MITECO_aporte_NAF_diario_estaciones.csv"
  ),
  na = ""
)


# 5. TRAZABILIDAD GEMAS + MITECO

traza <- readr::read_csv(
  ruta_traza_pm10,
  show_col_types = FALSE
) %>%
  mutate(
    fecha = as.Date(
      fecha
    ),
    fecha_muestra = as.Date(
      fecha_muestra
    ),
    sdi_pm10_ge50 = a_logico(
      sdi
    ),
    estacion_miteco = case_when(
      isla == "Tenerife" ~ "El_Rio",
      isla == "La Palma" ~ "Echedo_Valverde",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(
    cal %>%
      select(
        fecha,
        miteco_africano,
        episodio_miteco_id,
        estacion_4,
        regimen_estacional
      ),
    by = "fecha"
  ) %>%
  left_join(
    miteco_estaciones %>%
      select(
        fecha,
        estacion_miteco,
        pm10_miteco,
        perc40_miteco,
        aporte_naf_miteco,
        aporte_naf_usable
      ),
    by = c(
      "fecha",
      "estacion_miteco"
    )
  ) %>%
  mutate(
    miteco_africano = as.integer(
      miteco_africano
    ),
    aporte_naf_valido = !is.na(
      aporte_naf_usable
    ),
    miteco_evento_con_aporte_valido =
      miteco_africano == 1 &
        !is.na(
          aporte_naf_miteco
        )
  )

readr::write_csv(
  traza,
  file.path(
    ruta_tablas,
    "02_trazabilidad_MITECO_PM10_D1_D180.csv"
  ),
  na = ""
)


# 6. RESUMEN POR VENTANAS

resumir_ventana <- function(
  df,
  ndias
) {

  x <- df %>%
    filter(
      lag >= 1,
      lag <= ndias
    )

  n_miteco <- sum(
    x$miteco_africano == 1,
    na.rm = TRUE
  )

  n_miteco_valid_aporte <- sum(
    x$miteco_evento_con_aporte_valido,
    na.rm = TRUE
  )

  cobertura_evento <- if (
    n_miteco == 0
  ) {
    1
  } else {
    n_miteco_valid_aporte /
      n_miteco
  }

  aporte_ok <- cobertura_evento >=
    0.80

  aporte <- x$aporte_naf_usable

  list(
    miteco_any = as.integer(
      n_miteco > 0
    ),

    miteco_days = n_miteco,

    miteco_fraction = n_miteco /
      ndias,

    miteco_episodes = contar_episodios(
      x$episodio_miteco_id[
        x$miteco_africano == 1
      ]
    ),

    miteco_recency = if (
      n_miteco > 0
    ) {
      min(
        x$lag[
          x$miteco_africano == 1
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },

    naf_event_days_valid = n_miteco_valid_aporte,

    naf_event_coverage = cobertura_evento,

    naf_mean_calendar = if (
      aporte_ok
    ) {
      mean(
        aporte,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },

    naf_sum_calendar = if (
      aporte_ok
    ) {
      sum(
        aporte,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },

    naf_max = if (
      aporte_ok &&
        any(
          !is.na(
            aporte
          )
        )
    ) {
      max(
        aporte,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },

    naf_mean_event = if (
      aporte_ok &&
        n_miteco_valid_aporte >
          0
    ) {
      mean(
        x$aporte_naf_miteco[
          x$miteco_evento_con_aporte_valido
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    }
  )
}


ids <- unique(
  traza$id_muestra
)

out <- vector(
  "list",
  length(
    ids
  )
)

for (i in seq_along(ids)) {

  id <- ids[i]

  d <- traza %>%
    filter(
      id_muestra == id
    )

  w10 <- resumir_ventana(
    d,
    10
  )

  w30 <- resumir_ventana(
    d,
    30
  )

  w90 <- resumir_ventana(
    d,
    90
  )

  w180 <- resumir_ventana(
    d,
    180
  )

  out[[i]] <- tibble(
    id_muestra = id,
    fecha_muestra = first(
      d$fecha_muestra
    ),
    isla = first(
      d$isla
    ),
    estacion_miteco_proxy = first(
      d$estacion_miteco
    ),

    miteco_any_10d = w10$miteco_any,
    miteco_days_10d = w10$miteco_days,
    miteco_fraction_10d = w10$miteco_fraction,
    miteco_episodes_10d = w10$miteco_episodes,
    miteco_recency_10d = w10$miteco_recency,
    naf_event_coverage_10d = w10$naf_event_coverage,
    naf_mean_calendar_10d = w10$naf_mean_calendar,
    naf_sum_calendar_10d = w10$naf_sum_calendar,
    naf_max_10d = w10$naf_max,
    naf_mean_event_10d = w10$naf_mean_event,

    miteco_any_30d = w30$miteco_any,
    miteco_days_30d = w30$miteco_days,
    miteco_episodes_30d = w30$miteco_episodes,
    miteco_recency_30d = w30$miteco_recency,
    naf_event_coverage_30d = w30$naf_event_coverage,
    naf_mean_calendar_30d = w30$naf_mean_calendar,
    naf_sum_calendar_30d = w30$naf_sum_calendar,

    miteco_any_90d = w90$miteco_any,
    miteco_days_90d = w90$miteco_days,
    miteco_episodes_90d = w90$miteco_episodes,
    naf_event_coverage_90d = w90$naf_event_coverage,
    naf_mean_calendar_90d = w90$naf_mean_calendar,
    naf_sum_calendar_90d = w90$naf_sum_calendar,

    miteco_any_180d = w180$miteco_any,
    miteco_days_180d = w180$miteco_days,
    miteco_fraction_180d = w180$miteco_fraction,
    miteco_episodes_180d = w180$miteco_episodes,
    naf_event_coverage_180d = w180$naf_event_coverage,
    naf_mean_calendar_180d = w180$naf_mean_calendar,
    naf_sum_calendar_180d = w180$naf_sum_calendar
  )
}

exp_miteco <- bind_rows(
  out
)

# 7. COMPARAR CON DEFINICION OPERATIVA PM10 >=50
old <- readr::read_csv(
  ruta_caract_pm10,
  show_col_types = FALSE
) %>%
  mutate(
    fecha_muestra = as.Date(
      fecha_muestra
    )
  ) %>%
  select(
    id_muestra,
    sdi_exposicion_10d_episodios,
    n_dias_sdi_10d_episodios,
    pm10_media_10d_episodios,
    pm10_max_10d_episodios,
    auc50_10d_episodios,
    n_dias_sdi_180d_episodios,
    pm10_media_180d_episodios,
    auc50_180d_episodios
  )

master <- exp_miteco %>%
  left_join(
    old,
    by = "id_muestra"
  ) %>%
  mutate(
    categoria_MITECO_PM10_10d = case_when(
      miteco_any_10d == 0 &
        sdi_exposicion_10d_episodios == 0 ~
        "Ninguno",

      miteco_any_10d == 1 &
        sdi_exposicion_10d_episodios == 0 ~
        "MITECO_sin_PM10ge50",

      miteco_any_10d == 1 &
        sdi_exposicion_10d_episodios == 1 ~
        "MITECO_y_PM10ge50",

      miteco_any_10d == 0 &
        sdi_exposicion_10d_episodios == 1 ~
        "PM10ge50_sin_MITECO",

      TRUE ~ NA_character_
    ),

    mes_muestra = as.integer(
      format(
        fecha_muestra,
        "%m"
      )
    ),

    regimen_muestra = case_when(
      mes_muestra %in% c(
        11,
        12,
        1
      ) ~
        "Winter_proxy_NDJ",

      mes_muestra %in% c(
        7,
        8
      ) ~
        "Summer_proxy_JA",

      mes_muestra %in% 2:6 ~
        "Feb-Jun",

      mes_muestra %in% c(
        9,
        10
      ) ~
        "Sep-Oct",

      TRUE ~ NA_character_
    )
  )


# 8. QC COMPARATIVO

qc_muestra <- master %>%
  count(
    miteco_any_10d,
    sdi_exposicion_10d_episodios,
    name = "n"
  )

qc_categoria <- master %>%
  count(
    categoria_MITECO_PM10_10d,
    name = "n"
  )

qc_regimen <- master %>%
  count(
    regimen_muestra,
    categoria_MITECO_PM10_10d,
    name = "n"
  )

qc_isla <- master %>%
  count(
    isla,
    categoria_MITECO_PM10_10d,
    name = "n"
  )

qc_cobertura <- master %>%
  summarise(
    n = n(),
    n_naf10_disponible = sum(
      !is.na(
        naf_mean_calendar_10d
      )
    ),
    n_naf180_disponible = sum(
      !is.na(
        naf_mean_calendar_180d
      )
    ),
    mediana_cobertura_evento10 = median(
      naf_event_coverage_10d,
      na.rm = TRUE
    ),
    mediana_cobertura_evento180 = median(
      naf_event_coverage_180d,
      na.rm = TRUE
    )
  )


# 9. EXPORTAR


readr::write_csv(
  qc_anual,
  file.path(
    ruta_tablas,
    "03_QC_calendario_MITECO_anual.csv"
  ),
  na = ""
)

readr::write_csv(
  master,
  file.path(
    ruta_tablas,
    "04_exposicion_MITECO_PM10_por_muestra.csv"
  ),
  na = ""
)

readr::write_csv(
  qc_muestra,
  file.path(
    ruta_tablas,
    "05_QC_MITECO10_vs_PM10ge50.csv"
  ),
  na = ""
)

readr::write_csv(
  qc_categoria,
  file.path(
    ruta_tablas,
    "06_QC_categorias_combinadas.csv"
  ),
  na = ""
)

readr::write_csv(
  qc_regimen,
  file.path(
    ruta_tablas,
    "07_QC_categorias_por_regimen.csv"
  ),
  na = ""
)

readr::write_csv(
  qc_isla,
  file.path(
    ruta_tablas,
    "08_QC_categorias_por_isla.csv"
  ),
  na = ""
)

readr::write_csv(
  qc_cobertura,
  file.path(
    ruta_tablas,
    "09_QC_cobertura_aporte_NAF.csv"
  ),
  na = ""
)

notas <- tibble(
  aspecto = c(
    "MITECO_any",
    "PM10_ge50",
    "El_Rio",
    "Echedo",
    "Aporte_NAF",
    "Cobertura",
    "Interpretacion"
  ),
  detalle = c(
    "Presencia oficial regional de episodio africano en Canarias.",
    "Definicion operativa previa: PM10 local medio diario >=50 ug/m3.",
    "Estacion MITECO de fondo usada para Tenerife.",
    "Proxy MITECO para El Hierro, La Palma y La Gomera; no es una medicion residencial individual en La Palma.",
    "Campo Descuento oficial MITECO; cuantifica contribucion africana diaria estimada a PM10.",
    "Las metricas continuas NAF se calculan si >=80% de los dias MITECO de la ventana tienen aporte disponible.",
    "MITECO y PM10 local se mantienen como dimensiones complementarias; una no sustituye automaticamente a la otra."
  )
)

writexl::write_xlsx(
  list(
    Exposicion = master,
    QC_anual = qc_anual,
    QC_MITECO_PM10 = qc_muestra,
    Categorias = qc_categoria,
    Por_regimen = qc_regimen,
    Por_isla = qc_isla,
    Cobertura_NAF = qc_cobertura,
    Notas = notas
  ),
  file.path(
    ruta_tablas,
    "10_resumen_integracion_MITECO_PM10_GEMAS.xlsx"
  )
)

saveRDS(
  list(
    calendario = cal,
    miteco_estaciones = miteco_estaciones,
    trazabilidad = traza,
    exposicion_por_muestra = master
  ),
  file.path(
    ruta_objetos,
    "objetos_integracion_MITECO_PM10.rds"
  )
)


# 10. MENSAJE FINAL

cat(
  "\n============================================\n"
)

cat(
  "INTEGRACION MITECO + PM10 FINALIZADA\n"
)

cat(
  "============================================\n\n"
)

cat(
  "QC anual calendario MITECO:\n"
)

print(
  qc_anual
)

cat(
  "\nComparacion MITECO10 vs PM10>=50:\n"
)

print(
  qc_muestra
)

cat(
  "\nCategorias combinadas:\n"
)

print(
  qc_categoria
)

cat(
  "\nPor regimen estacional:\n"
)

print(
  qc_regimen
)

cat(
  "\nCobertura aporte NAF:\n"
)

print(
  qc_cobertura
)
