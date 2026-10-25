# RECONSTRUIR_EPISODIOS_PM10_GEMAS.R
# TFM GEMAS
#
# OBJETIVO
# Explotar la serie diaria completa de PM10 para construir,
# para cada muestra de saliva, caracteristicas de la exposicion que no
# quedan resumidas por el simple SDI10 si/no:
#
# - recencia del ultimo dia compatible con SDI;
# - duracion del episodio mas reciente;
# - intensidad maxima y media del episodio;
# - AUC50 del episodio;
# - numero de episodios en 10/30/90/180 dias;
# - racha maxima de dias consecutivos PM10 >= 50;
# - carga PM10/AUC50 en 10/30/90/180 dias;
# - una clasificacion descriptiva del patron reciente:
#     Sin_SDI_10d
#     Aguda_reciente
#     Sostenida_reciente
#     Recurrente_10d
#     Remota_D4_D10


# 0. PAQUETES

paquetes <- c(
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
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
    "\nInstalalos antes de ejecutar el script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(writexl)
})


# 1. PARAMETROS Y RUTAS

UMBRAL_SDI <- 50
VENTANA_MAX <- 180
MIN_VALIDOS_10 <- 7
MIN_VALIDOS_180 <- 153

ruta_cohorte <- file.path(
  "outputs",
  "06_cohorte_elegible_cruce_ambiental.csv"
)

ruta_aire <- file.path(
  "data_raw",
  "CALIMA",
  "auditoria_estaciones",
  "03_datos_diarios_normalizados.csv"
)

ruta_macro <- file.path(
  "outputs",
  "17_macrotabla_clinica_ambiental_por_muestra.csv"
)

ruta_salida <- file.path(
  "outputs",
  "analisis_microbiota",
  "patrones_exposicion_pm10"
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

for (f in c(ruta_cohorte, ruta_aire)) {
  if (!file.exists(f)) {
    stop("No existe: ", f)
  }
}

# 2. FUNCIONES

a_logico <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA,
    length(y)
  )

  out[
    y %in% c(
      "TRUE",
      "T",
      "1",
      "SI",
      "SÍ"
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

  as.logical(out)
}


a_num <- function(x) {

  if (is.numeric(x)) {
    return(x)
  }

  y <- trimws(
    as.character(x)
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
    as.numeric(y)
  )
}


normalizar_isla <- function(x) {

  y <- trimws(
    as.character(x)
  )

  y[
    y %in% c(
      "La_Palma",
      "La Palma"
    )
  ] <- "La Palma"

  y
}


suma_segura <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  sum(x, na.rm = TRUE)
}


media_segura <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}


max_seguro <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  max(x, na.rm = TRUE)
}


resumir_episodios <- function(
  d,
  fecha_muestra
) {

  # d debe contener fecha, pm10 y sdi, ordenado de antiguo -> reciente.
  d <- d %>%
    arrange(fecha)

  is_true <- !is.na(d$sdi) &
    d$sdi

  episodio_id <- rep(
    NA_integer_,
    nrow(d)
  )

  actual <- 0L
  anterior_true <- FALSE
  fecha_anterior <- as.Date(NA)

  for (i in seq_len(nrow(d))) {

    if (is_true[i]) {

      consecutivo <- anterior_true &&
        !is.na(fecha_anterior) &&
        d$fecha[i] ==
          fecha_anterior + 1

      if (!consecutivo) {
        actual <- actual + 1L
      }

      episodio_id[i] <- actual
      anterior_true <- TRUE
      fecha_anterior <- d$fecha[i]

    } else {

      anterior_true <- FALSE
      fecha_anterior <- d$fecha[i]
    }
  }

  if (all(is.na(episodio_id))) {

    return(
      list(
        episodios = tibble(),
        ultimo = NULL
      )
    )
  }

  dd <- d %>%
    mutate(
      episodio_id = episodio_id
    ) %>%
    filter(
      !is.na(episodio_id)
    )

  episodios <- dd %>%
    group_by(
      episodio_id
    ) %>%
    summarise(
      fecha_inicio = min(fecha),
      fecha_fin = max(fecha),
      duracion_dias = n(),
      pm10_media = mean(
        pm10,
        na.rm = TRUE
      ),
      pm10_max = max(
        pm10,
        na.rm = TRUE
      ),
      auc50 = sum(
        pmax(
          pm10 - UMBRAL_SDI,
          0
        ),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      lag_fin = as.integer(
        fecha_muestra -
          fecha_fin
      )
    ) %>%
    arrange(
      fecha_fin
    )

  ultimo <- episodios %>%
    slice_tail(
      n = 1
    )

  list(
    episodios = episodios,
    ultimo = ultimo
  )
}


resumen_ventana <- function(
  d,
  ndias
) {

  x <- d %>%
    filter(
      lag >= 1,
      lag <= ndias
    )

  validos <- sum(
    !is.na(
      x$pm10
    )
  )

  n_sdi <- sum(
    x$sdi %in% TRUE,
    na.rm = TRUE
  )

  list(
    n_validos = validos,
    pm10_media = media_segura(
      x$pm10
    ),
    pm10_max = max_seguro(
      x$pm10
    ),
    auc50 = suma_segura(
      pmax(
        x$pm10 - UMBRAL_SDI,
        0
      )
    ),
    n_sdi = n_sdi
  )
}


n_episodios_ventana <- function(
  episodios,
  ndias
) {

  if (nrow(episodios) == 0) {
    return(0L)
  }

  sum(
    episodios$lag_fin <= ndias
  )
}


racha_max_ventana <- function(
  d,
  ndias
) {

  x <- d %>%
    filter(
      lag >= 1,
      lag <= ndias
    ) %>%
    arrange(fecha)

  if (nrow(x) == 0) {
    return(NA_integer_)
  }

  y <- x$sdi %in% TRUE

  rr <- rle(y)

  if (!any(rr$values)) {
    return(0L)
  }

  max(
    rr$lengths[
      rr$values
    ]
  )
}


clasificar_patron_10d <- function(
  n_validos_10,
  n_sdi_10,
  n_episodios_10,
  racha_max_10,
  lag_ultimo_dia_sdi
) {

  if (
    is.na(n_validos_10) ||
      n_validos_10 <
        MIN_VALIDOS_10
  ) {
    return(NA_character_)
  }

  if (
    is.na(n_sdi_10) ||
      n_sdi_10 == 0
  ) {
    return("Sin_SDI_10d")
  }

  if (
    !is.na(n_episodios_10) &&
      n_episodios_10 >= 2
  ) {
    return("Recurrente_10d")
  }

  if (
    !is.na(lag_ultimo_dia_sdi) &&
      lag_ultimo_dia_sdi >= 4 &&
      lag_ultimo_dia_sdi <= 10
  ) {
    return("Remota_D4_D10")
  }

  if (
    !is.na(lag_ultimo_dia_sdi) &&
      lag_ultimo_dia_sdi <= 3 &&
      !is.na(racha_max_10) &&
      racha_max_10 >= 3
  ) {
    return("Sostenida_reciente")
  }

  if (
    !is.na(lag_ultimo_dia_sdi) &&
      lag_ultimo_dia_sdi <= 3
  ) {
    return("Aguda_reciente")
  }

  "Otro_patron"
}

# 3. CARGAR COHORTE Y SERIE DIARIA

cohorte <- readr::read_csv(
  ruta_cohorte,
  show_col_types = FALSE
) %>%
  mutate(
    id_muestra = as.character(
      id_muestra
    ),
    fecha_muestra = as.Date(
      fecha_muestra
    ),
    isla = normalizar_isla(
      isla
    )
  ) %>%
  filter(
    !is.na(
      fecha_muestra
    ),
    !is.na(
      isla
    )
  )

if (nrow(cohorte) != 205) {
  warning(
    "La cohorte ambiental tiene ",
    nrow(cohorte),
    " muestras; se esperaban 205."
  )
}

if (anyDuplicated(cohorte$id_muestra) > 0) {
  stop(
    "Hay id_muestra duplicados."
  )
}


aire <- readr::read_csv2(
  ruta_aire,
  show_col_types = FALSE
)

necesarias <- c(
  "isla",
  "estacion",
  "fecha",
  "n_horas_pm10",
  "pm10_media_diaria",
  "pm10_max_horario",
  "pm10_dia_valido"
)

faltan_aire <- setdiff(
  necesarias,
  names(aire)
)

if (length(faltan_aire) > 0) {
  stop(
    "Faltan columnas en serie diaria: ",
    paste(
      faltan_aire,
      collapse = ", "
    )
  )
}

aire <- aire %>%
  mutate(
    isla = normalizar_isla(
      isla
    ),
    estacion = trimws(
      as.character(
        estacion
      )
    ),
    fecha = as.Date(
      fecha
    ),
    pm10_dia_valido = a_logico(
      pm10_dia_valido
    ),
    pm10_media_diaria = a_num(
      pm10_media_diaria
    ),
    pm10_max_horario = a_num(
      pm10_max_horario
    ),
    n_horas_pm10 = a_num(
      n_horas_pm10
    )
  )


# 4. JERARQUIA PM10 

jerarquia <- tibble::tribble(
  ~isla,      ~estacion,                       ~prioridad, ~nivel,
  "Tenerife", "El Río",                        1L,         "primaria",
  "Tenerife", "Granadilla",                    2L,         "secundaria",
  "La Palma", "El Pilar-Sta Cruz de La Palma", 1L,         "primaria",
  "La Palma", "La Grama-Breña Alta",           2L,         "secundaria"
)

serie_pm <- aire %>%
  inner_join(
    jerarquia,
    by = c(
      "isla",
      "estacion"
    )
  ) %>%
  filter(
    pm10_dia_valido %in% TRUE,
    !is.na(
      pm10_media_diaria
    )
  ) %>%
  arrange(
    isla,
    fecha,
    prioridad
  ) %>%
  group_by(
    isla,
    fecha
  ) %>%
  slice(
    1
  ) %>%
  ungroup() %>%
  transmute(
    isla,
    fecha,
    estacion_origen = estacion,
    nivel_estacion = nivel,
    prioridad_estacion = prioridad,
    n_horas_pm10,
    pm10 = pm10_media_diaria,
    pm10_max_horario,
    sdi = pm10 >=
      UMBRAL_SDI
  )


# 5. CALENDARIO COMPLETO PARA NO CONFUNDIR MISSING CON NO-SDI

fecha_min <- min(
  cohorte$fecha_muestra,
  na.rm = TRUE
) - VENTANA_MAX

fecha_max <- max(
  cohorte$fecha_muestra,
  na.rm = TRUE
) - 1

calendario <- tidyr::crossing(
  isla = c(
    "La Palma",
    "Tenerife"
  ),
  fecha = seq(
    fecha_min,
    fecha_max,
    by = "day"
  )
) %>%
  left_join(
    serie_pm,
    by = c(
      "isla",
      "fecha"
    )
  ) %>%
  arrange(
    isla,
    fecha
  )

saveRDS(
  calendario,
  file.path(
    ruta_objetos,
    "serie_PM10_diaria_jerarquizada_completa.rds"
  )
)


# 6. TRAZABILIDAD 180 DIAS POR MUESTRA

trazabilidad <- cohorte %>%
  select(
    id_muestra,
    fecha_muestra,
    isla
  ) %>%
  tidyr::crossing(
    lag = 1:VENTANA_MAX
  ) %>%
  mutate(
    fecha = fecha_muestra -
      lag
  ) %>%
  left_join(
    calendario,
    by = c(
      "isla",
      "fecha"
    )
  ) %>%
  arrange(
    id_muestra,
    lag
  )

readr::write_csv(
  trazabilidad,
  file.path(
    ruta_tablas,
    "01_trazabilidad_PM10_D1_D180.csv"
  ),
  na = ""
)


# 7. CARACTERISTICAS POR MUESTRA

ids <- cohorte$id_muestra

res <- vector(
  "list",
  length(ids)
)

for (i in seq_along(ids)) {

  id <- ids[i]

  fila <- cohorte %>%
    filter(
      id_muestra == id
    )

  d <- trazabilidad %>%
    filter(
      id_muestra == id
    ) %>%
    arrange(
      fecha
    )

  epi <- resumir_episodios(
    d = d,
    fecha_muestra = fila$fecha_muestra[1]
  )

  w10 <- resumen_ventana(
    d,
    10
  )

  w30 <- resumen_ventana(
    d,
    30
  )

  w90 <- resumen_ventana(
    d,
    90
  )

  w180 <- resumen_ventana(
    d,
    180
  )

  lag_ultimo <- if (
    nrow(
      epi$episodios
    ) > 0
  ) {
    min(
      epi$episodios$lag_fin
    )
  } else {
    NA_integer_
  }

  ultimo <- epi$ultimo

  # Carga con ponderacion exponencial por recencia.
  # Semivida operacional = 3 dias.
  d10 <- d %>%
    filter(
      lag >= 1,
      lag <= 10
    )

  carga_ponderada_10d_hl3 <- if (
    all(
      is.na(
        d10$pm10
      )
    )
  ) {
    NA_real_
  } else {
    sum(
      pmax(
        d10$pm10 - UMBRAL_SDI,
        0
      ) *
        2^(
          -d10$lag / 3
        ),
      na.rm = TRUE
    )
  }

  racha10 <- racha_max_ventana(
    d,
    10
  )

  racha30 <- racha_max_ventana(
    d,
    30
  )

  racha180 <- racha_max_ventana(
    d,
    180
  )

  nepi10 <- n_episodios_ventana(
    epi$episodios,
    10
  )

  nepi30 <- n_episodios_ventana(
    epi$episodios,
    30
  )

  nepi90 <- n_episodios_ventana(
    epi$episodios,
    90
  )

  nepi180 <- n_episodios_ventana(
    epi$episodios,
    180
  )

  dias10 <- d %>%
    filter(
      lag >= 1,
      lag <= 10
    ) %>%
    select(
      lag,
      pm10
    ) %>%
    mutate(
      nombre = sprintf(
        "pm10_d%02d",
        lag
      )
    ) %>%
    select(
      nombre,
      pm10
    ) %>%
    pivot_wider(
      names_from = nombre,
      values_from = pm10
    )

  base <- tibble(
    id_muestra = id,
    fecha_muestra = fila$fecha_muestra[1],
    isla = fila$isla[1],

    n_validos_10d = w10$n_validos,
    n_validos_30d = w30$n_validos,
    n_validos_90d = w90$n_validos,
    n_validos_180d = w180$n_validos,

    pm10_media_10d_episodios = w10$pm10_media,
    pm10_max_10d_episodios = w10$pm10_max,
    auc50_10d_episodios = w10$auc50,
    n_dias_sdi_10d_episodios = w10$n_sdi,

    pm10_media_30d = w30$pm10_media,
    pm10_max_30d = w30$pm10_max,
    auc50_30d = w30$auc50,
    n_dias_sdi_30d = w30$n_sdi,

    pm10_media_90d = w90$pm10_media,
    pm10_max_90d = w90$pm10_max,
    auc50_90d = w90$auc50,
    n_dias_sdi_90d = w90$n_sdi,

    pm10_media_180d_episodios = w180$pm10_media,
    pm10_max_180d_episodios = w180$pm10_max,
    auc50_180d_episodios = w180$auc50,
    n_dias_sdi_180d_episodios = w180$n_sdi,

    sdi_exposicion_10d_episodios = ifelse(
      w10$n_validos >= MIN_VALIDOS_10,
      as.integer(
        w10$n_sdi > 0
      ),
      NA_integer_
    ),

    n_episodios_sdi_10d = nepi10,
    n_episodios_sdi_30d = nepi30,
    n_episodios_sdi_90d = nepi90,
    n_episodios_sdi_180d = nepi180,

    racha_max_sdi_10d = racha10,
    racha_max_sdi_30d = racha30,
    racha_max_sdi_180d = racha180,

    lag_ultimo_dia_sdi = lag_ultimo,

    ultimo_episodio_inicio = if (
      !is.null(ultimo)
    ) {
      ultimo$fecha_inicio[1]
    } else {
      as.Date(NA)
    },

    ultimo_episodio_fin = if (
      !is.null(ultimo)
    ) {
      ultimo$fecha_fin[1]
    } else {
      as.Date(NA)
    },

    duracion_ultimo_episodio = if (
      !is.null(ultimo)
    ) {
      ultimo$duracion_dias[1]
    } else {
      NA_integer_
    },

    pm10_media_ultimo_episodio = if (
      !is.null(ultimo)
    ) {
      ultimo$pm10_media[1]
    } else {
      NA_real_
    },

    pm10_max_ultimo_episodio = if (
      !is.null(ultimo)
    ) {
      ultimo$pm10_max[1]
    } else {
      NA_real_
    },

    auc50_ultimo_episodio = if (
      !is.null(ultimo)
    ) {
      ultimo$auc50[1]
    } else {
      NA_real_
    },

    carga_ponderada_10d_hl3 = carga_ponderada_10d_hl3
  )

  base$fenotipo_exposicion_10d <- clasificar_patron_10d(
    n_validos_10 = base$n_validos_10d,
    n_sdi_10 = base$n_dias_sdi_10d_episodios,
    n_episodios_10 = base$n_episodios_sdi_10d,
    racha_max_10 = base$racha_max_sdi_10d,
    lag_ultimo_dia_sdi = base$lag_ultimo_dia_sdi
  )

  res[[i]] <- bind_cols(
    base,
    dias10
  )
}

caracteristicas <- bind_rows(
  res
)


# 8. QC FRENTE A MACROTABLA 

qc_macro <- tibble()

if (file.exists(ruta_macro)) {

  macro <- readr::read_csv2(
    ruta_macro,
    show_col_types = FALSE
  )

  comp <- caracteristicas %>%
    select(
      id_muestra,
      sdi_exposicion_10d_episodios,
      n_dias_sdi_10d_episodios,
      auc50_10d_episodios
    ) %>%
    left_join(
      macro %>%
        select(
          id_muestra,
          sdi_exposicion_10d,
          n_dias_sdi_10d,
          pm10_auc50_10d
        ),
      by = "id_muestra"
    ) %>%
    mutate(
      coincide_sdi10 = sdi_exposicion_10d_episodios ==
        as.integer(
          sdi_exposicion_10d
        ),
      coincide_n_dias = n_dias_sdi_10d_episodios ==
        n_dias_sdi_10d,
      coincide_auc50 = abs(
        auc50_10d_episodios -
          pm10_auc50_10d
      ) <
        1e-8
    )

  qc_macro <- tibble(
    comprobacion = c(
      "sdi_exposicion_10d",
      "n_dias_sdi_10d",
      "pm10_auc50_10d"
    ),
    n_discrepancias = c(
      sum(
        !comp$coincide_sdi10,
        na.rm = TRUE
      ),
      sum(
        !comp$coincide_n_dias,
        na.rm = TRUE
      ),
      sum(
        !comp$coincide_auc50,
        na.rm = TRUE
      )
    )
  )

  readr::write_csv(
    comp,
    file.path(
      ruta_tablas,
      "02_QC_comparacion_macrotabla.csv"
    ),
    na = ""
  )
}


# 9. RESUMEN DESCRIPTIVO

resumen_fenotipos <- caracteristicas %>%
  count(
    isla,
    fenotipo_exposicion_10d,
    name = "n"
  ) %>%
  group_by(
    isla
  ) %>%
  mutate(
    porcentaje = 100 *
      n /
      sum(n)
  ) %>%
  ungroup()

resumen_global <- caracteristicas %>%
  summarise(
    n = n(),
    n_sdi10 = sum(
      sdi_exposicion_10d_episodios == 1,
      na.rm = TRUE
    ),
    n_no_sdi10 = sum(
      sdi_exposicion_10d_episodios == 0,
      na.rm = TRUE
    ),
    mediana_lag_ultimo_sdi = median(
      lag_ultimo_dia_sdi,
      na.rm = TRUE
    ),
    mediana_duracion_ultimo_episodio = median(
      duracion_ultimo_episodio,
      na.rm = TRUE
    ),
    mediana_auc_ultimo_episodio = median(
      auc50_ultimo_episodio,
      na.rm = TRUE
    ),
    max_racha_180d = max(
      racha_max_sdi_180d,
      na.rm = TRUE
    )
  )


# 10. EXPORTAR


readr::write_csv(
  caracteristicas,
  file.path(
    ruta_tablas,
    "03_caracteristicas_episodios_PM10_por_muestra.csv"
  ),
  na = ""
)

readr::write_csv(
  resumen_fenotipos,
  file.path(
    ruta_tablas,
    "04_resumen_fenotipos_exposicion_10d.csv"
  ),
  na = ""
)

writexl::write_xlsx(
  list(
    Caracteristicas = caracteristicas,
    Fenotipos = resumen_fenotipos,
    Resumen = resumen_global,
    QC_macrotabla = qc_macro
  ),
  file.path(
    ruta_tablas,
    "05_resumen_patrones_exposicion_PM10.xlsx"
  )
)

saveRDS(
  caracteristicas,
  file.path(
    ruta_objetos,
    "caracteristicas_episodios_PM10_por_muestra.rds"
  )
)

# 11. MENSAJE FINAL

cat(
  "\n============================================\n"
)

cat(
  "EPISODIOS PM10 RECONSTRUIDOS\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Muestras:",
  nrow(
    caracteristicas
  ),
  "\n\n"
)

print(
  resumen_global
)

cat(
  "\nFenotipos 10d:\n"
)

print(
  caracteristicas %>%
    count(
      fenotipo_exposicion_10d
    )
)

if (nrow(qc_macro) > 0) {
  cat(
    "\nQC frente a macrotabla:\n"
  )
  print(qc_macro)
}

cat(
  "\nSalida:\n",
  file.path(
    ruta_tablas,
    "05_resumen_patrones_exposicion_PM10.xlsx"
  ),
  "\n",
  sep = ""
)
