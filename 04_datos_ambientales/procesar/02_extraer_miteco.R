# Preparar las series ambientales públicas.
# No utiliza identificadores ni fechas de pacientes.

paquetes <- c("dplyr", "tidyr", "readr", "tibble", "stringr", "readxl")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) stop("Faltan paquetes: ", paste(faltan, collapse = ", "))
suppressPackageStartupMessages(invisible(lapply(paquetes, library, character.only = TRUE)))
ruta_miteco <- file.path("data_raw", "CALIMA", "CALIMA MITECO 2016-2020")
ruta_calendario <- file.path(ruta_miteco, "00_MITECO_Canarias_calendario_africanos_2016_2020.csv")
ruta_tablas <- file.path("outputs", "environmental", "miteco")
dir.create(ruta_tablas, recursive = TRUE, showWarnings = FALSE)
normalizar_texto <- function(x) {

  y <- as.character(x)

  y <- iconv(y, from = "", to = "ASCII//TRANSLIT")

  y <- tolower(trimws(y))

  y <- gsub("[^a-z0-9]+", " ", y)

  trimws(y)
}

a_num <- function(x) {

  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  y <- trimws(as.character(x))

  y[
    y == ""
  ] <- NA_character_

  y <- gsub(",", ".", y, fixed = TRUE)

  suppressWarnings(as.numeric(y))
}

a_logico <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(trimws(as.character(x)))

  out <- rep(NA, length(y))

  out[
    y %in% c("TRUE", "T", "1", "SI", "SÍ", "YES")
  ] <- TRUE

  out[
    y %in% c("FALSE", "F", "0", "NO")
  ] <- FALSE

  as.logical(out)
}

parse_fecha_excel <- function(x) {

  if (inherits(x, "Date")) {
    return(as.Date(x))
  }

  if (
    inherits(x, "POSIXct") ||
      inherits(x, "POSIXt")
  ) {
    return(as.Date(x))
  }

# Intentar fecha textual ISO / habitual.
  y <- trimws(as.character(x))

  out <- suppressWarnings(as.Date(y))

# Excel serial si no fue posible.
  idx <- is.na(out)

  if (any(idx)) {

    num <- suppressWarnings(as.numeric(y[ idx ]))

    out[
      idx
    ] <- as.Date(num, origin = "1899-12-30")
  }

  out
}

encontrar_bloque_metrica <- function(
  metricas,
  tipo = c("pm10", "perc40", "descuento")
) {

  tipo <- match.arg(tipo)

  if (
    tipo == "pm10"
  ) {
    idx <- which(grepl("^pm10$", metricas))
  }

  if (
    tipo == "perc40"
  ) {
    idx <- which(grepl("^perc ?40$", metricas) | grepl("percentil.*40", metricas))
  }

  if (
    tipo == "descuento"
  ) {
    idx <- which(grepl("descuento", metricas))
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

  pos <- which(grepl(patron_estacion, estaciones_pm10))

  if (
    length(pos) != 1
  ) {
    stop(
      "No se encontro una posicion unica para la estacion '",
      patron_estacion,
      "' dentro del bloque PM10. Posiciones: ",
      paste(pos, collapse = ", "),
      "\nEstaciones PM10 detectadas: ",
      paste(estaciones_pm10, collapse = " | ")
    )
  }

  pos
}

leer_canarias_miteco <- function(
  archivo
) {

  cat("Leyendo: ", basename(archivo), "\n", sep = "")

  raw <- readxl::read_excel(archivo, sheet = "CANARIAS", col_names = FALSE)

  if (nrow(raw) < 10) {
    stop("Hoja CANARIAS inesperadamente corta en ", basename(archivo))
  }

  estaciones <- normalizar_texto(unlist(raw[ 1, , drop = TRUE ]))

  metricas <- normalizar_texto(unlist(raw[ 2, , drop = TRUE ]))

# Las hojas antiguas 2016-2018 no siempre repiten el nombre
# de la estación sobre los bloques PERC40/Descuento.
# Por eso NO buscamos "El Rio + PERC40" en la misma columna.
# Detectamos primero los tres bloques de métricas y usamos la
# posición de cada estación dentro del bloque PM10 como mapa
# para PERC40 y Descuento.

  idx_pm10 <- encontrar_bloque_metrica(metricas, "pm10")

  idx_perc <- encontrar_bloque_metrica(metricas, "perc40")

  idx_desc <- encontrar_bloque_metrica(metricas, "descuento")

  if (
    length(idx_pm10) == 0
  ) {
    stop(
      "No se detecto el bloque PM10 en ",
      basename(archivo),
      ". Metricas leidas: ",
      paste(metricas, collapse = " | ")
    )
  }

# PERC40 es secundario para nuestro análisis. Si un año no lo
# contiene de forma reconocible, se conserva como NA sin abortar.
  tiene_perc <- length(idx_perc) > 0

  if (
    length(idx_desc) == 0
  ) {
    stop(
      "No se detecto el bloque Descuento en ",
      basename(archivo),
      ". Metricas leidas: ",
      paste(metricas, collapse = " | ")
    )
  }

  if (
    length(idx_desc) !=
      length(idx_pm10)
  ) {
    stop(
      "Los bloques PM10 y Descuento no tienen la misma longitud en ",
      basename(archivo),
      ": PM10=",
      length(idx_pm10),
      "; Descuento=",
      length(idx_desc)
    )
  }

  if (
    tiene_perc &&
      length(idx_perc) !=
        length(idx_pm10)
  ) {
    warning(
      "Bloque PERC40 con longitud distinta en ",
      basename(archivo),
      ". Se ignorara PERC40 para ese año."
    )

    tiene_perc <- FALSE
  }

  pos_elrio <- localizar_estacion_en_bloque_pm10(estaciones, idx_pm10, "el rio")

  pos_echedo <- localizar_estacion_en_bloque_pm10(estaciones, idx_pm10, "echedo")

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
    length(idx_pm10),
    "; PERC40=",
    length(idx_perc),
    "; Descuento=",
    length(idx_desc),
    "
",
    sep = ""
  )

  cat(
    "  El Rio columnas: PM10=",
    i_elrio_pm10,
    "; PERC40=",
    ifelse(is.na(i_elrio_perc), "NA", i_elrio_perc),
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
    ifelse(is.na(i_echedo_perc), "NA", i_echedo_perc),
    "; Descuento=",
    i_echedo_desc,
    "
",
    sep = ""
  )

  fecha <- parse_fecha_excel(raw[[1]][ -c(1, 2) ])

  extraer_estacion <- function(
    nombre,
    idx_pm10,
    idx_perc,
    idx_desc
  ) {

    pm10 <- a_num(raw[[idx_pm10]][ -c(1, 2) ])

    perc40 <- if (
      is.na(idx_perc)
    ) {
      rep(NA_real_, length(fecha))
    } else {
      a_num(raw[[idx_perc]][ -c(1, 2) ])
    }

    aporte <- a_num(raw[[idx_desc]][ -c(1, 2) ])

# QC físico:
# el aporte natural no puede exceder el PM10 observado.
    aporte[
      !is.na(aporte) &
        (is.na(pm10) | aporte < -0.001 | aporte > pm10 + 0.5)
    ] <- NA_real_

    tibble(
      fecha = fecha,
      estacion_miteco = nombre,
      pm10_miteco = pm10,
      perc40_miteco = perc40,
      aporte_naf_miteco = aporte
    ) %>%
      filter(!is.na(fecha))
  }

  bind_rows(
    extraer_estacion("El_Rio", i_elrio_pm10, i_elrio_perc, i_elrio_desc),
    extraer_estacion("Echedo_Valverde", i_echedo_pm10, i_echedo_perc, i_echedo_desc)
  )
}

contar_episodios <- function(
  ids
) {

  ids <- ids[
    !is.na(ids)
  ]

  if (length(ids) == 0) {
    return(0L)
  }

  length(unique(ids))
}

# Calendario oficial MITECO.

cal <- readr::read_csv(ruta_calendario, show_col_types = FALSE) %>%
  mutate(fecha = as.Date(fecha), miteco_africano = as.integer(miteco_africano))

if (
  min(cal$fecha) != as.Date("2016-01-01") ||
    max(cal$fecha) != as.Date("2020-12-31")
) {
  stop("El calendario MITECO no cubre exactamente 2016-2020.")
}

qc_anual <- cal %>%
  group_by(anio) %>%
  summarise(
    dias_miteco = sum(miteco_africano, na.rm = TRUE),
    dias_anio = n(),
    porcentaje = 100 *
      dias_miteco /
      dias_anio,
    .groups = "drop"
  )

esperado <- tibble(
  anio = 2016:2020,
  dias_esperados = c(154, 201, 129, 149, 171),
  pct_reportado = c(42, 55, 35, 41, 47)
)

qc_anual <- qc_anual %>%
  left_join(esperado, by = "anio") %>%
  mutate(coincide_dias = dias_miteco == dias_esperados)

if (!all(qc_anual$coincide_dias)) {
  stop("El calendario MITECO no reproduce los dias anuales esperados.")
}

# Aportes naf de las hojas MITECO.

archivos_descuentos <- list.files(
  ruta_miteco,
  pattern = "^pm10-descuentos_todas_estaciones_20(16|17|18|19|20)_.*\\.xls[x]?$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(archivos_descuentos) != 5) {
  stop(
    "Se esperaban 5 ficheros pm10-descuentos 2016-2020 y se encontraron ",
    length(archivos_descuentos)
  )
}

miteco_estaciones <- bind_rows(lapply(sort(archivos_descuentos), leer_canarias_miteco)) %>%
  left_join(cal %>% select(fecha, miteco_africano, episodio_miteco_id), by = "fecha") %>%
  mutate(
# Fuera de episodio oficial el aporte NAF se conoce como 0.
# Dentro de episodio, si la estacion carece de dato se mantiene NA.
    aporte_naf_usable = case_when(
      miteco_africano == 0 ~ 0,
      miteco_africano == 1 &
        !is.na(aporte_naf_miteco) ~ aporte_naf_miteco,
      TRUE ~ NA_real_
    )
  )

readr::write_csv(
  miteco_estaciones,
  file.path(ruta_tablas, "01_MITECO_aporte_NAF_diario_estaciones.csv"),
  na = ""
)

readr::write_csv(qc_anual, file.path(ruta_tablas, "QC_calendario_anual.csv"))
message("Series MITECO por estacion-fecha exportadas. No se han utilizado pacientes.")
