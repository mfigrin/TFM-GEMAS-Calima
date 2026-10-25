# Preparar las series ambientales públicas.
# No utiliza identificadores ni fechas de pacientes.

paquetes <- c("dplyr", "tidyr", "readr", "readxl", "stringr", "tibble")

faltan <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltan) > 0) {
  stop(
    "Faltan paquetes: ",
    paste(faltan, collapse = ", "),
    ". Ejecuta primero scripts/utilidades/01_instalar_dependencias.R"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
})

ruta_calima <- file.path("data_raw", "CALIMA")
ruta_salida <- file.path(ruta_calima, "auditoria_estaciones")
dir.create(ruta_salida, recursive = TRUE, showWarnings = FALSE)

archivos <- file.path(ruta_calima, paste0("Datos ", 2016:2020, ".xlsx"))

faltan_archivos <- archivos[!file.exists(archivos)]

if (length(faltan_archivos) > 0) {
  stop("Faltan los Excel horarios:\n- ", paste(faltan_archivos, collapse = "\n- "))
}

normalizar_nombre <- function(x) {
  y <- iconv(trimws(as.character(x)), from = "", to = "ASCII//TRANSLIT")
  y <- tolower(y)
  y <- gsub("[^a-z0-9]+", " ", y)
  trimws(y)
}

normalizar_isla <- function(estacion) {
  e <- normalizar_nombre(estacion)

  if (grepl("^el pilar", e) || grepl("^la grama", e)) {
    return("La Palma")
  }

  "Tenerife"
}

buscar_columna <- function(nombres, patrones) {
  n <- normalizar_nombre(nombres)

  idx <- which(
    vapply(
      n,
      function(x) {
        any(vapply(patrones, function(p) grepl(p, x), logical(1)))
      },
      logical(1)
    )
  )

  if (length(idx) == 0) {
    return(NA_integer_)
  }

  idx[1]
}

parsear_fecha <- function(x, anio_esperado = NULL) {
  if (inherits(x, "Date")) {
    out <- as.Date(x)
  } else if (inherits(x, "POSIXt")) {
    out <- as.Date(x)
  } else if (is.numeric(x)) {
    out <- as.Date(x, origin = "1899-12-30")
  } else {
    y <- trimws(as.character(x))
    out <- rep(as.Date(NA), length(y))
    formatos <- c("^\\d{4}-\\d{2}-\\d{2}$" = "%Y-%m-%d",
                  "^\\d{1,2}/\\d{1,2}/\\d{4}$" = "%d/%m/%Y",
                  "^\\d{1,2}/\\d{1,2}/\\d{2}$" = "%d/%m/%y",
                  "^\\d{1,2}-\\d{1,2}-\\d{2}$" = "%d-%m-%y")
    for (patron in names(formatos)) {
      idx <- !is.na(y) & grepl(patron, y)
      out[idx] <- as.Date(y[idx], format = formatos[[patron]])
    }
    serial <- !is.na(y) & grepl("^[0-9]{4,5}([.]0+)?$", y)
    out[serial] <- as.Date(as.numeric(y[serial]), origin = "1899-12-30")
    invalido <- !is.na(y) & nzchar(y) & is.na(out)
    if (any(invalido)) stop("Hay fechas no interpretables; no se eliminan silenciosamente.")
  }
  if (!is.null(anio_esperado) && any(!is.na(out) & as.integer(format(out, "%Y")) != anio_esperado)) {
    stop("El año de una fecha no coincide con el libro de origen.")
  }
  out
}

a_num <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "N", "NA", "NaN", "-")] <- NA_character_
  y <- gsub(",", ".", y, fixed = TRUE)
  suppressWarnings(as.numeric(y))
}

leer_estacion <- function(archivo, hoja, anio) {
  raw <- readxl::read_excel(
    archivo,
    sheet = hoja,
    skip = 1,
    col_types = "text",  # No inferir tipos a partir de las primeras horas sin dato.
    na = c("", "N", "NA", "-")
  )

  if (nrow(raw) == 0) {
    return(NULL)
  }

  nombres <- names(raw)

  i_fecha <- buscar_columna(nombres, c("^fecha$", "fecha"))
  i_hora <- buscar_columna(nombres, c("^hora$", "hora"))
  i_pm10 <- buscar_columna(nombres, c("^pm10$", "pm10"))
  i_pm25 <- buscar_columna(nombres, c("^pm2 5$", "^pm25$", "pm2 5"))
  i_no2 <- buscar_columna(nombres, c("^no2$", "no2"))
  i_o3 <- buscar_columna(nombres, c("^o3$", "o3"))

  if (is.na(i_fecha) || is.na(i_hora)) {
    stop("No se identificaron fecha/hora en ", basename(archivo), " - ", hoja)
  }

  obtener <- function(i) {
    if (is.na(i)) {
      return(rep(NA_real_, nrow(raw)))
    }
    a_num(raw[[i]])
  }

  tibble(
    anio = anio,
    isla = normalizar_isla(hoja),
    estacion = trimws(hoja),
    estacion_id = normalizar_nombre(hoja),
    fecha = parsear_fecha(raw[[i_fecha]], anio),
    hora = obtener(i_hora),
    pm10 = obtener(i_pm10),
    pm25 = obtener(i_pm25),
    no2 = obtener(i_no2),
    o3 = obtener(i_o3),
    archivo_origen = basename(archivo)
  ) %>%
    filter(!is.na(fecha))
}

horarios_lista <- list()
estructura_lista <- list()
contador <- 1L

for (archivo in archivos) {
  anio <- as.integer(stringr::str_extract(basename(archivo), "\\d{4}"))

  hojas <- readxl::excel_sheets(archivo)
  hojas <- hojas[normalizar_nombre(hojas) != "leyenda"]

  for (hoja in hojas) {
    x <- leer_estacion(archivo, hoja, anio)

    if (is.null(x)) {
      next
    }

    horarios_lista[[contador]] <- x

    estructura_lista[[contador]] <- tibble(
      anio = anio,
      isla = first(x$isla),
      estacion = first(x$estacion),
      archivo_origen = basename(archivo),
      n_registros_horarios = nrow(x),
      fecha_minima = min(x$fecha, na.rm = TRUE),
      fecha_maxima = max(x$fecha, na.rm = TRUE),
      pm10_disponible = any(!is.na(x$pm10)),
      pm25_disponible = any(!is.na(x$pm25)),
      no2_disponible = any(!is.na(x$no2)),
      o3_disponible = any(!is.na(x$o3))
    )

    contador <- contador + 1L
  }
}

horarios <- bind_rows(horarios_lista)
estructura <- bind_rows(estructura_lista)

duplicados <- horarios %>%
  count(anio, estacion_id, fecha, hora, name = "n") %>%
  filter(n > 1)

if (nrow(duplicados) > 0) {
  warning("Se detectaron ", nrow(duplicados), " combinaciones estación-fecha-hora duplicadas.")
}

horarios <- horarios %>%
  group_by(anio, isla, estacion, estacion_id, fecha, hora) %>%
  summarise(
    pm10 = first(pm10[!is.na(pm10)], default = NA_real_),
    pm25 = first(pm25[!is.na(pm25)], default = NA_real_),
    no2 = first(no2[!is.na(no2)], default = NA_real_),
    o3 = first(o3[!is.na(o3)], default = NA_real_),
    .groups = "drop"
  )

resumen_diario <- horarios %>%
  group_by(isla, estacion, estacion_id, fecha, anio) %>%
  summarise(
    n_horas_pm10 = sum(!is.na(pm10)),
    pm10_media_diaria = ifelse(n_horas_pm10 >= 18, mean(pm10, na.rm = TRUE), NA_real_),
    pm10_max_horario = ifelse(n_horas_pm10 >= 18, max(pm10, na.rm = TRUE), NA_real_),

    n_horas_pm25 = sum(!is.na(pm25)),
    pm25_media_diaria = ifelse(n_horas_pm25 >= 18, mean(pm25, na.rm = TRUE), NA_real_),
    pm25_max_horario = ifelse(n_horas_pm25 >= 18, max(pm25, na.rm = TRUE), NA_real_),

    n_horas_no2 = sum(!is.na(no2)),
    no2_media_diaria = ifelse(n_horas_no2 >= 18, mean(no2, na.rm = TRUE), NA_real_),

    n_horas_o3 = sum(!is.na(o3)),
    o3_media_diaria = ifelse(n_horas_o3 >= 18, mean(o3, na.rm = TRUE), NA_real_),

    .groups = "drop"
  ) %>%
  mutate(
    pm10_dia_valido = n_horas_pm10 >= 18,
    pm25_dia_valido = n_horas_pm25 >= 18,
    no2_dia_valido = n_horas_no2 >= 18,
    o3_dia_valido = n_horas_o3 >= 18,
    panel_pm_dia_valido = pm10_dia_valido & pm25_dia_valido,
    mes = as.Date(format(fecha, "%Y-%m-01"))
  ) %>%
  arrange(isla, estacion, fecha)

cobertura_anual <- resumen_diario %>%
  group_by(anio, isla, estacion) %>%
  summarise(
    dias = n(),
    dias_pm10_validos = sum(pm10_dia_valido),
    dias_pm25_validos = sum(pm25_dia_valido),
    dias_no2_validos = sum(no2_dia_valido),
    dias_o3_validos = sum(o3_dia_valido),
    .groups = "drop"
  )

cobertura_global <- resumen_diario %>%
  group_by(isla, estacion) %>%
  summarise(
    fecha_minima = min(fecha),
    fecha_maxima = max(fecha),
    dias = n(),
    porcentaje_pm10_valido = 100 * mean(pm10_dia_valido),
    porcentaje_pm25_valido = 100 * mean(pm25_dia_valido),
    porcentaje_no2_valido = 100 * mean(no2_dia_valido),
    porcentaje_o3_valido = 100 * mean(o3_dia_valido),
    .groups = "drop"
  )

cobertura_mensual <- resumen_diario %>%
  group_by(isla, estacion, mes) %>%
  summarise(
    dias = n(),
    porcentaje_pm10_valido = 100 * mean(pm10_dia_valido),
    porcentaje_pm25_valido = 100 * mean(pm25_dia_valido),
    porcentaje_no2_valido = 100 * mean(no2_dia_valido),
    porcentaje_o3_valido = 100 * mean(o3_dia_valido),
    .groups = "drop"
  )

write_csv2(estructura, file.path(ruta_salida, "01_estructura_hojas.csv"), na = "")

write_csv2(duplicados, file.path(ruta_salida, "02_duplicados_horarios.csv"), na = "")

write_csv2(resumen_diario, file.path(ruta_salida, "03_datos_diarios_normalizados.csv"), na = "")

write_csv2(cobertura_anual, file.path(ruta_salida, "04_cobertura_anual.csv"), na = "")

write_csv2(cobertura_global, file.path(ruta_salida, "05_cobertura_global.csv"), na = "")

write_csv2(cobertura_mensual, file.path(ruta_salida, "06_cobertura_mensual.csv"), na = "")

cat("\nSERIE DIARIA NORMALIZADA CREADA\n")
cat("Estaciones-año: ", nrow(estructura), "\n", sep = "")
cat("Duplicados horarios: ", nrow(duplicados), "\n", sep = "")
cat("Archivo principal: ", file.path(ruta_salida, "03_datos_diarios_normalizados.csv"), "\n", sep = "")
