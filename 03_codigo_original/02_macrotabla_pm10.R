# MACROTABLA AMBIENTAL, DICCIONARIO Y PAQUETE PARA ROCIO
# Proyecto: GEMAS Canarias - microbiota, asma y exposicion ambiental
#
# Este script:
#   1) Construye una serie diaria jerarquizada por contaminante.
#   2) Cruza cada muestra con D-1 a D-180.
#   3) Genera una macrotabla (una fila por muestra).
#   4) Calcula las variables preespecificadas de PM10/SDI
#      para D-10 a D-1.
#   5) Incluye PM2.5, NO2, O3 y la ventana de 180 dias como bloques secundarios/exploratorios.
#   6) Genera un diccionario completo de variables y una
#      carpeta preparada para compartir con Rocio.
#
# Convenciones:
#   - D0 = fecha de la muestra y NO se incluye en exposicion.
#   - D-1 es el dia inmediatamente previo a la muestra.
#   - Dia valido = >=18 registros horarios validos.
#   - SDI compatible / exposicion compatible con calima:
#       PM10 medio diario >=50 ug/m3.
#   - PM2.5 >=25 ug/m3 se conserva como indicador paralelo de
#     episodio de particulas.

# 0. PAQUETES

paquetes <- c("dplyr", "tidyr", "readr", "writexl")

faltantes <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltantes) > 0) {
  install.packages(faltantes)
}

invisible(lapply(paquetes, library, character.only = TRUE))

# 1. RUTAS

ruta_cohorte <- "outputs/06_cohorte_elegible_cruce_ambiental.csv"

ruta_aire <- paste0(
  "data_raw/CALIMA/auditoria_estaciones/",
  "03_datos_diarios_normalizados.csv"
)

ruta_phyloseq <- "data_raw/data_phyloseq_salivas_GEMAS_Canarias.rds"

ruta_outputs <- "outputs"
ruta_rocio <- file.path(ruta_outputs, "para_rocio")

if (!file.exists(ruta_cohorte)) {
  stop(
    "No se encuentra outputs/06_cohorte_elegible_cruce_ambiental.csv. ",
    "Ejecuta antes el script 02."
  )
}

if (!file.exists(ruta_aire)) {
  stop(
    "No se encuentra 03_datos_diarios_normalizados.csv en la ruta esperada."
  )
}

if (!dir.exists(ruta_outputs)) {
  dir.create(ruta_outputs, recursive = TRUE)
}

if (!dir.exists(ruta_rocio)) {
  dir.create(ruta_rocio, recursive = TRUE)
}

# 2. PARAMETROS PREESPECIFICADOS

dias_ventana_principal <- 10
dias_ventana_secundaria <- 180

min_dias_validos_10d <- 7
min_dias_validos_180d <- 153

# Definiciones operativas del proyecto.
umbral_sdi_pm10 <- 50
umbral_particulas_pm25 <- 25

# Puntos de corte de categorias de calidad del aire de la Red.
# Se aplican a las medias diarias disponibles.
categorias_calidad <- c(
  "Buena",
  "Razonablemente_buena",
  "Regular",
  "Desfavorable",
  "Muy_desfavorable",
  "Extremadamente_desfavorable"
)

# 3. FUNCIONES

a_logico <- function(x) {
  if (is.logical(x)) return(x)

  x <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(x))
  out[x %in% c("TRUE", "T", "1", "SI", "SÍ")] <- TRUE
  out[x %in% c("FALSE", "F", "0", "NO")] <- FALSE
  as.logical(out)
}

a_numero <- function(x) {
  if (is.numeric(x)) return(x)

  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

a_fecha <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  as.Date(as.character(x))
}

media_segura <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

suma_segura <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
}

maximo_seguro <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

sd_segura <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x)
}

p90_seguro <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = 0.90, names = FALSE))
}

auc_umbral <- function(x, umbral) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  sum(pmax(x - umbral, 0), na.rm = TRUE)
}

contar_verdaderos <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  as.integer(sum(x %in% TRUE, na.rm = TRUE))
}

algun_verdadero <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  as.integer(any(x %in% TRUE, na.rm = TRUE))
}

contar_categoria <- function(x, categoria) {
  if (all(is.na(x))) return(NA_integer_)
  as.integer(sum(x == categoria, na.rm = TRUE))
}

peor_nivel <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_integer_)
  as.integer(max(x))
}

peor_categoria <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  categorias_calidad[max(x)]
}

clasificar_calidad <- function(x, contaminante) {
  cortes <- switch(
    contaminante,
    pm10 = c(-Inf, 20, 40, 50, 100, 150, Inf),
    pm25 = c(-Inf, 10, 20, 25, 50, 75, Inf),
    no2  = c(-Inf, 40, 90, 120, 230, 340, Inf),
    o3   = c(-Inf, 50, 100, 130, 240, 380, Inf),
    stop("Contaminante no reconocido en clasificar_calidad().")
  )

  as.character(
    cut(
      x,
      breaks = cortes,
      labels = categorias_calidad,
      include.lowest = TRUE,
      right = TRUE
    )
  )
}

estacion_ano <- function(fecha) {
  mes <- as.integer(format(fecha, "%m"))

  ifelse(
    mes %in% c(12, 1, 2), "Invierno",
    ifelse(
      mes %in% c(3, 4, 5), "Primavera",
      ifelse(mes %in% c(6, 7, 8), "Verano", "Otono")
    )
  )
}

tipo_sdi_estacional <- function(fechas, sdi) {
  fechas <- as.Date(fechas)
  sdi <- sdi %in% TRUE

  if (!any(sdi, na.rm = TRUE)) return("No_expuesto")

  meses <- as.integer(format(fechas[sdi], "%m"))
  n_invernal <- sum(meses %in% c(10, 11, 12, 1, 2, 3))
  n_calido <- sum(meses %in% c(4, 5, 6, 7, 8, 9))

  if (n_invernal > n_calido) return("Invernal")
  if (n_calido > n_invernal) return("Calido")
  "Mixto_transicional"
}

grado_sdi <- function(n_dias) {
  dplyr::case_when(
    is.na(n_dias) ~ NA_character_,
    n_dias == 0 ~ "Sin_exposicion",
    n_dias >= 1 & n_dias <= 3 ~ "Baja_1_3_dias",
    n_dias >= 4 ~ "Alta_ge4_dias",
    TRUE ~ NA_character_
  )
}

patron_temporal_sdi <- function(temprana, tardia) {
  dplyr::case_when(
    is.na(temprana) & is.na(tardia) ~ NA_character_,
    temprana == 0 & tardia == 0 ~ "Sin_SDI",
    temprana == 1 & tardia == 0 ~ "Solo_temprana_D1_D3",
    temprana == 0 & tardia == 1 ~ "Solo_tardia_D4_D10",
    temprana == 1 & tardia == 1 ~ "Temprana_y_tardia",
    TRUE ~ NA_character_
  )
}

seleccionar_serie_pm <- function(aire, jerarquia) {
  # VERSION v5:
  # La seleccion de estacion para PM10 depende EXCLUSIVAMENTE
  # de la validez diaria de PM10. No se exige PM2.5 valido,
  # porque PM2.5 no forma parte del plan analitico final del TFM.
  aire %>%
    inner_join(jerarquia, by = c("isla", "estacion")) %>%
    filter(
      pm10_dia_valido %in% TRUE,
      !is.na(pm10_media_diaria)
    ) %>%
    arrange(isla, fecha, prioridad) %>%
    group_by(isla, fecha) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      isla,
      fecha,
      estacion_origen_pm = estacion,
      nivel_estacion_pm = nivel,
      prioridad_estacion_pm = prioridad,
      n_horas_pm10,
      pm10_media_diaria,
      pm10_max_horario,
      n_horas_pm25,
      pm25_media_diaria,
      pm25_max_horario,
      dia_pm_valido = TRUE
    )
}

seleccionar_serie_contaminante <- function(
    aire,
    jerarquia,
    contaminante
) {
  col_valido <- paste0(contaminante, "_dia_valido")
  col_media <- paste0(contaminante, "_media_diaria")
  col_horas <- paste0("n_horas_", contaminante)

  aire %>%
    inner_join(jerarquia, by = c("isla", "estacion")) %>%
    filter(
      .data[[col_valido]] %in% TRUE,
      !is.na(.data[[col_media]])
    ) %>%
    arrange(isla, fecha, prioridad) %>%
    group_by(isla, fecha) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      isla,
      fecha,
      estacion_origen = estacion,
      nivel_estacion = nivel,
      prioridad_estacion = prioridad,
      n_horas = .data[[col_horas]],
      media_diaria = .data[[col_media]],
      dia_valido = TRUE
    ) %>%
    rename(
      !!paste0("estacion_origen_", contaminante) := estacion_origen,
      !!paste0("nivel_estacion_", contaminante) := nivel_estacion,
      !!paste0("prioridad_estacion_", contaminante) := prioridad_estacion,
      !!paste0("n_horas_", contaminante) := n_horas,
      !!paste0(contaminante, "_media_diaria") := media_diaria,
      !!paste0("dia_", contaminante, "_valido") := dia_valido
    )
}

# 4. COHORTE CLINICA
cohorte <- readr::read_csv(ruta_cohorte, show_col_types = FALSE)

columnas_minimas <- c("id_muestra", "fecha_muestra", "isla")
faltantes_cohorte <- setdiff(columnas_minimas, names(cohorte))

if (length(faltantes_cohorte) > 0) {
  stop(
    "Faltan columnas esenciales en la cohorte: ",
    paste(faltantes_cohorte, collapse = ", ")
  )
}

cohorte <- cohorte %>%
  mutate(
    id_muestra = as.character(id_muestra),
    fecha_muestra = a_fecha(fecha_muestra),
    isla = trimws(as.character(isla))
  ) %>%
  filter(!is.na(fecha_muestra), !is.na(isla))

if (anyDuplicated(cohorte$id_muestra) > 0) {
  stop("Hay IDs de muestra duplicados en la cohorte.")
}

if (!"adulto" %in% names(cohorte) && "edad" %in% names(cohorte)) {
  cohorte$adulto <- ifelse(is.na(cohorte$edad), NA, cohorte$edad >= 18)
}

cohorte <- cohorte %>%
  mutate(
    anio_muestra = as.integer(format(fecha_muestra, "%Y")),
    mes_muestra = as.integer(format(fecha_muestra, "%m")),
    estacion_muestra = estacion_ano(fecha_muestra)
  )

cat("\nCohorte elegible cargada: ", nrow(cohorte), " muestras.\n", sep = "")
print(table(cohorte$isla))


# 5. DATOS AMBIENTALES DIARIOS

aire <- readr::read_csv2(ruta_aire, show_col_types = FALSE)

columnas_aire <- c(
  "isla", "estacion", "fecha",
  "n_horas_pm10", "pm10_media_diaria", "pm10_max_horario",
  "n_horas_pm25", "pm25_media_diaria", "pm25_max_horario",
  "n_horas_no2", "no2_media_diaria",
  "n_horas_o3", "o3_media_diaria",
  "pm10_dia_valido", "pm25_dia_valido", "panel_pm_dia_valido",
  "no2_dia_valido", "o3_dia_valido"
)

faltantes_aire <- setdiff(columnas_aire, names(aire))

if (length(faltantes_aire) > 0) {
  stop(
    "Faltan columnas en los datos ambientales diarios: ",
    paste(faltantes_aire, collapse = ", ")
  )
}

aire <- aire %>%
  mutate(
    isla = trimws(as.character(isla)),
    estacion = trimws(as.character(estacion)),
    fecha = a_fecha(fecha)
  )

for (columna in c(
  "pm10_dia_valido", "pm25_dia_valido", "panel_pm_dia_valido",
  "no2_dia_valido", "o3_dia_valido"
)) {
  aire[[columna]] <- a_logico(aire[[columna]])
}

for (columna in c(
  "n_horas_pm10", "pm10_media_diaria", "pm10_max_horario",
  "n_horas_pm25", "pm25_media_diaria", "pm25_max_horario",
  "n_horas_no2", "no2_media_diaria",
  "n_horas_o3", "o3_media_diaria"
)) {
  aire[[columna]] <- a_numero(aire[[columna]])
}


# 6. JERARQUIA DE ESTACIONES
#
# PM10/PM2.5: jerarquia validada en la auditoria previa.
# NO2/O3: bloque exploratorio; se deja documentada la
# jerarquia utilizada y se evaluó su cobertura antes de decidir su incorporación al analisis final.


jerarquia_pm <- tibble::tribble(
  ~isla,       ~estacion,                         ~prioridad, ~nivel,
  "Tenerife",  "El Río",                           1L,         "primaria",
  "Tenerife",  "Granadilla",                       2L,         "secundaria",
  "La Palma",  "El Pilar-Sta Cruz de La Palma",    1L,         "primaria",
  "La Palma",  "La Grama-Breña Alta",              2L,         "secundaria"
)

jerarquia_no2 <- tibble::tribble(
  ~isla,       ~estacion,                         ~prioridad, ~nivel,
  "Tenerife",  "El Río",                           1L,         "primaria",
  "Tenerife",  "Granadilla",                       2L,         "secundaria",
  "La Palma",  "El Pilar-Sta Cruz de La Palma",    1L,         "primaria",
  "La Palma",  "La Grama-Breña Alta",              2L,         "secundaria"
)

jerarquia_o3 <- tibble::tribble(
  ~isla,       ~estacion,                         ~prioridad, ~nivel,
  "Tenerife",  "El Río",                           1L,         "primaria",
  "Tenerife",  "Buzanada",                         2L,         "secundaria",
  "La Palma",  "El Pilar-Sta Cruz de La Palma",    1L,         "primaria",
  "La Palma",  "La Grama-Breña Alta",              2L,         "secundaria"
)

jerarquia_estaciones <- bind_rows(
  mutate(jerarquia_pm, contaminante = "PM10_y_PM2.5"),
  mutate(jerarquia_no2, contaminante = "NO2"),
  mutate(jerarquia_o3, contaminante = "O3")
) %>%
  select(contaminante, isla, estacion, prioridad, nivel)

estaciones_configuradas <- unique(jerarquia_estaciones$estacion)
estaciones_no_encontradas <- setdiff(estaciones_configuradas, unique(aire$estacion))

if (length(estaciones_no_encontradas) > 0) {
  warning(
    "No se encontraron estas estaciones configuradas: ",
    paste(estaciones_no_encontradas, collapse = ", "),
    ". Revisa los nombres exactos en unique(aire$estacion)."
  )
}

# 7. SERIES DIARIAS JERARQUIZADAS Y CLASIFICACION

serie_pm <- seleccionar_serie_pm(aire, jerarquia_pm)
serie_no2 <- seleccionar_serie_contaminante(aire, jerarquia_no2, "no2")
serie_o3 <- seleccionar_serie_contaminante(aire, jerarquia_o3, "o3")

serie_ambiental <- serie_pm %>%
  full_join(serie_no2, by = c("isla", "fecha")) %>%
  full_join(serie_o3, by = c("isla", "fecha")) %>%
  arrange(isla, fecha) %>%
  mutate(
    pm10_calidad_diaria = clasificar_calidad(pm10_media_diaria, "pm10"),
    pm25_calidad_diaria = clasificar_calidad(pm25_media_diaria, "pm25"),
    no2_calidad_diaria = clasificar_calidad(no2_media_diaria, "no2"),
    o3_calidad_diaria = clasificar_calidad(o3_media_diaria, "o3"),

    pm10_nivel_calidad = match(pm10_calidad_diaria, categorias_calidad),
    pm25_nivel_calidad = match(pm25_calidad_diaria, categorias_calidad),
    no2_nivel_calidad = match(no2_calidad_diaria, categorias_calidad),
    o3_nivel_calidad = match(o3_calidad_diaria, categorias_calidad),

    # Definicion principal del TFM para SDI compatible / calima operativa.
    sdi_pm10_ge50 = ifelse(
      is.na(pm10_media_diaria), NA, pm10_media_diaria >= umbral_sdi_pm10
    ),

    pm25_ge25 = ifelse(
      is.na(pm25_media_diaria), NA, pm25_media_diaria >= umbral_particulas_pm25
    ),

    episodio_particulas_cualquiera = case_when(
      is.na(sdi_pm10_ge50) & is.na(pm25_ge25) ~ NA,
      sdi_pm10_ge50 %in% TRUE | pm25_ge25 %in% TRUE ~ TRUE,
      TRUE ~ FALSE
    ),

    episodio_particulas_concordante = case_when(
      is.na(sdi_pm10_ge50) | is.na(pm25_ge25) ~ NA,
      sdi_pm10_ge50 %in% TRUE & pm25_ge25 %in% TRUE ~ TRUE,
      TRUE ~ FALSE
    ),

    pm10_desfavorable_o_peor = ifelse(
      is.na(pm10_nivel_calidad), NA, pm10_nivel_calidad >= 4
    ),
    pm25_desfavorable_o_peor = ifelse(
      is.na(pm25_nivel_calidad), NA, pm25_nivel_calidad >= 4
    ),
    no2_desfavorable_o_peor = ifelse(
      is.na(no2_nivel_calidad), NA, no2_nivel_calidad >= 4
    ),
    o3_desfavorable_o_peor = ifelse(
      is.na(o3_nivel_calidad), NA, o3_nivel_calidad >= 4
    )
  )

readr::write_csv2(
  serie_ambiental,
  file.path(ruta_outputs, "15_serie_diaria_multicontaminante_jerarquizada.csv"),
  na = ""
)


# 8. TRAZABILIDAD DIA A DIA POR MUESTRA

trazabilidad <- cohorte %>%
  tidyr::crossing(lag_dias = 1:dias_ventana_secundaria) %>%
  mutate(
    fecha_ambiental = fecha_muestra - lag_dias,
    pertenece_10d = lag_dias <= dias_ventana_principal
  ) %>%
  left_join(
    serie_ambiental,
    by = c("isla", "fecha_ambiental" = "fecha")
  )

for (columna in c("dia_pm_valido", "dia_no2_valido", "dia_o3_valido")) {
  trazabilidad[[columna]][is.na(trazabilidad[[columna]])] <- FALSE
}

readr::write_csv2(
  trazabilidad %>%
    select(
      id_muestra, fecha_muestra, isla, lag_dias, fecha_ambiental, pertenece_10d,
      dia_pm_valido, estacion_origen_pm, nivel_estacion_pm, prioridad_estacion_pm,
      n_horas_pm10, pm10_media_diaria, pm10_max_horario,
      pm10_calidad_diaria, pm10_nivel_calidad, sdi_pm10_ge50,
      pm10_desfavorable_o_peor,
      n_horas_pm25, pm25_media_diaria, pm25_max_horario,
      pm25_calidad_diaria, pm25_nivel_calidad, pm25_ge25,
      pm25_desfavorable_o_peor,
      episodio_particulas_cualquiera, episodio_particulas_concordante,
      dia_no2_valido, estacion_origen_no2, nivel_estacion_no2,
      prioridad_estacion_no2, n_horas_no2, no2_media_diaria,
      no2_calidad_diaria, no2_nivel_calidad, no2_desfavorable_o_peor,
      dia_o3_valido, estacion_origen_o3, nivel_estacion_o3,
      prioridad_estacion_o3, n_horas_o3, o3_media_diaria,
      o3_calidad_diaria, o3_nivel_calidad, o3_desfavorable_o_peor
    ),
  file.path(ruta_outputs, "16_trazabilidad_diaria_multicontaminante_por_muestra.csv"),
  na = ""
)

# 9. RESUMEN DE VENTANA


resumir_ventana <- function(datos, sufijo, min_dias_pm) {

  salida <- datos %>%
    group_by(id_muestra) %>%
    summarise(
      fecha_inicio = min(fecha_ambiental),
      fecha_fin = max(fecha_ambiental),
      n_dias_esperados = n(),

      # Cobertura
      n_dias_validos_pm = sum(dia_pm_valido, na.rm = TRUE),
      n_dias_no_validos_pm = sum(!dia_pm_valido, na.rm = TRUE),
      porcentaje_dias_validos_pm = round(100 * mean(dia_pm_valido), 2),

      n_dias_validos_no2 = sum(dia_no2_valido, na.rm = TRUE),
      n_dias_no_validos_no2 = sum(!dia_no2_valido, na.rm = TRUE),
      porcentaje_dias_validos_no2 = round(100 * mean(dia_no2_valido), 2),

      n_dias_validos_o3 = sum(dia_o3_valido, na.rm = TRUE),
      n_dias_no_validos_o3 = sum(!dia_o3_valido, na.rm = TRUE),
      porcentaje_dias_validos_o3 = round(100 * mean(dia_o3_valido), 2),

      # Trazabilidad de estaciones
      n_dias_pm_primaria = sum(nivel_estacion_pm == "primaria", na.rm = TRUE),
      n_dias_pm_secundaria = sum(nivel_estacion_pm == "secundaria", na.rm = TRUE),
      n_dias_no2_primaria = sum(nivel_estacion_no2 == "primaria", na.rm = TRUE),
      n_dias_no2_secundaria = sum(nivel_estacion_no2 == "secundaria", na.rm = TRUE),
      n_dias_o3_primaria = sum(nivel_estacion_o3 == "primaria", na.rm = TRUE),
      n_dias_o3_secundaria = sum(nivel_estacion_o3 == "secundaria", na.rm = TRUE),

      # PM10
      pm10_media = media_segura(pm10_media_diaria),
      pm10_maximo = maximo_seguro(pm10_media_diaria),
      pm10_acumulado = suma_segura(pm10_media_diaria),
      pm10_sd = sd_segura(pm10_media_diaria),
      pm10_p90 = p90_seguro(pm10_media_diaria),
      pm10_auc50 = auc_umbral(pm10_media_diaria, umbral_sdi_pm10),

      n_dias_sdi = contar_verdaderos(sdi_pm10_ge50),
      sdi_exposicion = algun_verdadero(sdi_pm10_ge50),

      n_dias_pm10_buena = contar_categoria(pm10_calidad_diaria, "Buena"),
      n_dias_pm10_razonablemente_buena = contar_categoria(pm10_calidad_diaria, "Razonablemente_buena"),
      n_dias_pm10_regular = contar_categoria(pm10_calidad_diaria, "Regular"),
      n_dias_pm10_desfavorable = contar_categoria(pm10_calidad_diaria, "Desfavorable"),
      n_dias_pm10_muy_desfavorable = contar_categoria(pm10_calidad_diaria, "Muy_desfavorable"),
      n_dias_pm10_extremadamente_desfavorable = contar_categoria(pm10_calidad_diaria, "Extremadamente_desfavorable"),
      n_dias_pm10_desfavorable_o_peor = contar_verdaderos(pm10_desfavorable_o_peor),
      pm10_exposicion_desfavorable_o_peor = algun_verdadero(pm10_desfavorable_o_peor),
      peor_nivel_calidad_pm10 = peor_nivel(pm10_nivel_calidad),

      # PM2.5
      pm25_media = media_segura(pm25_media_diaria),
      pm25_maximo = maximo_seguro(pm25_media_diaria),
      pm25_acumulado = suma_segura(pm25_media_diaria),
      pm25_sd = sd_segura(pm25_media_diaria),
      pm25_p90 = p90_seguro(pm25_media_diaria),
      pm25_auc25 = auc_umbral(pm25_media_diaria, umbral_particulas_pm25),

      n_dias_pm25_ge25 = contar_verdaderos(pm25_ge25),
      pm25_exposicion_ge25 = algun_verdadero(pm25_ge25),

      n_dias_pm25_buena = contar_categoria(pm25_calidad_diaria, "Buena"),
      n_dias_pm25_razonablemente_buena = contar_categoria(pm25_calidad_diaria, "Razonablemente_buena"),
      n_dias_pm25_regular = contar_categoria(pm25_calidad_diaria, "Regular"),
      n_dias_pm25_desfavorable = contar_categoria(pm25_calidad_diaria, "Desfavorable"),
      n_dias_pm25_muy_desfavorable = contar_categoria(pm25_calidad_diaria, "Muy_desfavorable"),
      n_dias_pm25_extremadamente_desfavorable = contar_categoria(pm25_calidad_diaria, "Extremadamente_desfavorable"),
      n_dias_pm25_desfavorable_o_peor = contar_verdaderos(pm25_desfavorable_o_peor),
      pm25_exposicion_desfavorable_o_peor = algun_verdadero(pm25_desfavorable_o_peor),
      peor_nivel_calidad_pm25 = peor_nivel(pm25_nivel_calidad),

      # Episodios combinados de particulas
      n_dias_episodio_particulas_cualquiera = contar_verdaderos(episodio_particulas_cualquiera),
      episodio_particulas_cualquiera_exposicion = algun_verdadero(episodio_particulas_cualquiera),
      n_dias_episodio_particulas_concordante = contar_verdaderos(episodio_particulas_concordante),
      episodio_particulas_concordante_exposicion = algun_verdadero(episodio_particulas_concordante),

      # NO2
      no2_media = media_segura(no2_media_diaria),
      no2_maximo = maximo_seguro(no2_media_diaria),
      no2_acumulado = suma_segura(no2_media_diaria),
      no2_sd = sd_segura(no2_media_diaria),
      no2_p90 = p90_seguro(no2_media_diaria),
      no2_auc120 = auc_umbral(no2_media_diaria, 120),

      n_dias_no2_buena = contar_categoria(no2_calidad_diaria, "Buena"),
      n_dias_no2_razonablemente_buena = contar_categoria(no2_calidad_diaria, "Razonablemente_buena"),
      n_dias_no2_regular = contar_categoria(no2_calidad_diaria, "Regular"),
      n_dias_no2_desfavorable = contar_categoria(no2_calidad_diaria, "Desfavorable"),
      n_dias_no2_muy_desfavorable = contar_categoria(no2_calidad_diaria, "Muy_desfavorable"),
      n_dias_no2_extremadamente_desfavorable = contar_categoria(no2_calidad_diaria, "Extremadamente_desfavorable"),
      n_dias_no2_desfavorable_o_peor = contar_verdaderos(no2_desfavorable_o_peor),
      no2_exposicion_desfavorable_o_peor = algun_verdadero(no2_desfavorable_o_peor),
      peor_nivel_calidad_no2 = peor_nivel(no2_nivel_calidad),

      # O3
      o3_media = media_segura(o3_media_diaria),
      o3_maximo = maximo_seguro(o3_media_diaria),
      o3_acumulado = suma_segura(o3_media_diaria),
      o3_sd = sd_segura(o3_media_diaria),
      o3_p90 = p90_seguro(o3_media_diaria),
      o3_auc130 = auc_umbral(o3_media_diaria, 130),

      n_dias_o3_buena = contar_categoria(o3_calidad_diaria, "Buena"),
      n_dias_o3_razonablemente_buena = contar_categoria(o3_calidad_diaria, "Razonablemente_buena"),
      n_dias_o3_regular = contar_categoria(o3_calidad_diaria, "Regular"),
      n_dias_o3_desfavorable = contar_categoria(o3_calidad_diaria, "Desfavorable"),
      n_dias_o3_muy_desfavorable = contar_categoria(o3_calidad_diaria, "Muy_desfavorable"),
      n_dias_o3_extremadamente_desfavorable = contar_categoria(o3_calidad_diaria, "Extremadamente_desfavorable"),
      n_dias_o3_desfavorable_o_peor = contar_verdaderos(o3_desfavorable_o_peor),
      o3_exposicion_desfavorable_o_peor = algun_verdadero(o3_desfavorable_o_peor),
      peor_nivel_calidad_o3 = peor_nivel(o3_nivel_calidad),

      .groups = "drop"
    ) %>%
    mutate(
      ventana_pm_valida = n_dias_validos_pm >= min_dias_pm,
      ventana_no2_valida = n_dias_validos_no2 >= min_dias_pm,
      ventana_o3_valida = n_dias_validos_o3 >= min_dias_pm,
      # La función peor_categoria() recibe aquí un vector de niveles,
      # por lo que debe aplicarse elemento a elemento. De otro modo
      # devolvería el máximo global de toda la cohorte y lo reciclaría.
      peor_calidad_pm10 = vapply(
        peor_nivel_calidad_pm10,
        peor_categoria,
        character(1)
      ),
      peor_calidad_pm25 = vapply(
        peor_nivel_calidad_pm25,
        peor_categoria,
        character(1)
      ),
      peor_calidad_no2 = vapply(
        peor_nivel_calidad_no2,
        peor_categoria,
        character(1)
      ),
      peor_calidad_o3 = vapply(
        peor_nivel_calidad_o3,
        peor_categoria,
        character(1)
      )
    )

  nombres <- names(salida)
  nombres[nombres != "id_muestra"] <- paste0(
    nombres[nombres != "id_muestra"], "_", sufijo
  )
  names(salida) <- nombres

  salida
}


# 10. RESUMEN D-10 A D-1 Y D-180 A D-1

resumen_10d <- resumir_ventana(
  trazabilidad %>% filter(lag_dias <= dias_ventana_principal),
  sufijo = "10d",
  min_dias_pm = min_dias_validos_10d
)

resumen_180d <- resumir_ventana(
  trazabilidad,
  sufijo = "180d",
  min_dias_pm = min_dias_validos_180d
)

# Variables temporales de SDI especificadas en el protocolo: solo D-10 a D-1.
temporal_sdi_10d <- trazabilidad %>%
  filter(lag_dias <= dias_ventana_principal) %>%
  group_by(id_muestra) %>%
  summarise(
    sdi_temprana_10d = algun_verdadero(
      sdi_pm10_ge50[lag_dias >= 1 & lag_dias <= 3]
    ),
    sdi_tardia_10d = algun_verdadero(
      sdi_pm10_ge50[lag_dias >= 4 & lag_dias <= 10]
    ),
    tipo_sdi_estacional_10d = tipo_sdi_estacional(
      fecha_ambiental, sdi_pm10_ge50
    ),
    .groups = "drop"
  )

temporal_sdi_10d <- temporal_sdi_10d %>%
  left_join(
    resumen_10d %>% select(id_muestra, n_dias_sdi_10d),
    by = "id_muestra"
  ) %>%
  mutate(
    grado_sdi_10d = grado_sdi(n_dias_sdi_10d),
    patron_temporal_sdi_10d = patron_temporal_sdi(
      sdi_temprana_10d, sdi_tardia_10d
    )
  ) %>%
  select(-n_dias_sdi_10d)


# 11. EXPOSICIONES DIARIAS D-1 A D-10 

diarias_10d_ancho <- trazabilidad %>%
  filter(lag_dias <= dias_ventana_principal) %>%
  mutate(lag = sprintf("d%02d", lag_dias)) %>%
  select(
    id_muestra, lag,
    pm10_media_diaria, pm25_media_diaria, no2_media_diaria, o3_media_diaria,
    sdi_pm10_ge50, pm25_ge25,
    episodio_particulas_cualquiera, episodio_particulas_concordante,
    pm10_calidad_diaria, pm25_calidad_diaria, no2_calidad_diaria, o3_calidad_diaria
  ) %>%
  rename(
    pm10 = pm10_media_diaria,
    pm25 = pm25_media_diaria,
    no2 = no2_media_diaria,
    o3 = o3_media_diaria
  ) %>%
  pivot_wider(
    names_from = lag,
    values_from = c(
      pm10, pm25, no2, o3,
      sdi_pm10_ge50, pm25_ge25,
      episodio_particulas_cualquiera, episodio_particulas_concordante,
      pm10_calidad_diaria, pm25_calidad_diaria,
      no2_calidad_diaria, o3_calidad_diaria
    ),
    names_glue = "{.value}_{lag}"
  )


# 12. MACROTABLA FINAL POR MUESTRA

macrotabla <- cohorte %>%
  left_join(resumen_10d, by = "id_muestra") %>%
  left_join(temporal_sdi_10d, by = "id_muestra") %>%
  left_join(resumen_180d, by = "id_muestra") %>%
  left_join(diarias_10d_ancho, by = "id_muestra") %>%
  mutate(
    estado_ventana_pm_10d = ifelse(
      ventana_pm_valida_10d %in% TRUE,
      "Valida_ge7_dias",
      "Cobertura_insuficiente"
    ),
    estado_ventana_pm_180d = ifelse(
      ventana_pm_valida_180d %in% TRUE,
      "Valida_ge153_dias",
      "Cobertura_insuficiente"
    ),
    estado_ventana_no2_10d = ifelse(
      ventana_no2_valida_10d %in% TRUE,
      "Valida_ge7_dias",
      "Cobertura_insuficiente"
    ),
    estado_ventana_no2_180d = ifelse(
      ventana_no2_valida_180d %in% TRUE,
      "Valida_ge153_dias",
      "Cobertura_insuficiente"
    ),
    estado_ventana_o3_10d = ifelse(
      ventana_o3_valida_10d %in% TRUE,
      "Valida_ge7_dias",
      "Cobertura_insuficiente"
    ),
    estado_ventana_o3_180d = ifelse(
      ventana_o3_valida_180d %in% TRUE,
      "Valida_ge153_dias",
      "Cobertura_insuficiente"
    )
  ) %>%
  arrange(isla, fecha_muestra, id_muestra)


# 13. RESUMEN DE COBERTURA

resumen_cobertura <- macrotabla %>%
  group_by(isla) %>%
  summarise(
    n_muestras = n(),
    pm_10d_validas = sum(ventana_pm_valida_10d %in% TRUE),
    pm_180d_validas = sum(ventana_pm_valida_180d %in% TRUE),
    no2_10d_validas = sum(ventana_no2_valida_10d %in% TRUE),
    no2_180d_validas = sum(ventana_no2_valida_180d %in% TRUE),
    o3_10d_validas = sum(ventana_o3_valida_10d %in% TRUE),
    o3_180d_validas = sum(ventana_o3_valida_180d %in% TRUE),
    media_dias_validos_pm_10d = round(mean(n_dias_validos_pm_10d, na.rm = TRUE), 2),
    media_dias_validos_pm_180d = round(mean(n_dias_validos_pm_180d, na.rm = TRUE), 2),
    media_dias_validos_no2_10d = round(mean(n_dias_validos_no2_10d, na.rm = TRUE), 2),
    media_dias_validos_no2_180d = round(mean(n_dias_validos_no2_180d, na.rm = TRUE), 2),
    media_dias_validos_o3_10d = round(mean(n_dias_validos_o3_10d, na.rm = TRUE), 2),
    media_dias_validos_o3_180d = round(mean(n_dias_validos_o3_180d, na.rm = TRUE), 2),
    .groups = "drop"
  )

# 14. DICCIONARIO DE VARIABLES

definiciones_calidad <- tibble::tribble(
  ~contaminante, ~categoria, ~rango_ug_m3, ~desfavorable_o_peor,
  "PM10",  "Buena",                         "0-20",     FALSE,
  "PM10",  "Razonablemente_buena",          "21-40",    FALSE,
  "PM10",  "Regular",                       "41-50",    FALSE,
  "PM10",  "Desfavorable",                  "51-100",   TRUE,
  "PM10",  "Muy_desfavorable",              "101-150",  TRUE,
  "PM10",  "Extremadamente_desfavorable",   ">=151",    TRUE,
  "PM2.5", "Buena",                         "0-10",     FALSE,
  "PM2.5", "Razonablemente_buena",          "11-20",    FALSE,
  "PM2.5", "Regular",                       "21-25",    FALSE,
  "PM2.5", "Desfavorable",                  "26-50",    TRUE,
  "PM2.5", "Muy_desfavorable",              "51-75",    TRUE,
  "PM2.5", "Extremadamente_desfavorable",   ">=76",     TRUE,
  "NO2",   "Buena",                         "0-40",     FALSE,
  "NO2",   "Razonablemente_buena",          "41-90",    FALSE,
  "NO2",   "Regular",                       "91-120",   FALSE,
  "NO2",   "Desfavorable",                  "121-230",  TRUE,
  "NO2",   "Muy_desfavorable",              "231-340",  TRUE,
  "NO2",   "Extremadamente_desfavorable",   ">=341",    TRUE,
  "O3",    "Buena",                         "0-50",     FALSE,
  "O3",    "Razonablemente_buena",          "51-100",   FALSE,
  "O3",    "Regular",                       "101-130",  FALSE,
  "O3",    "Desfavorable",                  "131-240",  TRUE,
  "O3",    "Muy_desfavorable",              "241-380",  TRUE,
  "O3",    "Extremadamente_desfavorable",   ">=381",    TRUE
)

unidad_variable <- function(v) {
  if (grepl("fecha", v)) return("fecha ISO (AAAA-MM-DD)")
  if (grepl("pm10_|pm25_|no2_|o3_", v) && grepl("_d[0-9]{2}$|_media_|_maximo_|_p90_", v)) return("ug/m3")
  if (grepl("_acumulado_|_auc", v)) return("ug*día/m3")
  if (grepl("^n_dias_", v) || grepl("dias_validos|dias_no_validos|dias_primaria|dias_secundaria", v)) return("días")
  if (grepl("porcentaje", v)) return("%")
  if (grepl("exposicion|sdi_temprana|sdi_tardia|sdi_pm10_ge50|episodio_particulas|ventana_.*valida|adulto|_ge25$", v)) return("0/1")
  if (grepl("edad", v)) return("años")
  if (grepl("anio", v)) return("año")
  if (grepl("mes", v)) return("1-12")
  if (grepl("nivel", v)) return("ordinal 1-6")
  if (grepl("calidad|grado|patron|estacion|isla|sexo|estado", v)) return("categórica")
  "según variable"
}

ventana_variable <- function(v) {
  if (grepl("_d[0-9]{2}$", v)) {
    lag <- sub("^.*_d", "D-", v)
    return(lag)
  }
  if (grepl("_10d$", v)) return("D-10 a D-1")
  if (grepl("_180d$", v)) return("D-180 a D-1")
  if (v %in% c("fecha_muestra", "anio_muestra", "mes_muestra", "estacion_muestra")) return("D0 / muestra basal")
  "No aplica"
}

origen_variable <- function(v) {
  if (v %in% names(cohorte)) return("Metadatos clinicos/técnicos del objeto phyloseq GEMAS")
  if (grepl("pm10|pm25|no2|o3|sdi|episodio_particulas|calidad|ventana|dias_validos|estacion_origen|nivel_estacion|prioridad_estacion|fecha_inicio|fecha_fin", v)) {
    return("Derivada de datos diarios de la Red Canaria de Calidad del Aire")
  }
  "Derivada en este script"
}

descripcion_variable <- function(v) {
  # Identificación y covariables.
  fijas <- c(
    id_muestra = "Identificador pseudonimizado de la muestra de saliva.",
    fecha_muestra = "Fecha de toma de la muestra basal de saliva (D0). D0 no entra en los calculos ambientales.",
    isla = "Isla asociada a la muestra: Tenerife o La Palma.",
    edad = "Edad del participante, en años.",
    sexo = "Sexo registrado en metadatos GEMAS.",
    sexo_codigo = "Codigo numerico original de sexo, si está disponible en la cohorte.",
    exacerbacion_6m = "Exacerbaciones en los seis meses previos: 0 = no exacerbador; 1 = exacerbador.",
    exacerbaciones_6m_codigo = "Codigo original de exacerbaciones a 6 meses, previo a la recodificacion.",
    exacerbaciones_12m_codigo = "Codigo original de exacerbaciones a 12 meses.",
    antibiotico_2m = "Uso de antibioticos en los dos meses previos, segun la codificacion original GEMAS.",
    pool = "Pool/lote tecnico de procesamiento, segun metadatos GEMAS.",
    caries = "Variable de caries registrada en metadatos GEMAS.",
    chicle_30min = "Consumo de chicle en los 30 minutos previos, segun metadatos GEMAS.",
    adulto = "Indicador de edad >=18 años.",
    anio_muestra = "Año de la muestra basal.",
    mes_muestra = "Mes de la muestra basal.",
    estacion_muestra = "Estacion del año derivada de la fecha de la muestra."
  )
  if (v %in% names(fijas)) return(unname(fijas[v]))

  # Valores diarios D-1...D-10.
  if (grepl("^pm10_d[0-9]{2}$", v)) {
    lag <- sub("^pm10_d", "D-", v)
    return(paste0("Media diaria de PM10 correspondiente a ", lag, " respecto a la muestra."))
  }
  if (grepl("^pm25_d[0-9]{2}$", v)) {
    lag <- sub("^pm25_d", "D-", v)
    return(paste0("Media diaria de PM2.5 correspondiente a ", lag, " respecto a la muestra."))
  }
  if (grepl("^no2_d[0-9]{2}$", v)) {
    lag <- sub("^no2_d", "D-", v)
    return(paste0("Media diaria de NO2 correspondiente a ", lag, " respecto a la muestra."))
  }
  if (grepl("^o3_d[0-9]{2}$", v)) {
    lag <- sub("^o3_d", "D-", v)
    return(paste0("Media diaria de O3 correspondiente a ", lag, " respecto a la muestra."))
  }

  # PM10.
  if (grepl("^pm10_media_", v)) return("Media de las concentraciones medias diarias de PM10 dentro de la ventana.")
  if (grepl("^pm10_maximo_", v)) return("Mayor media diaria de PM10 observada dentro de la ventana.")
  if (grepl("^pm10_acumulado_", v)) return("Suma de las medias diarias de PM10 dentro de la ventana.")
  if (grepl("^pm10_sd_", v)) return("Desviacion estandar de las medias diarias de PM10 dentro de la ventana.")
  if (grepl("^pm10_p90_", v)) return("Percentil 90 de las medias diarias de PM10 dentro de la ventana.")
  if (grepl("^pm10_auc50_", v)) return("Carga acumulada de PM10 por encima de 50 ug/m3: suma de max(PM10 diario - 50, 0).")
  if (grepl("^n_dias_sdi_", v)) return("Numero de dias con PM10 medio diario >=50 ug/m3; definicion operativa principal de SDI/exposicion compatible con calima.")
  if (grepl("^sdi_exposicion_", v)) return("Indicador binario: 1 si hubo al menos un dia con PM10 >=50 ug/m3; 0 si no hubo ninguno.")
  if (grepl("^grado_sdi_", v)) return("Grado de exposicion SDI en 10 dias: Sin_exposicion (0), Baja_1_3_dias (1-3) o Alta_ge4_dias (>=4).")
  if (grepl("^sdi_temprana_", v)) return("Indicador binario de al menos un dia SDI en D-1 a D-3.")
  if (grepl("^sdi_tardia_", v)) return("Indicador binario de al menos un dia SDI en D-4 a D-10.")
  if (grepl("^patron_temporal_sdi_", v)) return("Patron temporal de SDI: Sin_SDI, Solo_temprana_D1_D3, Solo_tardia_D4_D10 o Temprana_y_tardia.")
  if (grepl("^tipo_sdi_estacional_", v)) return("Patron estacional predominante entre los dias SDI: Invernal (octubre-marzo), Calido (abril-septiembre) o Mixto_transicional.")
  if (grepl("^n_dias_pm10_", v) && grepl("buena|regular|desfavorable", v)) return("Numero de dias de la ventana clasificados en la categoria indicada de calidad del aire para PM10.")
  if (grepl("^pm10_exposicion_desfavorable_o_peor_", v)) return("Indicador binario de al menos un dia PM10 en categoria Desfavorable, Muy desfavorable o Extremadamente desfavorable.")
  if (grepl("^n_dias_pm10_desfavorable_o_peor_", v)) return("Numero de dias PM10 en categoria Desfavorable, Muy desfavorable o Extremadamente desfavorable.")
  if (grepl("^peor_calidad_pm10_", v)) return("Peor categoria diaria de calidad del aire por PM10 observada en la ventana.")
  if (grepl("^peor_nivel_calidad_pm10_", v)) return("Nivel ordinal de la peor categoria PM10: 1=Buena a 6=Extremadamente_desfavorable.")

  # PM2.5.
  if (grepl("^pm25_media_", v)) return("Media de las concentraciones medias diarias de PM2.5 dentro de la ventana.")
  if (grepl("^pm25_maximo_", v)) return("Mayor media diaria de PM2.5 observada dentro de la ventana.")
  if (grepl("^pm25_acumulado_", v)) return("Suma de las medias diarias de PM2.5 dentro de la ventana.")
  if (grepl("^pm25_sd_", v)) return("Desviacion estandar de las medias diarias de PM2.5 dentro de la ventana.")
  if (grepl("^pm25_p90_", v)) return("Percentil 90 de las medias diarias de PM2.5 dentro de la ventana.")
  if (grepl("^pm25_auc25_", v)) return("Carga acumulada de PM2.5 por encima de 25 ug/m3: suma de max(PM2.5 diario - 25, 0).")
  if (grepl("^n_dias_pm25_ge25_", v)) return("Numero de dias con PM2.5 medio diario >=25 ug/m3.")
  if (grepl("^pm25_exposicion_ge25_", v)) return("Indicador binario de al menos un dia con PM2.5 >=25 ug/m3.")
  if (grepl("^n_dias_pm25_", v) && grepl("buena|regular|desfavorable", v)) return("Numero de dias de la ventana clasificados en la categoria indicada de calidad del aire para PM2.5.")
  if (grepl("^pm25_exposicion_desfavorable_o_peor_", v)) return("Indicador binario de al menos un dia PM2.5 en categoria Desfavorable, Muy desfavorable o Extremadamente desfavorable.")
  if (grepl("^n_dias_pm25_desfavorable_o_peor_", v)) return("Numero de dias PM2.5 en categoria Desfavorable, Muy desfavorable o Extremadamente desfavorable.")
  if (grepl("^peor_calidad_pm25_", v)) return("Peor categoria diaria de calidad del aire por PM2.5 observada en la ventana.")
  if (grepl("^peor_nivel_calidad_pm25_", v)) return("Nivel ordinal de la peor categoria PM2.5: 1=Buena a 6=Extremadamente_desfavorable.")

  # NO2 y O3.
  for (cont in c("no2", "o3")) {
    nombre <- ifelse(cont == "no2", "NO2", "O3")
    if (grepl(paste0("^", cont, "_media_"), v)) return(paste0("Media de las concentraciones medias diarias de ", nombre, " dentro de la ventana."))
    if (grepl(paste0("^", cont, "_maximo_"), v)) return(paste0("Mayor media diaria de ", nombre, " observada dentro de la ventana."))
    if (grepl(paste0("^", cont, "_acumulado_"), v)) return(paste0("Suma de las medias diarias de ", nombre, " dentro de la ventana."))
    if (grepl(paste0("^", cont, "_sd_"), v)) return(paste0("Desviacion estandar de las medias diarias de ", nombre, " dentro de la ventana."))
    if (grepl(paste0("^", cont, "_p90_"), v)) return(paste0("Percentil 90 de las medias diarias de ", nombre, " dentro de la ventana."))
    if (grepl(paste0("^", cont, "_auc"), v)) {
      u <- ifelse(cont == "no2", "120", "130")
      return(paste0("Carga acumulada de ", nombre, " por encima de ", u, " ug/m3: suma de max(valor diario - umbral, 0)."))
    }
    if (grepl(paste0("^n_dias_", cont, "_"), v) && grepl("buena|regular|desfavorable", v)) return(paste0("Numero de dias de la ventana en la categoria indicada de calidad del aire para ", nombre, "."))
    if (grepl(paste0("^", cont, "_exposicion_desfavorable_o_peor_"), v)) return(paste0("Indicador binario de al menos un dia ", nombre, " en categoria Desfavorable, Muy desfavorable o Extremadamente desfavorable."))
    if (grepl(paste0("^n_dias_", cont, "_desfavorable_o_peor_"), v)) return(paste0("Numero de dias ", nombre, " en categoria Desfavorable, Muy desfavorable o Extremadamente desfavorable."))
    if (grepl(paste0("^peor_calidad_", cont, "_"), v)) return(paste0("Peor categoria diaria de calidad del aire por ", nombre, " observada en la ventana."))
    if (grepl(paste0("^peor_nivel_calidad_", cont, "_"), v)) return(paste0("Nivel ordinal de la peor categoria ", nombre, ": 1=Buena a 6=Extremadamente_desfavorable."))
  }

  # Episodios combinados.
  if (grepl("^n_dias_episodio_particulas_cualquiera_", v)) return("Numero de dias con PM10 >=50 ug/m3 o PM2.5 >=25 ug/m3.")
  if (grepl("^episodio_particulas_cualquiera_exposicion_", v)) return("Indicador binario de al menos un dia con PM10 >=50 ug/m3 o PM2.5 >=25 ug/m3.")
  if (grepl("^n_dias_episodio_particulas_concordante_", v)) return("Numero de dias con PM10 >=50 ug/m3 y PM2.5 >=25 ug/m3 simultaneamente.")
  if (grepl("^episodio_particulas_concordante_exposicion_", v)) return("Indicador binario de al menos un dia con PM10 >=50 ug/m3 y PM2.5 >=25 ug/m3 simultaneamente.")

  # Cobertura y trazabilidad.
  if (grepl("^fecha_inicio_", v)) return("Primer dia incluido en la ventana ambiental.")
  if (grepl("^fecha_fin_", v)) return("Ultimo dia incluido en la ventana ambiental (D-1).")
  if (grepl("^n_dias_esperados_", v)) return("Numero total de dias teoricos de la ventana.")
  if (grepl("^n_dias_validos_", v)) return("Numero de dias con datos validos para el contaminante o panel correspondiente (>=18 horas validas).")
  if (grepl("^n_dias_no_validos_", v)) return("Numero de dias sin datos validos para el contaminante o panel correspondiente.")
  if (grepl("^porcentaje_dias_validos_", v)) return("Porcentaje de dias validos dentro de la ventana.")
  if (grepl("^ventana_.*_valida_", v)) return("Indicador de cobertura minima: >=7/10 dias en ventana de 10 dias o >=153/180 dias en ventana de 180 dias.")
  if (grepl("^estado_ventana_", v)) return("Estado legible de validez de la ventana ambiental segun el criterio de cobertura.")
  if (grepl("^n_dias_.*_(primaria|secundaria)_", v)) return("Numero de dias en los que se utilizo la estacion primaria o secundaria de la jerarquia correspondiente.")

  # Variables diarias auxiliares.
  if (grepl("^(sdi_pm10_ge50|pm25_ge25|episodio_particulas|.*_calidad_diaria)_d", v)) {
    return("Indicador o categoria diaria correspondiente al lag especificado (D-1 a D-10).")
  }

  "Variable derivada o conservada en la macrotabla; consultar script de generacion."
}

grupo_variable <- function(v) {
  if (v %in% c("id_muestra", "fecha_muestra", "isla")) return("Identificacion")
  if (v %in% names(cohorte) && !grepl("fecha_muestra|isla", v)) return("Covariables_clinicas_y_tecnicas")
  if (grepl("anio_muestra|mes_muestra|estacion_muestra", v)) return("Temporalidad_muestra")
  if (grepl("^pm10|^sdi_|^grado_sdi|^patron_temporal_sdi|^tipo_sdi|^n_dias_sdi", v)) return("PM10_y_SDI_principal")
  if (grepl("^pm25|^episodio_particulas", v)) return("PM2.5_y_episodios_particulas")
  if (grepl("^no2", v)) return("NO2_exploratorio")
  if (grepl("^o3", v)) return("O3_exploratorio")
  if (grepl("ventana|dias_validos|dias_no_validos|porcentaje_dias_validos|estado_ventana", v)) return("Cobertura_y_calidad_dato")
  if (grepl("primaria|secundaria|estacion_origen|nivel_estacion|prioridad_estacion", v)) return("Trazabilidad_estacion")
  "Otros"
}

diccionario <- tibble(
  orden = seq_along(names(macrotabla)),
  variable = names(macrotabla),
  tipo_R = vapply(macrotabla, function(x) class(x)[1], character(1)),
  grupo = vapply(names(macrotabla), grupo_variable, character(1)),
  ventana = vapply(names(macrotabla), ventana_variable, character(1)),
  unidad_o_codificacion = vapply(names(macrotabla), unidad_variable, character(1)),
  definicion_operativa = vapply(names(macrotabla), descripcion_variable, character(1)),
  origen = vapply(names(macrotabla), origen_variable, character(1)),
  observaciones = ifelse(
    grepl("sdi|episodio_particulas", names(macrotabla)),
    "Variable operativa de exposicion compatible con calima/SDI; no equivale por si sola a confirmacion meteorologica de intrusión sahariana.",
    ""
  )
)

# 15. EXPORTAR RESULTADOS

ruta_macro_csv <- file.path(
  ruta_outputs, "17_macrotabla_clinica_ambiental_por_muestra.csv"
)
ruta_diccionario_csv <- file.path(
  ruta_outputs, "21_diccionario_macrotabla_completo.csv"
)
ruta_cobertura_csv <- file.path(
  ruta_outputs, "18_resumen_cobertura_multicontaminante.csv"
)
ruta_jerarquia_csv <- file.path(
  ruta_outputs, "22_jerarquia_estaciones_contaminantes.csv"
)
ruta_definiciones_csv <- file.path(
  ruta_outputs, "19_definiciones_calidad_aire.csv"
)

readr::write_csv2(macrotabla, ruta_macro_csv, na = "")
readr::write_csv2(diccionario, ruta_diccionario_csv, na = "")
readr::write_csv2(resumen_cobertura, ruta_cobertura_csv, na = "")
readr::write_csv2(jerarquia_estaciones, ruta_jerarquia_csv, na = "")
readr::write_csv2(definiciones_calidad, ruta_definiciones_csv, na = "")

# Copias para la carpeta de entrega a Rocio.
readr::write_csv2(
  macrotabla,
  file.path(ruta_rocio, "01_macrotabla_GEMAS_microbiota_ambiente.csv"),
  na = ""
)
readr::write_csv2(
  diccionario,
  file.path(ruta_rocio, "02_diccionario_variables_GEMAS.csv"),
  na = ""
)
readr::write_csv2(
  definiciones_calidad,
  file.path(ruta_rocio, "03_definiciones_calidad_aire.csv"),
  na = ""
)
readr::write_csv2(
  jerarquia_estaciones,
  file.path(ruta_rocio, "04_jerarquia_estaciones_contaminantes.csv"),
  na = ""
)
readr::write_csv2(
  resumen_cobertura,
  file.path(ruta_rocio, "05_resumen_cobertura_multicontaminante.csv"),
  na = ""
)

writexl::write_xlsx(
  list(
    Macrotabla = macrotabla,
    Diccionario = diccionario,
    Cobertura = resumen_cobertura,
    Calidad_aire = definiciones_calidad,
    Jerarquia_estaciones = jerarquia_estaciones
  ),
  file.path(ruta_rocio, "01_paquete_GEMAS_ambiente_y_diccionario.xlsx")
)

if (file.exists(ruta_phyloseq)) {
  file.copy(
    from = ruta_phyloseq,
    to = file.path(ruta_rocio, basename(ruta_phyloseq)),
    overwrite = TRUE
  )
}

readme <- c(
  "PAQUETE GEMAS CANARIAS: MICROBIOTA Y EXPOSICION AMBIENTAL",
  "",
  paste0("Fecha de generacion: ", Sys.Date()),
  "",
  "CONTENIDO",
  "1. 01_macrotabla_GEMAS_microbiota_ambiente.csv",
  "   Una fila por muestra con covariables clinicas y exposiciones ambientales.",
  "2. 02_diccionario_variables_GEMAS.csv",
  "   Diccionario completo: definicion, ventana, unidad/codificacion y origen.",
  "3. 03_definiciones_calidad_aire.csv",
  "   Categorias y umbrales aplicados a las medias diarias.",
  "4. 04_jerarquia_estaciones_contaminantes.csv",
  "   Estaciones primaria/secundaria usadas por contaminante e isla.",
  "5. 05_resumen_cobertura_multicontaminante.csv",
  "   Completitud de ventanas por isla y contaminante.",
  "6. data_phyloseq_salivas_GEMAS_Canarias.rds",
  "   Objeto phyloseq, cuando existe en data_raw/.",
  "",
  "VENTANAS",
  "- Principal: D-10 a D-1.",
  "- Secundaria/exploratoria: D-180 a D-1.",
  "- D0 (dia de la muestra) queda excluido.",
  "",
  "VARIABLE PRINCIPAL SDI / EXPOSICION COMPATIBLE CON CALIMA",
  "- PM10 medio diario >=50 ug/m3.",
  "- Es una definicion operativa del TFM basada en PM10; no confirma por si sola origen sahariano meteorologico.",
  "",
  "CONFIDENCIALIDAD",
  "- Contiene IDs pseudonimizados, fechas y variables clinicas.",
  "- No subir a GitHub ni a repositorios publicos."
)

writeLines(readme, file.path(ruta_rocio, "README_entrega_Rocio.txt"))


# 16. MENSAJE FINAL

cat("\n=================================================\n")
cat("MACROTABLA, DICCIONARIO Y PAQUETE PARA ROCIO CREADOS\n")
cat("=================================================\n")
cat("Dimensiones de la macrotabla: ", nrow(macrotabla), " x ", ncol(macrotabla), "\n", sep = "")
cat("\nCobertura por isla:\n")
print(resumen_cobertura, width = Inf)
cat("\nCarpeta de entrega: ", ruta_rocio, "\n", sep = "")
cat("\nSiguiente paso: ejecutar 05_cribado_exposiciones_ambientales_GEMAS.R\n")
