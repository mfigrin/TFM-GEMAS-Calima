# Integrar MITECO con los objetos analíticos.
# Se comprueban fechas e identificadores antes de derivar las ventanas anteriores a D0.

ps_amb_f <- readRDS(ruta_ps_amb_f)

ps_amb_u <- readRDS(ruta_ps_amb_u)

ps_saliva_u <- readRDS(ruta_ps_saliva_u)

miteco <- readr::read_csv(ruta_miteco_muestra, show_col_types = FALSE) %>%
  mutate(fecha_muestra = as.Date(fecha_muestra))

traza <- readr::read_csv(ruta_miteco_traza, show_col_types = FALSE) %>%
  mutate(
    fecha = as.Date(fecha),
    fecha_muestra = as.Date(fecha_muestra),
    miteco_africano = as.integer(tidyr::replace_na(miteco_africano, 0))
  )

# Ventanas temporales oficiales MITECO.

ventanas <- traza %>%
  group_by(id_muestra) %>%
  summarise(
    miteco_days_D1_3 = sum(miteco_africano[ between(lag, 1, 3) ], na.rm = TRUE),

    miteco_days_D4_10 = sum(miteco_africano[ between(lag, 4, 10) ], na.rm = TRUE),

    miteco_days_D1_10_check = sum(miteco_africano[ between(lag, 1, 10) ], na.rm = TRUE),

    miteco_days_D11_30 = sum(miteco_africano[ between(lag, 11, 30) ], na.rm = TRUE),

    miteco_days_D31_90 = sum(miteco_africano[ between(lag, 31, 90) ], na.rm = TRUE),

    miteco_days_D1_30 = sum(miteco_africano[ between(lag, 1, 30) ], na.rm = TRUE),

    miteco_days_D1_90 = sum(miteco_africano[ between(lag, 1, 90) ], na.rm = TRUE),

    miteco_recency_D1_10 = ifelse(
      any(miteco_africano[ between(lag, 1, 10) ] == 1, na.rm = TRUE),
      min(lag[ between(lag, 1, 10) & miteco_africano == 1 ], na.rm = TRUE),
      NA_real_
    ),

    .groups = "drop"
  )

miteco <- miteco %>%
  left_join(ventanas, by = "id_muestra") %>%
  mutate(
# Auditoria: la suma diaria debe reproducir script 30.
    qc_dias10_coincide = miteco_days_D1_10_check ==
      miteco_days_10d,

    miteco10 = as.numeric(miteco_any_10d),

    miteco10_factor = factor(
      ifelse(miteco10 == 1, "MITECO_si", "MITECO_no"),
      levels = c("MITECO_no", "MITECO_si")
    ),

    miteco_days10_z = zscore(miteco_days_10d),

    naf10_z = zscore(naf_mean_calendar_10d),

    miteco_days180_z = zscore(miteco_days_180d),

    recent_D1_10_z = zscore(miteco_days_10d),

    residual_D11_30_z = zscore(miteco_days_D11_30),

    prior_D31_90_z = zscore(miteco_days_D31_90),

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
      levels = c("0_sin_MITECO", "1_3_dias", "4plus_dias")
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
      levels = c("Sin_MITECO_10d", "Solo_D1_3", "Solo_D4_10", "D1_3_y_D4_10")
    )
  )

if (!all(
  miteco$qc_dias10_coincide,
  na.rm = TRUE
)) {
  stop("La reconstruccion MITECO D1-D10 no reproduce miteco_days_10d del script 30.")
}

registrar(
  "Integracion MITECO",
  "OK",
  paste0("n=", nrow(miteco), "; MITECO10+=", sum(miteco$miteco10 == 1, na.rm = TRUE))
)

# Metadata analitica.

meta_base <- sample_data_df(ps_amb_f)

# El objeto phyloseq analitico del script 07 YA contiene una columna
# fecha_muestra procedente de la macrotabla PM10. Al unir de nuevo la
# tabla MITECO, dplyr crearia fecha_muestra.x / fecha_muestra.y y
# desapareceria el nombre exacto fecha_muestra.
# Conservamos la fecha antigua solo para auditoria y dejamos que
# fecha_muestra sea la fecha de la integracion MITECO.
if (
  "fecha_muestra" %in% names(meta_base)
) {
  names(meta_base)[
    names(meta_base) == "fecha_muestra"
  ] <- "fecha_muestra_phyloseq"
}

meta_base <- meta_base %>%
  mutate(
    age = a_num(Age),

    sex_num = a_binario(Sex.letter),

    sex = factor(
      ifelse(sex_num == 1, "Male", ifelse(sex_num == 0, "Female", NA_character_)),
      levels = c("Female", "Male")
    ),

    island = factor(normalizar_isla(Island), levels = c("La Palma", "Tenerife")),

    pool = factor(Pool),

    antibiotic = factor(a_binario(Antibiotic.2mo), levels = c(0, 1)),

    cavities = factor(a_binario(Cavities), levels = c(0, 1)),

    gum = factor(a_binario(Chewing.gum.30min), levels = c(0, 1)),

    exac6m_num = a_num(Exacerbations.6M),

    exac6m_num = ifelse(exac6m_num %in% c(0, 1), exac6m_num, NA_real_),

    exac6m = factor(exac6m_num, levels = c(0, 1))
  ) %>%
  left_join(miteco, by = "id_muestra")

# QC DE LA UNION DE FECHAS

if (
  !"fecha_muestra" %in% names(meta_base)
) {
  stop(
    "Tras unir MITECO no existe fecha_muestra. ",
    "Columnas disponibles: ",
    paste(names(meta_base), collapse = ", ")
  )
}

meta_base$fecha_muestra <- as.Date(meta_base$fecha_muestra)

if (
  "fecha_muestra_phyloseq" %in% names(meta_base)
) {
  meta_base$fecha_muestra_phyloseq <- as.Date(meta_base$fecha_muestra_phyloseq)

  n_mismatch_fecha_phyloseq <- sum(
    !is.na(meta_base$fecha_muestra) &
      !is.na(meta_base$fecha_muestra_phyloseq) &
      meta_base$fecha_muestra !=
        meta_base$fecha_muestra_phyloseq
  )
} else {
  n_mismatch_fecha_phyloseq <- NA_integer_
}

fecha_original_phyloseq <- as.Date(meta_base$Date)

n_mismatch_fecha_original <- sum(
  !is.na(meta_base$fecha_muestra) &
    !is.na(fecha_original_phyloseq) &
    meta_base$fecha_muestra !=
      fecha_original_phyloseq
)

n_sin_fecha_miteco <- sum(is.na(meta_base$fecha_muestra))

qc_integracion_fechas <- tibble(
  n_muestras = nrow(meta_base),
  n_sin_fecha_MITECO = n_sin_fecha_miteco,
  discrepancias_fecha_vs_phyloseq_ambiental =
    n_mismatch_fecha_phyloseq,
  discrepancias_fecha_vs_Date_original =
    n_mismatch_fecha_original
)

print(qc_integracion_fechas)

if (
  n_sin_fecha_miteco > 0
) {
  stop("Hay ", n_sin_fecha_miteco, " muestras sin fecha tras la union MITECO.")
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
  !is.na(n_mismatch_fecha_phyloseq) &&
    n_mismatch_fecha_phyloseq > 0
) {
  warning(
    "La fecha MITECO discrepa de la antigua fecha_muestra ambiental en ",
    n_mismatch_fecha_phyloseq,
    " muestras, aunque coincide con Date original."
  )
}

meta_base$lecturas_totales <- as.numeric(sample_sums(ps_amb_f)[ meta_base$id_muestra ])

meta_base$lecturas_log10 <- log10(meta_base$lecturas_totales + 1)

meta_base <- meta_base %>%
  mutate(
    fecha = as.Date(fecha_muestra),

    year = factor(format(fecha, "%Y")),

    doy = as.integer(format(fecha, "%j")),

    sin_doy = sin(2 * pi * doy / 365.25),

    cos_doy = cos(2 * pi * doy / 365.25),

    cluster_isla_fecha = factor(paste(as.character(island), fecha, sep = "__")),

    regimen_strict = factor(
      case_when(
        format(fecha, "%m") %in% c("11", "12", "01") ~
          "Winter_proxy_NDJ",

        format(fecha, "%m") %in% c("07", "08") ~
          "Summer_proxy_JA",

        TRUE ~ NA_character_
      ),
      levels = c("Winter_proxy_NDJ", "Summer_proxy_JA")
    ),

    miteco_sex_int = miteco10 *
      sex_num,

    miteco_exac_int = miteco10 *
      exac6m_num
  )

# Actualizar phyloseq con metadata MITECO.

actualizar_ps <- function(
  ps_obj,
  meta
) {
  ids <- intersect(sample_names(ps_obj), meta$id_muestra)

  p <- prune_samples(ids, ps_obj)

  m <- meta %>%
    filter(id_muestra %in% ids) %>%
    as.data.frame(stringsAsFactors = FALSE)

  rownames(m) <- m$id_muestra

  m <- m[
    sample_names(p),
    ,
    drop = FALSE
  ]

  sample_data(p) <- sample_data(m)

  p
}

ps_amb_f_m <- actualizar_ps(ps_amb_f, meta_base)

ps_amb_u_m <- actualizar_ps(ps_amb_u, meta_base)

saveRDS(ps_amb_f_m, file.path(ruta_objetos, "ps_ambiental_205_filtrado_con_MITECO.rds"))

saveRDS(ps_amb_u_m, file.path(ruta_objetos, "ps_ambiental_205_sin_filtrar_con_MITECO.rds"))

# Qc descriptivo exposicion.

qc_exposicion <- meta_base %>%
  summarise(
    n = n(),

    miteco10_si = sum(miteco10 == 1, na.rm = TRUE),

    miteco10_no = sum(miteco10 == 0, na.rm = TRUE),

    miteco30_si = sum(miteco_days_D1_30 > 0, na.rm = TRUE),

    miteco30_no = sum(miteco_days_D1_30 == 0, na.rm = TRUE),

    miteco90_si = sum(miteco_days_D1_90 > 0, na.rm = TRUE),

    miteco90_no = sum(miteco_days_D1_90 == 0, na.rm = TRUE),

    mediana_dias_MITECO_D1_10 = median(miteco_days_10d, na.rm = TRUE),

    mediana_dias_MITECO_D11_30 = median(miteco_days_D11_30, na.rm = TRUE),

    mediana_dias_MITECO_D31_90 = median(miteco_days_D31_90, na.rm = TRUE)
  )

qc_grado <- meta_base %>%
  count(grado_miteco10, name = "n")

qc_patron <- meta_base %>%
  count(patron_miteco10, name = "n")

qc_sexo <- meta_base %>%
  count(sex, miteco10_factor, name = "n")

qc_isla <- meta_base %>%
  count(island, miteco10_factor, name = "n")

qc_regimen <- meta_base %>%
  filter(!is.na(regimen_strict)) %>%
  count(regimen_strict, miteco10_factor, name = "n")

write_csv(qc_integracion_fechas, file.path(ruta_tablas, "00_QC_integracion_fechas_MITECO.csv"), na = "")

write_csv(qc_exposicion, file.path(ruta_tablas, "00_QC_exposicion_MITECO.csv"), na = "")

write_csv(qc_grado, file.path(ruta_tablas, "00b_QC_grado_dias_MITECO10.csv"), na = "")

write_csv(qc_patron, file.path(ruta_tablas, "00c_QC_patron_temporal_MITECO10.csv"), na = "")

write_csv(qc_regimen, file.path(ruta_tablas, "00d_QC_winter_summer_MITECO.csv"), na = "")
