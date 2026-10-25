# Caracterización ambiental.

ventanas <- tibble(
  ventana = factor(
    c("D-1 a D-10", "D-1 a D-30", "D-1 a D-90", "D-1 a D-180"),
    levels = c("D-1 a D-10", "D-1 a D-30", "D-1 a D-90", "D-1 a D-180")
  ),
  variable_any = c("miteco10", "miteco_any_30d", "miteco_any_90d", "miteco_any_180d"),
  variable_days = c("miteco_days_10d", "miteco_days_30d", "miteco_days_90d", "miteco_days_180d")
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
  mutate(porcentaje = 100 * n / sum(n)) %>%
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
  stop("Alguna exposición ambiental no es constante dentro de un cluster isla-fecha.")
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
  mutate(island = factor(island, levels = c("La Palma", "Tenerife")))

calcular_correlacion <- function(data, x, y, estrato = "Total") {
  dd <- data %>%
    filter(!is.na(.data[[x]]), !is.na(.data[[y]]))

  if (nrow(dd) < 4) {
    return(tibble(estrato = estrato, x = x, y = y, n = nrow(dd), rho_spearman = NA_real_, p = NA_real_))
  }

  ct <- suppressWarnings(cor.test(dd[[x]], dd[[y]], method = "spearman", exact = FALSE))

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
  calcular_correlacion(env_cluster, "miteco_days_10d", "pm10_media_10d", "Total"),
  calcular_correlacion(env_cluster, "miteco_days_10d", "pm10_auc50_10d", "Total"),
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

# Distancias microbiológicas.

ps_rel <- transform_sample_counts(ps_amb, function(x) x / sum(x))

dist_bray <- phyloseq::distance(ps_rel, method = "bray")

distancias <- list(Bray_Curtis = dist_bray)

if (!is.na(ruta_wunifrac) && file.exists(ruta_wunifrac)) {
  distancias$Weighted_UniFrac <- readRDS(ruta_wunifrac)
}

if (!is.na(ruta_uunifrac) && file.exists(ruta_uunifrac)) {
  distancias$Unweighted_UniFrac_rarefactado <- readRDS(ruta_uunifrac)
}

for (nm in names(distancias)) {
  n_comunes <- length(intersect(labels(distancias[[nm]]), meta$id_muestra))

  if (n_comunes < 200) {
    warning(
      "La distancia ", nm, " solo contiene ", n_comunes,
      " muestras comunes con la tabla maestra."
    )
  }
}

registrar("Distancias microbiológicas", "OK", paste(names(distancias), collapse = ", "))

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

write_csv(qc_distancias, file.path(TAB_DIR, "00_QC_alineacion_distancias_metadata.csv"), na = "")

msg("QC de alineación de distancias:")
print(qc_distancias)

