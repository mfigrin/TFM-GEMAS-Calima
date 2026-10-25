# Síntesis ANCOM-BC2 frente a DESeq2.

ancom <- leer_csv(ruta_ancom)

deseq_lista <- lapply(
  names(rutas_deseq),
  function(nm) {
    leer_csv(rutas_deseq[[nm]]) %>%
      mutate(analisis = nm)
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
      clasificacion = clasificar_senal(q_ANCOM, robusto_ANCOM, FDR_DESeq2, direccion_concordante),
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

comparacion_metodos_completa <- bind_rows(lapply(names(rutas_deseq), integrar_metodos))

comparacion_metodos_seleccion <- comparacion_metodos_completa %>%
  filter(seleccionado) %>%
  arrange(
    factor(analisis, levels = names(rutas_deseq)),
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
    n_concordantes_entre_comunes = sum(direccion_concordante %in% TRUE, na.rm = TRUE),
    n_comunes_con_direccion = sum(!is.na(direccion_concordante)),
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
  paste0(nrow(comparacion_metodos_seleccion), " filas seleccionadas mediante reglas predefinidas.")
)

# Composición taxonómica y microbioma núcleo.

composicion_genus <- leer_csv(ruta_comp_genus)
prevalencia_genus <- leer_csv(ruta_prev_genus)

top15_genus <- composicion_genus %>%
  arrange(desc(abundancia_relativa_media)) %>%
  slice_head(n = 15) %>%
  mutate(taxon = factor(taxon, levels = rev(taxon)))

prevalencia_plot <- prevalencia_genus %>%
  mutate(
    core = case_when(
      prevalencia_prop >= 0.95 ~ "Núcleo ≥95 %",
      prevalencia_prop >= 0.90 ~ "Núcleo 90–<95 %",
      TRUE ~ "Fuera del núcleo"
    ),
    core = factor(core, levels = c("Fuera del núcleo", "Núcleo 90–<95 %", "Núcleo ≥95 %")),
    etiquetar = taxon %in% c(
      head(composicion_genus %>% arrange(desc(abundancia_relativa_media)) %>% pull(taxon), 8),
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

# Figuras.

# Diversidad ajustada por sexo.
predicciones_plot <- predicciones_sexo %>%
  mutate(
    sexo = factor(sexo, levels = c("Mujeres", "Hombres")),
    exposicion = factor(exposicion, levels = c("Sin intrusión oficial", "Con intrusión oficial"))
  )

crear_fig_sexo <- function(outcome_name, ylab) {
  ggplot(
    predicciones_plot %>% filter(.data$outcome == outcome_name),
    aes(x = exposicion, y = estimacion_ajustada, color = sexo, group = sexo)
  ) +
    geom_line(linewidth = 0.75, position = position_dodge(width = 0.12)) +
    geom_point(size = 2.8, position = position_dodge(width = 0.12)) +
    geom_errorbar(
      aes(ymin = IC95_inf, ymax = IC95_sup),
      width = 0.08,
      linewidth = 0.65,
      position = position_dodge(width = 0.12)
    ) +
    scale_color_manual(values = c("Mujeres" = COLOR_GRIS, "Hombres" = COLOR_ROJO)) +
    labs(title = outcome_name, x = NULL, y = ylab, color = "Sexo") +
    tema_tfm() +
    theme(axis.text.x = element_text(angle = 12, hjust = 1))
}

fig_sexo_shannon <- crear_fig_sexo("Shannon", "Índice de Shannon ajustado")

fig_sexo_simpson <- crear_fig_sexo("Simpson", "Índice de Simpson ajustado")

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

# Tabla integrada de indicadores ambientales.

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
beta_miteco_existente <- bind_rows(
  beta_miteco_existente,
  beta_final %>% transmute(distancia, exposure, n, R2, F, p, formula,
                          p_BH_distancias = q_BH_distancias))

etiquetas_exposicion <- c(
  miteco10 = "Intrusión oficial D-1 a D-10",
  miteco_days10_z = "Número de días de intrusión D-1 a D-10, por 1 DE",
  naf10_z = "Aporte natural africano D-1 a D-10, por 1 DE",
  pm10_media_10d_z = "PM10 superficial media D-1 a D-10, por 1 DE",
  pm10_auc50_10d_z = "Carga superficial sobre 50 µg/m³ D-1 a D-10, por 1 DE"
)

alpha_integrada <- bind_rows(
  alpha_miteco_existente %>%
    select(outcome, exposure, n, beta, IC95_inf, IC95_sup, p, q = p_BH_outcomes),
  alpha_pm10_10d %>%
    select(outcome, exposure, n, beta, IC95_inf, IC95_sup, p, q = q_BH_outcomes)
) %>%
  mutate(
    indicador = unname(etiquetas_exposicion[exposure]),
    tipo_resultado = "Diversidad alfa",
    metrica = outcome
  )

beta_integrada <- bind_rows(
  beta_miteco_existente %>%
    select(distancia, exposure, n, R2, F, p, q = p_BH_distancias),
  beta_exposiciones %>%
    filter(exposure %in% c("pm10_media_10d_z", "pm10_auc50_10d_z")) %>%
    select(distancia, exposure, n, R2, F, p, q = q_BH_distancias)
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

