# Mantener la secuencia aleatoria del código histórico sin repetir sus ajustes.
# Este módulo no produce resultados. La prueba principal y PERMDISP se calculan en 06.
# La equivalencia del estado aleatorio se debe comprobar con 01_analisis/utilidades/verificar_aleatoriedad.R.

betadisper_test <- function(
  dist_obj,
  distancia,
  grupo,
  data = meta_base
) {
  ids <- labels(dist_obj)

  dd <- data[
    match(ids, data$id_muestra),
    ,
    drop = FALSE
  ]

  rownames(dd) <- dd$id_muestra

  cc <- !is.na(dd[[grupo]])

  dd <- droplevels(dd[ cc, , drop = FALSE ])

  d2 <- subdist(dist_obj, rownames(dd))

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  g <- factor(dd[[grupo]])

  if (
    nlevels(droplevels(g)) < 2
  ) {
    return(NULL)
  }

  tfm_avanzar_permutaciones(N_PERM, nrow(dd))
  invisible(NULL)
}

meta_base$miteco_sex_group <- factor(
  ifelse(
    is.na(meta_base$miteco10) |
      is.na(meta_base$sex),
    NA_character_,
    paste0(as.character(meta_base$sex), "_", ifelse(meta_base$miteco10 == 1, "MITECO+", "MITECO-"))
  )
)

for (nm in names(distancias)) {
  betadisper_test(distancias[[nm]], nm, "miteco10_factor")
  betadisper_test(distancias[[nm]], nm, "miteco_sex_group")
}

# Beta: sensibilidad cluster isla-fecha, bray.

otu_rel <- otu_taxa_rows(ps_rel)

clusters <- as.character(meta_base$cluster_isla_fecha)

names(clusters) <- meta_base$id_muestra

ids_otu <- colnames(otu_rel)

clusters_otu <- clusters[
  ids_otu
]

cluster_levels <- unique(clusters_otu)

qc_cluster_constancia <- meta_base %>%
  group_by(cluster_isla_fecha) %>%
  summarise(
    n_miteco10 = n_distinct(miteco10[ !is.na(miteco10) ]),
    n_days10 = n_distinct(miteco_days_10d[ !is.na(miteco_days_10d) ]),
    .groups = "drop"
  )

if (
  any(qc_cluster_constancia$n_miteco10 > 1, na.rm = TRUE) ||
    any(qc_cluster_constancia$n_days10 > 1, na.rm = TRUE)
) {
  stop("La exposicion MITECO no es constante dentro de algun cluster isla-fecha.")
}

meta_cluster <- meta_base %>%
  group_by(cluster_isla_fecha) %>%
  summarise(
    n_muestras = n(),

    island = first(as.character(island)),

    fecha = first(fecha),

    miteco10 = first(miteco10),

    miteco_days10_z = first(miteco_days10_z),

    naf10_z = first(naf10_z),

    recent_D1_10_z = first(recent_D1_10_z),

    residual_D11_30_z = first(residual_D11_30_z),

    prior_D31_90_z = first(prior_D31_90_z),

    year = first(as.character(year)),

    sin_doy = first(sin_doy),

    cos_doy = first(cos_doy),

    age_mean = mean(age, na.rm = TRUE),

    prop_male = mean(sex_num, na.rm = TRUE),

    prop_antibiotic = mean(as.numeric(as.character(antibiotic)), na.rm = TRUE),

    prop_cavities = mean(as.numeric(as.character(cavities)), na.rm = TRUE),

    prop_gum = mean(as.numeric(as.character(gum)), na.rm = TRUE),

    prop_exac6m = mean(exac6m_num, na.rm = TRUE),

    prop_pool2 = mean(as.character(pool) == "2", na.rm = TRUE),

    prop_pool3 = mean(as.character(pool) == "3", na.rm = TRUE),

    prop_pool6 = mean(as.character(pool) == "6", na.rm = TRUE),

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

permanova_cluster <- function(
  exposure
) {
  ids <- cluster_levels

  dd <- meta_cluster[
    ids,
    ,
    drop = FALSE
  ]

  vars <- unique(c(exposure, covars_cluster))

  cc <- complete.cases(dd[ , vars, drop = FALSE ])

  dd <- droplevels(dd[ cc, , drop = FALSE ])

  cov2 <- covars_activas(dd, covars_cluster)

  ftxt <- paste("d2 ~", paste(c(exposure, cov2), collapse = " + "))

  tfm_avanzar_permutaciones(N_PERM, nrow(dd))
  invisible(NULL)
}

invisible(lapply(c("miteco10", "miteco_days10_z", "naf10_z"), permanova_cluster))
