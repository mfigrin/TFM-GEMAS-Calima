# Asociación con el fenotipo exacerbador retrospectivo.
# Se comparan los ajustes con y sin antibióticos; no se modelan episodios futuros.

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
    cov_exac <- c(cov_exac, "antibiotic")
  }

  vars <- unique(c(outcome, "exac6m", cov_exac, "cluster_isla_fecha"))

  dd <- metricas[
    complete.cases(metricas[ , vars, drop = FALSE ]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  cov2 <- covars_activas(dd, cov_exac)

  fit <- lm(as.formula(paste(outcome, "~ exac6m +", paste(cov2, collapse = " + "))), data = dd)

  term <- names(coef(fit))

  term <- term[
    grepl("^exac6m", term)
  ]

  if (
    length(term) != 1
  ) {
    return(NULL)
  }

  z <- extraer_coef_robusto(fit, dd$cluster_isla_fecha, term)

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    mutate(
      outcome = outcome,
      modelo = ifelse(incluir_antibiotico, "Con_antibiotico", "Sin_antibiotico"),
      n = nrow(dd),
      n_clusters = n_distinct(dd$cluster_isla_fecha),
      .before = 1
    )
}

exac_alpha <- bind_rows(
  lapply(
    alpha_outcomes,
    function(o) {
      bind_rows(ajustar_exac_alpha(o, TRUE), ajustar_exac_alpha(o, FALSE))
    }
  )
) %>%
  group_by(modelo) %>%
  mutate(p_BH_outcomes = p.adjust(p, method = "BH")) %>%
  ungroup()

# Carga ambiental miteco180 -> exacerbacion 6m..
ajustar_exac_outcome <- function(
  incluir_antibiotico = TRUE
) {
  cov_exac <- c("age", "sex", "island", "pool", "cavities", "gum", "year", "sin_doy", "cos_doy")

  if (
    incluir_antibiotico
  ) {
    cov_exac <- c(cov_exac, "antibiotic")
  }

  vars <- unique(c("exac6m", "miteco_days180_z", cov_exac, "cluster_isla_fecha"))

  dd <- metricas[
    complete.cases(metricas[ , vars, drop = FALSE ]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  cov2 <- covars_activas(dd, cov_exac)

  fit <- glm(
    as.formula(paste("exac6m ~ miteco_days180_z +", paste(cov2, collapse = " + "))),
    family = binomial(),
    data = dd
  )

  z <- extraer_coef_robusto(fit, dd$cluster_isla_fecha, "miteco_days180_z")

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    transmute(
      modelo = ifelse(incluir_antibiotico, "Con_antibiotico", "Sin_antibiotico"),
      n = nrow(dd),
      n_clusters = n_distinct(dd$cluster_isla_fecha),
      beta = beta,
      se = se,
      OR_por_1DE = exp(beta),
      OR_IC95_inf = exp(IC95_inf),
      OR_IC95_sup = exp(IC95_sup),
      p = p
    )
}

exac_miteco180 <- bind_rows(ajustar_exac_outcome(TRUE), ajustar_exac_outcome(FALSE))
