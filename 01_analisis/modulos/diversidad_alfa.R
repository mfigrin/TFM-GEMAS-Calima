# Diversidad alfa, Aitchison, interacciones y ventanas.
# Cada familia conserva su corrección de Benjamini–Hochberg.

alpha <- estimate_richness(ps_amb_f_m, measures = c("Observed", "Shannon", "Simpson"))

alpha$id_muestra <- rownames(alpha)

# Aitchison a centro de isla:
# misma definicion usada en el analisis exploratorio previo.
psg_ait <- tax_glom(ps_amb_u_m, taxrank = "Genus", NArm = TRUE)

otu_g <- otu_taxa_rows(psg_ait)

prev_g <- rowMeans(otu_g > 0)

keep_g <- names(prev_g)[
  prev_g >= 0.10
]

psg_ait <- prune_taxa(keep_g, psg_ait)

otu_g <- otu_taxa_rows(psg_ait)

x_clr <- t(otu_g + PSEUDOCOUNT_CLR)

clr <- log(x_clr) -
  rowMeans(log(x_clr))

meta_ait <- meta_base[
  match(rownames(clr), meta_base$id_muestra),
  ,
  drop = FALSE
]

centros_isla <- lapply(
  levels(droplevels(meta_ait$island)),
  function(isl) {
    ids <- meta_ait$id_muestra[
      meta_ait$island == isl
    ]

    colMeans(clr[ ids, , drop = FALSE ])
  }
)

names(centros_isla) <- levels(droplevels(meta_ait$island))

aitchison_isla <- vapply(
  meta_ait$id_muestra,
  function(id) {
    isl <- as.character(meta_ait$island[ meta_ait$id_muestra == id ][1])

    sqrt(sum((clr[ id, , drop = TRUE ] - centros_isla[[isl]])^2))
  },
  numeric(1)
)

metricas <- meta_base %>%
  left_join(alpha, by = "id_muestra") %>%
  mutate(aitchison_centro_isla = aitchison_isla[ id_muestra ])

alpha_outcomes <- c("Observed", "Shannon", "Simpson", "aitchison_centro_isla")

covars_alpha <- c(
  "age",
  "sex",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m",
  "lecturas_log10",
  "year",
  "sin_doy",
  "cos_doy"
)

ajustar_alpha_cont <- function(
  outcome,
  exposure,
  data = metricas,
  covars = covars_alpha
) {
  vars <- unique(c(outcome, exposure, covars, "cluster_isla_fecha"))

  dd <- data[
    complete.cases(data[ , vars, drop = FALSE ]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  cov2 <- covars_activas(dd, covars)

  if (
    n_niveles(dd[[exposure]]) < 2
  ) {
    return(NULL)
  }

  fit <- lm(as.formula(paste(outcome, "~", paste(c(exposure, cov2), collapse = " + "))), data = dd)

  z <- extraer_coef_robusto(fit, dd$cluster_isla_fecha, exposure)

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    mutate(
      outcome = outcome,
      exposure = exposure,
      n = nrow(dd),
      n_clusters = n_distinct(dd$cluster_isla_fecha),
      .before = 1
    )
}

exposiciones_alpha <- c("miteco10", "miteco_days10_z", "naf10_z")

alpha_global <- bind_rows(
  lapply(
    exposiciones_alpha,
    function(e) {
      bind_rows(
        lapply(
          alpha_outcomes,
          function(o) {
            ajustar_alpha_cont(o, e)
          }
        )
      )
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(p_BH_outcomes = p.adjust(p, method = "BH")) %>%
  ungroup()

# Alfa: grado / patron temporal.

ajustar_alpha_factor <- function(
  outcome,
  exposure,
  data = metricas,
  covars = covars_alpha
) {
  vars <- unique(c(outcome, exposure, covars, "cluster_isla_fecha"))

  dd <- data[
    complete.cases(data[ , vars, drop = FALSE ]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  if (
    n_niveles(dd[[exposure]]) < 2
  ) {
    return(NULL)
  }

  cov2 <- covars_activas(dd, covars)

  fit <- lm(as.formula(paste(outcome, "~", paste(c(exposure, cov2), collapse = " + "))), data = dd)

  terms <- names(coef(fit))

  terms <- terms[
    grepl(paste0("^", exposure), terms)
  ]

  w <- wald_terms(fit, dd$cluster_isla_fecha, terms)

  tibble(
    outcome = outcome,
    exposure = exposure,
    n = nrow(dd),
    n_clusters = n_distinct(dd$cluster_isla_fecha),
    chi2 = w$chi2,
    df = w$df,
    p = w$p
  )
}

alpha_factor <- bind_rows(
  lapply(
    c("grado_miteco10", "patron_miteco10"),
    function(e) {
      bind_rows(
        lapply(
          alpha_outcomes,
          function(o) {
            ajustar_alpha_factor(o, e)
          }
        )
      )
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(p_BH_outcomes = p.adjust(p, method = "BH")) %>%
  ungroup()

# Alfa: interacciones.

ajustar_alpha_interaccion <- function(
  outcome,
  exposure,
  modifier,
  data = metricas,
  covars = covars_alpha
) {
  cov2 <- setdiff(covars, modifier)

  vars <- unique(c(outcome, exposure, modifier, cov2, "cluster_isla_fecha"))

  dd <- data[
    complete.cases(data[ , vars, drop = FALSE ]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  if (
    n_niveles(dd[[exposure]]) < 2 ||
      n_niveles(dd[[modifier]]) < 2
  ) {
    return(NULL)
  }

  cov2 <- covars_activas(dd, cov2)

  fit <- lm(
    as.formula(
      paste(
        outcome,
        "~",
        paste0(exposure, " * ", modifier),
        if (
          length(cov2) > 0
        ) {
          paste0(" + ", paste(cov2, collapse = " + "))
        } else {
          ""
        }
      )
    ),
    data = dd
  )

  terms <- buscar_terminos_interaccion(fit, exposure, modifier)

  w <- wald_terms(fit, dd$cluster_isla_fecha, terms)

  tibble(
    outcome = outcome,
    exposure = exposure,
    modifier = modifier,
    n = nrow(dd),
    n_clusters = n_distinct(dd$cluster_isla_fecha),
    chi2 = w$chi2,
    df = w$df,
    p_interaccion = w$p
  )
}

alpha_int_sexo <- bind_rows(
  lapply(
    c("miteco10", "miteco_days10_z", "naf10_z"),
    function(e) {
      bind_rows(
        lapply(
          alpha_outcomes,
          function(o) {
            ajustar_alpha_interaccion(o, e, "sex")
          }
        )
      )
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(p_BH_outcomes = p.adjust(p_interaccion, method = "BH")) %>%
  ungroup()

alpha_int_isla <- bind_rows(
  lapply(
    alpha_outcomes,
    function(o) {
      ajustar_alpha_interaccion(o, "miteco10", "island")
    }
  )
) %>%
  mutate(p_BH_outcomes = p.adjust(p_interaccion, method = "BH"))

alpha_int_exac <- bind_rows(
  lapply(
    alpha_outcomes,
    function(o) {
      ajustar_alpha_interaccion(o, "miteco10", "exac6m")
    }
  )
) %>%
  mutate(p_BH_outcomes = p.adjust(p_interaccion, method = "BH"))

metricas_regimen <- metricas %>%
  filter(!is.na(regimen_strict)) %>%
  droplevels()

alpha_int_regimen <- bind_rows(
  lapply(
    alpha_outcomes,
    function(o) {
      ajustar_alpha_interaccion(o, "miteco10", "regimen_strict", data = metricas_regimen)
    }
  )
) %>%
  mutate(p_BH_outcomes = p.adjust(p_interaccion, method = "BH"))

# Alfa: estratificado por sexo para interpretacion.

alpha_sexo_estrat <- bind_rows(
  lapply(
    c("Female", "Male"),
    function(s) {
      d <- metricas %>%
        filter(sex == s) %>%
        droplevels()

      cov_s <- setdiff(covars_alpha, "sex")

      bind_rows(
        lapply(
          c("miteco10", "miteco_days10_z"),
          function(e) {
            bind_rows(
              lapply(
                alpha_outcomes,
                function(o) {
                  z <- ajustar_alpha_cont(o, e, data = d, covars = cov_s)

                  if (is.null(z)) {
                    return(NULL)
                  }

                  z %>%
                    mutate(sex = s, .before = 1)
                }
              )
            )
          }
        )
      )
    }
  )
)

# Ventanas temporales simultaneas.

windows <- c("recent_D1_10_z", "residual_D11_30_z", "prior_D31_90_z")

ajustar_alpha_tres_ventanas <- function(
  outcome
) {
  vars <- unique(c(outcome, windows, covars_alpha, "cluster_isla_fecha"))

  dd <- metricas[
    complete.cases(metricas[ , vars, drop = FALSE ]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  cov2 <- covars_activas(dd, covars_alpha)

  fit <- lm(as.formula(paste(outcome, "~", paste(c(windows, cov2), collapse = " + "))), data = dd)

  bind_rows(
    lapply(
      windows,
      function(w) {
        z <- extraer_coef_robusto(fit, dd$cluster_isla_fecha, w)

        if (is.null(z)) {
          return(NULL)
        }

        z %>%
          mutate(
            outcome = outcome,
            window = w,
            n = nrow(dd),
            n_clusters = n_distinct(dd$cluster_isla_fecha),
            .before = 1
          )
      }
    )
  )
}

alpha_ventanas <- bind_rows(lapply(alpha_outcomes, ajustar_alpha_tres_ventanas)) %>%
  group_by(window) %>%
  mutate(p_BH_outcomes = p.adjust(p, method = "BH")) %>%
  ungroup()

ajustar_alpha_ventanas_sexo <- function(
  outcome
) {
  cov2 <- setdiff(covars_alpha, "sex")

  vars <- unique(c(outcome, windows, "sex", cov2, "cluster_isla_fecha"))

  dd <- metricas[
    complete.cases(metricas[ , vars, drop = FALSE ]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  cov2 <- covars_activas(dd, cov2)

  fit <- lm(
    as.formula(
      paste(
        outcome,
        "~ (",
        paste(windows, collapse = " + "),
        ") * sex",
        if (
          length(cov2) > 0
        ) {
          paste0(" + ", paste(cov2, collapse = " + "))
        } else {
          ""
        }
      )
    ),
    data = dd
  )

  bind_rows(
    lapply(
      windows,
      function(w) {
        terms <- buscar_terminos_interaccion(fit, w, "sex")

        wt <- wald_terms(fit, dd$cluster_isla_fecha, terms)

        tibble(
          outcome = outcome,
          window = w,
          n = nrow(dd),
          n_clusters = n_distinct(dd$cluster_isla_fecha),
          chi2 = wt$chi2,
          df = wt$df,
          p_interaccion = wt$p
        )
      }
    )
  )
}

alpha_ventanas_sexo <- bind_rows(lapply(alpha_outcomes, ajustar_alpha_ventanas_sexo)) %>%
  group_by(window) %>%
  mutate(p_BH_outcomes = p.adjust(p_interaccion, method = "BH")) %>%
  ungroup()

# Dosis-respuesta no lineal.

analizar_spline <- function(
  outcome,
  exposure_raw,
  data = metricas
) {
  vars <- unique(c(outcome, exposure_raw, covars_alpha))

  dd <- data[
    complete.cases(data[ , vars, drop = FALSE ]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  cov2 <- covars_activas(dd, covars_alpha)

  rhs_cov <- if (
    length(cov2) > 0
  ) {
    paste(cov2, collapse = " + ")
  } else {
    "1"
  }

  m0 <- lm(as.formula(paste(outcome, "~", rhs_cov)), data = dd)

  m1 <- lm(as.formula(paste(outcome, "~", exposure_raw, "+", rhs_cov)), data = dd)

  m2 <- lm(as.formula(paste(outcome, "~ splines::ns(", exposure_raw, ", df=3) +", rhs_cov)), data = dd)

  a01 <- anova(m0, m1)

  a12 <- anova(m1, m2)

  tibble(
    outcome = outcome,
    exposure = exposure_raw,
    n = nrow(dd),
    p_global_lineal = a01$`Pr(>F)`[2],
    p_no_lineal_sobre_lineal = a12$`Pr(>F)`[2],
    AIC_lineal = AIC(m1),
    AIC_spline = AIC(m2)
  )
}

spline_res <- bind_rows(
  lapply(
    c("miteco_days_10d", "naf_mean_calendar_10d"),
    function(e) {
      bind_rows(
        lapply(
          c("Observed", "Shannon", "Simpson"),
          function(o) {
            analizar_spline(o, e)
          }
        )
      )
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(q_no_lineal = p.adjust(p_no_lineal_sobre_lineal, method = "BH")) %>%
  ungroup()
