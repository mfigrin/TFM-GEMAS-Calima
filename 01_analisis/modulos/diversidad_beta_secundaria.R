# Contrastes beta secundarios: carga, interacciones y ventanas.
# Se mantienen 999 permutaciones; la prueba principal se calcula aparte con 9.999.
# El avance aleatorio de compatibilidad no ajusta ni presenta modelos descartados.

ps_rel <- transform_sample_counts(
  ps_amb_f_m,
  function(x) {
    x /
      sum(x)
  }
)

dist_bray <- phyloseq::distance(ps_rel, method = "bray")

distancias <- list(Bray_Curtis = dist_bray)

if (
  EJECUTAR_UNIFRAC_SI_EXISTE &&
    !is.na(ruta_wunifrac) &&
    file.exists(ruta_wunifrac)
) {
  distancias$Weighted_UniFrac <- readRDS(ruta_wunifrac)

  registrar("Weighted UniFrac", "OK", ruta_wunifrac)
} else {
  registrar("Weighted UniFrac", "OMITIDO", "No existe el objeto 08B.")
}

if (
  EJECUTAR_UNIFRAC_SI_EXISTE &&
    !is.na(ruta_uunifrac) &&
    file.exists(ruta_uunifrac)
) {
  distancias$Unweighted_UniFrac_rarefactado <- readRDS(ruta_uunifrac)

  registrar("Unweighted UniFrac", "OK", ruta_uunifrac)
} else {
  registrar("Unweighted UniFrac", "OMITIDO", "No existe el objeto 08B.")
}

# PERMANOVA general.

covars_beta <- c(
  "age",
  "sex",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m",
  "year",
  "sin_doy",
  "cos_doy"
)

permanova_cont <- function(
  dist_obj,
  distancia,
  exposure,
  data = meta_base,
  covars = covars_beta
) {
  ids <- intersect(labels(dist_obj), data$id_muestra)

  dd <- data[
    match(ids, data$id_muestra),
    ,
    drop = FALSE
  ]

  rownames(dd) <- dd$id_muestra

  vars <- unique(c(exposure, covars))

  cc <- complete.cases(dd[ , vars, drop = FALSE ])

  dd <- dd[
    cc,
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)

  d2 <- subdist(dist_obj, rownames(dd))

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  cov2 <- covars_activas(dd, covars)

  if (
    n_niveles(dd[[exposure]]) < 2
  ) {
    return(NULL)
  }

  ftxt <- paste("d2 ~", paste(c(exposure, cov2), collapse = " + "))

  if (exposure == "miteco10") {
    tfm_avanzar_permutaciones(N_PERM, nrow(dd))
    return(NULL)
  }

  fit <- vegan::adonis2(as.formula(ftxt), data = dd, permutations = N_PERM, by = "margin")

  tab <- as.data.frame(fit)

  tab$term <- rownames(tab)

  fila <- tab[
    tab$term == exposure,
    ,
    drop = FALSE
  ]

  if (
    nrow(fila) != 1
  ) {
    return(NULL)
  }

  tibble(
    distancia = distancia,
    exposure = exposure,
    n = nrow(dd),
    R2 = fila$R2,
    F = fila$F,
    p = fila$`Pr(>F)`,
    formula = ftxt
  )
}

beta_global <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      d <- distancias[[nm]]

      bind_rows(
        lapply(
          c("miteco10", "miteco_days10_z", "naf10_z"),
          function(e) {
            permanova_cont(d, nm, e)
          }
        )
      )
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(p_BH_distancias = p.adjust(p, method = "BH")) %>%
  ungroup()

# PERMANOVA interacciones.

permanova_interaccion <- function(
  dist_obj,
  distancia,
  exposure,
  modifier,
  data = meta_base,
  covars = covars_beta
) {
  ids <- intersect(labels(dist_obj), data$id_muestra)

  dd <- data[
    match(ids, data$id_muestra),
    ,
    drop = FALSE
  ]

  rownames(dd) <- dd$id_muestra

  cov2 <- setdiff(covars, modifier)

  vars <- unique(c(exposure, modifier, cov2))

  cc <- complete.cases(dd[ , vars, drop = FALSE ])

  dd <- droplevels(dd[ cc, , drop = FALSE ])

  d2 <- subdist(dist_obj, rownames(dd))

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  cov2 <- covars_activas(dd, cov2)

  if (
    n_niveles(dd[[exposure]]) < 2 ||
      n_niveles(dd[[modifier]]) < 2
  ) {
    return(NULL)
  }

  ftxt <- paste(
    "d2 ~",
    paste0(exposure, " * ", modifier),
    if (
      length(cov2) > 0
    ) {
      paste0(" + ", paste(cov2, collapse = " + "))
    } else {
      ""
    }
  )

  fit <- vegan::adonis2(as.formula(ftxt), data = dd, permutations = N_PERM, by = "margin")

  tab <- as.data.frame(fit)

  tab$term <- rownames(tab)

  idx <- grepl(exposure, tab$term, fixed = TRUE) &
    grepl(modifier, tab$term, fixed = TRUE) &
    grepl(":", tab$term, fixed = TRUE)

  fila <- tab[
    idx,
    ,
    drop = FALSE
  ]

  if (
    nrow(fila) != 1
  ) {
    return(
      tibble(
        distancia = distancia,
        exposure = exposure,
        modifier = modifier,
        n = nrow(dd),
        R2 = NA_real_,
        F = NA_real_,
        p = NA_real_,
        formula = ftxt
      )
    )
  }

  tibble(
    distancia = distancia,
    exposure = exposure,
    modifier = modifier,
    n = nrow(dd),
    R2 = fila$R2,
    F = fila$F,
    p = fila$`Pr(>F)`,
    formula = ftxt
  )
}

beta_int_sexo <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_interaccion(distancias[[nm]], nm, "miteco10", "sex")
    }
  )
) %>%
  mutate(p_BH_distancias = p.adjust(p, method = "BH"))

beta_int_isla <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_interaccion(distancias[[nm]], nm, "miteco10", "island")
    }
  )
) %>%
  mutate(p_BH_distancias = p.adjust(p, method = "BH"))

beta_int_exac <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_interaccion(distancias[[nm]], nm, "miteco10", "exac6m")
    }
  )
) %>%
  mutate(p_BH_distancias = p.adjust(p, method = "BH"))

beta_int_regimen <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_interaccion(
        distancias[[nm]],
        nm,
        "miteco10",
        "regimen_strict",
        data = meta_base %>%
          filter(!is.na(regimen_strict)) %>%
          droplevels()
      )
    }
  )
) %>%
  mutate(p_BH_distancias = p.adjust(p, method = "BH"))

# Beta: ventanas temporales simultaneas.

permanova_ventanas <- function(
  dist_obj,
  distancia
) {
  ids <- labels(dist_obj)

  dd <- meta_base[
    match(ids, meta_base$id_muestra),
    ,
    drop = FALSE
  ]

  rownames(dd) <- dd$id_muestra

  vars <- unique(c(windows, covars_beta))

  cc <- complete.cases(dd[ , vars, drop = FALSE ])

  dd <- droplevels(dd[ cc, , drop = FALSE ])

  d2 <- subdist(dist_obj, rownames(dd))

  dd <- dd[
    labels(d2),
    ,
    drop = FALSE
  ]

  cov2 <- covars_activas(dd, covars_beta)

  ftxt <- paste("d2 ~", paste(c(windows, cov2), collapse = " + "))

  fit <- vegan::adonis2(as.formula(ftxt), data = dd, permutations = N_PERM, by = "margin")

  tab <- as.data.frame(fit)

  tab$term <- rownames(tab)

  bind_rows(
    lapply(
      windows,
      function(w) {
        fila <- tab[
          tab$term == w,
          ,
          drop = FALSE
        ]

        if (
          nrow(fila) != 1
        ) {
          return(NULL)
        }

        tibble(
          distancia = distancia,
          window = w,
          n = nrow(dd),
          R2 = fila$R2,
          F = fila$F,
          p = fila$`Pr(>F)`,
          formula = ftxt
        )
      }
    )
  )
}

beta_ventanas <- bind_rows(
  lapply(
    names(distancias),
    function(nm) {
      permanova_ventanas(distancias[[nm]], nm)
    }
  )
) %>%
  group_by(window) %>%
  mutate(p_BH_distancias = p.adjust(p, method = "BH")) %>%
  ungroup()
