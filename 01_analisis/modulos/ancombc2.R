# Abundancia diferencial con ANCOM-BC2.
# Se mantienen los filtros, las covariables y la sensibilidad a pseudoconteos.

prep_meta_ancom <- function(
  ps_obj
) {
  m <- sample_data_df(ps_obj)

# Las variables ya estan integradas porque ps_amb_u_m se actualizo.
  m$age <- a_num(m$age)

  m$sex_num <- a_num(m$sex_num)

  m$miteco10 <- a_num(m$miteco10)

  m$miteco_days10_z <- a_num(m$miteco_days10_z)

  m$naf10_z <- a_num(m$naf10_z)

  m$exac6m_num <- a_num(m$exac6m_num)

  m$miteco_sex_int <- m$miteco10 *
    m$sex_num

  m$miteco_exac_int <- m$miteco10 *
    m$exac6m_num

  m$island <- factor(m$island)

  m$pool <- factor(m$pool)

  m$antibiotic <- factor(m$antibiotic)

  m$cavities <- factor(m$cavities)

  m$gum <- factor(m$gum)

  m$year <- factor(m$year)

  m
}

preparar_nivel_ancom <- function(
  ps_obj,
  nivel
) {
  if (
    nivel == "ASV"
  ) {
    p <- ps_obj
  } else {
    if (
      !nivel %in%
        rank_names(ps_obj)
    ) {
      stop("No existe el rango ", nivel)
    }

    p <- tax_glom(ps_obj, taxrank = nivel, NArm = TRUE)
  }

  o <- otu_taxa_rows(p)

  nmin <- ceiling(PREV_ANCOM * ncol(o))

  keep <- rownames(o)[
    rowSums(o > 0) >=
      nmin
  ]

  p <- prune_taxa(keep, p)

  p <- prune_taxa(taxa_sums(p) > 0, p)

  p
}

ejecutar_ancom <- function(
  ps_obj,
  nivel,
  analisis,
  main_var,
  covars
) {
  meta <- prep_meta_ancom(ps_obj)

  vars <- c(main_var, covars)

  vars <- vars[
    vars %in% names(meta)
  ]

  cc <- complete.cases(meta[ , vars, drop = FALSE ])

  meta2 <- droplevels(meta[ cc, , drop = FALSE ])

  cov2 <- covars_activas(meta2, covars)

  vars_final <- c(main_var, cov2)

  cc <- complete.cases(meta[ , vars_final, drop = FALSE ])

  ids <- meta$id_muestra[
    cc
  ]

  p <- prune_samples(ids, ps_obj)

  p <- prune_taxa(taxa_sums(p) > 0, p)

  p <- preparar_nivel_ancom(p, nivel)

  meta3 <- prep_meta_ancom(p)

  rownames(meta3) <- meta3$id_muestra

  sample_data(p) <- sample_data(meta3)

  fix_formula <- paste(c(main_var, cov2), collapse = " + ")

  fit <- ANCOMBC::ancombc2(
    data = p,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = fix_formula,
    rand_formula = NULL,
    p_adj_method = "BH",
    pseudo_sens = TRUE,
    prv_cut = 0,
    lib_cut = 0,
    s0_perc = 0.05,
    group = NULL,
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = ALPHA_FDR,
    n_cl = 1,
    verbose = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE
  )

  saveRDS(fit, file.path(ruta_objetos, paste0("ANCOMBC2_", analisis, "_", nivel, ".rds")))

  rr <- as.data.frame(fit$res, stringsAsFactors = FALSE, check.names = FALSE)

  if (
    "taxon" %in%
      names(rr)
  ) {
    rr$taxon_id <- as.character(rr$taxon)
  } else {
    rr$taxon_id <- rownames(rr)
  }

  lfc_cols <- grep("^lfc_", names(rr), value = TRUE)

  lc <- lfc_cols[
    grepl(main_var, lfc_cols, fixed = TRUE)
  ]

  if (
    length(lc) != 1
  ) {
    stop(
      "No se identifica de forma unica el coeficiente ",
      main_var,
      " en ",
      analisis,
      " ",
      nivel,
      ". Candidatos: ",
      paste(lc, collapse = ", ")
    )
  }

  suf <- sub("^lfc_", "", lc)

  getcol <- function(
    pref
  ) {
    z <- paste0(pref, "_", suf)

    if (
      z %in%
        names(rr)
    ) {
      z
    } else {
      NA_character_
    }
  }

  cq <- getcol("q")

  cp <- getcol("p")

  cd <- getcol("diff")

  cpass <- getcol("passed_ss")

  crob <- getcol("diff_robust")

  if (
    is.na(cq) ||
      is.na(cp)
  ) {
    stop("No se encuentran columnas p/q en ", analisis, " ", nivel)
  }

  out <- data.frame(
    analisis = analisis,
    nivel = nivel,
    n_muestras = nsamples(p),
    n_taxa = ntaxa(p),
    taxon_id = rr$taxon_id,
    beta = rr[[lc]],
    p = rr[[cp]],
    q = rr[[cq]],
    diff_q = if (
      !is.na(cd)
    ) {
      rr[[cd]]
    } else {
      NA
    },
    passed_sensitivity = if (
      !is.na(cpass)
    ) {
      rr[[cpass]]
    } else {
      NA
    },
    diff_robust = if (
      !is.na(crob)
    ) {
      rr[[crob]]
    } else {
      NA
    },
    stringsAsFactors = FALSE
  )

  tt <- as.data.frame(tax_table(p), stringsAsFactors = FALSE, check.names = FALSE)

  tt$taxon_id <- rownames(tt)

  out <- left_join(out, tt, by = "taxon_id")

  out$significativo_q <- !is.na(out$q) &
    out$q <
      ALPHA_FDR

  out$robusto <- ifelse(
    !is.na(out$diff_robust),
    as.logical(out$diff_robust),
    out$significativo_q &
      (is.na(out$passed_sensitivity) | as.logical(out$passed_sensitivity))
  )

  out
}

cov_ancom_global <- c(
  "age",
  "sex_num",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m_num",
  "year",
  "sin_doy",
  "cos_doy"
)

cov_ancom_sexint <- c(
  "miteco10",
  "sex_num",
  "age",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "exac6m_num",
  "year",
  "sin_doy",
  "cos_doy"
)

cov_ancom_exacint <- c(
  "miteco10",
  "exac6m_num",
  "age",
  "sex_num",
  "island",
  "pool",
  "antibiotic",
  "cavities",
  "gum",
  "year",
  "sin_doy",
  "cos_doy"
)

# ANCOM-BC2 genero.

ancom_genus <- list()

ancom_genus[[1]] <- ejecutar_ancom(ps_amb_u_m, "Genus", "MITECO10", "miteco10", cov_ancom_global)

ancom_genus[[2]] <- ejecutar_ancom(
  ps_amb_u_m,
  "Genus",
  "MITECO_dias10",
  "miteco_days10_z",
  cov_ancom_global
)

ancom_genus[[3]] <- ejecutar_ancom(ps_amb_u_m, "Genus", "NAF10", "naf10_z", cov_ancom_global)

ancom_genus[[4]] <- ejecutar_ancom(
  ps_amb_u_m,
  "Genus",
  "MITECO10_x_sexo",
  "miteco_sex_int",
  cov_ancom_sexint
)

ancom_genus[[5]] <- ejecutar_ancom(
  ps_amb_u_m,
  "Genus",
  "MITECO10_x_exac6m",
  "miteco_exac_int",
  cov_ancom_exacint
)

ancom_genus_all <- bind_rows(ancom_genus)

ancom_genus_resumen <- ancom_genus_all %>%
  group_by(analisis, nivel) %>%
  summarise(
    n_muestras = first(n_muestras),

    n_taxa = first(n_taxa),

    n_q005 = sum(significativo_q, na.rm = TRUE),

    n_robustos = sum(robusto, na.rm = TRUE),

    .groups = "drop"
  )

registrar("ANCOM-BC2 Genus", "OK", paste0("Modelos=", n_distinct(ancom_genus_all$analisis)))

# ANCOM-BC2 sensibilidad family / genus / asv.

niveles_ancom <- c("Family", "Genus")

if (
  EJECUTAR_ANCOM_ASV
) {
  niveles_ancom <- c(niveles_ancom, "ASV")
}

ancom_sens <- list()
kk <- 1

for (niv in niveles_ancom) {
  ancom_sens[[kk]] <- ejecutar_ancom(
    ps_amb_u_m,
    niv,
    "MITECO10_sensibilidad",
    "miteco10",
    cov_ancom_global
  )

  kk <- kk +
    1

  ancom_sens[[kk]] <- ejecutar_ancom(
    ps_amb_u_m,
    niv,
    "MITECO10_x_sexo_sensibilidad",
    "miteco_sex_int",
    cov_ancom_sexint
  )

  kk <- kk +
    1
}

ancom_sens_all <- bind_rows(ancom_sens)

ancom_sens_resumen <- ancom_sens_all %>%
  group_by(analisis, nivel) %>%
  summarise(
    n_muestras = first(n_muestras),

    n_taxa = first(n_taxa),

    n_q005 = sum(significativo_q, na.rm = TRUE),

    n_robustos = sum(robusto, na.rm = TRUE),

    .groups = "drop"
  )

registrar("ANCOM-BC2 sensibilidad taxonomica", "OK", paste(niveles_ancom, collapse = "/"))
