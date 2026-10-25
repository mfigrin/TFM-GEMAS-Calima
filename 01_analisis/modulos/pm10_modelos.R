# Diversidad alfa: PM10 armonizado.

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

ajustar_alpha <- function(
  outcome,
  exposure,
  data = meta,
  covars = covars_alpha
) {
  vars <- unique(c(outcome, exposure, covars, "cluster_isla_fecha"))

  dd <- data[
    complete.cases(data[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)
  cov2 <- covariables_activas(dd, covars)

  if (n_niveles(dd[[exposure]]) < 2) {
    return(NULL)
  }

  formula_txt <- paste(outcome, "~", paste(c(exposure, cov2), collapse = " + "))

  fit <- lm(as.formula(formula_txt), data = dd)

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
      formula = formula_txt,
      .before = 1
    )
}

exposiciones_pm10_10d <- c("pm10_media_10d_z", "pm10_auc50_10d_z")

alpha_pm10_10d <- bind_rows(
  lapply(
    exposiciones_pm10_10d,
    function(e) {
      bind_rows(lapply(alpha_outcomes, function(o) ajustar_alpha(o, e)))
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(q_BH_outcomes = p.adjust(p, method = "BH")) %>%
  ungroup()

exposiciones_pm10_180d <- c("pm10_media_180d_z", "pm10_auc50_180d_z")

alpha_pm10_180d <- bind_rows(
  lapply(
    exposiciones_pm10_180d,
    function(e) {
      bind_rows(lapply(alpha_outcomes, function(o) ajustar_alpha(o, e)))
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(q_BH_outcomes = p.adjust(p, method = "BH")) %>%
  ungroup()

# Interacción PM10 x sexo.

ajustar_alpha_interaccion_sexo <- function(
  outcome,
  exposure,
  data = meta
) {
  covars <- setdiff(covars_alpha, "sex")

  vars <- unique(c(outcome, exposure, "sex", covars, "cluster_isla_fecha"))

  dd <- data[
    complete.cases(data[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)
  cov2 <- covariables_activas(dd, covars)

  formula_txt <- paste(
    outcome,
    "~",
    paste0(exposure, " * sex"),
    if (length(cov2) > 0) {
      paste0(" + ", paste(cov2, collapse = " + "))
    } else {
      ""
    }
  )

  fit <- lm(as.formula(formula_txt), data = dd)

  terminos_int <- names(coef(fit))[
    grepl(":", names(coef(fit)), fixed = TRUE) &
      grepl(exposure, names(coef(fit)), fixed = TRUE) &
      grepl("sex", names(coef(fit)), fixed = TRUE)
  ]

  w <- wald_robusto(fit, dd$cluster_isla_fecha, terminos_int)

  tibble(
    outcome = outcome,
    exposure = exposure,
    modifier = "sex",
    n = nrow(dd),
    n_clusters = n_distinct(dd$cluster_isla_fecha),
    chi2 = w$chi2,
    df = w$df,
    p_interaccion = w$p,
    formula = formula_txt
  )
}

alpha_pm10_x_sexo <- bind_rows(
  lapply(
    exposiciones_pm10_10d,
    function(e) {
      bind_rows(lapply(c("Shannon", "Simpson"), function(o) ajustar_alpha_interaccion_sexo(o, e)))
    }
  )
) %>%
  group_by(exposure) %>%
  mutate(q_BH_outcomes = p.adjust(p_interaccion, method = "BH")) %>%
  ungroup()

ajustar_alpha_estratificado_sexo <- function(
  outcome,
  exposure,
  sexo,
  data = meta
) {
  covars <- setdiff(covars_alpha, "sex")

  dd0 <- data %>%
    filter(sex == sexo) %>%
    droplevels()

  z <- ajustar_alpha(outcome = outcome, exposure = exposure, data = dd0, covars = covars)

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    mutate(sex = sexo, .before = 1)
}

alpha_pm10_sexo_estrat <- bind_rows(
  lapply(
    exposiciones_pm10_10d,
    function(e) {
      bind_rows(
        lapply(
          c("Female", "Male"),
          function(s) {
            bind_rows(
              lapply(
                c("Shannon", "Simpson"),
                function(o) {
                  ajustar_alpha_estratificado_sexo(o, e, s)
                }
              )
            )
          }
        )
      )
    }
  )
)

# Predicciones ajustadas: intrusión x sexo.

predicciones_marginales <- function(outcome) {
  covars <- setdiff(covars_alpha, "sex")

  vars <- unique(c(outcome, "miteco10", "sex", covars, "cluster_isla_fecha"))

  dd <- meta[
    complete.cases(meta[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)
  cov2 <- covariables_activas(dd, covars)

  formula_txt <- paste(
    outcome,
    "~ miteco10 * sex",
    if (length(cov2) > 0) {
      paste0(" + ", paste(cov2, collapse = " + "))
    } else {
      ""
    }
  )

  fit <- lm(as.formula(formula_txt), data = dd)

  V <- vcov_cluster(fit, dd$cluster_isla_fecha)

  b <- coef(fit)
  validos <- !is.na(b)
  b <- b[validos]
  V <- V[names(b), names(b), drop = FALSE]

  combos <- expand.grid(miteco10 = c(0, 1), sex = c("Female", "Male"), stringsAsFactors = FALSE)

  n_clusters <- n_distinct(dd$cluster_isla_fecha)
  crit <- qt(0.975, df = max(n_clusters - 1, 1))

  bind_rows(
    lapply(
      seq_len(nrow(combos)),
      function(i) {
        nd <- dd
        nd$miteco10 <- combos$miteco10[i]
        nd$sex <- factor(combos$sex[i], levels = levels(dd$sex))

        X <- model.matrix(delete.response(terms(fit)), data = nd)

        faltan_X <- setdiff(names(b), colnames(X))

        if (length(faltan_X) > 0) {
          extra <- matrix(0, nrow = nrow(X), ncol = length(faltan_X), dimnames = list(NULL, faltan_X))
          X <- cbind(X, extra)
        }

        X <- X[, names(b), drop = FALSE]
        L <- colMeans(X)

        estimacion <- as.numeric(sum(L * b))
        se <- sqrt(as.numeric(t(L) %*% V %*% L))

        tibble(
          outcome = outcome,
          miteco10 = combos$miteco10[i],
          exposicion = ifelse(
            combos$miteco10[i] == 1,
            "Con intrusión oficial",
            "Sin intrusión oficial"
          ),
          sex = combos$sex[i],
          sexo = ifelse(combos$sex[i] == "Male", "Hombres", "Mujeres"),
          n_modelo = nrow(dd),
          n_clusters = n_clusters,
          estimacion_ajustada = estimacion,
          se = se,
          IC95_inf = estimacion - crit * se,
          IC95_sup = estimacion + crit * se
        )
      }
    )
  )
}

predicciones_sexo <- bind_rows(predicciones_marginales("Shannon"), predicciones_marginales("Simpson"))

# PERMANOVA completa.

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

permanova_completa <- function(
  dist_obj,
  distancia,
  exposure,
  data = meta,
  covars = covars_beta,
  nperm = N_PERM
) {
  labs_dist <- dist_labels(dist_obj)

  ids <- intersect(labs_dist, as.character(data$id_muestra))

  if (length(ids) < 4) {
    stop(
      "PERMANOVA: no se pudieron alinear las muestras para ",
      distancia, " / ", exposure, ". ",
      "IDs comunes encontrados: ", length(ids), "."
    )
  }

  idx <- match(ids, data$id_muestra)
  dd <- data[idx, , drop = FALSE]
  dd <- as.data.frame(dd, stringsAsFactors = FALSE)

  vars <- unique(c(exposure, covars))
  faltan_vars <- setdiff(vars, names(dd))

  if (length(faltan_vars) > 0) {
    stop(
      "PERMANOVA ", distancia, " / ", exposure,
      ": faltan variables: ",
      paste(faltan_vars, collapse = ", ")
    )
  }

  cc <- complete.cases(dd[, vars, drop = FALSE])
  dd <- dd[cc, , drop = FALSE]
  dd <- droplevels(dd)

  if (nrow(dd) < 4) {
    stop(
      "PERMANOVA ", distancia, " / ", exposure,
      ": quedan menos de 4 casos completos (n = ", nrow(dd), ")."
    )
  }

# Submuestreamos la matriz de distancia con los IDs reales, no con
# rownames de un tibble.
  d2 <- subdist(dist_obj, dd$id_muestra)

# Reordenamos metadata exactamente como la matriz de distancia.
  labs_d2 <- dist_labels(d2)
  dd <- dd[
    match(labs_d2, dd$id_muestra),
    ,
    drop = FALSE
  ]
  rownames(dd) <- dd$id_muestra

  if (!identical(rownames(dd), labs_d2)) {
    stop(
      "PERMANOVA ", distancia, " / ", exposure,
      ": fallo de alineación entre metadata y matriz de distancia."
    )
  }

  cov2 <- covariables_activas(dd, covars)

  if (n_niveles(dd[[exposure]]) < 2) {
    warning(
      "PERMANOVA omitida para ", distancia, " / ", exposure,
      ": la exposición tiene menos de dos valores distintos tras casos completos."
    )
    return(NULL)
  }

  formula_txt <- paste("d2 ~", paste(c(exposure, cov2), collapse = " + "))

  set.seed(SEED)

  fit <- vegan::adonis2(as.formula(formula_txt), data = dd, permutations = nperm, by = "margin")

  tab <- as.data.frame(fit)
  tab$term <- rownames(tab)

# vegan puede cambiar ligeramente los nombres de columnas entre versiones.
  col_p <- intersect(c("Pr(>F)", "Pr..F."), names(tab))
  if (length(col_p) == 0) {
    stop(
      "PERMANOVA ", distancia, " / ", exposure,
      ": no se encontró la columna de p en adonis2. Columnas disponibles: ",
      paste(names(tab), collapse = ", ")
    )
  }
  col_p <- col_p[1]

  out <- tibble(
    distancia = rep(distancia, nrow(tab)),
    exposure = rep(exposure, nrow(tab)),
    n = rep(nrow(dd), nrow(tab)),
    term = as.character(tab$term),
    Df = as.numeric(tab$Df),
    SumOfSqs = as.numeric(tab$SumOfSqs),
    R2 = as.numeric(tab$R2),
    F = as.numeric(tab$F),
    p = as.numeric(tab[[col_p]]),
    formula = rep(formula_txt, nrow(tab))
  )

  idx_pruebas <- !out$term %in% c("Residual", "Total") &
    !is.na(out$p)

  out$q_BH_modelo <- NA_real_
  out$q_BH_modelo[idx_pruebas] <- p.adjust(out$p[idx_pruebas], method = "BH")

  out
}

exposiciones_beta <- c("pm10_media_10d_z", "pm10_auc50_10d_z")

resultados_permanova <- list()

for (nm in names(distancias)) {
  for (e in exposiciones_beta) {
    msg("PERMANOVA: %s ~ %s", nm, e)

    tmp <- tryCatch(
      permanova_completa(distancias[[nm]], nm, e),
      error = function(err) {
        registrar(paste0("PERMANOVA ", nm, " / ", e), "ERROR", conditionMessage(err))
        stop(err)
      }
    )

    if (!is.null(tmp) && nrow(tmp) > 0) {
      resultados_permanova[[paste(nm, e, sep = "__")]] <- tmp
    }
  }
}

permanova_todas <- dplyr::bind_rows(resultados_permanova)

if (
  nrow(permanova_todas) == 0 ||
  !all(c("term", "exposure", "distancia", "p") %in% names(permanova_todas))
) {
  stop(
    "No se generó una tabla PERMANOVA válida. ",
    "Columnas disponibles: ",
    paste(names(permanova_todas), collapse = ", "),
    ". Revisa los mensajes inmediatamente anteriores."
  )
}

beta_exposiciones <- permanova_todas %>%
  dplyr::filter(.data$term == .data$exposure) %>%
  dplyr::group_by(.data$exposure) %>%
  dplyr::mutate(q_BH_distancias = p.adjust(.data$p, method = "BH")) %>%
  dplyr::ungroup()

if (nrow(beta_exposiciones) == 0) {
  stop(
    "Las PERMANOVA se ejecutaron, pero no se localizaron las filas ",
    "correspondientes a las exposiciones. Términos disponibles: ",
    paste(unique(permanova_todas$term), collapse = ", ")
  )
}

# Incorporar el resultado principal ya calculado con 9.999 permutaciones.
carpeta_final <- file.path(PROJECT_DIR, "outputs/23_PERMANOVA_FINAL_ESTABLE")
beta_final <- read_csv(file.path(carpeta_final, "01_PERMANOVA_FINAL_9999.csv"), show_col_types = FALSE)
beta_final <- beta_final %>% transmute(
  distancia, exposure = "miteco10", n, term = "miteco10", R2,
  F = pseudo_F, p, formula, q_BH_distancias = q_BH_3_distancias,
  permutaciones)
beta_exposiciones$permutaciones <- N_PERM
beta_exposiciones <- bind_rows(beta_final, beta_exposiciones)
permanova_bray_miteco <- read_csv(
  file.path(carpeta_final, "02_PERMANOVA_BRAY_COVARIABLES_9999.csv"),
  show_col_types = FALSE)

permanova_bray_pm10 <- permanova_todas %>%
  filter(distancia == "Bray_Curtis", exposure == "pm10_media_10d_z")

# Reutilizar las comprobaciones finales, sin repetirlas con 999 permutaciones.
bd_final <- read_csv(file.path(carpeta_final, "03_PERMDISP_BRAY_9999.csv"), show_col_types = FALSE)
betadisper_miteco <- bd_final %>% transmute(
  distancia = "Bray_Curtis", grupo = "Intrusión oficial D-1 a D-10",
  n, F, p_betadisper = p, permutaciones)
cluster_final <- read_csv(file.path(carpeta_final, "04_SENSIBILIDAD_CLUSTER_ISLA_FECHA_9999.csv"), show_col_types = FALSE)
beta_cluster <- cluster_final %>% transmute(
  distancia = "Bray_cluster_isla_fecha", exposure = "miteco10", n_clusters,
  R2, F = pseudo_F, p, formula, permutaciones)

# PM10 a 180 días y fenotipo exacerbador.

ajustar_exacerbacion_outcome <- function(
  exposure,
  incluir_antibiotico = TRUE
) {
  covars <- c("age", "sex", "island", "pool", "cavities", "gum", "year", "sin_doy", "cos_doy")

  if (incluir_antibiotico) {
    covars <- c(covars, "antibiotic")
  }

  vars <- unique(c("exac6m_num", exposure, covars, "cluster_isla_fecha"))

  dd <- meta[
    complete.cases(meta[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  dd <- droplevels(dd)
  cov2 <- covariables_activas(dd, covars)

  formula_txt <- paste("exac6m_num ~", paste(c(exposure, cov2), collapse = " + "))

  fit <- glm(as.formula(formula_txt), family = binomial(), data = dd)

  z <- extraer_coef_robusto(fit, dd$cluster_isla_fecha, exposure)

  if (is.null(z)) {
    return(NULL)
  }

  z %>%
    transmute(
      exposure = exposure,
      modelo = ifelse(
        incluir_antibiotico,
        "Con ajuste por antibióticos",
        "Sin ajuste por antibióticos"
      ),
      n = nrow(dd),
      n_clusters = n_distinct(dd$cluster_isla_fecha),
      beta = beta,
      se = se,
      OR_por_1_DE = exp(beta),
      OR_IC95_inf = exp(IC95_inf),
      OR_IC95_sup = exp(IC95_sup),
      p = p,
      formula = formula_txt
    )
}

exac_pm10_180d <- bind_rows(
  lapply(
    exposiciones_pm10_180d,
    function(e) {
      bind_rows(ajustar_exacerbacion_outcome(e, TRUE), ajustar_exacerbacion_outcome(e, FALSE))
    }
  )
) %>%
  group_by(modelo) %>%
  mutate(q_BH_exposiciones = p.adjust(p, method = "BH")) %>%
  ungroup()

