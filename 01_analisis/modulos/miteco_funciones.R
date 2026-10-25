# Funciones de armonización y modelos ajustados.
# Se mantienen las codificaciones, los casos completos y los errores robustos originales.

a_num <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)), fixed = TRUE)))
}

a_binario <- function(x) {
  if (is.logical(x)) {
    return(as.numeric(x))
  }

  y <- tolower(trimws(as.character(x)))

  out <- rep(NA_real_, length(y))

  out[
    y %in% c(
      "1",
      "true",
      "t",
      "si",
      "sí",
      "yes",
      "y",
      "expuesto",
      "exposed",
      "exacerbador",
      "male",
      "m",
      "masculino",
      "hombre"
    )
  ] <- 1

  out[
    y %in% c(
      "0",
      "false",
      "f",
      "no",
      "n",
      "no_expuesto",
      "non-exposed",
      "nonexposed",
      "no_exacerbador",
      "female",
      "f",
      "femenino",
      "mujer"
    )
  ] <- 0

  suppressWarnings({
    z <- as.numeric(y)
  })

  idx <- is.na(out) &
    !is.na(z) &
    z %in% c(0, 1)

  out[idx] <- z[idx]

  out
}

limpiar_texto <- function(x) {
  y <- trimws(as.character(x))

  y[
    y %in% c("", "NA", "NaN", "NULL", "null")
  ] <- NA_character_

  y
}

normalizar_isla <- function(x) {
  y <- limpiar_texto(x)

  y[
    y %in% c("La_Palma", "La palma", "LA PALMA")
  ] <- "La Palma"

  y[
    y %in% c("TENERIFE", "tenerife")
  ] <- "Tenerife"

  y
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)

  if (
    is.na(s) ||
      s == 0
  ) {
    return(rep(NA_real_, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) /
    s
}

sample_data_df <- function(ps_obj) {
  sd_obj <- sample_data(ps_obj)

  out <- data.frame(
    setNames(sd_obj@.Data, names(sd_obj)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(out) <- sample_names(ps_obj)
  out$id_muestra <- rownames(out)

  out
}

otu_taxa_rows <- function(ps_obj) {
  m <- as(otu_table(ps_obj), "matrix")

  if (!taxa_are_rows(ps_obj)) {
    m <- t(m)
  }

  m
}

tax_df <- function(ps_obj) {
  tt <- as.data.frame(tax_table(ps_obj), stringsAsFactors = FALSE, check.names = FALSE)

  for (v in names(tt)) {
    tt[[v]] <- limpiar_texto(tt[[v]])
  }

  tt$taxon_id <- rownames(tt)
  tt
}

n_niveles <- function(x) {
  if (is.factor(x)) {
    return(nlevels(droplevels(x)))
  }

  length(unique(x[ !is.na(x) ]))
}

covars_activas <- function(
  d,
  vars
) {
  vars[
    vapply(
      vars,
      function(v) {
        v %in% names(d) &&
          n_niveles(d[[v]]) >= 2
      },
      logical(1)
    )
  ]
}

subdist <- function(
  dist_obj,
  ids
) {
  m <- as.matrix(dist_obj)

  ids <- intersect(ids, rownames(m))

  as.dist(m[ ids, ids, drop = FALSE ])
}

vcov_cluster <- function(
  fit,
  cluster
) {
  sandwich::vcovCL(fit, cluster = cluster, type = "HC1")
}

extraer_coef_robusto <- function(
  fit,
  cluster,
  termino
) {
  ct <- lmtest::coeftest(fit, vcov. = vcov_cluster(fit, cluster))

  if (!termino %in% rownames(ct)) {
    return(NULL)
  }

  beta <- ct[
    termino,
    "Estimate"
  ]

  se <- ct[
    termino,
    "Std. Error"
  ]

  tibble(
    termino = termino,
    beta = beta,
    se = se,
    IC95_inf = beta -
      1.96 *
      se,
    IC95_sup = beta +
      1.96 *
      se,
    p = ct[
      termino,
      ncol(ct)
    ]
  )
}

wald_terms <- function(
  fit,
  cluster,
  terms
) {
  b <- coef(fit)

  V <- vcov_cluster(fit, cluster)

  terms <- intersect(terms, names(b))

  terms <- terms[
    !is.na(b[ terms ])
  ]

  if (length(terms) == 0) {
    return(tibble(chi2 = NA_real_, df = 0L, p = NA_real_))
  }

  bb <- b[
    terms
  ]

  VV <- V[
    terms,
    terms,
    drop = FALSE
  ]

  inv <- tryCatch(
    solve(VV),
    error = function(e) {
      qr.solve(VV)
    }
  )

  W <- as.numeric(t(bb) %*% inv %*% bb)

  tibble(chi2 = W, df = length(bb), p = pchisq(W, df = length(bb), lower.tail = FALSE))
}

buscar_terminos_interaccion <- function(
  fit,
  token1,
  token2
) {
  nm <- names(coef(fit))

  nm[
    grepl(token1, nm, fixed = TRUE) &
      grepl(token2, nm, fixed = TRUE) &
      grepl(":", nm, fixed = TRUE)
  ]
}

manifest <- tibble(bloque = character(), estado = character(), detalle = character())

registrar <- function(
  bloque,
  estado,
  detalle = ""
) {
  manifest <<- bind_rows(manifest, tibble(bloque = bloque, estado = estado, detalle = detalle))
}
