# Agrupar taxones por linaje y conservar las lecturas.
# Estas funciones son las del script descriptivo original.

limpiar_texto <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  y
}

a_num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)), fixed = TRUE)))
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

etiquetas_rango <- function(ps_obj, rango) {
  tt <- tax_df(ps_obj)

  if (!rango %in% names(tt)) {
    stop("El rango ", rango, " no existe en tax_table.")
  }

  orden_rangos <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  orden_rangos <- orden_rangos[
    orden_rangos %in% names(tt)
  ]

  pos <- match(rango, orden_rangos)

  if (is.na(pos)) {
    stop("No se pudo localizar el rango ", rango, " en la jerarquia taxonomica.")
  }

# CLAVE TAXONOMICA
# La agregacion se hace usando TODA la linea taxonomica hasta
# el rango solicitado. De esta forma dos taxones cuyo Genus sea
# literalmente "uncultured" NO se fusionan si pertenecen a
# familias/ordenes diferentes.

  rangos_clave <- orden_rangos[
    seq_len(pos)
  ]

  clave_df <- tt[
    ,
    rangos_clave,
    drop = FALSE
  ]

  for (v in names(clave_df)) {
    clave_df[[v]][
      is.na(clave_df[[v]])
    ] <- "<NA>"
  }

  clave <- apply(clave_df, 1, paste, collapse = " | ")

# ETIQUETA HUMANA

  etiqueta <- tt[[rango]]

  if (rango == "Phylum") {
    etiqueta[
      is.na(etiqueta)
    ] <- "Unclassified_Phylum"
  }

  if (rango == "Family") {

    order <- if (
      "Order" %in% names(tt)
    ) {
      tt$Order
    } else {
      rep(NA_character_, nrow(tt))
    }

    idx <- is.na(etiqueta) &
      !is.na(order)

    etiqueta[idx] <- paste0("Unclassified_", order[idx])

    etiqueta[
      is.na(etiqueta)
    ] <- "Unclassified_Family"
  }

  if (rango == "Genus") {

    family <- if (
      "Family" %in% names(tt)
    ) {
      tt$Family
    } else {
      rep(NA_character_, nrow(tt))
    }

    order <- if (
      "Order" %in% names(tt)
    ) {
      tt$Order
    } else {
      rep(NA_character_, nrow(tt))
    }

    idx <- is.na(etiqueta) &
      !is.na(family)

    etiqueta[idx] <- paste0("Unclassified_", family[idx])

    idx <- is.na(etiqueta) &
      !is.na(order)

    etiqueta[idx] <- paste0("Unclassified_", order[idx])

    etiqueta[
      is.na(etiqueta)
    ] <- "Unclassified_Genus"
  }

# Solo se eliminan caracteres de CONTROL.
# IMPORTANTE: no usar "[\\r\\n\\t]" con doble escape, porque
# puede interpretarse como un conjunto que incluye las letras
# r, n y t. La clase POSIX evita ese problema.
  etiqueta <- gsub("[[:cntrl:]]", " ", etiqueta)

  etiqueta <- trimws(etiqueta)

# MAPA CLAVE -> ETIQUETA

  mapa <- data.frame(clave = clave, etiqueta = etiqueta, stringsAsFactors = FALSE) %>%
    distinct(clave, .keep_all = TRUE)

# Si la misma etiqueta visible aparece en varias lineas
# taxonomicas, se añade contexto taxonomico para que no exista
# ninguna fusion ni ambiguedad en tablas/figuras.
  duplicada <- duplicated(mapa$etiqueta) |
    duplicated(mapa$etiqueta, fromLast = TRUE)

  if (
    any(duplicada)
  ) {

    tax_context <- tt %>%
      mutate(.clave = clave) %>%
      distinct(.clave, .keep_all = TRUE)

    mapa <- mapa %>%
      left_join(
        tax_context %>%
          select(.clave, any_of(c("Order", "Family"))),
        by = c("clave" = ".clave")
      )

    idx <- duplicada

    if (
      "Family" %in% names(mapa)
    ) {
      idx_family <- idx &
        !is.na(mapa$Family)

      mapa$etiqueta[idx_family] <- paste0(mapa$etiqueta[idx_family], " [", mapa$Family[idx_family], "]")
    }

# Si aun quedan etiquetas duplicadas, añadir Order.
    duplicada2 <- duplicated(mapa$etiqueta) |
      duplicated(mapa$etiqueta, fromLast = TRUE)

    if (
      any(duplicada2) &&
      "Order" %in% names(mapa)
    ) {

      idx_order <- duplicada2 &
        !is.na(mapa$Order)

      mapa$etiqueta[idx_order] <- paste0(mapa$etiqueta[idx_order], " [", mapa$Order[idx_order], "]")
    }
  }

# Ultima garantia de unicidad para visualizacion.
  mapa$etiqueta <- make.unique(mapa$etiqueta, sep = "_")

  asignacion_clave <- setNames(mapa$etiqueta, mapa$clave)

  etiqueta_asv <- unname(asignacion_clave[ clave ])

  names(etiqueta_asv) <- tt$taxon_id

  attr(etiqueta_asv, "clave_taxonomica") <- setNames(clave, tt$taxon_id)

  attr(etiqueta_asv, "mapa_clave_etiqueta") <- mapa

  etiqueta_asv
}

agregar_rango <- function(ps_obj, rango) {

  otu <- otu_taxa_rows(ps_obj)

  etiquetas_obj <- etiquetas_rango(ps_obj, rango)

  etiquetas <- etiquetas_obj[
    rownames(otu)
  ]

  claves <- attr(etiquetas_obj, "clave_taxonomica")

  claves <- claves[
    rownames(otu)
  ]

  if (
    any(is.na(etiquetas)) ||
    any(is.na(claves))
  ) {
    stop("Hay ASVs sin etiqueta/clave agregable en rango ", rango, ".")
  }

# Agregar por CLAVE TAXONOMICA, no por nombre visible.
  agregado <- rowsum(otu, group = claves, reorder = FALSE)

  mapa_asv <- data.frame(clave = claves, etiqueta = etiquetas, stringsAsFactors = FALSE) %>%
    distinct(clave, .keep_all = TRUE)

  etiquetas_finales <- mapa_asv$etiqueta[
    match(rownames(agregado), mapa_asv$clave)
  ]

  if (
    any(is.na(etiquetas_finales))
  ) {
    stop("No se pudo asignar etiqueta final a alguna clave de ", rango, ".")
  }

  if (
    anyDuplicated(etiquetas_finales) > 0
  ) {
    stop("Persisten etiquetas duplicadas tras desambiguar ", rango, ".")
  }

  rownames(agregado) <- etiquetas_finales

  agregado <- agregado[
    rowSums(agregado) > 0,
    ,
    drop = FALSE
  ]

# QC esencial: la agregacion NO puede cambiar las lecturas
# de ninguna muestra.
  original_lib <- colSums(otu)

  agregado_lib <- colSums(agregado)

  if (
    !identical(names(original_lib), names(agregado_lib))
  ) {
    stop("Las muestras cambiaron de orden/nombre al agregar ", rango, ".")
  }

  if (
    any(original_lib != agregado_lib)
  ) {
    stop("La agregacion en ", rango, " no conserva exactamente el numero de lecturas.")
  }

  agregado
}

a_relativa <- function(counts) {
  libs <- colSums(counts)

  if (any(libs <= 0)) {
    stop("Hay muestras con suma de lecturas <=0.")
  }

  sweep(counts, 2, libs, "/")
}

resumen_taxones <- function(counts) {
  rel <- a_relativa(counts)

  data.frame(
    taxon = rownames(counts),
    abundancia_total = rowSums(counts),
    prevalencia_n = rowSums(counts > 0),
    prevalencia_prop = rowMeans(counts > 0),
    abundancia_relativa_media = rowMeans(rel),
    abundancia_relativa_mediana = apply(rel, 1, median),
    abundancia_relativa_P25 = apply(rel, 1, quantile, probs = 0.25, names = FALSE),
    abundancia_relativa_P75 = apply(rel, 1, quantile, probs = 0.75, names = FALSE),
    abundancia_relativa_max = apply(rel, 1, max),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(abundancia_relativa_media))
}
