# COMPOSICION_TAXONOMICA_DESCRIPTIVA_GEMAS.R
# TFM GEMAS - Composición taxonómica descriptiva
#
# Objetivo:
#   Caracterizar la microbiota de vía aérea superior de la cohorte
#   GEMAS y generar tablas/figuras descriptivas definitivas.

# Cohortes:
#   - Saliva completa: 206 muestras
#     -> descripción global, isla, exacerbación 6M y sexo.
#   - Ambiental: 205 muestras
#     -> descripción según exposición SDI/calima 10 días.
#
# Rango taxonómico:
#   - Phylum
#   - Family
#   - Genus
#
# Salida:
# outputs/analisis_microbiota/composicion_taxonomica/

# 0. PAQUETES Y CONFIGURACIÓN
paquetes <- c(
  "phyloseq",
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "ggplot2",
  "writexl"
)

faltan <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltan) > 0) {
  stop(
    "Faltan paquetes: ",
    paste(faltan, collapse = ", "),
    ". Instálalos antes de ejecutar el script."
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(writexl)
})

set.seed(20260822)

top_phylum <- 8
top_family <- 15
top_genus <- 20

top_phylum_grupos <- 8
top_genus_grupos <- 12


# 1. RUTAS

ruta_objetos <- file.path(
  "outputs", "analisis_microbiota", "objetos_phyloseq"
)

ruta_ps_saliva <- file.path(
  ruta_objetos,
  "ps_saliva_206_sin_filtrado.rds"
)

ruta_ps_ambiental <- file.path(
  ruta_objetos,
  "ps_ambiental_205_sin_filtrado.rds"
)

if (!file.exists(ruta_ps_saliva)) {
  stop("No existe: ", ruta_ps_saliva)
}

if (!file.exists(ruta_ps_ambiental)) {
  stop("No existe: ", ruta_ps_ambiental)
}

ruta_salida <- file.path(
  "outputs", "analisis_microbiota", "composicion_taxonomica_v3"
)

ruta_tablas <- file.path(ruta_salida, "tablas")
ruta_figuras <- file.path(ruta_salida, "figuras")
ruta_objetos_salida <- file.path(ruta_salida, "objetos")

dir.create(ruta_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_objetos_salida, recursive = TRUE, showWarnings = FALSE)


# 2. FUNCIONES AUXILIARES

limpiar_texto <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  y
}

a_num <- function(x) {
  suppressWarnings(
    as.numeric(
      gsub(",", ".", trimws(as.character(x)), fixed = TRUE)
    )
  )
}

a_binario <- function(x) {
  y <- tolower(trimws(as.character(x)))

  out <- rep(NA_real_, length(y))

  out[y %in% c(
    "1", "true", "t", "si", "sí", "yes", "y",
    "expuesto", "exposed",
    "exacerbador"
  )] <- 1

  out[y %in% c(
    "0", "false", "f", "no", "n",
    "no_expuesto", "non-exposed", "nonexposed",
    "no_exacerbador"
  )] <- 0

  suppressWarnings({
    z <- as.numeric(y)
  })

  idx <- is.na(out) & !is.na(z) & z %in% c(0, 1)
  out[idx] <- z[idx]

  out
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

  out %>% relocate(id_muestra)
}

preparar_metadata <- function(ps_obj) {
  meta <- sample_data_df(ps_obj)

  if ("Age" %in% names(meta)) {
    meta$edad <- a_num(meta$Age)
  }

  if ("Island" %in% names(meta)) {
    isla <- limpiar_texto(meta$Island)

    isla[isla %in% c("La_Palma", "La Palma")] <- "La Palma"
    isla[isla %in% c("Tenerife")] <- "Tenerife"

    meta$isla_factor <- factor(
      isla,
      levels = c("La Palma", "Tenerife")
    )
  }

  if ("Sex.letter" %in% names(meta)) {
    sx <- tolower(trimws(as.character(meta$Sex.letter)))

    sexo <- ifelse(
      sx %in% c("female", "f", "femenino", "mujer"),
      "Mujer",
      ifelse(
        sx %in% c("male", "m", "masculino", "hombre"),
        "Hombre",
        NA_character_
      )
    )

    meta$sexo_factor <- factor(
      sexo,
      levels = c("Mujer", "Hombre")
    )
  }

  if ("Exacerbations.6M" %in% names(meta)) {
    ex <- a_num(meta$Exacerbations.6M)
    ex[!ex %in% c(0, 1)] <- NA_real_

    meta$exac6m_factor <- factor(
      ifelse(
        is.na(ex),
        NA,
        ifelse(ex == 1, "Exacerbador", "No exacerbador")
      ),
      levels = c("No exacerbador", "Exacerbador")
    )
  }

  if ("sdi_exposicion_10d" %in% names(meta)) {
    sdi <- a_binario(meta$sdi_exposicion_10d)

    meta$sdi10_factor <- factor(
      ifelse(
        is.na(sdi),
        NA,
        ifelse(sdi == 1, "Expuesto", "No expuesto")
      ),
      levels = c("No expuesto", "Expuesto")
    )
  }

  meta
}

otu_taxa_rows <- function(ps_obj) {
  m <- as(otu_table(ps_obj), "matrix")

  if (!taxa_are_rows(ps_obj)) {
    m <- t(m)
  }

  m
}

tax_df <- function(ps_obj) {
  tt <- as.data.frame(
    tax_table(ps_obj),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

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

  orden_rangos <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )

  orden_rangos <- orden_rangos[
    orden_rangos %in% names(tt)
  ]

  pos <- match(
    rango,
    orden_rangos
  )

  if (is.na(pos)) {
    stop(
      "No se pudo localizar el rango ",
      rango,
      " en la jerarquia taxonomica."
    )
  }

  # ----------------------------------------------------------
  # CLAVE TAXONOMICA

  # La agregacion se hace usando TODA la linea taxonomica hasta el rango solicitado. De esta forma dos taxones cuyo Genus sea
  # literalmente "uncultured" NO se fusionan si pertenecen a familias/ordenes diferentes.

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

  clave <- apply(
    clave_df,
    1,
    paste,
    collapse = " | "
  )


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
      rep(
        NA_character_,
        nrow(tt)
      )
    }

    idx <- is.na(etiqueta) &
      !is.na(order)

    etiqueta[idx] <- paste0(
      "Unclassified_",
      order[idx]
    )

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
      rep(
        NA_character_,
        nrow(tt)
      )
    }

    order <- if (
      "Order" %in% names(tt)
    ) {
      tt$Order
    } else {
      rep(
        NA_character_,
        nrow(tt)
      )
    }

    idx <- is.na(etiqueta) &
      !is.na(family)

    etiqueta[idx] <- paste0(
      "Unclassified_",
      family[idx]
    )

    idx <- is.na(etiqueta) &
      !is.na(order)

    etiqueta[idx] <- paste0(
      "Unclassified_",
      order[idx]
    )

    etiqueta[
      is.na(etiqueta)
    ] <- "Unclassified_Genus"
  }

  etiqueta <- gsub(
    "[[:cntrl:]]",
    " ",
    etiqueta
  )

  etiqueta <- trimws(
    etiqueta
  )

  # MAPA CLAVE -> ETIQUETA

  mapa <- data.frame(
    clave = clave,
    etiqueta = etiqueta,
    stringsAsFactors = FALSE
  ) %>%
    distinct(
      clave,
      .keep_all = TRUE
    )

  # Si la misma etiqueta visible aparece en varias lineas
  # taxonomicas, se añade contexto taxonomico para que no exista
  # ninguna fusion ni ambiguedad en tablas/figuras.
  duplicada <- duplicated(
    mapa$etiqueta
  ) |
    duplicated(
      mapa$etiqueta,
      fromLast = TRUE
    )

  if (
    any(
      duplicada
    )
  ) {

    tax_context <- tt %>%
      mutate(
        .clave = clave
      ) %>%
      distinct(
        .clave,
        .keep_all = TRUE
      )

    mapa <- mapa %>%
      left_join(
        tax_context %>%
          select(
            .clave,
            any_of(
              c(
                "Order",
                "Family"
              )
            )
          ),
        by = c(
          "clave" = ".clave"
        )
      )

    idx <- duplicada

    if (
      "Family" %in% names(mapa)
    ) {
      idx_family <- idx &
        !is.na(mapa$Family)

      mapa$etiqueta[idx_family] <- paste0(
        mapa$etiqueta[idx_family],
        " [",
        mapa$Family[idx_family],
        "]"
      )
    }

    # Si aun quedan etiquetas duplicadas, añadir Order.
    duplicada2 <- duplicated(
      mapa$etiqueta
    ) |
      duplicated(
        mapa$etiqueta,
        fromLast = TRUE
      )

    if (
      any(
        duplicada2
      ) &&
      "Order" %in% names(mapa)
    ) {

      idx_order <- duplicada2 &
        !is.na(mapa$Order)

      mapa$etiqueta[idx_order] <- paste0(
        mapa$etiqueta[idx_order],
        " [",
        mapa$Order[idx_order],
        "]"
      )
    }
  }

  mapa$etiqueta <- make.unique(
    mapa$etiqueta,
    sep = "_"
  )

  asignacion_clave <- setNames(
    mapa$etiqueta,
    mapa$clave
  )

  etiqueta_asv <- unname(
    asignacion_clave[
      clave
    ]
  )

  names(
    etiqueta_asv
  ) <- tt$taxon_id

  attr(
    etiqueta_asv,
    "clave_taxonomica"
  ) <- setNames(
    clave,
    tt$taxon_id
  )

  attr(
    etiqueta_asv,
    "mapa_clave_etiqueta"
  ) <- mapa

  etiqueta_asv
}


agregar_rango <- function(ps_obj, rango) {

  otu <- otu_taxa_rows(
    ps_obj
  )

  etiquetas_obj <- etiquetas_rango(
    ps_obj,
    rango
  )

  etiquetas <- etiquetas_obj[
    rownames(otu)
  ]

  claves <- attr(
    etiquetas_obj,
    "clave_taxonomica"
  )

  claves <- claves[
    rownames(otu)
  ]

  if (
    any(
      is.na(etiquetas)
    ) ||
    any(
      is.na(claves)
    )
  ) {
    stop(
      "Hay ASVs sin etiqueta/clave agregable en rango ",
      rango,
      "."
    )
  }

  # Agregar por CLAVE TAXONOMICA, no por nombre visible.
  agregado <- rowsum(
    otu,
    group = claves,
    reorder = FALSE
  )

  mapa_asv <- data.frame(
    clave = claves,
    etiqueta = etiquetas,
    stringsAsFactors = FALSE
  ) %>%
    distinct(
      clave,
      .keep_all = TRUE
    )

  etiquetas_finales <- mapa_asv$etiqueta[
    match(
      rownames(agregado),
      mapa_asv$clave
    )
  ]

  if (
    any(
      is.na(
        etiquetas_finales
      )
    )
  ) {
    stop(
      "No se pudo asignar etiqueta final a alguna clave de ",
      rango,
      "."
    )
  }

  if (
    anyDuplicated(
      etiquetas_finales
    ) > 0
  ) {
    stop(
      "Persisten etiquetas duplicadas tras desambiguar ",
      rango,
      "."
    )
  }

  rownames(
    agregado
  ) <- etiquetas_finales

  agregado <- agregado[
    rowSums(
      agregado
    ) > 0,
    ,
    drop = FALSE
  ]

  original_lib <- colSums(
    otu
  )

  agregado_lib <- colSums(
    agregado
  )

  if (
    !identical(
      names(original_lib),
      names(agregado_lib)
    )
  ) {
    stop(
      "Las muestras cambiaron de orden/nombre al agregar ",
      rango,
      "."
    )
  }

  if (
    any(
      original_lib !=
        agregado_lib
    )
  ) {
    stop(
      "La agregacion en ",
      rango,
      " no conserva exactamente el numero de lecturas."
    )
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
    abundancia_relativa_P25 = apply(
      rel,
      1,
      quantile,
      probs = 0.25,
      names = FALSE
    ),
    abundancia_relativa_P75 = apply(
      rel,
      1,
      quantile,
      probs = 0.75,
      names = FALSE
    ),
    abundancia_relativa_max = apply(rel, 1, max),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(abundancia_relativa_media))
}

matriz_relativa_df <- function(rel) {
  as.data.frame(
    t(rel),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    rownames_to_column("id_muestra")
}

crear_top_barplot <- function(resumen, top_n, rango) {
  dd <- resumen %>%
    slice_head(n = top_n) %>%
    mutate(
      taxon_plot = factor(
        taxon,
        levels = rev(taxon)
      )
    )

  p <- ggplot(
    dd,
    aes(
      x = abundancia_relativa_media,
      y = taxon_plot
    )
  ) +
    geom_col() +
    labs(
      title = paste0(
        "Taxones más abundantes - ",
        rango
      ),
      subtitle = "Media de abundancia relativa en 206 muestras de saliva",
      x = "Abundancia relativa media",
      y = NULL
    ) +
    scale_x_continuous(
      labels = function(x) paste0(round(100 * x, 1), "%")
    ) +
    theme_minimal(base_size = 12)

  archivo <- file.path(
    ruta_figuras,
    paste0(
      "01_top_",
      tolower(rango),
      "_global.png"
    )
  )

  ggsave(
    archivo,
    p,
    width = 8,
    height = max(5, 0.35 * nrow(dd) + 2),
    dpi = 300
  )

  archivo
}

crear_prevalencia_abundancia <- function(resumen_genero) {
  dd <- resumen_genero %>%
    mutate(
      abundancia_pct = 100 * abundancia_relativa_media,
      prevalencia_pct = 100 * prevalencia_prop
    )

  p <- ggplot(
    dd,
    aes(
      x = prevalencia_pct,
      y = abundancia_pct
    )
  ) +
    geom_point(alpha = 0.65) +
    labs(
      title = "Prevalencia y abundancia media de los géneros",
      subtitle = "Cohorte completa de saliva (n=206)",
      x = "Prevalencia (%)",
      y = "Abundancia relativa media (%)"
    ) +
    theme_minimal(base_size = 12)

  archivo <- file.path(
    ruta_figuras,
    "04_prevalencia_vs_abundancia_genero.png"
  )

  ggsave(
    archivo,
    p,
    width = 8,
    height = 6,
    dpi = 300
  )

  archivo
}

resumen_por_grupo <- function(
  rel,
  meta,
  variable_grupo,
  top_taxa,
  nombre_grupo
) {
  ids_comunes <- intersect(
    colnames(rel),
    meta$id_muestra
  )

  rel2 <- rel[, ids_comunes, drop = FALSE]

  meta2 <- meta %>%
    filter(id_muestra %in% ids_comunes)

  meta2 <- meta2[
    match(colnames(rel2), meta2$id_muestra),
    ,
    drop = FALSE
  ]

  grupo <- meta2[[variable_grupo]]

  keep <- !is.na(grupo)

  rel2 <- rel2[, keep, drop = FALSE]
  grupo <- droplevels(factor(grupo[keep]))

  taxa_presentes <- intersect(top_taxa, rownames(rel2))

  top_mat <- rel2[taxa_presentes, , drop = FALSE]

  resto <- pmax(
    0,
    1 - colSums(top_mat)
  )

  plot_mat <- rbind(
    top_mat,
    Other = resto
  )

  dd <- as.data.frame(
    t(plot_mat),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    mutate(
      grupo = as.character(grupo)
    ) %>%
    pivot_longer(
      cols = -grupo,
      names_to = "taxon",
      values_to = "abundancia_relativa"
    ) %>%
    group_by(grupo, taxon) %>%
    summarise(
      abundancia_relativa_media = mean(abundancia_relativa),
      abundancia_relativa_mediana = median(abundancia_relativa),
      .groups = "drop"
    )

  orden_taxa <- c(taxa_presentes, "Other")

  dd$taxon <- factor(
    dd$taxon,
    levels = orden_taxa
  )

  dd$grupo <- factor(
    dd$grupo,
    levels = levels(grupo)
  )

  dd$tipo_grupo <- nombre_grupo

  dd
}

crear_stacked_grupo <- function(
  dd,
  titulo,
  archivo_nombre
) {
  p <- ggplot(
    dd,
    aes(
      x = grupo,
      y = abundancia_relativa_media,
      fill = taxon
    )
  ) +
    geom_col(
      position = "fill",
      width = 0.72
    ) +
    labs(
      title = titulo,
      subtitle = "Composición media descriptiva",
      x = NULL,
      y = "Abundancia relativa",
      fill = "Taxón"
    ) +
    scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      expand = c(0, 0)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right"
    )

  archivo <- file.path(
    ruta_figuras,
    archivo_nombre
  )

  ggsave(
    archivo,
    p,
    width = 9,
    height = 6,
    dpi = 300
  )

  archivo
}

resumen_cohorte <- function(meta_saliva, meta_ambiental) {
  data.frame(
    indicador = c(
      "Muestras saliva",
      "Muestras ambientales",
      "Edad media saliva",
      "Edad SD saliva",
      "Mujeres",
      "Hombres",
      "La Palma",
      "Tenerife",
      "No exacerbadores 6M",
      "Exacerbadores 6M",
      "No expuestos SDI10",
      "Expuestos SDI10"
    ),
    valor = c(
      nrow(meta_saliva),
      nrow(meta_ambiental),
      round(mean(meta_saliva$edad, na.rm = TRUE), 2),
      round(sd(meta_saliva$edad, na.rm = TRUE), 2),
      sum(meta_saliva$sexo_factor == "Mujer", na.rm = TRUE),
      sum(meta_saliva$sexo_factor == "Hombre", na.rm = TRUE),
      sum(meta_saliva$isla_factor == "La Palma", na.rm = TRUE),
      sum(meta_saliva$isla_factor == "Tenerife", na.rm = TRUE),
      sum(meta_saliva$exac6m_factor == "No exacerbador", na.rm = TRUE),
      sum(meta_saliva$exac6m_factor == "Exacerbador", na.rm = TRUE),
      sum(meta_ambiental$sdi10_factor == "No expuesto", na.rm = TRUE),
      sum(meta_ambiental$sdi10_factor == "Expuesto", na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

# 3. CARGAR OBJETOS


ps_saliva <- readRDS(ruta_ps_saliva)
ps_ambiental <- readRDS(ruta_ps_ambiental)

meta_saliva <- preparar_metadata(ps_saliva)
meta_ambiental <- preparar_metadata(ps_ambiental)

cat("\n============================================\n")
cat("OBJETOS CARGADOS\n")
cat("============================================\n")
cat(
  "Saliva:",
  nsamples(ps_saliva),
  "muestras;",
  ntaxa(ps_saliva),
  "ASVs\n"
)
cat(
  "Ambiental:",
  nsamples(ps_ambiental),
  "muestras;",
  ntaxa(ps_ambiental),
  "ASVs\n"
)


# ------------------------------------------------------------
# 4. AGREGACIÓN TAXONÓMICA GLOBAL
# ------------------------------------------------------------

counts_phylum_saliva <- agregar_rango(
  ps_saliva,
  "Phylum"
)

counts_family_saliva <- agregar_rango(
  ps_saliva,
  "Family"
)

counts_genus_saliva <- agregar_rango(
  ps_saliva,
  "Genus"
)

rel_phylum_saliva <- a_relativa(counts_phylum_saliva)
rel_family_saliva <- a_relativa(counts_family_saliva)
rel_genus_saliva <- a_relativa(counts_genus_saliva)

res_phylum <- resumen_taxones(counts_phylum_saliva)
res_family <- resumen_taxones(counts_family_saliva)
res_genus <- resumen_taxones(counts_genus_saliva)

cat("\nTaxones agregados:\n")
cat("Phylum:", nrow(counts_phylum_saliva), "\n")
cat("Family:", nrow(counts_family_saliva), "\n")
cat("Genus:", nrow(counts_genus_saliva), "\n")


# 4B. QC DE AGREGACION TAXONOMICA 

otu_original_saliva <- otu_taxa_rows(
  ps_saliva
)

qc_agregacion <- data.frame(
  rango = c(
    "Phylum",
    "Family",
    "Genus"
  ),
  n_categorias = c(
    nrow(counts_phylum_saliva),
    nrow(counts_family_saliva),
    nrow(counts_genus_saliva)
  ),
  lecturas_originales = rep(
    sum(
      otu_original_saliva
    ),
    3
  ),
  lecturas_agregadas = c(
    sum(counts_phylum_saliva),
    sum(counts_family_saliva),
    sum(counts_genus_saliva)
  ),
  diferencia_lecturas = c(
    sum(counts_phylum_saliva) -
      sum(otu_original_saliva),
    sum(counts_family_saliva) -
      sum(otu_original_saliva),
    sum(counts_genus_saliva) -
      sum(otu_original_saliva)
  ),
  etiquetas_duplicadas = c(
    anyDuplicated(
      rownames(
        counts_phylum_saliva
      )
    ),
    anyDuplicated(
      rownames(
        counts_family_saliva
      )
    ),
    anyDuplicated(
      rownames(
        counts_genus_saliva
      )
    )
  ),
  stringsAsFactors = FALSE
)

if (
  any(
    qc_agregacion$diferencia_lecturas != 0
  )
) {
  stop(
    "QC v3: la agregacion taxonomica no conserva todas las lecturas."
  )
}

if (
  any(
    qc_agregacion$etiquetas_duplicadas != 0
  )
) {
  stop(
    "QC v3: persisten etiquetas taxonomicas duplicadas."
  )
}

# Auditoria explicita de etiquetas genericas repetidas en Genus.
tt_saliva_v3 <- tax_df(
  ps_saliva
)

auditoria_uncultured <- tt_saliva_v3 %>%
  filter(
    !is.na(Genus),
    tolower(
      trimws(
        Genus
      )
    ) == "uncultured"
  ) %>%
  select(
    taxon_id,
    any_of(
      c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus"
      )
    )
  )

write_csv(
  qc_agregacion,
  file.path(
    ruta_tablas,
    "00_QC_agregacion_taxonomica_v3.csv"
  ),
  na = ""
)

write_csv(
  auditoria_uncultured,
  file.path(
    ruta_tablas,
    "00b_auditoria_linajes_uncultured.csv"
  ),
  na = ""
)


# 5. EXPORTAR TABLAS 

write_csv(
  res_phylum,
  file.path(
    ruta_tablas,
    "01_resumen_taxonomico_phylum.csv"
  ),
  na = ""
)

write_csv(
  res_family,
  file.path(
    ruta_tablas,
    "02_resumen_taxonomico_family.csv"
  ),
  na = ""
)

write_csv(
  res_genus,
  file.path(
    ruta_tablas,
    "03_resumen_taxonomico_genus.csv"
  ),
  na = ""
)

write_csv(
  matriz_relativa_df(rel_genus_saliva),
  file.path(
    ruta_tablas,
    "04_matriz_abundancia_relativa_genero_saliva206.csv"
  ),
  na = ""
)

saveRDS(
  rel_phylum_saliva,
  file.path(
    ruta_objetos_salida,
    "matriz_relativa_phylum_saliva206.rds"
  )
)

saveRDS(
  rel_family_saliva,
  file.path(
    ruta_objetos_salida,
    "matriz_relativa_family_saliva206.rds"
  )
)

saveRDS(
  rel_genus_saliva,
  file.path(
    ruta_objetos_salida,
    "matriz_relativa_genus_saliva206.rds"
  )
)

# 6. FIGURAS 

fig_phylum <- crear_top_barplot(
  res_phylum,
  top_phylum,
  "Phylum"
)

fig_family <- crear_top_barplot(
  res_family,
  top_family,
  "Family"
)

fig_genus <- crear_top_barplot(
  res_genus,
  top_genus,
  "Genus"
)

fig_prev <- crear_prevalencia_abundancia(
  res_genus
)


# 7. TOP TAXA FIJOS PARA COMPARACIONES DESCRIPTIVAS

top_phylum_lista <- res_phylum %>%
  slice_head(n = top_phylum_grupos) %>%
  pull(taxon)

top_genus_lista <- res_genus %>%
  slice_head(n = top_genus_grupos) %>%
  pull(taxon)


# 8. COMPOSICIÓN POR ISLA - 206 SALIVA

isla_phylum <- resumen_por_grupo(
  rel = rel_phylum_saliva,
  meta = meta_saliva,
  variable_grupo = "isla_factor",
  top_taxa = top_phylum_lista,
  nombre_grupo = "Isla"
)

isla_genus <- resumen_por_grupo(
  rel = rel_genus_saliva,
  meta = meta_saliva,
  variable_grupo = "isla_factor",
  top_taxa = top_genus_lista,
  nombre_grupo = "Isla"
)

write_csv(
  isla_phylum,
  file.path(
    ruta_tablas,
    "05_composicion_phylum_por_isla.csv"
  ),
  na = ""
)

write_csv(
  isla_genus,
  file.path(
    ruta_tablas,
    "06_composicion_genero_por_isla.csv"
  ),
  na = ""
)

fig_isla_phylum <- crear_stacked_grupo(
  isla_phylum,
  "Composición taxonómica por isla - Phylum",
  "05_phylum_por_isla.png"
)

fig_isla_genus <- crear_stacked_grupo(
  isla_genus,
  "Composición taxonómica por isla - Género",
  "06_genero_por_isla.png"
)


# 9. COMPOSICIÓN POR EXACERBACIÓN 6M - 206 SALIVA

exac_phylum <- resumen_por_grupo(
  rel = rel_phylum_saliva,
  meta = meta_saliva,
  variable_grupo = "exac6m_factor",
  top_taxa = top_phylum_lista,
  nombre_grupo = "Exacerbación 6M"
)

exac_genus <- resumen_por_grupo(
  rel = rel_genus_saliva,
  meta = meta_saliva,
  variable_grupo = "exac6m_factor",
  top_taxa = top_genus_lista,
  nombre_grupo = "Exacerbación 6M"
)

write_csv(
  exac_phylum,
  file.path(
    ruta_tablas,
    "07_composicion_phylum_por_exacerbacion6M.csv"
  ),
  na = ""
)

write_csv(
  exac_genus,
  file.path(
    ruta_tablas,
    "08_composicion_genero_por_exacerbacion6M.csv"
  ),
  na = ""
)

fig_exac_phylum <- crear_stacked_grupo(
  exac_phylum,
  "Composición taxonómica por exacerbación 6M - Phylum",
  "07_phylum_por_exacerbacion6M.png"
)

fig_exac_genus <- crear_stacked_grupo(
  exac_genus,
  "Composición taxonómica por exacerbación 6M - Género",
  "08_genero_por_exacerbacion6M.png"
)


# 10. COMPOSICIÓN POR SEXO - 206 SALIVA

sexo_phylum <- resumen_por_grupo(
  rel = rel_phylum_saliva,
  meta = meta_saliva,
  variable_grupo = "sexo_factor",
  top_taxa = top_phylum_lista,
  nombre_grupo = "Sexo"
)

sexo_genus <- resumen_por_grupo(
  rel = rel_genus_saliva,
  meta = meta_saliva,
  variable_grupo = "sexo_factor",
  top_taxa = top_genus_lista,
  nombre_grupo = "Sexo"
)

write_csv(
  sexo_phylum,
  file.path(
    ruta_tablas,
    "09_composicion_phylum_por_sexo.csv"
  ),
  na = ""
)

write_csv(
  sexo_genus,
  file.path(
    ruta_tablas,
    "10_composicion_genero_por_sexo.csv"
  ),
  na = ""
)

fig_sexo_phylum <- crear_stacked_grupo(
  sexo_phylum,
  "Composición taxonómica por sexo - Phylum",
  "09_phylum_por_sexo.png"
)

fig_sexo_genus <- crear_stacked_grupo(
  sexo_genus,
  "Composición taxonómica por sexo - Género",
  "10_genero_por_sexo.png"
)


# 11. COMPOSICIÓN SEGÚN SDI/CALIMA - 205 AMBIENTAL

counts_phylum_amb <- agregar_rango(
  ps_ambiental,
  "Phylum"
)

counts_genus_amb <- agregar_rango(
  ps_ambiental,
  "Genus"
)

rel_phylum_amb <- a_relativa(counts_phylum_amb)
rel_genus_amb <- a_relativa(counts_genus_amb)

sdi_phylum <- resumen_por_grupo(
  rel = rel_phylum_amb,
  meta = meta_ambiental,
  variable_grupo = "sdi10_factor",
  top_taxa = top_phylum_lista,
  nombre_grupo = "SDI/calima 10d"
)

sdi_genus <- resumen_por_grupo(
  rel = rel_genus_amb,
  meta = meta_ambiental,
  variable_grupo = "sdi10_factor",
  top_taxa = top_genus_lista,
  nombre_grupo = "SDI/calima 10d"
)

write_csv(
  sdi_phylum,
  file.path(
    ruta_tablas,
    "11_composicion_phylum_por_SDI10.csv"
  ),
  na = ""
)

write_csv(
  sdi_genus,
  file.path(
    ruta_tablas,
    "12_composicion_genero_por_SDI10.csv"
  ),
  na = ""
)

fig_sdi_phylum <- crear_stacked_grupo(
  sdi_phylum,
  "Composición taxonómica según exposición a SDI/calima - Phylum",
  "11_phylum_por_SDI10.png"
)

fig_sdi_genus <- crear_stacked_grupo(
  sdi_genus,
  "Composición taxonómica según exposición a SDI/calima - Género",
  "12_genero_por_SDI10.png"
)


# 12. TABLA DE PREVALENCIA DE GÉNEROS

prev_generos <- res_genus %>%
  mutate(
    presente_10pct = prevalencia_prop >= 0.10,
    presente_25pct = prevalencia_prop >= 0.25,
    presente_50pct = prevalencia_prop >= 0.50,
    presente_75pct = prevalencia_prop >= 0.75
  )

write_csv(
  prev_generos,
  file.path(
    ruta_tablas,
    "13_prevalencia_generos_saliva206.csv"
  ),
  na = ""
)

# 13. RESUMEN DE COHORTE

cohorte <- resumen_cohorte(
  meta_saliva,
  meta_ambiental
)

write_csv(
  cohorte,
  file.path(
    ruta_tablas,
    "14_resumen_cohorte_composicion.csv"
  ),
  na = ""
)


# 14. ÍNDICE DE FIGURAS

indice_figuras <- data.frame(
  figura = c(
    "Top Phylum global",
    "Top Family global",
    "Top Genus global",
    "Prevalencia vs abundancia género",
    "Phylum por isla",
    "Género por isla",
    "Phylum por exacerbación 6M",
    "Género por exacerbación 6M",
    "Phylum por sexo",
    "Género por sexo",
    "Phylum por SDI10",
    "Género por SDI10"
  ),
  archivo = c(
    fig_phylum,
    fig_family,
    fig_genus,
    fig_prev,
    fig_isla_phylum,
    fig_isla_genus,
    fig_exac_phylum,
    fig_exac_genus,
    fig_sexo_phylum,
    fig_sexo_genus,
    fig_sdi_phylum,
    fig_sdi_genus
  ),
  stringsAsFactors = FALSE
)

write_csv(
  indice_figuras,
  file.path(
    ruta_tablas,
    "15_indice_figuras_composicion.csv"
  ),
  na = ""
)

# 15. EXCEL RESUMEN

notas <- data.frame(
  aspecto = c(
    "Tipo de análisis",
    "Cohorte global",
    "Cohorte ambiental",
    "Datos",
    "Abundancia relativa",
    "Phylum",
    "Family",
    "Genus",
    "Taxones no clasificados",
    "Comparaciones de grupos",
    "Inferencia"
  ),
  detalle = c(
    "Descriptivo.",
    "206 muestras de saliva.",
    "205 muestras con cruce ambiental válido.",
    "Conteos originales del objeto phyloseq sin filtrado previo de ASVs.",
    "Se calcula dentro de cada muestra después de agregar al rango taxonómico.",
    "Descripción global y por grupos.",
    "Descripción global.",
    "Descripción global, prevalencia y por grupos.",
    "En Family/Genus se conservan como categorías Unclassified usando el rango superior cuando es posible.",
    "Isla, exacerbación 6M, sexo y SDI/calima 10 días.",
    "No se calculan p-valores ni se interpreta causalidad; ANCOM-BC2 v3 es el análisis inferencial de abundancia diferencial."
  ),
  stringsAsFactors = FALSE
)

writexl::write_xlsx(
  list(
    Cohorte = cohorte,
    QC_agregacion = qc_agregacion,
    Linajes_uncultured = auditoria_uncultured,
    Phylum = res_phylum,
    Family = res_family,
    Genus = res_genus,
    Prev_genus = prev_generos,
    Isla_genus = isla_genus,
    Exac6M_genus = exac_genus,
    Sexo_genus = sexo_genus,
    SDI10_genus = sdi_genus,
    Isla_phylum = isla_phylum,
    Exac6M_phylum = exac_phylum,
    Sexo_phylum = sexo_phylum,
    SDI10_phylum = sdi_phylum,
    Figuras = indice_figuras,
    Notas = notas
  ),
  path = file.path(
    ruta_tablas,
    "16_resumen_composicion_taxonomica_GEMAS_v3.xlsx"
  )
)


# 16. MENSAJE FINAL

cat("\n============================================\n")
cat("COMPOSICIÓN TAXONÓMICA v3 FINALIZADA\n")
cat("============================================\n\n")

cat("Cohorte saliva:", nsamples(ps_saliva), "muestras\n")
cat("Cohorte ambiental:", nsamples(ps_ambiental), "muestras\n\n")

cat("Número de categorías taxonómicas:\n")
cat("Phylum:", nrow(res_phylum), "\n")
cat("Family:", nrow(res_family), "\n")
cat("Genus:", nrow(res_genus), "\n\n")

cat("QC agregación v3:\n")
print(qc_agregacion)

cat("\nTop 10 géneros por abundancia relativa media:\n")
print(
  res_genus %>%
    select(
      taxon,
      prevalencia_prop,
      abundancia_relativa_media
    ) %>%
    slice_head(n = 10)
)

cat(
  "\nExcel resumen:\n",
  file.path(
    ruta_tablas,
    "16_resumen_composicion_taxonomica_GEMAS_v3.xlsx"
  ),
  "\n",
  sep = ""
)
