# Calcular UniFrac una sola vez.
# UniFrac ponderado usa conteos sin rarefacción.
# UniFrac no ponderado usa 10.000 lecturas, semilla 20260822 y muestreo sin reemplazo.

paquetes <- c("phyloseq", "ape", "vegan", "dplyr", "tibble", "readr", "ggplot2", "writexl")

faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]

if (length(faltan) > 0) {
  stop("Faltan paquetes: ", paste(faltan, collapse = ", "), ". Instálalos antes de ejecutar el script.")
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(ape)
  library(vegan)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(writexl)
})

set.seed(20260822)

profundidad_unweighted <- 10000

ruta_ps_ambiental <- file.path(
  "outputs", "analisis_microbiota", "objetos_phyloseq",
  "ps_ambiental_205_filtrado_prev2_total10.rds"
)

ruta_ps_saliva <- file.path(
  "outputs", "analisis_microbiota", "objetos_phyloseq",
  "ps_saliva_206_filtrado_prev2_total10.rds"
)

ruta_salida <- file.path("outputs", "analisis_microbiota", "unifrac_sensibilidad")

ruta_tablas <- file.path(ruta_salida, "tablas")
ruta_figuras <- file.path(ruta_salida, "figuras")
ruta_objetos <- file.path(ruta_salida, "objetos")

dir.create(ruta_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_objetos, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ruta_ps_ambiental)) {
  stop("No existe: ", ruta_ps_ambiental)
}

if (!file.exists(ruta_ps_saliva)) {
  stop("No existe: ", ruta_ps_saliva)
}

resolver_arbol <- function(ps_obj, nombre_dataset) {
  arbol_original <- phy_tree(ps_obj)

  qc_original <- data.frame(
    dataset = nombre_dataset,
    version_arbol = "original",
    n_tips = length(arbol_original$tip.label),
    n_nodos = arbol_original$Nnode,
    binario = ape::is.binary(arbol_original),
    enraizado = ape::is.rooted(arbol_original),
    tips_coinciden_asv = setequal(taxa_names(ps_obj), arbol_original$tip.label),
    tips_duplicados = anyDuplicated(arbol_original$tip.label),
    ramas_na = sum(is.na(arbol_original$edge.length)),
    ramas_cero = sum(arbol_original$edge.length == 0, na.rm = TRUE),
    ramas_negativas = sum(arbol_original$edge.length < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  arbol_binario <- ape::multi2di(arbol_original, random = FALSE)

  qc_binario <- data.frame(
    dataset = nombre_dataset,
    version_arbol = "binario_multi2di",
    n_tips = length(arbol_binario$tip.label),
    n_nodos = arbol_binario$Nnode,
    binario = ape::is.binary(arbol_binario),
    enraizado = ape::is.rooted(arbol_binario),
    tips_coinciden_asv = setequal(taxa_names(ps_obj), arbol_binario$tip.label),
    tips_duplicados = anyDuplicated(arbol_binario$tip.label),
    ramas_na = sum(is.na(arbol_binario$edge.length)),
    ramas_cero = sum(arbol_binario$edge.length == 0, na.rm = TRUE),
    ramas_negativas = sum(arbol_binario$edge.length < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (!ape::is.binary(arbol_binario)) {
    stop("El árbol no quedó binario en ", nombre_dataset)
  }

  if (!ape::is.rooted(arbol_binario)) {
    stop("El árbol no está enraizado en ", nombre_dataset)
  }

  if (!setequal(taxa_names(ps_obj), arbol_binario$tip.label)) {
    stop("ASVs y tips no coinciden en ", nombre_dataset)
  }

  if (any(is.na(arbol_binario$edge.length))) {
    stop("Hay longitudes de rama NA en ", nombre_dataset)
  }

  if (any(arbol_binario$edge.length < 0)) {
    stop("Hay longitudes de rama negativas en ", nombre_dataset)
  }

  ps_out <- ps_obj
  phy_tree(ps_out) <- arbol_binario

  list(ps = ps_out, qc = bind_rows(qc_original, qc_binario))
}

control_distancia <- function(dist_obj, dataset, distancia) {
  m <- as.matrix(dist_obj)

  data.frame(
    dataset = dataset,
    distancia = distancia,
    n_muestras = nrow(m),
    n_na = sum(is.na(m)),
    n_no_finitos = sum(!is.finite(m)),
    diagonal_max_abs = max(abs(diag(m))),
    asimetria_max_abs = max(abs(m - t(m)), na.rm = TRUE),
    distancia_min = min(m, na.rm = TRUE),
    distancia_max = max(m, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# Cargar objetos.

ps_ambiental <- readRDS(ruta_ps_ambiental)
ps_saliva <- readRDS(ruta_ps_saliva)

cat("OBJETOS CARGADOS\n")
cat("Ambiental:", nsamples(ps_ambiental), "muestras;", ntaxa(ps_ambiental), "ASVs\n")
cat("Saliva:", nsamples(ps_saliva), "muestras;", ntaxa(ps_saliva), "ASVs\n")

# Resolver politomias de forma determinista.

res_arbol_ambiental <- resolver_arbol(ps_ambiental, "ambiental_205")

res_arbol_saliva <- resolver_arbol(ps_saliva, "saliva_206")

ps_ambiental_bin <- res_arbol_ambiental$ps
ps_saliva_bin <- res_arbol_saliva$ps

qc_arbol <- bind_rows(res_arbol_ambiental$qc, res_arbol_saliva$qc)

readr::write_csv(qc_arbol, file.path(ruta_tablas, "01_qc_arboles_unifrac.csv"), na = "")

saveRDS(ps_ambiental_bin, file.path(ruta_objetos, "ps_ambiental_205_arbol_binario.rds"))

saveRDS(ps_saliva_bin, file.path(ruta_objetos, "ps_saliva_206_arbol_binario.rds"))

# Weighted UniFrac.

cat("\nCalculando Weighted UniFrac ambiental...\n")

dist_w_ambiental <- phyloseq::UniFrac(
  ps_ambiental_bin,
  weighted = TRUE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

cat("Weighted UniFrac ambiental: OK\n")

cat("\nCalculando Weighted UniFrac saliva...\n")

dist_w_saliva <- phyloseq::UniFrac(
  ps_saliva_bin,
  weighted = TRUE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

cat("Weighted UniFrac saliva: OK\n")

saveRDS(dist_w_ambiental, file.path(ruta_objetos, "dist_weighted_unifrac_ambiental_205.rds"))

saveRDS(dist_w_saliva, file.path(ruta_objetos, "dist_weighted_unifrac_saliva_206.rds"))

# Unweighted UniFrac con rarefaccion a 10.000.

min_ambiental <- min(sample_sums(ps_ambiental_bin))
min_saliva <- min(sample_sums(ps_saliva_bin))

if (min_ambiental < profundidad_unweighted) {
  stop(
    "Hay muestras ambientales con menos de ",
    profundidad_unweighted,
    " lecturas. Mínimo: ",
    min_ambiental
  )
}

if (min_saliva < profundidad_unweighted) {
  stop(
    "Hay muestras de saliva con menos de ",
    profundidad_unweighted,
    " lecturas. Mínimo: ",
    min_saliva
  )
}

cat(
  "\nRarefaccion para Unweighted UniFrac a ",
  profundidad_unweighted,
  " lecturas/muestra...\n",
  sep = ""
)

set.seed(20260822)

ps_ambiental_rare <- phyloseq::rarefy_even_depth(
  ps_ambiental_bin,
  sample.size = profundidad_unweighted,
  rngseed = 20260822,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

set.seed(20260822)

ps_saliva_rare <- phyloseq::rarefy_even_depth(
  ps_saliva_bin,
  sample.size = profundidad_unweighted,
  rngseed = 20260822,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

# La poda posterior a la rarefaccion puede reintroducir politomias.
ps_ambiental_rare <- resolver_arbol(ps_ambiental_rare, "ambiental_205_rarefactado")$ps

ps_saliva_rare <- resolver_arbol(ps_saliva_rare, "saliva_206_rarefactado")$ps

cat("\nCalculando Unweighted UniFrac ambiental...\n")

dist_u_ambiental <- phyloseq::UniFrac(
  ps_ambiental_rare,
  weighted = FALSE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

cat("Unweighted UniFrac ambiental: OK\n")

cat("\nCalculando Unweighted UniFrac saliva...\n")

dist_u_saliva <- phyloseq::UniFrac(
  ps_saliva_rare,
  weighted = FALSE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

cat("Unweighted UniFrac saliva: OK\n")

saveRDS(
  dist_u_ambiental,
  file.path(ruta_objetos, "dist_unweighted_unifrac_ambiental_205_rarefactado_10000.rds")
)

saveRDS(
  dist_u_saliva,
  file.path(ruta_objetos, "dist_unweighted_unifrac_saliva_206_rarefactado_10000.rds")
)

# Control de calidad de distancias.

qc_distancias <- bind_rows(
  control_distancia(dist_w_ambiental, "ambiental_205", "Weighted_UniFrac"),
  control_distancia(dist_u_ambiental, "ambiental_205", "Unweighted_UniFrac_rarefactado_10000"),
  control_distancia(dist_w_saliva, "saliva_206", "Weighted_UniFrac"),
  control_distancia(dist_u_saliva, "saliva_206", "Unweighted_UniFrac_rarefactado_10000")
)

readr::write_csv(qc_distancias, file.path(ruta_tablas, "02_qc_distancias_unifrac.csv"), na = "")
