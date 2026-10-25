# Revisar la estructura de entrada sin repetir los cálculos del análisis.
local({
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Falta phyloseq.")
  work <- file.path(getwd(), "02_ejemplo_sintetico/trabajo")
  if (!file.exists(file.path(work, "SYNTHETIC_ONLY.txt"))) stop("Falta el marcador sintético.")
  ps <- readRDS(file.path(work, "phyloseq_sintetico.rds"))
  stopifnot(inherits(ps, "phyloseq"), methods::validObject(ps))
  m <- as(phyloseq::sample_data(ps), "data.frame")
  if (!"synthetic" %in% names(m) || !isTRUE(all(as.logical(m$synthetic))) ||
      !all(grepl("^SYN_", rownames(m)))) stop("La entrada no está identificada como sintética.")

  # Los conteos, la taxonomía y el árbol deben referirse a las mismas muestras y ASV.
  ct <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) ct <- t(ct)
  tx <- as(phyloseq::tax_table(ps), "matrix")
  tree <- phyloseq::phy_tree(ps)
  stopifnot(all(is.finite(ct)), all(ct >= 0), all(ct == floor(ct)),
            identical(colnames(ct), rownames(m)),
            setequal(rownames(ct), rownames(tx)),
            setequal(rownames(ct), tree$tip.label),
            length(tree$edge.length) == nrow(tree$edge),
            all(is.finite(tree$edge.length)), all(tree$edge.length > 0))

  # La fixture conserva los tamaños de las ramas para probar las verificaciones originales.
  saliva <- m$sample.type == "Saliva"
  ambiental <- saliva & !is.na(m$Date)
  stopifnot(sum(saliva) == 206L, sum(ambiental) == 205L,
            all(colSums(ct)[ambiental] >= 10000))
  capture.output(sessionInfo(), file = file.path(work, "logs/sessionInfo_entrada.txt"))
  message("Entrada sintética comprobada. Las distancias se calcularán en la etapa 02.")
})
