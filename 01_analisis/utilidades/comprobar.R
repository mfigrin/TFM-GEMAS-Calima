# Revisar dependencias y sintaxis sin instalar ni actualizar paquetes.
paquetes <- c("phyloseq", "DESeq2", "ANCOMBC", "ape", "vegan", "permute", "dplyr", "tidyr",
              "readr", "tibble", "ggplot2", "writexl", "sandwich", "lmtest", "stringr",
              "readxl", "lubridate", "scales", "patchwork", "purrr")
instalados <- vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
version <- vapply(paquetes, function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
}, character(1))
print(data.frame(paquete = paquetes, instalado = instalados, version), row.names = FALSE)
if (any(!instalados)) stop("Faltan paquetes: ", paste(paquetes[!instalados], collapse = ", "))
archivos <- unlist(lapply(c("01_analisis", "02_ejemplo_sintetico", "04_datos_ambientales/procesar"),
  list.files, pattern = "[.]R$", recursive = TRUE, full.names = TRUE), use.names = FALSE)
archivos <- archivos[!grepl("/trabajo/", archivos)]
for (f in archivos) parse(f, encoding = "UTF-8")
message("Sintaxis comprobada en ", length(archivos), " archivos activos. No se han ajustado modelos.")
