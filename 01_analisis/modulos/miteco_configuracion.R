# Configurar el análisis MITECO y localizar sus entradas.
# La exposición principal es el calendario oficial; PM10 es complementario.

set.seed(20260826)

N_PERM <- 999
ALPHA_FDR <- 0.05
PREV_ANCOM <- 0.10
PSEUDOCOUNT_CLR <- 0.5

# Si ANCOM-BC2 a ASV resulta demasiado lento en el primer ensayo,
# se puede poner FALSE para obtener el borrador preliminar y
# reactivarlo antes de la entrega final.
EJECUTAR_ANCOM_ASV <- TRUE

# Si no existen los objetos UniFrac creados previamente, el script
# continuara sin detenerse y lo registrara en el manifest.
EJECUTAR_UNIFRAC_SI_EXISTE <- TRUE

# Paquetes.

paquetes <- c(
  "phyloseq",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "ggplot2",
  "vegan",
  "ape",
  "writexl",
  "sandwich",
  "lmtest",
  "ANCOMBC",
  "splines"
)

faltan <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltan) > 0) {
  stop(
    "Faltan paquetes: ",
    paste(faltan, collapse = ", "),
    "\nInstalalos antes de ejecutar el reanalisis."
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(vegan)
  library(ape)
  library(writexl)
  library(sandwich)
  library(lmtest)
  library(ANCOMBC)
})

# Rutas.

# Los scripts previos del proyecto no han usado siempre exactamente
# la misma subcarpeta para guardar los objetos. Para evitar depender
# de una ruta historica concreta, localizamos cada archivo por su
# nombre exacto dentro de outputs/ de forma recursiva.

localizar_archivo_unico <- function(
  nombre,
  raiz = "outputs",
  obligatorio = TRUE
) {
  encontrados <- list.files(
    raiz,
    pattern = paste0("^", gsub("([.(){}+*?^$|\\\\])", "\\\\\\1", nombre), "$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )

  if (length(encontrados) == 1) {
    cat("Encontrado: ", nombre, " -> ", encontrados, "\n", sep = "")

    return(encontrados)
  }

  if (length(encontrados) > 1) {
    stop(
      "Se encontraron varias copias de '",
      nombre,
      "'. No puedo elegir de forma segura:\n- ",
      paste(encontrados, collapse = "\n- ")
    )
  }

  if (obligatorio) {
    stop(
      "No se encontro '",
      nombre,
      "' dentro de ",
      raiz,
      "/.\n",
      "Ejecuta primero el script que genera ese objeto o comprueba su nombre real."
    )
  }

  NA_character_
}

ruta_ps_amb_f <- localizar_archivo_unico("ps_ambiental_205_filtrado_prev2_total10.rds")

ruta_ps_amb_u <- localizar_archivo_unico("ps_ambiental_205_sin_filtrado.rds")

ruta_ps_saliva_u <- localizar_archivo_unico("ps_saliva_206_sin_filtrado.rds")

ruta_miteco_muestra <- localizar_archivo_unico("04_exposicion_MITECO_PM10_por_muestra.csv")

ruta_miteco_traza <- localizar_archivo_unico("02_trazabilidad_MITECO_PM10_D1_D180.csv")

ruta_wunifrac <- localizar_archivo_unico("dist_weighted_unifrac_ambiental_205.rds", obligatorio = TRUE)

ruta_uunifrac <- localizar_archivo_unico(
  "dist_unweighted_unifrac_ambiental_205_rarefactado_10000.rds",
  obligatorio = TRUE
)

ruta_salida <- file.path("outputs", "analisis_microbiota", "reanalisis_completo_MITECO")

ruta_tablas <- file.path(ruta_salida, "tablas")

ruta_figuras <- file.path(ruta_salida, "figuras")

ruta_objetos <- file.path(ruta_salida, "objetos")

for (r in c(
  ruta_salida,
  ruta_tablas,
  ruta_figuras,
  ruta_objetos
)) {
  dir.create(r, recursive = TRUE, showWarnings = FALSE)
}
