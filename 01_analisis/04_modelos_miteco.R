# Ejecutar los modelos finales de MITECO por bloques.
# El orden mantiene la dependencia entre metadatos, métricas y modelos.

raiz_codigo <- Sys.getenv("TFM_CODE_ROOT")
if (!nzchar(raiz_codigo)) stop("Ejecutar desde 02_ejemplo_sintetico/ejecutar.R.")
source(file.path(raiz_codigo, "01_analisis/modulos/rng.R"))
modulos <- c(
  "miteco_configuracion.R", "miteco_funciones.R", "miteco_integracion.R",
  "diversidad_alfa.R", "diversidad_beta_secundaria.R", "compatibilidad_aleatoria.R",
  "nucleo_composicion.R", "ancombc2.R", "exacerbaciones.R", "exportar_miteco.R",
  "tabla_maestra.R")
for (archivo in modulos) {
  message("  ", archivo)
  sys.source(file.path(raiz_codigo, "01_analisis/modulos", archivo), envir = globalenv())
}
