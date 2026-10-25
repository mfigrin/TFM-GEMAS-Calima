# Modelos complementarios de PM10 y tablas de resultados.
# Reutiliza la prueba principal MITECO de 06; no vuelve a ajustarla.
raiz_codigo <- Sys.getenv("TFM_CODE_ROOT")
if (!nzchar(raiz_codigo)) stop("Ejecutar desde 02_ejemplo_sintetico/ejecutar.R.")
for (archivo in c("pm10_preparacion.R", "pm10_descripcion.R", "pm10_modelos.R",
                  "pm10_sintesis.R", "pm10_exportacion.R")) {
  message("  ", archivo)
  sys.source(file.path(raiz_codigo, "01_analisis/modulos", archivo), envir = globalenv())
}
