# Compatibilidad de la secuencia aleatoria al retirar ajustes redundantes.
# Solo genera las permutaciones que consumía el ajuste omitido. No calcula un contraste.
# La equivalencia depende de las versiones de vegan/permute; ejecutar 01_analisis/utilidades/verificar_aleatoriedad.R.

tfm_avanzar_permutaciones <- function(nperm, n) {
  invisible(permute::shuffleSet(n, control = permute::how(nperm = nperm)))
}
