# Comprobar en la versión instalada que el avance RNG sustituye al modelo omitido.
# Si falla, no se ejecuta el flujo depurado: hay que revisar la compatibilidad con vegan.
local({
  source("01_analisis/modulos/rng.R", local = TRUE)
  n <- 40L
  x <- matrix((seq_len(n * 7L) %% 19L) + 1, nrow = n)
  d <- vegan::vegdist(x, method = "bray")
  meta <- data.frame(grupo = factor(rep(c("A", "B"), each = n / 2)), edad = seq_len(n))
  set.seed(20260826)
  invisible(vegan::adonis2(d ~ grupo + edad, data = meta, permutations = 999, by = "margin"))
  estado_adonis <- .Random.seed
  set.seed(20260826)
  tfm_avanzar_permutaciones(999, n)
  if (!identical(estado_adonis, .Random.seed))
    stop("El avance RNG no coincide con adonis2 en esta instalación. No continuar.")
  bd <- vegan::betadisper(d, meta$grupo)
  set.seed(20260826)
  invisible(vegan::permutest(bd, permutations = 999))
  estado_disp <- .Random.seed
  set.seed(20260826)
  tfm_avanzar_permutaciones(999, n)
  if (!identical(estado_disp, .Random.seed))
    stop("El avance RNG no coincide con permutest.betadisper. No continuar.")
  message("Compatibilidad RNG comprobada para adonis2 y PERMDISP en esta instalación.")
})
