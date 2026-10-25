# PREPARACION DE OBJETOS PHYLOSEQ ANALITICOS
# GEMAS CANARIAS - MICROBIOTA, PM10/CALIMA, EXACERBACIONES Y SEXO
#
# Objetivo:
# - Crear objetos phyloseq limpios para los analisis posteriores.
# - Separar muestras clinicas de saliva y controles.
# - Unir metadatos ambientales derivados de PM10/calima.
# - Crear un dataset clinico completo de saliva (n esperado = 206).
# - Crear un dataset ambiental con fecha valida (n esperado = 205).
# - Crear versiones sin filtrado y con filtrado conservador de ASVs.
# - Exportar metadatos, control de calidad y resumen de objetos.


# 0. PAQUETES

paquetes_cran <- c(
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "writexl"
)

paquetes_faltantes <- paquetes_cran[
  !vapply(
    paquetes_cran,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(paquetes_faltantes) > 0) {
  install.packages(paquetes_faltantes)
}

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  stop(
    "No esta instalado phyloseq. Instala con: BiocManager::install('phyloseq')"
  )
}

library(phyloseq)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(writexl)


# 1. RUTAS

ruta_phyloseq <- "data_raw/data_phyloseq_salivas_GEMAS_Canarias.rds"
ruta_macrotabla <- "outputs/17_macrotabla_clinica_ambiental_por_muestra.csv"

ruta_salida <- "outputs/analisis_microbiota"
ruta_objetos <- file.path(ruta_salida, "objetos_phyloseq")
ruta_tablas <- file.path(ruta_salida, "tablas_qc")

if (!file.exists(ruta_phyloseq)) {
  stop(
    paste0(
      "No se encuentra el objeto phyloseq en: ",
      ruta_phyloseq,
      "\nComprueba getwd() y la ubicacion del archivo."
    )
  )
}

if (!file.exists(ruta_macrotabla)) {
  stop(
    paste0(
      "No se encuentra la macrotabla ambiental en: ",
      ruta_macrotabla,
      "\nEjecuta antes el script 04."
    )
  )
}

for (ruta in c(ruta_salida, ruta_objetos, ruta_tablas)) {
  if (!dir.exists(ruta)) {
    dir.create(ruta, recursive = TRUE)
  }
}


# 2. FUNCIONES 

extraer_sample_data <- function(ps_obj) {
  sd_obj <- sample_data(ps_obj)
  sd_df <- data.frame(
    setNames(
      sd_obj@.Data,
      names(sd_obj)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(sd_df) <- sample_names(ps_obj)
  sd_df$id_muestra <- rownames(sd_df)
  sd_df <- sd_df %>% relocate(id_muestra)
  sd_df
}

convertir_numero <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

convertir_logico <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  x <- toupper(trimws(as.character(x)))
  salida <- rep(NA, length(x))
  salida[x %in% c("TRUE", "T", "1", "SI", "SÍ", "YES", "Y", "EXPUESTO")] <- TRUE
  salida[x %in% c("FALSE", "F", "0", "NO", "N", "NO EXPUESTO", "NO_EXPUESTO")] <- FALSE
  as.logical(salida)
}

normalizar_isla <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x[x %in% c("LaPalma", "La palma", "LA PALMA")] <- "La Palma"
  x[x %in% c("TENERIFE", "tenerife")] <- "Tenerife"
  x
}

normalizar_sexo <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("F", "Female", "female", "Mujer", "mujer")] <- "Female"
  x[x %in% c("M", "Male", "male", "Hombre", "hombre")] <- "Male"
  x
}

resumen_phyloseq <- function(ps_obj, nombre_objeto) {
  profundidad <- sample_sums(ps_obj)
  otu <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) {
    prevalencia <- rowSums(otu > 0)
  } else {
    prevalencia <- colSums(otu > 0)
  }
  data.frame(
    objeto = nombre_objeto,
    n_muestras = nsamples(ps_obj),
    n_asv = ntaxa(ps_obj),
    lecturas_totales = sum(profundidad),
    lecturas_media = round(mean(profundidad), 1),
    lecturas_sd = round(sd(profundidad), 1),
    lecturas_mediana = round(median(profundidad), 1),
    lecturas_p25 = round(as.numeric(quantile(profundidad, 0.25)), 1),
    lecturas_p75 = round(as.numeric(quantile(profundidad, 0.75)), 1),
    lecturas_min = min(profundidad),
    lecturas_max = max(profundidad),
    prevalencia_asv_mediana = round(median(prevalencia), 1),
    prevalencia_asv_p75 = round(as.numeric(quantile(prevalencia, 0.75)), 1),
    stringsAsFactors = FALSE
  )
}

calcular_tabla_profundidad <- function(ps_obj, nombre_objeto) {
  profundidad <- sample_sums(ps_obj)
  otu <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) {
    asv_observadas <- colSums(otu > 0)
  } else {
    asv_observadas <- rowSums(otu > 0)
  }
  data.frame(
    objeto = nombre_objeto,
    id_muestra = names(profundidad),
    lecturas_totales = as.numeric(profundidad),
    asv_observadas = as.numeric(asv_observadas),
    stringsAsFactors = FALSE
  )
}

filtrar_asv_conservador <- function(
    ps_obj,
    min_prevalencia = 2,
    min_abundancia_total = 10
) {
  otu <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) {
    prevalencia <- rowSums(otu > 0)
    abundancia_total <- rowSums(otu)
  } else {
    prevalencia <- colSums(otu > 0)
    abundancia_total <- colSums(otu)
  }
  mantener <- prevalencia >= min_prevalencia &
    abundancia_total >= min_abundancia_total
  prune_taxa(mantener, ps_obj)
}

prevalencia_taxa <- function(ps_obj) {
  otu <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) {
    prevalencia <- rowSums(otu > 0)
    abundancia_total <- rowSums(otu)
  } else {
    prevalencia <- colSums(otu > 0)
    abundancia_total <- colSums(otu)
  }
  data.frame(
    ASV = names(prevalencia),
    prevalencia_muestras = as.numeric(prevalencia),
    abundancia_total = as.numeric(abundancia_total),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(prevalencia_muestras), desc(abundancia_total))
}

actualizar_sample_data <- function(ps_obj, metadata_df) {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  rownames(metadata_df) <- metadata_df$id_muestra
  metadata_df <- metadata_df[sample_names(ps_obj), , drop = FALSE]
  sample_data(ps_obj) <- sample_data(metadata_df)
  ps_obj
}

# 3. CARGAR DATOS

ps <- readRDS(ruta_phyloseq)

cat("\n============================================\n")
cat("OBJETO PHYLOSEQ ORIGINAL\n")
cat("============================================\n")
print(ps)

metadata_original <- extraer_sample_data(ps)

# Normalizar algunas variables del phyloseq
metadata_original <- metadata_original %>%
  mutate(
    Age = convertir_numero(Age),
    Sex.number = convertir_numero(Sex.number),
    Sex.letter = normalizar_sexo(Sex.letter),
    Exacerbations.12M = convertir_numero(Exacerbations.12M),
    Antibiotic.2mo = convertir_numero(Antibiotic.2mo),
    Pool = as.character(Pool),
    Cavities = convertir_numero(Cavities),
    Chewing.gum.30min = convertir_numero(Chewing.gum.30min),
    Exacerbations.6M = convertir_numero(Exacerbations.6M),
    Date = as.Date(Date),
    Island = normalizar_isla(Island)
  )

ps <- actualizar_sample_data(ps, metadata_original)

macrotabla <- read_csv2(
  ruta_macrotabla,
  show_col_types = FALSE
)

# Asegurar id y nombres basicos
if (!"id_muestra" %in% names(macrotabla)) {
  stop("La macrotabla no contiene id_muestra.")
}

# Seleccionar solo variables ambientales y algunas variables de enlace.
# Evitamos duplicar edad/sexo/isla clinicas, que ya proceden del phyloseq.
variables_ambientales_interes <- c(
  "id_muestra",
  "fecha_muestra",
  "adulto",
  "pm10_d01",
  "pm10_d02",
  "pm10_d03",
  "pm10_d04",
  "pm10_d05",
  "pm10_d06",
  "pm10_d07",
  "pm10_d08",
  "pm10_d09",
  "pm10_d10",
  "pm10_media_10d",
  "pm10_maximo_10d",
  "pm10_acumulado_10d",
  "pm10_auc50_10d",
  "n_dias_sdi_10d",
  "sdi_exposicion_10d",
  "grado_sdi_10d",
  "sdi_temprana_10d",
  "sdi_tardia_10d",
  "patron_temporal_sdi_10d",
  "tipo_sdi_estacional_10d",
  "pm10_media_180d",
  "pm10_maximo_180d",
  "pm10_acumulado_180d",
  "pm10_auc50_180d",
  "n_dias_sdi_180d",
  "sdi_exposicion_180d",
  "grado_sdi_180d"
)

variables_ambientales_presentes <- intersect(
  variables_ambientales_interes,
  names(macrotabla)
)

macro_amb <- macrotabla %>%
  select(all_of(variables_ambientales_presentes))

# Normalizar variables ambientales clave si existen
for (v in intersect(
  c(
    "pm10_d01", "pm10_d02", "pm10_d03", "pm10_d04", "pm10_d05",
    "pm10_d06", "pm10_d07", "pm10_d08", "pm10_d09", "pm10_d10",
    "pm10_media_10d", "pm10_maximo_10d", "pm10_acumulado_10d",
    "pm10_auc50_10d", "n_dias_sdi_10d", "pm10_media_180d",
    "pm10_maximo_180d", "pm10_acumulado_180d", "pm10_auc50_180d",
    "n_dias_sdi_180d"
  ),
  names(macro_amb)
)) {
  macro_amb[[v]] <- convertir_numero(macro_amb[[v]])
}

for (v in intersect(
  c(
    "sdi_exposicion_10d", "sdi_temprana_10d", "sdi_tardia_10d",
    "sdi_exposicion_180d", "adulto"
  ),
  names(macro_amb)
)) {
  macro_amb[[v]] <- convertir_logico(macro_amb[[v]])
}

if ("fecha_muestra" %in% names(macro_amb)) {
  macro_amb$fecha_muestra <- as.Date(macro_amb$fecha_muestra)
}

# 4. SEPARAR CONTROLES Y SALIVA

if (!"sample.type" %in% sample_variables(ps)) {
  stop("No existe sample.type en sample_data(ps). No se pueden separar saliva y controles.")
}

ps_saliva_206 <- subset_samples(ps, sample.type == "Saliva")
ps_controles <- subset_samples(ps, sample.type != "Saliva")

ps_saliva_206 <- prune_taxa(taxa_sums(ps_saliva_206) > 0, ps_saliva_206)
ps_controles <- prune_taxa(taxa_sums(ps_controles) > 0, ps_controles)

cat("\nMuestras saliva: ", nsamples(ps_saliva_206), "\n", sep = "")
cat("Muestras controles: ", nsamples(ps_controles), "\n", sep = "")


# 5. UNIR METADATOS AMBIENTALES AL OBJETO DE SALIVA

metadata_saliva <- extraer_sample_data(ps_saliva_206)

metadata_saliva <- metadata_saliva %>%
  mutate(
    Age = convertir_numero(Age),
    Sex.number = convertir_numero(Sex.number),
    Sex.letter = normalizar_sexo(Sex.letter),
    Exacerbations.12M = convertir_numero(Exacerbations.12M),
    Antibiotic.2mo = convertir_numero(Antibiotic.2mo),
    Pool = as.character(Pool),
    Cavities = convertir_numero(Cavities),
    Chewing.gum.30min = convertir_numero(Chewing.gum.30min),
    Exacerbations.6M = convertir_numero(Exacerbations.6M),
    Date = as.Date(Date),
    Island = normalizar_isla(Island)
  )

metadata_saliva_amb <- metadata_saliva %>%
  left_join(
    macro_amb,
    by = "id_muestra"
  )

metadata_saliva_amb <- metadata_saliva_amb %>%
  mutate(
    tiene_cruce_ambiental = !is.na(sdi_exposicion_10d)
  )

ids_sin_cruce_ambiental <- metadata_saliva_amb %>%
  filter(!tiene_cruce_ambiental) %>%
  select(
    id_muestra,
    Age,
    Sex.letter,
    Exacerbations.6M,
    Exacerbations.12M,
    Date,
    Island
  )

ps_saliva_206 <- actualizar_sample_data(ps_saliva_206, metadata_saliva_amb)

ps_ambiental_205 <- subset_samples(
  ps_saliva_206,
  tiene_cruce_ambiental == TRUE
)

ps_ambiental_205 <- prune_taxa(
  taxa_sums(ps_ambiental_205) > 0,
  ps_ambiental_205
)

cat("\nMuestras saliva con cruce ambiental: ", nsamples(ps_ambiental_205), "\n", sep = "")
cat("Muestras saliva sin cruce ambiental: ", nrow(ids_sin_cruce_ambiental), "\n", sep = "")


# 6. FILTRADO CONSERVADOR

min_prevalencia_asv <- 2
min_abundancia_total_asv <- 10

ps_saliva_206_filtrado <- filtrar_asv_conservador(
  ps_saliva_206,
  min_prevalencia = min_prevalencia_asv,
  min_abundancia_total = min_abundancia_total_asv
)

ps_ambiental_205_filtrado <- filtrar_asv_conservador(
  ps_ambiental_205,
  min_prevalencia = min_prevalencia_asv,
  min_abundancia_total = min_abundancia_total_asv
)

ps_saliva_206_filtrado <- prune_taxa(
  taxa_sums(ps_saliva_206_filtrado) > 0,
  ps_saliva_206_filtrado
)

ps_ambiental_205_filtrado <- prune_taxa(
  taxa_sums(ps_ambiental_205_filtrado) > 0,
  ps_ambiental_205_filtrado
)

# 7. GUARDAR OBJETOS

saveRDS(
  ps_saliva_206,
  file.path(ruta_objetos, "ps_saliva_206_sin_filtrado.rds")
)

saveRDS(
  ps_ambiental_205,
  file.path(ruta_objetos, "ps_ambiental_205_sin_filtrado.rds")
)

saveRDS(
  ps_saliva_206_filtrado,
  file.path(
    ruta_objetos,
    "ps_saliva_206_filtrado_prev2_total10.rds"
  )
)

saveRDS(
  ps_ambiental_205_filtrado,
  file.path(
    ruta_objetos,
    "ps_ambiental_205_filtrado_prev2_total10.rds"
  )
)

saveRDS(
  ps_controles,
  file.path(ruta_objetos, "ps_controles_mock_ntc.rds")
)

# 8. TABLAS DE CONTROL DE CALIDAD

resumen_objetos <- bind_rows(
  resumen_phyloseq(ps, "ps_original_218"),
  resumen_phyloseq(ps_saliva_206, "ps_saliva_206_sin_filtrado"),
  resumen_phyloseq(ps_ambiental_205, "ps_ambiental_205_sin_filtrado"),
  resumen_phyloseq(ps_saliva_206_filtrado, "ps_saliva_206_filtrado_prev2_total10"),
  resumen_phyloseq(ps_ambiental_205_filtrado, "ps_ambiental_205_filtrado_prev2_total10"),
  resumen_phyloseq(ps_controles, "ps_controles_mock_ntc")
)

write_csv2(
  resumen_objetos,
  file.path(ruta_tablas, "01_resumen_objetos_phyloseq_analiticos.csv"),
  na = ""
)

profundidad_todos <- bind_rows(
  calcular_tabla_profundidad(ps_saliva_206, "ps_saliva_206_sin_filtrado"),
  calcular_tabla_profundidad(ps_ambiental_205, "ps_ambiental_205_sin_filtrado"),
  calcular_tabla_profundidad(ps_saliva_206_filtrado, "ps_saliva_206_filtrado_prev2_total10"),
  calcular_tabla_profundidad(ps_ambiental_205_filtrado, "ps_ambiental_205_filtrado_prev2_total10"),
  calcular_tabla_profundidad(ps_controles, "ps_controles_mock_ntc")
)

write_csv2(
  profundidad_todos,
  file.path(ruta_tablas, "02_profundidad_por_muestra_objetos.csv"),
  na = ""
)

metadata_saliva_final <- extraer_sample_data(ps_saliva_206)
metadata_ambiental_final <- extraer_sample_data(ps_ambiental_205)

write_csv2(
  metadata_saliva_final,
  file.path(ruta_tablas, "03_metadata_saliva_206_con_ambiental.csv"),
  na = ""
)

write_csv2(
  metadata_ambiental_final,
  file.path(ruta_tablas, "04_metadata_ambiental_205.csv"),
  na = ""
)

write_csv2(
  ids_sin_cruce_ambiental,
  file.path(ruta_tablas, "05_muestras_saliva_sin_cruce_ambiental.csv"),
  na = ""
)

prevalencia_asv_saliva <- prevalencia_taxa(ps_saliva_206)
prevalencia_asv_ambiental <- prevalencia_taxa(ps_ambiental_205)

write_csv2(
  prevalencia_asv_saliva,
  file.path(ruta_tablas, "06_prevalencia_asv_saliva_206.csv"),
  na = ""
)

write_csv2(
  prevalencia_asv_ambiental,
  file.path(ruta_tablas, "07_prevalencia_asv_ambiental_205.csv"),
  na = ""
)

# Resumen de variables clave en los objetos analiticos
resumen_variables_clave <- metadata_saliva_final %>%
  summarise(
    n_saliva = n(),
    n_con_cruce_ambiental = sum(tiene_cruce_ambiental %in% TRUE, na.rm = TRUE),
    n_sin_cruce_ambiental = sum(!tiene_cruce_ambiental %in% TRUE, na.rm = TRUE),
    edad_media = round(mean(Age, na.rm = TRUE), 1),
    edad_sd = round(sd(Age, na.rm = TRUE), 1),
    n_female = sum(Sex.letter == "Female", na.rm = TRUE),
    n_male = sum(Sex.letter == "Male", na.rm = TRUE),
    n_exacerbador_6m = sum(Exacerbations.6M == 1, na.rm = TRUE),
    n_no_exacerbador_6m = sum(Exacerbations.6M == 0, na.rm = TRUE),
    n_exacerbacion_12m_codigo_1 = sum(Exacerbations.12M == 1, na.rm = TRUE),
    n_exacerbacion_12m_codigo_2 = sum(Exacerbations.12M == 2, na.rm = TRUE),
    n_tenerife = sum(Island == "Tenerife", na.rm = TRUE),
    n_la_palma = sum(Island == "La Palma", na.rm = TRUE),
    n_antibiotico_2m = sum(Antibiotic.2mo == 1, na.rm = TRUE),
    n_chicle_30min = sum(Chewing.gum.30min == 1, na.rm = TRUE),
    n_caries = sum(Cavities == 1, na.rm = TRUE)
  )

write_csv2(
  resumen_variables_clave,
  file.path(ruta_tablas, "08_resumen_variables_clave_saliva.csv"),
  na = ""
)

# Tablas de frecuencia utiles
frecuencias_clave <- bind_rows(
  metadata_saliva_final %>% count(Sex.letter, name = "n") %>% mutate(variable = "Sex.letter", valor = as.character(Sex.letter)) %>% select(variable, valor, n),
  metadata_saliva_final %>% count(Island, name = "n") %>% mutate(variable = "Island", valor = as.character(Island)) %>% select(variable, valor, n),
  metadata_saliva_final %>% count(Exacerbations.6M, name = "n") %>% mutate(variable = "Exacerbations.6M", valor = as.character(Exacerbations.6M)) %>% select(variable, valor, n),
  metadata_saliva_final %>% count(Exacerbations.12M, name = "n") %>% mutate(variable = "Exacerbations.12M", valor = as.character(Exacerbations.12M)) %>% select(variable, valor, n),
  metadata_saliva_final %>% count(Antibiotic.2mo, name = "n") %>% mutate(variable = "Antibiotic.2mo", valor = as.character(Antibiotic.2mo)) %>% select(variable, valor, n),
  metadata_saliva_final %>% count(Chewing.gum.30min, name = "n") %>% mutate(variable = "Chewing.gum.30min", valor = as.character(Chewing.gum.30min)) %>% select(variable, valor, n),
  metadata_saliva_final %>% count(Cavities, name = "n") %>% mutate(variable = "Cavities", valor = as.character(Cavities)) %>% select(variable, valor, n),
  metadata_saliva_final %>% count(Pool, name = "n") %>% mutate(variable = "Pool", valor = as.character(Pool)) %>% select(variable, valor, n),
  metadata_saliva_final %>% count(tiene_cruce_ambiental, name = "n") %>% mutate(variable = "tiene_cruce_ambiental", valor = as.character(tiene_cruce_ambiental)) %>% select(variable, valor, n)
) %>%
  group_by(variable) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  ungroup()

write_csv2(
  frecuencias_clave,
  file.path(ruta_tablas, "09_frecuencias_variables_clave_saliva.csv"),
  na = ""
)

# Decision sobre analisis funcional inferido
estado_picrust2 <- data.frame(
  elemento = c(
    "refseq_en_phyloseq",
    "archivo_fasta_o_qza_en_proyecto",
    "decision"
  ),
  valor = c(
    !is.null(refseq(ps, errorIfNULL = FALSE)),
    FALSE,
    "No incluido en el plan analitico final: no se dispone de secuencias representativas de ASV."
  ),
  stringsAsFactors = FALSE
)

write_csv2(
  estado_picrust2,
  file.path(ruta_tablas, "10_decision_analisis_funcional.csv"),
  na = ""
)

write_xlsx(
  list(
    Resumen_objetos = resumen_objetos,
    Variables_clave = resumen_variables_clave,
    Frecuencias_clave = frecuencias_clave,
    Sin_cruce_ambiental = ids_sin_cruce_ambiental,
    Decision_funcional = estado_picrust2
  ),
  file.path(ruta_tablas, "11_resumen_preparacion_phyloseq_analitico.xlsx")
)

# 9. RESUMEN  Y FIN
cat("\n============================================\n")
cat("PREPARACION PHYLOSEQ ANALITICO FINALIZADA\n")
cat("============================================\n")

cat("\nResumen de objetos:\n")
print(resumen_objetos)

cat("\nMuestras sin cruce ambiental:\n")
print(ids_sin_cruce_ambiental)

cat("\nVariables clave:\n")
print(resumen_variables_clave)

cat("\nFrecuencias clave:\n")
print(frecuencias_clave)

cat("\nEstado PICRUSt2:\n")
print(estado_picrust2)

cat("\nCarpeta de salida:\n")
cat(ruta_salida, "\n")
