# Comprobar e instalar phyloseq 

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install(
    "phyloseq",
    ask = FALSE,
    update = FALSE
  )
}

library(phyloseq)
# Localizar el objeto phyloseq actualizado
list.files("data_raw", full.names = TRUE)
# Cargar el objeto actualizado
archivo_microbiota <- "data_raw/data_phyloseq_salivas_GEMAS_Canarias.rds"

ps_gemas <- readRDS(archivo_microbiota)

class(ps_gemas)
print(ps_gemas)
# Extraer los metadatos de todas las muestras
metadata_gemas <- as.data.frame(sample_data(ps_gemas))

# Conservar el identificador original de cada muestra
metadata_gemas$id_muestra_phyloseq <- rownames(metadata_gemas)

# Ver los nombres de todas las variables disponibles
names(metadata_gemas)
# Revisar únicamente las muestras clínicas de saliva
metadata_saliva <- metadata_gemas[
  metadata_gemas$sample.type == "Saliva",
  ,
  drop = FALSE
]

cat("Número de muestras de saliva:", nrow(metadata_saliva), "\n\n")

# Revisar formato y ejemplos de las nuevas variables
str(metadata_saliva[, c(
  "id_muestra_phyloseq",
  "Exacerbations.6M",
  "Exacerbations.12M",
  "Date",
  "Island"
)])

cat("\nDistribución de isla:\n")
print(table(metadata_saliva$Island, useNA = "ifany"))

cat("\nExacerbaciones en 6 meses:\n")
print(table(metadata_saliva$Exacerbations.6M, useNA = "ifany"))

cat("\nExacerbaciones en 12 meses:\n")
print(table(metadata_saliva$Exacerbations.12M, useNA = "ifany"))

cat("\nRelación entre exacerbaciones a 6 y 12 meses:\n")
print(table(
  Exacerbaciones_6m = metadata_saliva$Exacerbations.6M,
  Exacerbaciones_12m = metadata_saliva$Exacerbations.12M,
  useNA = "ifany"
))

cat("\nPrimeras 15 muestras clínicas:\n")
print(
  head(
    metadata_saliva[, c(
      "id_muestra_phyloseq",
      "Age",
      "Sex.letter",
      "Exacerbations.6M",
      "Exacerbations.12M",
      "Date",
      "Island"
    )],
    15
  )
)
############################################################
# PREPARACIÓN DE LA COHORTE DE MICROBIOTA GEMAS CANARIAS
#
# Objetivo:
# - Cargar el objeto phyloseq actualizado.
# - Mantener solo las muestras clínicas de saliva.
# - Crear una tabla de metadatos limpia para el cruce ambiental.
# - Confirmar y recodificar exacerbaciones a 6 meses:
#       0 = no exacerbador
#       1 = exacerbador
# - Marcar muestras sin fecha para excluirlas únicamente del
#   cruce ambiental.


# 0. PAQUETES Y COMPROBACIONES INICIALES

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  stop(
    "No está instalado el paquete phyloseq. ",
    "Ejecuta primero: BiocManager::install('phyloseq')"
  )
}

library(phyloseq)

# Crear carpetas 
if (!dir.exists("outputs")) {
  dir.create("outputs")
}

if (!dir.exists("data_raw")) {
  stop(
    "No existe la carpeta 'data_raw'. ",
    "Comprueba que estás trabajando dentro del proyecto TFM_GEMAS."
  )
}

# 1. CARGAR EL OBJETO PHYLOSEQ 

archivo_microbiota <- "data_raw/data_phyloseq_salivas_GEMAS_Canarias.rds"

if (!file.exists(archivo_microbiota)) {
  stop(
    paste0(
      "No se encuentra el archivo:\n",
      archivo_microbiota,
      "\n\nComprueba el nombre del archivo dentro de data_raw."
    )
  )
}

ps_gemas <- readRDS(archivo_microbiota)

if (!inherits(ps_gemas, "phyloseq")) {
  stop("El archivo cargado no contiene un objeto phyloseq.")
}

cat("\n============================================\n")
cat("OBJETO PHYLOSEQ CARGADO CORRECTAMENTE\n")
cat("============================================\n")

print(ps_gemas)


# 2. EXTRAER LOS METADATOS COMO DATA.FRAME BASE

# sample_data() conserva una clase especial de phyloseq.
# Esta conversión evita errores posteriores con dplyr/tidyverse.

metadata_completo <- data.frame(
  setNames(
    sample_data(ps_gemas)@.Data,
    names(sample_data(ps_gemas))
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(metadata_completo) <- sample_names(ps_gemas)

# Comprobar que existen las columnas necesarias
columnas_requeridas <- c(
  "Age",
  "Sex.number",
  "Sex.letter",
  "Exacerbations.6M",
  "Exacerbations.12M",
  "Antibiotic.2mo",
  "Pool",
  "Cavities",
  "Chewing.gum.30min",
  "sample.type",
  "Date",
  "Island"
)

columnas_ausentes <- setdiff(
  columnas_requeridas,
  names(metadata_completo)
)

if (length(columnas_ausentes) > 0) {
  stop(
    paste0(
      "Faltan las siguientes variables en el objeto actualizado:\n",
      paste(columnas_ausentes, collapse = ", ")
    )
  )
}

cat("\n============================================\n")
cat("VARIABLES DISPONIBLES EN EL OBJETO\n")
cat("============================================\n")

print(names(metadata_completo))

# 3. SELECCIONAR MUESTRAS CLÍNICAS DE SALIVA

metadata_saliva <- metadata_completo[
  metadata_completo$sample.type == "Saliva",
  ,
  drop = FALSE
]

cat("\n============================================\n")
cat("TIPOS DE MUESTRA\n")
cat("============================================\n")

print(table(metadata_completo$sample.type, useNA = "ifany"))

cat("\nNúmero de muestras clínicas de saliva: ",
    nrow(metadata_saliva),
    "\n",
    sep = "")


# 4. FUNCIONES DE LIMPIEZA

convertir_blanco_a_na <- function(x) {
  
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  
  return(x)
}

normalizar_isla <- function(x) {
  
  x <- trimws(as.character(x))
  
  resultado <- ifelse(
    x == "Tenerife",
    "Tenerife",
    ifelse(
      x %in% c("La_Palma", "La Palma"),
      "La Palma",
      NA_character_
    )
  )
  
  return(resultado)
}

# 5. CREAR LA COHORTE PARA EL CRUCE AMBIENTAL

cohorte_microbiota <- data.frame(
  
  id_muestra = rownames(metadata_saliva),
  
  fecha_muestra = as.Date(metadata_saliva$Date),
  
  isla = normalizar_isla(metadata_saliva$Island),
  
  edad = suppressWarnings(
    as.numeric(as.character(metadata_saliva$Age))
  ),
  
  sexo = convertir_blanco_a_na(
    metadata_saliva$Sex.letter
  ),
  
  sexo_codigo = convertir_blanco_a_na(
    metadata_saliva$Sex.number
  ),
  
  # Codificación confirmada por equipo GEMAS:
  # 0 = no exacerbador en los 6 meses previos
  # 1 = exacerbador en los 6 meses previos
  
  exacerbaciones_6m_codigo = convertir_blanco_a_na(
    metadata_saliva$Exacerbations.6M
  ),
  
  exacerbacion_6m = ifelse(
    convertir_blanco_a_na(metadata_saliva$Exacerbations.6M) == "1",
    1,
    ifelse(
      convertir_blanco_a_na(metadata_saliva$Exacerbations.6M) == "0",
      0,
      NA
    )
  ),
  
  exacerbaciones_12m_codigo = convertir_blanco_a_na(
    metadata_saliva$Exacerbations.12M
  ),
  
  antibiotico_2m = convertir_blanco_a_na(
    metadata_saliva$Antibiotic.2mo
  ),
  
  pool = convertir_blanco_a_na(
    metadata_saliva$Pool
  ),
  
  caries = convertir_blanco_a_na(
    metadata_saliva$Cavities
  ),
  
  chicle_30min = convertir_blanco_a_na(
    metadata_saliva$Chewing.gum.30min
  ),
  
  stringsAsFactors = FALSE
)

# Variables de elegibilidad y trazabilidad

cohorte_microbiota$adulto <- ifelse(
  is.na(cohorte_microbiota$edad),
  NA,
  cohorte_microbiota$edad >= 18
)

cohorte_microbiota$fecha_valida <- !is.na(
  cohorte_microbiota$fecha_muestra
)

cohorte_microbiota$isla_valida <- !is.na(
  cohorte_microbiota$isla
)

# Puede entrar al cruce ambiental: requiere fecha e isla, pero NO excluye por edad.
cohorte_microbiota$incluible_cruce_ambiental <- 
  cohorte_microbiota$fecha_valida &
  cohorte_microbiota$isla_valida

# Motivo de no inclusión, para trazabilidad.
cohorte_microbiota$motivo_no_cruce_ambiental <- ifelse(
  cohorte_microbiota$incluible_cruce_ambiental,
  NA_character_,
  ifelse(
    !cohorte_microbiota$fecha_valida &
      !cohorte_microbiota$isla_valida,
    "Fecha e isla ausentes",
    ifelse(
      !cohorte_microbiota$fecha_valida,
      "Fecha de muestra ausente",
      "Isla de muestra ausente"
    )
  )
)

# Ordenar para facilitar auditoría manual
cohorte_microbiota <- cohorte_microbiota[
  order(
    cohorte_microbiota$isla,
    cohorte_microbiota$fecha_muestra,
    cohorte_microbiota$id_muestra,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

# 6. AUDITORÍA DE LA COHORTE

cat("\n============================================\n")
cat("RESUMEN DE LA COHORTE MICROBIOTA\n")
cat("============================================\n")

resumen_cohorte <- data.frame(
  
  total_muestras_saliva = nrow(cohorte_microbiota),
  
  ids_unicos = length(
    unique(cohorte_microbiota$id_muestra)
  ),
  
  fecha_minima = min(
    cohorte_microbiota$fecha_muestra,
    na.rm = TRUE
  ),
  
  fecha_maxima = max(
    cohorte_microbiota$fecha_muestra,
    na.rm = TRUE
  ),
  
  fecha_ausente = sum(
    is.na(cohorte_microbiota$fecha_muestra)
  ),
  
  isla_ausente = sum(
    is.na(cohorte_microbiota$isla)
  ),
  
  muestras_incluibles_cruce_ambiental = sum(
    cohorte_microbiota$incluible_cruce_ambiental,
    na.rm = TRUE
  ),
  
  muestras_no_incluibles_cruce_ambiental = sum(
    !cohorte_microbiota$incluible_cruce_ambiental,
    na.rm = TRUE
  ),
  
  adultos = sum(
    cohorte_microbiota$adulto,
    na.rm = TRUE
  ),
  
  menores_18 = sum(
    !cohorte_microbiota$adulto,
    na.rm = TRUE
  ),
  
  edad_ausente = sum(
    is.na(cohorte_microbiota$edad)
  ),
  
  exacerbadores_6m = sum(
    cohorte_microbiota$exacerbacion_6m == 1,
    na.rm = TRUE
  ),
  
  no_exacerbadores_6m = sum(
    cohorte_microbiota$exacerbacion_6m == 0,
    na.rm = TRUE
  ),
  
  exacerbacion_6m_ausente = sum(
    is.na(cohorte_microbiota$exacerbacion_6m)
  )
)

print(resumen_cohorte)


cat("\n============================================\n")
cat("MUESTRAS POR ISLA\n")
cat("============================================\n")

print(
  table(
    cohorte_microbiota$isla,
    useNA = "ifany"
  )
)


cat("\n============================================\n")
cat("EXACERBACIONES EN LOS 6 MESES PREVIOS\n")
cat("============================================\n")

print(
  table(
    cohorte_microbiota$exacerbacion_6m,
    useNA = "ifany"
  )
)


cat("\n============================================\n")
cat("ADULTOS Y MENORES DE 18 AÑOS\n")
cat("============================================\n")

print(
  table(
    cohorte_microbiota$adulto,
    useNA = "ifany"
  )
)


cat("\n============================================\n")
cat("MUESTRAS SIN FECHA O ISLA PARA CRUCE AMBIENTAL\n")
cat("============================================\n")

muestras_revision <- cohorte_microbiota[
  !cohorte_microbiota$incluible_cruce_ambiental,
  c(
    "id_muestra",
    "fecha_muestra",
    "isla",
    "edad",
    "sexo",
    "exacerbacion_6m",
    "motivo_no_cruce_ambiental"
  ),
  drop = FALSE
]

print(muestras_revision)


cat("\n============================================\n")
cat("MENORES DE 18 AÑOS\n")
cat("============================================\n")

muestras_menores <- cohorte_microbiota[
  !is.na(cohorte_microbiota$edad) &
    cohorte_microbiota$edad < 18,
  c(
    "id_muestra",
    "edad",
    "sexo",
    "fecha_muestra",
    "isla",
    "exacerbacion_6m"
  ),
  drop = FALSE
]

print(muestras_menores)


cat("\n============================================\n")
cat("VALORES AUSENTES POR VARIABLE\n")
cat("============================================\n")

valores_ausentes <- data.frame(
  variable = names(cohorte_microbiota),
  n_ausentes = colSums(is.na(cohorte_microbiota)),
  porcentaje_ausente = round(
    100 * colSums(is.na(cohorte_microbiota)) /
      nrow(cohorte_microbiota),
    2
  ),
  row.names = NULL
)

print(valores_ausentes)

# 7. EXPORTAR RESULTADOS

write.csv(
  cohorte_microbiota,
  "outputs/01_cohorte_microbiota_para_exposicion_ambiental.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  resumen_cohorte,
  "outputs/02_resumen_cohorte_microbiota.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  muestras_revision,
  "outputs/03_muestras_no_incluibles_cruce_ambiental.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  muestras_menores,
  "outputs/04_muestras_menores_18.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  valores_ausentes,
  "outputs/05_valores_ausentes_cohorte.csv",
  row.names = FALSE,
  na = ""
)

# Se exporta una versión exclusivamente apta para el cruce ambiental: fecha e isla obligatorias.
# MUESTRA_SIN_FECHA no aparecerá aquí mientras siga sin fecha.

cohorte_cruce_ambiental <- cohorte_microbiota[
  cohorte_microbiota$incluible_cruce_ambiental,
  ,
  drop = FALSE
]

write.csv(
  cohorte_cruce_ambiental,
  "outputs/06_cohorte_elegible_cruce_ambiental.csv",
  row.names = FALSE,
  na = ""
)

# 8. MENSAJE FINAL

cat("\n============================================\n")
cat("PREPARACIÓN DE COHORTE FINALIZADA\n")
cat("============================================\n")

cat(
  "\nArchivos creados en la carpeta outputs:\n",
  "- 01_cohorte_microbiota_para_exposicion_ambiental.csv\n",
  "- 02_resumen_cohorte_microbiota.csv\n",
  "- 03_muestras_no_incluibles_cruce_ambiental.csv\n",
  "- 04_muestras_menores_18.csv\n",
  "- 05_valores_ausentes_cohorte.csv\n",
  "- 06_cohorte_elegible_cruce_ambiental.csv\n",
  sep = ""
)

cat(
  "\n\nIMPORTANTE: estos archivos contienen IDs pseudonimizados, ",
  "fechas e información clínica. Deben mantenerse fuera de GitHub.\n",
  sep = ""
)