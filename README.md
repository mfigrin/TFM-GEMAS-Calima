# GEMAS–Calima | Análisis ambiental y microbiota salival
Código asociado al TFM sobre intrusiones de polvo sahariano y microbiota salival en pacientes con asma de Canarias.

## Carpetas
 #`analisis_final/`  Siete etapas numeradas del análisis final; `modulos/` reúne las funciones auxiliares. Se usa desde la demostración. |
 #`demo/` Tablas y árbol **sintéticos**, para generar `phyloseq`. No se incluyen datos de participantes reales. 
 # `codigo_original_tfm/`  Once scripts **No ejecutarlos como un pipeline**
 #`ambiente/`  Importación de la Red Canaria y registros MITECO, sin muestras de participantes. 
 #`data/environmental/` Fuentes ambientales originales 2016–2020 
 #`herramientas/` y `tests/`  Dependencias de ejecución y control técnico de la demostración. 


## Ejecutar en RStudio

# Comprobar paquetes y sintaxis, sin instalar nada automáticamente.
source("herramientas/comprobar.R")

# Comprobar reconstrucción de los datos ambientales públicos.
source("ambiente/ejecutar.R")

# Demostración completa SOLO con datos sintéticos.
source("demo/ejecutar.R")
```

La demostración genera `demo/trabajo/` y registros por etapa. 

El código ejecuta: preparación de objetos, diversidad beta/UniFrac, taxonomía, modelos de intrusión MITECO, DESeq2, PERMANOVA principal y complementos de PM10/tablas. La PERMANOVA principal usa **9.999 permutaciones**, semilla **20260831**; UniFrac no ponderado requiere rarefacción a **10.000 lecturas**, semilla **20260822**. 
Esta demostración comprueba el funcionamiento técnico; **no reproduce los resultados clínicos publicados en el TFM**. Para hacerlo serían necesarios los datos originales de GEMAS y los metadatos clínicos obtenidos con autorización. 

Los once scripts conservan internamente las fórmulas, revisiones y comentarios de la investigación original. Se adaptaron previamente rutas locales y dos ejemplos de identificadores de muestras para evitar difusión de información individual. Esta versión **solo les cambia el nombre de archivo**, de la siguiente forma:

## Fuentes, privacidad y difusión

- **Fuentes MITECO:** [Fuentes naturales e informes históricos](https://www.miteco.gob.es/es/calidad-y-evaluacion-ambiental/temas/atmosfera-y-calidad-del-aire/evaluacion-y-datos-de-calidad-del-aire/fuentes-naturales/anuales.html).
- **Red Canaria:** [Red de Control y Vigilancia de la Calidad del Aire](https://www3.gobiernodecanarias.org/medioambiente/calidaddelaire/inicio.do).
- **Lecturas microbiológicas originales:** BioProject [PRJNA878647](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA878647), custodiadas por SRA

