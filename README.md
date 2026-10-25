# GEMAS–Calima | Intrusiones saharianas y microbiota salival en asma

Repositorio de código del Trabajo Fin de Máster (TFM) sobre la asociación entre intrusiones de polvo sahariano y microbiota salival en personas con asma de Tenerife y La Palma. Incluye un flujo de análisis en R, datos ambientales de fuentes públicas, una demostración **completamente sintética** y el código histórico para consulta docente.

## Organización del proyecto

| Carpeta | Contenido |
| --- | --- |
| `01_analisis/` | Siete etapas numeradas del análisis final; `modulos/` y `utilidades/` contienen funciones auxiliares. |
| `02_ejemplo_sintetico/` | Tablas, árbol y metadatos inventados; preparación del objeto `phyloseq` y ejecución de la demostración. |
| `03_codigo_original/` | Once scripts históricos del trabajo, conservados con fines docentes. **No deben ejecutarse secuencialmente como un pipeline**. |
| `04_datos_ambientales/` | Datos ambientales originales y procesados; `procesar/` contiene los lectores de MITECO y Red Canaria. |

Los archivos `TFM-GEMAS-Calima.Rproj`, `CITATION.cff` y `LICENSE` están en la raíz. Los resultados temporales y objetos generados en cada ejecución no se incluyen en el control de versiones.

## Ejecución en RStudio

Requiere R y los paquetes enumerados en `01_analisis/utilidades/comprobar.R`. El flujo se probó en Windows 11 con R 4.5.2 y las versiones de los paquetes documentadas en el TFM. Abrir `TFM-GEMAS-Calima.Rproj` y ejecutar desde la raíz:

```r
# Comprobar dependencias y sintaxis.
source("01_analisis/utilidades/comprobar.R")

# Reconstruir y comprobar las series ambientales públicas.
source("04_datos_ambientales/procesar/ejecutar.R")

# Ejecutar las siete etapas con datos enteramente sintéticos.
source("02_ejemplo_sintetico/ejecutar.R")
```

La demostración construye su objeto `phyloseq` exclusivamente a partir de los archivos ficticios de `02_ejemplo_sintetico/datos/` y registra las salidas en `02_ejemplo_sintetico/trabajo/`. Para repetir una ejecución, conservar una copia de las salidas que se quieran guardar y renombrar o eliminar **solo ese directorio de trabajo sintético**. No es necesario disponer de datos clínicos reales.

Las etapas comprenden preparación de objetos, distancias ecológicas, composición taxonómica, modelos de exposición MITECO, DESeq2, PERMANOVA principal y análisis complementarios de PM10. La PERMANOVA principal utiliza 9.999 permutaciones (semilla 20260831); los contrastes beta secundarios conservan 999 permutaciones. UniFrac no ponderado utiliza rarefacción a 10.000 lecturas por muestra (semilla 20260822).

**Alcance de la reproducibilidad:** la demostración permite comprobar técnicamente la ejecución, pero **no reproduce los resultados ni las estimaciones clínicas del TFM**. Para ello se necesita acceso autorizado al objeto microbiológico procesado y a los metadatos individuales de GEMAS. No se distribuyen esos datos ni el objeto `phyloseq` real.

## Fuentes de datos y difusión

- **MITECO:** [Fuentes naturales e informes históricos](https://www.miteco.gob.es/es/calidad-y-evaluacion-ambiental/temas/atmosfera-y-calidad-del-aire/evaluacion-y-datos-de-calidad-del-aire/fuentes-naturales/anuales.html).
- **Red Canaria de Calidad del Aire:** [Portal oficial](https://www3.gobiernodecanarias.org/medioambiente/calidaddelaire/inicio.do).
- **Lecturas microbiológicas originales de GEMAS:** [BioProject PRJNA878647](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA878647). No se reprocesan FASTQ en el análisis secundario del TFM.

La licencia MIT del repositorio se aplica al código propio. Los datos y documentos de terceros mantienen sus condiciones de uso y atribución de origen.

## Citar este repositorio

Usar la referencia indicada en `CITATION.cff`. Una vez depositada una versión en Zenodo, el DOI específico de esa versión se añadirá aquí.
