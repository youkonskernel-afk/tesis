# tesis — sRNAs no descritos vía PU learning

## Regla de ubicación (aplica siempre, sin preguntar)

**El código va a GitHub. La data va a Google Drive.** Nada de data pesada entra
al repo, ni siquiera "temporalmente".

| Va a git | Va a Drive | No va a ninguna parte |
| :-- | :-- | :-- |
| Scripts del pipeline, código del modelo | BAMs, salidas YASMA, QC | `.sra` (público, `prefetch` lo recupera) |
| `config.sh`, `environment.yml` | Matrices de features, **genomas** | Índices bowtie (derivados del genoma) |
| `organismos.tsv`, manifiestos | Modelos entrenados, checkpoints | Temporales de `fasterq-dump` |
| Vitácoras, docs, figuras vectoriales | Figuras raster pesadas | `__pycache__`, logs |

Criterio: si es texto, chico y necesario para **reproducir**, va a git. Si es
salida binaria y pesada, va a Drive. Si es derivado y se re-genera en minutos,
no se respalda.

**Los genomas sí se respaldan** (`70_genomas/`), aunque sean públicos. Es una
excepción deliberada al criterio anterior: un ensamblado se puede retirar,
renumerar o reemplazar por una versión nueva, y entonces el alineamiento deja
de ser reproducible. Guardar el FASTA exacto con su `sha256` es la única forma
de poder decir, dentro de dos años, contra qué se alineó. Los índices bowtie no
se respaldan: se reconstruyen del FASTA en minutos.

Drive es `Mi unidad/tesis/` en `seb.ugazm@gmail.com`. Los IDs de carpeta están
en `data/DRIVE.md` — esa es la única fuente de verdad de la ruta.

## Qué es este proyecto

Objetivo: **descubrir sRNAs no descritos** con un clasificador **PU learning**
(positive-unlabeled). Los sRNAs anotados en miRBase/Rfam son el conjunto
positivo; el resto de loci anotados por YASMA son *unlabeled*, no negativos —
de ahí PU y no clasificación binaria clásica. Asumir que lo no anotado es
negativo sesga el modelo justo contra lo que buscamos.

El pipeline bioinformático es el **upstream** que produce los loci candidatos,
no el aporte de la tesis:

1. Descarga ~425 corridas sRNA-seq de 9 organismos (~7.4 G spots)
2. `fastp` (15-50 nt) → `bowtie1` estilo ShortStack3 → `samtools`
3. Anotación de loci con YASMA v1.1.1
4. Features por locus → PU learning → candidatos priorizados

## Dataset

**9 organismos, 3 por reino, cada uno con BioProject primario y duplicado.**
`data/organismos.tsv` es la especificación (entrada de `gen_manifest.sh`).
El manifiesto resuelto de esta ronda **todavía no existe**: hay que generarlo
con `./gen_manifest.sh` contra la ENA. `data/srr_manifest_r1.tsv` es el de la
ronda anterior, conservado como referencia.

El primario es el set de **descubrimiento**; el duplicado es un experimento
independiente para **validar** candidatos. Un locus que el modelo prioriza en el
primario y que reaparece en el duplicado es mucho más defendible que uno que
solo existe en un experimento. Esa separación hay que respetarla: **el duplicado
no entra al entrenamiento**, o la validación deja de ser independiente.

| org | especie | reino | primario | duplicado |
| :-- | :-- | :-- | :-- | :-- |
| rhirr | *Rhizophagus irregularis* | Fungi (Glomeromycota) | PRJEB29180 | PRJNA722321 |
| sclsc | *Sclerotinia sclerotiorum* | Fungi (Ascomycota) | PRJNA477286 | PRJNA985401 |
| cloro | *Clonostachys rosea* | Fungi (Ascomycota) | PRJEB43636 | PRJEB51338 |
| phypa | *Physcomitrium patens* | Plantae (briofita) | PRJNA222997 | PRJNA277372 |
| prupe | *Prunus persica* | Plantae (rosácea) | PRJNA929031 | PRJNA780811 |
| maldo | *Malus domestica* | Plantae (rosácea) | PRJNA681626 | PRJNA784097 |
| gadmo | *Gadus morhua* | Animalia (pez) | PRJNA284846 | PRJNA328800 |
| galga | *Gallus gallus* | Animalia (ave) | PRJEB12164 | PRJNA694114 |
| maggi | *Magallana gigas* | Animalia (molusco) | PRJNA154615 + PRJNA232734 | PRJNA1254880 |

Escala declarada: ~195 corridas / ~3.0 G spots en los primarios, ~230 / ~4.4 G
en los duplicados. **~425 corridas y ~7.4 G spots en total, contra 169 y 2.18 G
del set anterior** — del orden de 3.4× más. Ver "Consecuencias del cambio".

`gadmo`, `galga` y los dos primarios de `maggi` están **declarados por
MirGeneDB**, o sea que su set de miRNAs positivos es curado. Para PU learning
eso importa más que la profundidad: la calidad del conjunto positivo define el
techo del modelo.

### Consecuencias del cambio de dataset (pendientes)

- **Se caen `arath`, `danre` y `nemve`.** Sus `.sra`, BAMs e índices son ahora
  peso muerto. No borrar sin confirmar: el caché es re-descargable, pero los
  BAMs de `danre` costaron horas de alineamiento.
- **Faltan 6 ensamblados** (`cloro`, `prupe`, `maldo`, `gadmo`, `galga`,
  `maggi`). Están en `data/genomas.tsv` con estado `candidato`: son sugerencias
  **sin verificar**, escritas de memoria y no comprobadas contra ninguna base.
  `fetch_genomes.sh` se niega a bajarlas hasta que una persona las confirme con
  `./scripts/fetch_genomes.sh resolve` y cambie el estado a `verificado`. El
  paso manual es a propósito: un ensamblado equivocado no falla ruidosamente,
  alinea peor y contamina la anotación.
- **`Magallana gigas` = `Crassostrea gigas`.** El género se renombró; Ensembl
  Metazoa y buena parte de las bases todavía usan *Crassostrea*. Buscar el
  genoma por el nombre nuevo no va a encontrarlo.
- **`maldo` duplicado (PRJNA784097) trae reads de 151 nt.** Son librerías sin
  recortar: hay que revisar el pre-trim antes de que `fastp` las vea, o la
  ventana 15-50 nt las descarta enteras.
- **Dos proyectos no son miRNA-Seq**: `sclsc` duplicado es RNA-Seq y `cloro`
  primario es ncRNA-Seq. Igual que `phypa`, hay que declararlo en métodos.
- **Re-estimar tiempo y espacio.** Con 3.4× más reads, las 8-12 h de
  alineamiento pasan al orden de 30-40 h, y los ~100 GB de BAMs al orden de
  340 GB. Entra sin problema en 1.6 TB, pero el cronograma cambia.

## Trampas conocidas — no re-introducir

- **`fastp --disable_length_filtering`**: desactiva el filtro de longitud, no el
  de Ns (eso es `--n_base_limit`). Con ese flag entraba mRNA fragmentado de
  50-150 nt y distorsionaba las *size class*. Ya se eliminó; no volver a ponerlo.
- **Ventana 15-50 nt**, no 15-40: con 40 se truncaban tRFs (30-40 nt) y los
  piRNAs (24-32 nt) quedaban sin margen.
- **`-m 50` en bowtie se mantiene.** Medido en `danre` (ya fuera del set):
  descartaba el 76% de los reads, pero el 98.8% de lo que se recuperaría al
  subirlo eran fragmentos de tRNA (tRF-5) de 34 nt exactos, en cientos de copias
  génicas. Subirlo inunda la anotación con una sola especie repetitiva.
  El hallazgo es del organismo viejo, pero el fenómeno no: los tRF multimapean
  en cualquier genoma, y `gadmo`, `galga` y `maggi` son animales igual que
  `danre`. Si alguno muestra una fracción alineada anormalmente baja, revisar
  la distribución de longitudes antes de tocar `-m` — y declararlo en métodos,
  porque un revisor lo va a preguntar.
- **YASMA v1.1.1 escribe en `annotations/<nombre>/loci.gff3`**, no en
  `annotation/`. Un chequeo contra la ruta vieja hacía abortar el pipeline tras
  el primer organismo.
- **`fasterq-dump` sin `-t`** crea temporales en el CWD; llegaron a 109 GB.
- **Al matar el pipeline**: matar también `bowtie-align-s`, `fastp`,
  `fasterq-dump` y `samtools`, o quedan huérfanos escribiendo el mismo BAM.

## Entorno

micromamba, entorno `srna2` (se llama así porque el primer intento usó Python
3.11 y YASMA exige >= 3.12). Pin en `environment.yml` + `environment_pip.txt`.
`./check_env.sh` verifica que cada herramienta **arranque** de verdad, incluido
`import RNA` de ViennaRNA — el fallo clásico.
