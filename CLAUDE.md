# tesis — sRNAs no descritos vía PU learning

## Regla de ubicación (aplica siempre, sin preguntar)

**El código va a GitHub. La data va a Google Drive.** Nada de data pesada entra
al repo, ni siquiera "temporalmente".

| Va a git | Va a Drive | No va a ninguna parte |
| :-- | :-- | :-- |
| Scripts del pipeline, código del modelo | BAMs, salidas YASMA, QC | `.sra` (público, `prefetch` lo recupera) |
| `config.sh`, `environment.yml` | Matrices de features | Genomas e índices (re-descargables) |
| `srr_manifest.tsv` | Modelos entrenados, checkpoints | Temporales de `fasterq-dump` |
| Vitácoras, docs, figuras vectoriales | Figuras raster pesadas | `__pycache__`, logs |

Criterio: si es texto, chico y necesario para **reproducir**, va a git. Si es
salida binaria y pesada, va a Drive. Si se puede re-descargar de una base
pública, no se respalda.

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

1. Descarga 169 corridas sRNA-seq de 6 organismos (2.18 G reads)
2. `fastp` (15-50 nt) → `bowtie1` estilo ShortStack3 → `samtools`
3. Anotación de loci con YASMA v1.1.1
4. Features por locus → PU learning → candidatos priorizados

## Dataset

169 corridas, 6 organismos, 4 reinos. `data/srr_manifest.tsv` es la fuente de
verdad (regenerable con `gen_manifest.sh` desde la ENA Portal API).

| org | especie | reino | corridas |
| :-- | :-- | :-- | --: |
| rhirr | *Rhizophagus irregularis* | Fungi (Glomeromycota) | 9 |
| sclsc | *Sclerotinia sclerotiorum* | Fungi (Ascomycota) | 18 |
| arath | *Arabidopsis thaliana* | Plantae (rósida) | 3 |
| phypa | *Physcomitrium patens* | Plantae (briofita) | 30 |
| danre | *Danio rerio* | Animalia (vertebrado) | 42 |
| nemve | *Nematostella vectensis* | Animalia (cnidario) | 67 |

## Trampas conocidas — no re-introducir

- **`fastp --disable_length_filtering`**: desactiva el filtro de longitud, no el
  de Ns (eso es `--n_base_limit`). Con ese flag entraba mRNA fragmentado de
  50-150 nt y distorsionaba las *size class*. Ya se eliminó; no volver a ponerlo.
- **Ventana 15-50 nt**, no 15-40: con 40 se truncaban tRFs (30-40 nt) y los
  piRNAs (24-32 nt) quedaban sin margen.
- **`-m 50` en bowtie se mantiene.** En danre descarta el 76% de los reads, pero
  el 98.8% de lo que se recuperaría son fragmentos de tRNA (tRF-5) de 34 nt
  exactos, en cientos de copias génicas. Subirlo inundaría la anotación con una
  sola especie repetitiva. Esto hay que declararlo en métodos: un revisor va a
  preguntar por qué la fracción alineada de pez cebra es tan baja.
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
