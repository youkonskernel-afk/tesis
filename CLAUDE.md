# tesis — sRNAs no descritos vía PU learning

## Regla de ubicación (aplica siempre, sin preguntar)

**El código va a GitHub. La data va a Google Drive.** Nada de data pesada entra
al repo, ni siquiera "temporalmente".

| Va a git | Va a Drive | No va a ninguna parte |
| :-- | :-- | :-- |
| Scripts del pipeline, código del modelo | BAMs, salidas YASMA, QC, **`.sra`** | Índices bowtie (derivados) |
| `config.sh`, `environment.yml` | Matrices de features, **genomas** | Temporales de `fasterq-dump` |
| `organismos.tsv`, manifiestos, checksums | Modelos entrenados, checkpoints | `__pycache__`, logs |
| Vitácoras, docs, notebooks | Figuras raster pesadas | |

Criterio: si es texto, chico y necesario para **reproducir**, va a git. Si es
salida binaria y pesada, va a Drive. Si es derivado y se re-genera en minutos,
no se respalda.

**Los genomas y los `.sra` sí se respaldan** (`70_genomas/`, `80_sra/`), aunque
sean públicos. Es una excepción deliberada al criterio anterior: un ensamblado
se puede retirar,
renumerar o reemplazar por una versión nueva, y entonces el alineamiento deja
de ser reproducible. Guardar el FASTA exacto con su `sha256` es la única forma
de poder decir, dentro de dos años, contra qué se alineó. Los `.sra` se guardan
por otro motivo: **Colab los baja y las sesiones de Colab se mueren**, así que si
no persistieran se perdería la descarga; además desacopla el alineamiento local
de la descarga. Los índices bowtie no se respaldan: se reconstruyen en minutos.

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
`data/organismos.tsv` es la especificación: 19 accessions en 18 slots (el
primario de `maggi` son dos BioProjects combinados). El manifiesto resuelto de
esta ronda **todavía no existe**: se genera con
`./scripts/fetch_runs.sh manifest` contra la ENA, y después
`./scripts/fetch_runs.sh prefetch` baja los `.sra`. En la práctica esto corre
desde Colab (`notebooks/10_descarga_runs.ipynb`), no en la máquina local.
`data/srr_manifest_r1.tsv` es el de la ronda anterior, como referencia.

`scripts/fetch_runs.sh` **reemplaza al `gen_manifest.sh` de R1**, que asumía un
proyecto por organismo. Filtra a datos de RNA con `library_source =
TRANSCRIPTOMIC` —el filtro duro que excluye las corridas GENOMIC que varios
BioProjects mezclan— y exige `SINGLE` para RNA-Seq, porque el PAIRED de un
proyecto de RNA-Seq no es sRNA-seq.

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

### Entrenamiento vs aplicación

**El modelo se entrena solo en `gadmo`, `galga` y `maggi`** — los tres con
positivos curados por MirGeneDB. Los otros seis son conjunto de **aplicación**:
se predice sobre ellos, no se entrena. La columna `set_modelo` de
`data/organismos.tsv` lo hace explícito. En PU learning un falso positivo es
peor que un *unlabeled*: el método asume la clase positiva limpia, y las
entradas dudosas de miRBase invertirían ese supuesto sin vuelta atrás.

**El costo de esto es que los tres de entrenamiento son los tres animales**, y
los seis de aplicación son hongos y plantas. Es transferencia entre reinos y es
la parte más frágil del diseño: los precursores de plantas son más largos y
heterogéneos, las clases de tamaño vegetales son 21 y 24 nt contra el pico
animal de ~22, y en hongos predominan milRNA y siRNA Dicer-dependiente sobre el
miRNA canónico. Un modelo que aprendió la firma animal puede no encontrar nada
en plantas por buscar la forma equivocada, y eso se confunde fácil con "no hay
nada que encontrar".

Dos cosas sostienen el diseño, y están en `docs/positivos.md`:

1. **miRBase como evaluación, nunca como entrenamiento.** En los seis de
   aplicación, cuántos miRNAs ya descritos recupera el modelo es la medida
   directa de si la transferencia funciona. No contamina nada.
2. **Validación cruzada dejando un organismo afuera entre los tres curados.**
   `maggi` es molusco y `gadmo`/`galga` vertebrados. Si el modelo no transfiere
   de pez a molusco, no va a transferir de pez a musgo — y conviene saberlo
   antes de correr los nueve organismos. Implementado en `scripts/loo_cv.py`
   (`--self-test` verifica la lógica sin datos).

**Elkan-Noto no cambia el ranking** respecto de tratar los no etiquetados como
negativos: divide por una constante. Sirve para calibrar, no para reordenar.
Para que el orden cambie hace falta bagging PU o nnPU. Verificado en el
self-test; declararlo en métodos.

**El prior de clase π se estima por organismo.** El π de animales no vale en
plantas ni hongos; si no se puede estimar bien, reportar ranking en vez de
probabilidad calibrada.

### Consecuencias del cambio de dataset (pendientes)

- **Se caen `arath`, `danre` y `nemve`.** Sus `.sra`, BAMs e índices son ahora
  peso muerto. No borrar sin confirmar: el caché es re-descargable, pero los
  BAMs de `danre` costaron horas de alineamiento.
- **Faltan 6 ensamblados.** Están en `data/genomas.tsv` como `candidato`, con
  una columna `confianza` que dice cuánto pesa cada propuesta. Cinco tienen
  candidato concreto (`prupe`, `maldo`, `gadmo`, `galga`, `maggi`) y `cloro`
  ninguno: depende de la cepa, y `PRJEB43636` son mutantes Dicer-like del grupo
  de Karlsson (SLU), casi seguro IK726. Ninguno está comprobado contra NCBI.
  `fetch_genomes.sh` se niega a bajarlos hasta que una persona corra
  `./scripts/fetch_genomes.sh resolve` y ponga `verificado`. Alternativa que
  evita el disco local: `notebooks/descarga_genomas.ipynb` en Colab verifica y
  baja directo a `70_genomas/` — ver `docs/colab.md`. El paso manual es
  a propósito: un ensamblado equivocado no falla ruidosamente, alinea peor y
  contamina la anotación. Ya no hace falta que coincida con el ensamblado de
  MirGeneDB — con el etiquetado por secuencia se elige por contigüidad y
  completitud.
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
- **Los positivos se etiquetan por secuencia, no por coordenada.** Decidido;
  ver `docs/positivos.md`. Un locus es positivo si su RNA mayoritario coincide
  con un sRNA descrito, no si su intervalo se solapa con una anotación. Esto
  desacopla el proyecto del ensamblado: si alineáramos contra un genoma y la
  base publicara sobre otro, los miRNAs conocidos no caerían sobre nuestros
  loci, quedarían como *unlabeled*, y el modelo aprendería que un miRNA real es
  un candidato novedoso — el modo de falla exacto que PU learning evita, y sin
  ningún error visible.
  Dos cosas que se siguen de esto y es fácil hacer mal:
  **no** comparar la secuencia genómica completa del locus contra el maduro (el
  locus mide cientos de nt y el maduro ~22, la identidad global no significa
  nada); y **no** usar un umbral de identidad plano, porque el 5' define la
  semilla y es preciso mientras el 3' varía de rutina por isomiRs. El criterio
  va anclado en 5' con holgura en 3'.
- **Ensembl y NCBI no nombran los cromosomas igual** (`1` vs `NC_006088.5`).
  Los 3 organismos heredados vienen de EnsemblGenomes y los 6 nuevos de NCBI.
  Con el etiquetado por secuencia esto dejó de afectar a los positivos, pero
  sigue valiendo para cualquier otro cruce con coordenadas externas.
- **YASMA v1.1.1 escribe en `annotations/<nombre>/loci.gff3`**, no en
  `annotation/`. Un chequeo contra la ruta vieja hacía abortar el pipeline tras
  el primer organismo.
- **`fasterq-dump` sin `-t`** crea temporales en el CWD; llegaron a 109 GB.
- **Al matar el pipeline**: matar también `bowtie-align-s`, `fastp`,
  `fasterq-dump` y `samtools`, o quedan huérfanos escribiendo el mismo BAM.

## Red

Esta sesión cloud tiene **bloqueada por política** la salida a NCBI, Ensembl,
la ENA, miRBase y MirGeneDB: el gateway responde 403 al `CONNECT`. Solo pasan
repositorios de paquetes (npm, PyPI, crates) y GitHub. O sea que desde acá no
se puede resolver un manifiesto ni bajar un genoma.

Lo que necesita red va en otro lado:

| tarea | dónde | con qué |
| :-- | :-- | :-- |
| setup de la sesión | Colab | `notebooks/00_setup.ipynb` |
| genomas | Colab | `notebooks/descarga_genomas.ipynb` |
| manifiesto y `.sra` | Colab | `notebooks/10_descarga_runs.ipynb` |
| ver qué falta | Colab | `notebooks/90_estado.ipynb` |
| traer `.sra` para alinear | máquina local | `scripts/drive_pull.sh sra <org> --go` |
| alineamiento y YASMA | máquina local | `orchestrate.sh` |
| BAMs a Drive | máquina local | `scripts/drive_push.sh` |

**Colab es el administrador de datos**: baja, valida y escribe a Drive sin pasar
por el disco local. Con Free el cómputo largo (30-40 h de bowtie) se queda en la
máquina local, que trae los `.sra` de a un organismo con `drive_pull.sh` y los
purga después. Ver `docs/colab.md` y `docs/plan_datos_colab.md`.

## Entorno

micromamba, entorno `srna2` (se llama así porque el primer intento usó Python
3.11 y YASMA exige >= 3.12). Pin en `environment.yml` + `environment_pip.txt`.
`./check_env.sh` verifica que cada herramienta **arranque** de verdad, incluido
`import RNA` de ViennaRNA — el fallo clásico.
