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
primario de `maggi` son dos BioProjects combinados). `data/srr_manifest.tsv` es
el manifiesto resuelto: **415 corridas**, generado con
`./scripts/fetch_runs.sh manifest` contra la ENA desde Colab
(`notebooks/10_descarga_runs.ipynb`).

**415 de 415 bajadas en `80_sra/`**, con los 415 md5 en `data/sra_md5.tsv` y
reconciliado contra el manifiesto en las dos direcciones. Eran 416 hasta que
`SRR23277331` salió por no ser sRNA-seq — ver más abajo.

**Falta una corrida, y el conteo no lo delata.** El manifiesto resuelve 18 de
los 19 BioProjects de la spec: el duplicado nuevo de `sclsc` (`PRJNA1135930`)
entró a `organismos.tsv` y el manifiesto todavía no se regeneró contra la ENA,
así que `415 de 415` está completo **respecto de la spec anterior**. Se ve con
`./scripts/check_docs.py`, que compara los dos ficheros en vez de creerle al
total; a ojo no se ve, porque el manifiesto al día también da 416 —una corrida
menos por la exclusión, una más por `sclsc`—.
Las dos últimas —`SRR317135` y `SRR1066790`, del primario de `maggi`— costaron
una ronda entera porque **solo existen en formato SRA Lite** y el script buscaba
únicamente `.sra`: `prefetch` las bajaba bien, salía con código 0, y el
`.sralite` se quedaba en el staging sin que nadie lo mirara. Ver la trampa de
`prefetch` más abajo, que es donde está lo que hay que declarar en métodos.
`data/srr_manifest_r1.tsv` es el de la ronda anterior, como referencia.

**Los genomas están cerrados** (9 de 9 con `sha256` versionado) y los `.sra`
casi: falta la corrida de `sclsc` de arriba. Se comprueba con
`./scripts/fetch_runs.sh estado`, `./scripts/fetch_genomes.sh verificar` —que
recalcula los sha256 contra el ledger en vez de confiar en el tamaño— y
`./scripts/check_docs.py`, que cruza los docs contra `data/`. Lo que sigue es
el alineamiento.

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
| sclsc | *Sclerotinia sclerotiorum* | Fungi (Ascomycota) | PRJNA477286 | PRJNA1135930 |
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
- **Los 9 ensamblados fijados y verificados contra NCBI.** Bajados a
  `70_genomas/` con `sha256` en `data/genomas.sha256`; `verificar` recalcula los
  checksums contra ese ledger, así que la integridad está comprobada y no
  inferida del tamaño. El porqué de cada elección está en la columna `nota` de
  `data/genomas.tsv`.

  **Los 3 heredados de R1 ya tienen accession**, que era el hueco más serio:
  hasta ahora no había registro de contra qué FASTA se había alineado. Se
  resolvieron **sin `config.sh`**, desde el nombre del assembly. Los dos de
  Ensembl eran el mismo malentendido —NCBI está en una versión posterior del
  mismo linaje: `ASM14694v1` → `ASM14694v2` y `Phypa_V3` → `Phypa V5`, mismo
  accession base— y `rhirr` tenía el suyo **retirado** (`ASM43914v3`), pero su
  reemplazo `GCF_026210795.1` es la referencia vigente, es DAOM 197198 igual que
  el viejo, y es 17× mejor en N50.

  **Tres ensamblados estaban `suppressed`** y dos de ellos son los más citados
  de su especie. Ver la trampa correspondiente abajo: es la justificación
  concreta de por qué existe el paso de verificación.

  **`cloro` usa la cepa de los datos, no la referencia de la especie.**
  `GCA_902827195.2` (`C_rosea_IK726`), porque `PRJEB43636` son mutantes
  Dicer-like del grupo de Karlsson (SLU) sobre IK726, y la coincidencia de cepa
  no se recupera de ninguna otra forma. `resolve` dice `DIFIERE` y está bien que
  lo diga. Son 70.7 Mb contra una mediana de 55.2 Mb en las otras 20 cepas y
  ~58 Mb que declara la publicación; la anomalía quedó sin explicar, así que
  **hay que chequear después del alineamiento si aparecen loci duplicados con el
  mismo RNA mayoritario**, que es cómo se vería contenido duplicado.

  El paso manual sigue siendo a propósito: un ensamblado equivocado no falla
  ruidosamente, alinea peor y contamina la anotación. Ya no hace falta que
  coincida con el ensamblado de MirGeneDB — con el etiquetado por secuencia se
  elige por contigüidad y completitud.
- **`Magallana gigas` = `Crassostrea gigas`.** El género se renombró; Ensembl
  Metazoa y buena parte de las bases todavía usan *Crassostrea*. Buscar el
  genoma por el nombre nuevo no va a encontrarlo.
- **`SRR23277331` se sacó del manifiesto: era mRNA, no sRNA-seq.** Medido, no
  supuesto: `perfil` encontró **0 de 40 000 reads con adaptador 3'**. Es
  concluyente porque la lista incluye `AGATCGGAAGAGC`, el universal de Illumina,
  así que un inserto corto habría dado read-through con cualquier kit — cero
  significa que todos los insertos pasan los 150 nt. **El manifiesto queda en
  415**, `prupe` primario en 8 corridas, y el ledger reconcilia en cero por
  ambos lados. El `.sra` de 3 GB sigue en Drive como peso muerto; `90_estado` lo
  va a listar como "sobra", que es correcto.

  **La lección importante no es esa corrida, es que la etiqueta engañó.** Decía
  `miRNA-Seq` y pasó el filtro. `avg_len` tampoco distingue: un sRNA de 22 nt
  corrido en 2×150 da 273 nt igual que un mRNA. Lo único que separa los dos
  casos es dónde empieza el adaptador, y por eso existe
  `./scripts/fetch_runs.sh perfil`.
- **Los 19 proyectos verificados: todos son sRNA-seq.** `perfil --proyectos`
  mide dónde empieza el adaptador 3' en una corrida de cada uno; de los 18 que
  tenía el manifiesto, 17 dieron `PARECE sRNA-seq` y `cloro PRJEB43636`
  `YA RECORTADA`. El 19º, `PRJNA1135930`, se perfiló aparte antes de adoptarlo
  —por eso `perfil` acepta un `PRJ*` suelto— y da `PARECE sRNA-seq`.
  Los tres primarios etiquetados `RNA-Seq` que preocupaban son sRNA-seq de
  verdad: la etiqueta de la ENA está mal, la librería no. **El único dato que no era sRNA-seq en todo
  el dataset fue `SRR23277331`**, ya fuera del manifiesto.

  Tres cosas del perfilado que conviene tener a mano cuando se mire el QC de
  `fastp`, para no perseguir un bug que no existe:

  - **`cloro PRJEB43636` viene ya recortado** (reads de 35 nt, 0% de adaptador
    porque el read *es* el inserto). No necesita recorte.
  - **`gadmo PRJNA328800`: el 37% de los reads tienen inserto de 10 nt** y
    `phypa PRJNA222997` un 12% de inserto 0 — dímeros de adaptador. `fastp` los
    descarta por el piso de 15 nt, así que esos dos proyectos van a mostrar una
    retención baja **a propósito**, no por un fallo.
  - **`maldo PRJNA784097` (151 nt) está resuelto**: 100% de adaptador con
    inserto modal 24 nt, o sea que el recorte normal de `fastp` lo recupera. No
    hace falta pre-trim especial.
- **Reads de más de 50 nt: 143 de 415, y ninguna es un problema.** Verificado
  con `perfil`: 51 nt es una librería de 50 ciclos sin recortar, 65-75 nt una de
  75, y los 151 nt de `maldo` traen el adaptador al nt 24. `fastp` las resuelve
  todas recortando. El resumen del manifiesto agrupa por proyecto en vez de
  listar corrida por corrida, porque cortaba en 20 de 144 y enterraba lo que
  había que ver.
- **`sclsc` ya tiene duplicado: `PRJNA1135930`.** Reemplaza a `PRJNA985401`,
  que era RNA-Seq PAIRED y caía entero en el filtro. Verificado con `perfil`
  antes de adoptarlo: 98% de adaptador, inserto modal 22 nt, 97% dentro de la
  ventana. Con esto **los 9 organismos tienen primario y duplicado**.

  **Es UNA sola corrida (12 M spots), y eso va declarado en métodos**: la
  validación de `sclsc` es más débil que la de los otros 8. Alcanza para lo que
  el duplicado existe —preguntar si un locus priorizado reaparece en un
  experimento independiente— pero no da para nada estadístico.

  El otro candidato, `PRJNA379694` (6 corridas, 758 M spots, etiquetado
  `miRNA-Seq`), **se descartó**: `perfil` dio 1% de adaptador con reads de 100 nt
  e insertos de 71-87 nt. Es mRNA. Los 126 M spots por corrida ya lo hacían
  sospechoso y perfilarlo antes de adoptarlo evitó cambiar un duplicado que no
  servía por otro que tampoco.
- **`cloro` primario es `ncRNA-Seq`**, no miRNA-Seq. Igual que `phypa`, hay que
  declararlo en métodos. (Que venga ya recortado está anotado arriba.)
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
- **Un ensamblado se retira y nadie avisa.** Tres de los nueve tenían
  ensamblados que NCBI ya había marcado `suppressed`: GDDH13 v1.1
  (`GCF_002114115.1`, `maldo`), cgigas_uk_roslin_v1 (`GCF_902806645.1`, `maggi`)
  y ASM43914v3 (`GCF_000439145.1`, `rhirr`). Los dos primeros son de los más
  citados de su especie; el tercero **estaba en uso desde R1**. Descargarlos
  habría funcionado sin error visible.
  Es la justificación concreta de por qué `fetch_genomes.sh` no baja nada en
  estado `candidato` y por qué `resolve` ahora grita `RETIRADO por NCBI`. Los
  reemplazos son `GCF_042453785.1` (GDT2T_hap1, Golden Delicious T2T, el mismo
  cultivar del que deriva GDDH13) y `GCF_963853765.1` (xbMagGiga1.1, de Darwin
  Tree of Life). El de ostra importa más porque `maggi` es de entrenamiento: el
  de Roslin retuvo haplotigos —la misma región dos veces— y eso infla el
  multimapeo justo bajo `bowtie -m 50`.
- **Sacar una corrida editando el manifiesto a mano no sirve**: la próxima
  `fetch_runs.sh manifest` la vuelve a agregar, en silencio. Por eso existe
  `data/excluidas.tsv`, que el generador aplica y reporta. Una exclusión solo se
  agrega **con evidencia medida** —la columna `motivo` dice qué se midió— y hoy
  tiene una sola fila, `SRR23277331`.
- **Un 0% de adaptador no significa lo mismo con reads cortos que con largos.**
  `perfil` daba `NO PARECE sRNA-seq` a cualquier cosa sin adaptador, y marcó las
  34 corridas de `cloro PRJEB43636` — que vienen **ya recortadas**, con reads de
  30-34 nt donde el read es el inserto. Casi se tiran datos buenos por un defecto
  de la herramienta. Ahora reporta la longitud mediana del read y separa
  `YA RECORTADA` (sin adaptador, reads cortos) de `NO PARECE` (sin adaptador,
  reads largos, o sea inserto mayor que el read).
  Dos defectos más del mismo veredicto, corregidos a la vez: usaba una ventana
  propia de 18-30 nt **que contradecía la de `fastp`** —15-50, elegida para no
  truncar tRFs ni piRNAs—, y su lista de adaptadores era solo moderna.
  Ese último resultó ser el peor: con los adaptadores de 2011 agregados,
  `maggi PRJNA154615` pasó de **0% a 92%** y `phypa PRJNA277372` de **0% a 97%**.
  Los dos daban `NO PARECE sRNA-seq` y son sRNA-seq perfectamente normales —
  y uno es el primario de un organismo de entrenamiento. **Los cuatro proyectos
  que la herramienta marcó eran defectos de la herramienta, no de los datos.**
- **Un chequeo que mira el mensaje y no el dato no chequea nada.** El banco de
  `perfil` verificaba la línea `>>> YA RECORTADA` que se le imprime a la
  persona, pero no la variable `ver` que va a la tabla del resumen y **decide el
  exit code**. Son dos strings distintos en el mismo `awk`: se puede cambiar el
  veredicto de la tabla —regresionando el arreglo de `cloro`, el que casi tiró
  34 corridas buenas— y los 176 chequeos seguían verdes. Se descubrió mutando el
  código a propósito y viendo qué banco *no* se ponía rojo, que es la única
  forma de saber si un banco sirve. Ahora hay una fila `YA RECORTADA` en la
  tabla de `--proyectos` y otra en el caso "todos buenos", que es la que cubre
  el exit code. Al agregar un chequeo: afirmar sobre el valor que el programa
  *usa*, no sobre el texto que imprime al lado.

- **Los docs derivan igual que los datos, y el total no lo delata.** La tabla
  del dataset decía que el duplicado de `sclsc` era `PRJNA985401` cuando la spec
  ya decía `PRJNA1135930`, y la prosa correcta estaba 130 líneas más abajo —
  quien lee la tabla no llega ahí. El `README` era peor: mandaba a correr cinco
  comandos y cuatro no existen en este repo. Y el conteo no avisa: el manifiesto
  viejo y el al día dan los dos 416, una corrida menos por la exclusión y una
  más por `sclsc`. Por eso existe `./scripts/check_docs.py`, que recalcula desde
  `data/` en vez de creerle a lo escrito a mano, y falla si un doc cita un
  fichero que no está. Los seis del pipeline de alineamiento (`config.sh`,
  `orchestrate.sh`, `verify.sh`, `check_env.sh`, `environment.yml`,
  `environment_pip.txt`) están declarados ahí como deuda: el día que lleguen del
  `main` local se borra la línea y el chequeo vuelve a exigirlos.

- **Un campo que no se lee se ve igual que un campo vacío.** `resolve` mostraba
  `GCF_026210795.1` (`rhirr`) sin cepa, y parecía un ensamblado sin aislado
  declarado. Era que el formateador leía la cepa **solo** de
  `.organism.infraspecific_names.strain`, y NCBI también la publica en
  `assembly_info.biosample.attributes[]` como `strain` o `isolate`. Ahora cae en
  cascada por los tres y, cuando de verdad no hay ninguno, imprime `cepa=?` en
  vez de omitir la columna — porque una columna ausente se lee como un hecho.
  Mismo patrón que la paginación truncada de abajo: el límite de la herramienta
  disfrazado de dato.
- **Un listado paginado que se trunca miente en silencio.** `listar_cepas` pedía
  `page_size=20` y para `cloro` devolvió exactamente 20 — el límite— e imprimía
  el parcial como si fuera todo. *Clonostachys rosea* tiene **78** ensamblados, y
  el de IK726 era uno de los 58 que no se veían. Se estuvo a un paso de elegir
  otra cepa por un defecto de la herramienta, no por los datos. Ahora pagina
  hasta agotar y avisa si `total_count` no coincide con lo que trajo. Lo mismo
  vale para cualquier otra consulta paginada que se agregue.
- **Ensembl y NCBI no nombran los cromosomas igual** (`1` vs `NC_006088.5`).
  Los 3 organismos heredados vienen de EnsemblGenomes y los 6 nuevos de NCBI.
  Con el etiquetado por secuencia esto dejó de afectar a los positivos, pero
  sigue valiendo para cualquier otro cruce con coordenadas externas.
- **YASMA v1.1.1 escribe en `annotations/<nombre>/loci.gff3`**, no en
  `annotation/`. Un chequeo contra la ruta vieja hacía abortar el pipeline tras
  el primer organismo.
- **El BAM que va a YASMA tiene que traer `@RG` por corrida.** `yasma tradeoff`
  hace `header['RG']` sin `.get()`: un BAM sin read groups no degrada, tira
  `KeyError: 'RG'`. Y no es cosmético — agrega profundidad por read group, así
  que sin `@RG` no podría separar librerías aunque no se cayera. Medido; ver
  `docs/yasma.md`.
- **`yasma align` no es nuestro camino.** Envuelve a `ShortStack` y pediría las
  librerías ya recortadas; usarlo reemplazaría nuestro bowtie y con él el
  `-m 50`, que está medido. Lo que usamos es `yasma tradeoff -a <BAM>`, que
  consume **nuestro** alineamiento.
- **`yasma adapter` no reemplaza a `perfil`: su criterio es más débil.** Marca
  `PRE-TRIMMED` con `read_length_freq < 0.8 and best_perc < 0.10`, o sea que usa
  la **dispersión** del largo del read, no su magnitud. Corrido cabeza a cabeza
  contra `perfil` sobre las mismas lecturas, coinciden en `cloro` y en el caso
  mRNA, pero **una librería recortada a un solo largo le sale `None`, igual que
  el mRNA** — el mismo bug que `perfil` tuvo y se arregló. Ninguno de los 19
  proyectos está en ese estado, así que acierta en los 19, pero por suerte: el
  único pre-recortado es `cloro`, y sus largos de 30-34 nt son justo la
  variación que su criterio necesita. El que decide sigue siendo `perfil`.
- **Clonar YASMA por defecto no te da la versión pineada.** El branch por
  defecto del repo es `library-scaling`, no `main`, y deja 1.1.0. Peor: el
  `pyproject.toml` **dentro del tag `v1.1.1` declara `version = "1.1.0"`**, así
  que los metadatos del paquete instalado no sirven para verificar qué versión
  es. Se pinea por ref —`git+https://github.com/NateyJay/YASMA@v1.1.1`— y no por
  número. También: `import RNA` de ViennaRNA está en el nivel superior de
  `hairpin`, que `__init__` importa, así que sin ViennaRNA **no arranca ningún
  subcomando**, ni los que no tocan estructura.
- **`prefetch` puede salir con código 0 sin dejar un `.sra`.** Le pasa a las
  corridas que no tienen *SRA Normalized Format*: baja un `.sralite` y sale
  contento. También sale 0 cuando el resolver no encuentra nada. Por eso
  `fetch_runs.sh` no confía en el código de salida, busca el archivo, y cuando
  no está imprime el log de `prefetch` y lista el staging en vez de fallar
  mudo. No volver a mandar la salida de `prefetch` a `/dev/null`: un fallo que
  no se puede leer cuesta una corrida entera para diagnosticarse.
- **Dos corridas son SRA Lite y eso va en métodos.** `SRR317135` y
  `SRR1066790` (primario de `maggi`) solo existen en ese formato — `prefetch`
  dice explícitamente que prefiere el normalizado y que cae a lite *due to
  current file availability*, o sea que no hay versión full que pedir. Se
  aceptan a propósito, pero **SRA Lite guarda una sola calidad sintética para
  todas las bases**: en esas 2 corridas el filtro de calidad de `fastp` no
  descarta nada y en las otras 414 sí. La columna `formato` de
  `data/sra_md5.tsv` dice cuáles son; el archivo en Drive no, porque se guarda
  como `<RUN>.sra` igual que los demás. **Pendiente de confirmar contra
  `config.sh`**: si `bowtie` corre en modo `-v` (conteo de mismatches, el
  estilo ShortStack3) la calidad se ignora y el asunto se agota en `fastp`; si
  corriera `-n`/`-e` (suma de calidades del seed), las calidades sintéticas
  cambiarían el alineamiento de esas 2 corridas.
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

## Chequeos

Todo corre sin red y en segundos. Antes de cada push:

```bash
./tests/run_all.sh              # 9 bancos, 196 chequeos, binarios falsos en el PATH
./tests/mutar.py                # rompe el codigo y exige que algun banco grite
./scripts/check_docs.py         # lo que afirman los docs contra data/
./scripts/validate_notebooks.py # los .ipynb parsean y no hay duplicados
```

Los bancos encontraron **diez bugs** que ninguna lectura del código había visto,
y `check_docs.py` dos más.

**Un banco que pasa no prueba nada.** Prueba algo el día que se rompe lo que
cubre y el banco se queja, y la única forma de saberlo es romper el código a
propósito: eso es `tests/mutar.py`, 14 mutaciones que tienen que dar todas
`[OK]`. Un `[HUECO]` es un chequeo que falta; un `[VIEJA]` es una mutación cuyo
patrón ya no existe, que tampoco prueba nada. Así aparecieron los dos huecos que
ninguna otra cosa mostró — el veredicto de `perfil` que iba a la tabla sin estar
cubierto, y `estado` sin banco.

`check_docs.py` falla hoy a propósito, con `PRJNA1135930` sin corridas en el
manifiesto. Es el pendiente de `sclsc`, no un falso positivo: se cierra
regenerando el manifiesto desde Colab.

## Entorno

micromamba, entorno `srna2` (se llama así porque el primer intento usó Python
3.11 y YASMA exige >= 3.12). Pin en `environment.yml` + `environment_pip.txt`.
YASMA se pinea por ref de git, no por número de versión — ver `docs/yasma.md`,
que documenta qué consume, qué decide y dónde no hay que confiarle.
`./check_env.sh` verifica que cada herramienta **arranque** de verdad, incluido
`import RNA` de ViennaRNA — el fallo clásico.
