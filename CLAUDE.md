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

1. Descarga 417 corridas sRNA-seq de 9 organismos (7.90 G spots)
2. `yasma trim` (cutadapt, 15-50 nt) → `yasma align` (bowtie1 `-v 1 -m 50`,
   pesado estilo ShortStack3)
3. Anotación de loci con YASMA v1.1.1
4. Features por locus → PU learning → candidatos priorizados

## Dataset

**9 organismos, 3 por reino, cada uno con BioProject primario y duplicado.**
`data/organismos.tsv` es la especificación: 19 accessions en 18 slots (el
primario de `maggi` son dos BioProjects combinados). `data/srr_manifest.tsv` es
el manifiesto resuelto: **417 corridas**, generado con
`./scripts/fetch_runs.sh manifest` contra la ENA desde Colab
(`notebooks/10_descarga_runs.ipynb`).

**417 de 417 bajadas en `80_sra/`**, con los 417 md5 en `data/sra_md5.tsv` y
reconciliado contra el manifiesto en las dos direcciones. `./scripts/check_docs.py`
lo comprueba y hoy pasa en verde. Eran 416 hasta que `SRR23277331` salió por no
ser sRNA-seq, y 415 hasta que el duplicado nuevo de `sclsc` sumó 2 — ver abajo.

**El conteo no delata una spec desactualizada**, y por eso el chequeo compara
los ficheros en vez de los totales: mientras faltaba `sclsc`, el manifiesto
viejo y el al día daban los dos 416 —una corrida menos por la exclusión, una más
por `sclsc`— así que a ojo eran indistinguibles.
Las dos últimas —`SRR317135` y `SRR1066790`, del primario de `maggi`— costaron
una ronda entera porque **solo existen en formato SRA Lite** y el script buscaba
únicamente `.sra`: `prefetch` las bajaba bien, salía con código 0, y el
`.sralite` se quedaba en el staging sin que nadie lo mirara. Ver la trampa de
`prefetch` más abajo, que es donde está lo que hay que declarar en métodos.
`data/srr_manifest_r1.tsv` es el de la ronda anterior, como referencia.

**El upstream está cerrado**: 417 `.sra` y 9 genomas en Drive, los dos con
checksum versionado, y los 19 BioProjects de la spec resueltos. Se comprueba con
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

Escala **medida sobre el manifiesto resuelto**, no la declarada por la spec:
193 corridas / 3.92 G spots en los primarios, 224 / 3.98 G en los duplicados.
**417 corridas y 7.90 G spots en total, contra 169 y 2.18 G del set anterior** —
del orden de 3.6× más. Ver "Consecuencias del cambio".

Hasta acá se venía declarando "~425 corridas y ~7.4 G spots", que salía de sumar
la columna `spots_M` de `data/organismos.tsv`. Esa columna es de cuando se
eligió cada proyecto y es **estimación**: el total real es medio G de spots más
alto. Es la misma equivocación de altura que la columna `runs` —ver la trampa
correspondiente— y por eso los números de escala salen ahora de
`data/srr_manifest.tsv`, que es lo que la ENA devolvió.

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
  417 sin ella**, `prupe` primario en 8 corridas, y el ledger reconcilia en
  cero por ambos lados. El `.sra` de 3 GB sigue en Drive como peso muerto; `90_estado` lo
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

  **Y con qué adaptador cada uno**, que es lo que el recorte necesita:
  14 de los 19 con RA3, 2 con el universal de Illumina (`maggi PRJNA1254880` y
  `sclsc PRJNA1135930`), 2 con el small-RNA de 2011 (`maggi PRJNA154615` y
  `phypa PRJNA277372`) y `cloro PRJEB43636` sin ninguno porque ya viene
  recortado. Está en `data/adaptadores.tsv`, emitido por
  `perfil --proyectos --tsv`.

  **Que cada prefijo esté inmediatamente 3' del inserto se verificó con los
  datos, no se supuso.** Si un prefijo estuviera desplazado aguas abajo del
  inicio real del adaptador, el "inserto" medido saldría inflado por ese
  desplazamiento. Dentro de `maggi` las **tres** familias dan el mismo inserto
  modal (22/22/22), y en `phypa` las dos dan 21. Coinciden, así que los tres
  anclan bien como `-a` de cutadapt.

  Cuatro cosas del perfilado que conviene tener a mano cuando se mire el QC del
  recorte, para no perseguir un bug que no existe:

  - **La retención es la columna `retencion_est`, no `adapt_pct`.** cutadapt
    descarta dos cosas: los reads sin adaptador (`--trimmed-only`) **y** los que
    quedan fuera de 15-50. La retención es la intersección, y la brecha puede
    ser enorme: `galga PRJEB12164` tiene **95% de adaptador y retiene 52%**,
    porque el 20% de sus insertos mide 6-7 nt y muere en el piso de 15. Los tres
    más bajos son `gadmo PRJNA328800` (51%), `galga PRJEB12164` (52%) y
    `cloro PRJEB51338` (57%).
  - **`cloro PRJEB43636` retiene 100% por un motivo distinto del resto**: es
    `PRE-TRIMMED`, así que YASMA la pasa de largo en el recorte y **no le aplica
    ningún filtro**, ni el de longitud. Es la única de las 19 que sale del
    recorte sin pasar por la ventana 15-50. Dónde se filtra entonces depende del
    alineador: `yasma align` aplica 15-50 y descarta los reads con N **al
    alinear** (`XY:Z:F`), así que ahí se empareja con las otras 18; un bowtie
    propio no lo haría y habría que filtrarla aparte.
  - **`cloro PRJEB43636` viene ya recortado** (reads de 36 nt, 0% de adaptador
    porque el read *es* el inserto). Va `PRE-TRIMMED`: YASMA lo pasa de largo.
  - **`gadmo PRJNA328800`: el 28% de los reads tienen inserto de 10 nt** y
    `galga PRJEB12164` un 20% entre 6 y 7 nt — dímeros de adaptador. El piso de
    15 nt los descarta, así que esos dos pierden **el doble**: lo que no tiene
    adaptador y lo que queda demasiado corto.
  - **`maldo PRJNA784097` (151 nt) está resuelto**: 99% de adaptador con inserto
    modal 24 nt. No hace falta pre-trim especial.
- **Reads de más de 50 nt: 143 de 417, y ninguna es un problema.** Verificado
  con `perfil`: 51 nt es una librería de 50 ciclos sin recortar, 65-75 nt una de
  75, y los 151 nt de `maldo` traen el adaptador al nt 24. `fastp` las resuelve
  todas recortando. El resumen del manifiesto agrupa por proyecto en vez de
  listar corrida por corrida, porque cortaba en 20 de 144 y enterraba lo que
  había que ver.
- **`sclsc` ya tiene duplicado: `PRJNA1135930`.** Reemplaza a `PRJNA985401`,
  que era RNA-Seq PAIRED y caía entero en el filtro. Verificado con `perfil`
  antes de adoptarlo: 98% de adaptador, inserto modal 22 nt, 97% dentro de la
  ventana. Con esto **los 9 organismos tienen primario y duplicado**.

  **Son 2 corridas y 32.4 M spots** (`SRR31851668` 12.4 M, `SRR31851669`
  20.0 M). Al elegirlo se había anotado **1 corrida de 12 M**, y eso se escribió
  acá como un hecho de métodos —"la validación de `sclsc` es más débil"— cuando
  venía de la columna `runs` de la spec, que es una estimación. El manifiesto
  resuelto contra la ENA devolvió 2. **Sigue siendo el duplicado más chico de
  los 9** —la mediana está en 23 corridas— así que la validación de `sclsc`
  efectivamente es la más débil del set y eso va en métodos; pero no tanto como
  se había declarado.

  **Y es más débil todavía de lo que dicen esos 32.4 M spots: alineado contra
  `GCF_000146945.2`, el 67.1% de los reads no alinea en ninguna parte.**
  Primer proyecto alineado del set, y lo midió `align.sh verificar`: 4.1%
  único, 12.9% multimapeado con peso, 15.6% por encima de `--max_random 3`,
  0.1% por encima de `-m 50`, 67.1% sin ningún alineamiento válido. O sea que
  **~33% de la librería es de este hongo** y la profundidad efectiva de la
  validación de `sclsc` es de ~11 M reads, no 32.
  No es el genoma: contra un ensamblado equivocado bowtie da ~0% alineado, no
  33%, y el recorte tampoco —99% de retención, inserto modal 22 nt, cero reads
  descartados por el filtro de longitud (`XY:Z:F` = 0). La hipótesis a
  comprobar es que **`PRJNA1135930` sea un experimento *in planta***:
  *S. sclerotiorum* es necrótrofo y su sRNA se estudia sobre todo por RNAi
  entre reinos, así que una librería de tejido infectado sería mayoritariamente
  del huésped. **La comprobación barata es alinear `sclsc_primario`**
  (`PRJNA477286`, que debería ser cultivo puro) y comparar la fracción: si da
  80-90%, el 67% es el experimento y no el pipeline. Hasta entonces no se
  declara nada en métodos.
  Lo otro que llama la atención y **sí es esperable**: el multimapeado (12.9%)
  triplica al único (4.1%). Los sRNA de *Sclerotinia* que se describen como
  efectores son derivados de retrotransposones, o sea multicopia por
  definición. Es señal, no ruido — pero refuerza lo de `-m 50`: mirar la
  distribución de longitudes antes de tocarlo.

  El otro candidato, `PRJNA379694` (6 corridas, 758 M spots, etiquetado
  `miRNA-Seq`), **se descartó**: `perfil` dio 1% de adaptador con reads de 100 nt
  e insertos de 71-87 nt. Es mRNA. Los 126 M spots por corrida ya lo hacían
  sospechoso y perfilarlo antes de adoptarlo evitó cambiar un duplicado que no
  servía por otro que tampoco.
- **`cloro` primario es `ncRNA-Seq`**, no miRNA-Seq. Igual que `phypa`, hay que
  declararlo en métodos. (Que venga ya recortado está anotado arriba.)
- **Re-estimado con datos, no con reglas de tres.** Con 3.4× más reads, las
  8-12 h de alineamiento pasan al orden de **~100 h** (medidos: ~58 s por
  millón de reads recortados). Los BAMs **no** escalan igual: se habían
  declarado en ~340 GB escalando los ~100 GB del set anterior por 3.4, pero
  los ~100 de partida salían de estimar 45 B por alineamiento y lo medido es
  14.1, así que los 18 proyectos dan **~100 GB en total**. Entra sin problema en
  1.6 TB; lo que cambia es el cronograma, no el espacio.

## Trampas conocidas — no re-introducir

- **El recorte es `yasma trim`, que envuelve cutadapt — no `fastp`.** Sus
  defaults `--min_length 15 --max_length 50` coinciden con la ventana del
  proyecto, pero `scripts/trim.sh` los pasa explícitos igual: un default que
  cambie de versión no avisa. Tres consecuencias que van en métodos, y que
  **no** son las de `fastp`: `--trimmed-only` **descarta los reads sin
  adaptador**, así que la retención esperada es el `adapt_pct` de
  `data/adaptadores.tsv` y no ~100%; `--max-n 0` tira cualquier read con una N;
  y **no hay filtro de calidad** (cutadapt va sin `-q`), lo que cierra la mitad
  del asunto de SRA Lite — en el recorte la calidad sintética no distorsiona
  nada porque no se mira. Detalle que cuesta una corrida: `yasma trim` **no
  tiene `--override`**, aunque `yasma adapter` sí.
- **`yasma trim` pisa `trimmed_libraries` en vez de acumularlo, y eso rompe
  cualquier recorte por tandas.** Hace `ic.inputs['trimmed_libraries'] = []` al
  entrar y `= <lo de esta llamada>` al salir. Medido: tras recortar RUNA el json
  dice `['trim/RUNA.t.fq.gz']`, y tras recortar RUNB dice
  `['trim/RUNB.t.fq.gz']` — RUNA **sigue en disco y desapareció del registro**.
  Como `trim.sh` lee de ahí para saber qué está hecho, cada tanda desmentía a la
  anterior: `estado` reporta la corrida como faltante y `correr` la vuelve a
  volcar y recortar, para siempre. Por eso hay un ledger propio
  (`recortadas.tsv`, acumulativo) y por eso después de cada tanda se re-escribe
  `inputs.json` con la lista completa, que es lo que los comandos YASMA de aguas
  abajo leen. Pasarle todas las librerías en cada llamada tampoco sirve: no
  saltea las que ya tienen salida, las re-recorta.
  Del mismo palo, medidos a la vez: **`trim/log.txt` se trunca en cada llamada**
  (`Logger` lo abre con `"w"`), así que las estadísticas de cutadapt de las
  tandas previas se pierden y cada tanda se loguea aparte; y **`--cleanup` no se
  puede usar** —itera `ic.inputs['srrs']`, que en nuestro json es `None` →
  `TypeError`, y además vaciaría `untrimmed_libraries`, que para una librería
  `PRE-TRIMMED` es donde está la salida.
- **Volcar los 417 `.sra` a fastq antes de recortar son ~1.1 TB.** Calculado
  desde `base_count` del manifiesto (`2*bases + 35*reads`): `galga` solo son
  376 GB —245 el duplicado y 131 el primario— y eso no entra en el disco junto
  con los `.sra` y los ~100 GB de BAMs. `trim.sh` trabaja en **tandas acotadas
  por `PRESUPUESTO_GB`** (40 por defecto), borra el fastq sin recortar apenas la
  tanda termina, y vuelca la tanda siguiente mientras recorta la actual
  (`SOLAPAR=0` lo apaga). El pico es ~2× el presupuesto.
  **La excepción es `cloro PRJEB43636`**: al ser `PRE-TRIMMED`, YASMA anota el
  fastq sin recortar como su propia salida, así que ese fichero no se puede
  borrar — se guarda comprimido (53 GB → ~14).
- **El primario y el duplicado son proyectos YASMA distintos:
  `trim/<org>_<rol>/`.** 18 directorios, no 9. El duplicado es la validación
  independiente, y aguas abajo el directorio de proyecto **es** la unidad de
  `yasma tradeoff`: si comparten directorio, la separación queda en "acordarse
  de filtrar por la columna `rol`" y la validación deja de ser independiente el
  día que alguien no se acuerde. `maggi_primario` junta `PRJNA154615` y
  `PRJNA232734` en un solo directorio y está bien: el `adapters` de
  `inputs.json` es un dict **por fichero**, así que cada uno lleva su secuencia.
- **Recortar con la secuencia equivocada no falla: deja el `.t.fq.gz` casi
  vacío.** cutadapt corre con `--trimmed-only`, así que lo que no matchea se
  descarta y el pipeline sigue en verde — `estado` diría "recortadas, 0 faltan".
  Para eso existe `./scripts/trim.sh verificar`, que compara la retención
  **medida** (de los conteos de cutadapt, guardados en el ledger) contra la
  `retencion_est` de `data/adaptadores.tsv`, y grita `VACIA` por debajo del 5% y
  `DESVIADA` más allá de 15 puntos. Es el único chequeo del recorte que atrapa
  el caso de `maggi`: dos BioProjects de 2011 y 2014 con kits distintos, donde
  una sola fila de adaptador habría vaciado el proyecto que no corresponde.
- **Un BioProject puede mezclar kits, y `perfil --proyectos` no lo puede ver.**
  `gadmo/duplicado` (`PRJNA328800`, 12 corridas) recortó con la fila del
  proyecto y **6 de 12 retuvieron 0.6-2.0% contra el 51% esperado**; las otras 6
  dieron 38-58%, o sea bien. `perfil --proyectos` mide **una** corrida por
  proyecto —`!(($1 FS $3) in v)`, la primera— y le tocó una de las buenas, así
  que las otras 11 nunca se midieron. El comentario de `secuencia_de` decía
  *"se busca por proyecto porque el adaptador es del kit, no de la corrida"*: el
  kit sí es del proyecto, lo que no es cierto es que un BioProject use un solo
  kit.
  Tres piezas, y las tres hacían falta:
  `data/adaptadores.tsv` lleva **columna `run`** —`-` vale para todo el proyecto
  y una fila con el RUN exacto le gana—, existe
  **`perfil --corridas <org>/<rol>`** que mide todas una por una y emite esas
  filas con la clave llena, y **`trim.sh rehacer <org>/<rol> RUN...`** saca esas
  corridas del registro y borra su `.t.fq.gz`, porque el recorte es idempotente
  y si no las saltea para siempre.
  Lo que **sí** funcionó es el guardia: `trim.sh verificar` las marcó `VACIA`,
  el `assert` de §3 cortó, y nada llegó al alineamiento. Es el mismo caso de
  `maggi` un nivel más abajo — ahí eran dos BioProjects con kits distintos, acá
  es uno solo.
  **Cortar no alcanzaba: el paso siguiente eran seis accessions transcritas a
  mano** desde una tabla formateada, que es donde se cuela el error que después
  no falla ruidosamente. Por eso §3b del notebook saca `MALAS` de la columna
  `CORRIDA` de la propia tabla de `verificar` —con `correr(..., capturar=True)`—
  y §3c pega las filas y llama a `rehacer` con esa lista. Tienen banco propio
  (`tests/test_celda_trim.py`, 8 escenarios sobre la tabla real de `gadmo`) y
  cinco mutaciones, entre ellas que **si la tabla cambia de forma la celda
  revienta** en vez de terminar en un `rehacer` vacío que parece éxito.
- **`maggi` primario son dos BioProjects de eras distintas, y eso es una fila
  de adaptador cada uno.** `vdb-dump --info` da la fecha de carga:
  `SRR317135` (`PRJNA154615`) es de **julio de 2011** y `SRR1066790`
  (`PRJNA232734`) de **enero de 2014**. Justo el rango donde cambió el kit — es
  el proyecto que pasó de 0% a 92% de adaptador al agregar los prefijos viejos.
  `data/adaptadores.tsv` lleva una fila por proyecto, no por organismo, así que
  esto queda cubierto; pero llenarla con un solo adaptador para "maggi primario"
  sería un error, y no fallaría ruidosamente: `--trimmed-only` simplemente
  descartaría casi todo el proyecto cuya secuencia no corresponde.
- **La secuencia universal va con el prefijo común de 21 nt, no con la de
  TruSeq.** `AGATCGGAAGAGCACACGTCT` lo comparten TruSeq y el 3' SR de NEBNext
  Small RNA, y divergen después del nt 21. Un adaptador inmediatamente 3' de un
  inserto de 22 nt en una librería de sRNA es NEBNext antes que TruSeq; darle la
  versión larga de TruSeq a una librería NEBNext obliga a cutadapt a rechazar el
  alineamiento largo (6 diferencias en 27 nt = 22%, sobre el 10% por defecto) y
  recién después aceptar el corto. Funciona, pero dependiendo del umbral de
  error. Con 21 nt no hay mismatch posible para ninguno de los dos kits.
- **Una comilla simple adentro de un programa `awk` cierra la cadena del
  shell.** Todo el `awk` va entre comillas simples, así que un `'-'` en un
  **comentario** trunca el programa. El síntoma engaña: `awk` reporta un error
  de sintaxis apuntando a una línea de comentario, que es donde se quedó sin
  texto, no donde está el problema.
- **Un porcentaje sobre 5 reads no es una medición, y el mismo comando no puede
  dar dos respuestas.** `cloro PRJEB43636` salió con
  `adaptador: RA3 (100% de los que tienen)` — de **5 reads de 20 000**, o sea
  coincidencias al azar. La tabla para leer mostraba `RA3` y `19 nt` mientras el
  bloque `--tsv` del **mismo comando** ponía `-`: se había parcheado solo el
  segundo. Ahora el umbral está en un lugar —`hay_ad = con >= 0.2*total`— y es
  **el mismo que usa el veredicto** para decir que no hay adaptador, así que no
  pueden contradecirse. Una familia secundaria que redondea a 0% tampoco se
  lista.
- **El veredicto `PARECE sRNA-seq` pasa con el 30% de los insertos en la
  ventana, así que no puede afirmar que "el inserto cae dentro".**
  `gadmo PRJNA328800` pasa con el **51%** y su inserto modal es de **10 nt** —
  dímeros de adaptador—, o sea fuera de 15-50. El mensaje decía "el inserto cae
  dentro de la ventana" al lado de un `INSERTO 10 nt` en la tabla, que se lee
  como una contradicción. Ahora da la fracción medida y avisa aparte cuando el
  modal queda fuera, porque es justo lo que el recorte va a descartar.
- **`perfil --proyectos --tsv` emite las filas de `data/adaptadores.tsv`.**
  La tabla de `--proyectos` es para leer; pasar 19 filas de ahí a mano es
  exactamente donde se cuela un error que después no falla ruidosamente, solo
  recorta con la secuencia equivocada. El `--tsv` traduce familia → secuencia
  completa con el mapa `SECUENCIAS` del script, pone `PRE-TRIMMED` en las
  `YA RECORTADA`, y deja `-` en las que no se pueden recortar.
- **El adaptador no se detecta en tiempo de corrida: sale de
  `data/adaptadores.tsv`.** `yasma trim` ante un adaptador `"None"` **descarta
  la librería** —no la suma a `trimmed_libraries` y desaparece del pipeline sin
  error— y `yasma adapter` da `None` tanto para una ya recortada a largo fijo
  como para mRNA. Una librería mal clasificada se perdería en silencio. Un
  proyecto que no esté en la tabla hace fallar `trim.sh`, que es mejor que
  recortar con una secuencia adivinada.
- **El prefijo que detecta el adaptador no es la secuencia con la que se
  recorta.** Cotejando los 6 prefijos de `perfil` contra los 161 adaptadores de
  YASMA: `TGGAATTCTCGGG` es RA3 y lo comparten **50** entradas de la familia RPI
  (identifica la familia, no un adaptador); `GATCGTCGGACTG` es el
  `RNA_Adapter_(RA5)`, o sea un adaptador **5'** — encontrarlo es dímero o
  quimera, no read-through, y como `-a` de cutadapt sería incorrecto; y
  `CGCCTTGGCCGT` **no aparece en ninguno de los 161**. Por eso esos dos llevan
  `5p:` y `??:` en la lista, `perfil` avisa, y `trim.sh` se niega a recortar con
  ellos. (`ATCTCGTATGCCG` y `TCGTATGCCGTCTTCTGCTTG` sí son el adaptador de 2011;
  sus 110 coincidencias son constructos modernos que lo contienen aguas abajo.)
- **Los nombres de salida de `yasma trim` no se pueden adivinar.** Es
  `<RUN>.t.fq.gz`, no `<RUN>.tfq.gz` —hace `'.t' + library_format` y ese formato
  ya trae el punto— y una librería `PRE-TRIMMED` **no produce fichero nuevo**:
  anota la ruta original. Con el nombre adivinado nada contaba como recortado, o
  sea ni idempotencia ni `estado`. `trim.sh` lee `trimmed_libraries` de
  `inputs.json`, que es el registro que YASMA deja de lo que hizo.
- **La columna `runs` de `organismos.tsv` es una estimación, no un hecho** —
  difiere del manifiesto resuelto en **7 de los 19** proyectos (`prupe
  PRJNA780811` declara 5 y trae 7; `galga PRJNA694114` declara 96 y trae 95;
  los dos de `maggi` declaran `NA`). Eso está bien: el manifiesto es la verdad
  y la columna es de cuando se eligió el proyecto. **Lo que no está bien es
  convertirla en una afirmación de métodos**, que es lo que pasó con `sclsc`:
  declaraba 1 corrida, se escribió tres veces en este fichero que su validación
  era "más débil" por eso, y la ENA devolvió **2**.
- **Lo que le falta a git son dos mitades, no una.** §4 reportaba solo el
  ledger de md5, así que las filas del manifiesto había que sacarlas con un
  `grep` a mano — y commitear un md5 cuya corrida no está en el manifiesto hace
  fallar a `check_docs.py` con *"en el ledger y no en el manifiesto"*. Ahora
  reporta las dos. Y compara contra **lo que tiene git** (`git show HEAD:...`),
  no contra el working tree: §1 pisa `data/srr_manifest.tsv` del clon con la
  copia de Drive, así que leer ese fichero da el manifiesto de Drive y el diff
  sale vacío siempre. Tiene banco (`tests/test_celda4.py`, 6 escenarios).
- **Si la celda ya sabe que algo está viejo, que lo arregle ella.** §1 del
  notebook de descargas estuvo mal **tres veces**, y las tres fueron la misma
  equivocación de altura: primero reusaba la copia de Drive siempre (cambiar
  `organismos.tsv` o `excluidas.tsv` no tenía efecto visible); después detectaba
  la copia vieja pero solo la **imprimía** y seguía, así que §2-§4 corrían sobre
  el manifiesto viejo y §4 ofrecía la fila excluida para commitear; después
  **cortaba** con `RuntimeError`, que frena el daño pero deja a la persona
  moviendo un flag a mano para algo que la celda ya sabía hacer. Ahora
  **regenera sola** cuando el chequeo falla, y solo corta si la ENA sigue sin
  coincidir con la spec **después** de consultarla — que ya no es una copia
  vieja sino un problema de datos, y el mensaje manda a `fetch_runs.sh buscar`
  en vez de a mover un flag. `REGENERAR = True` quedó como forzado, no como
  requisito. Tiene banco propio (`tests/test_celda1.py`, 7 escenarios) porque el
  código de un notebook se rompe igual que el de `scripts/`.
- **Un `grep` cuyo patrón empieza con `-` lo toma como opción, y en un
  `notiene()` eso hace que el chequeo pase SIEMPRE.** grep aborta, el exit
  nonzero cae en la rama de éxito, y se imprime `ok`. Estaba en los **14**
  helpers de los bancos; todos llevan `--` ahora. Al escribir un helper que
  reciba un patrón como dato: `grep -qF -- "$2"`.
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
- **`yasma align` NO envuelve a ShortStack, y esta línea decía que sí.** El
  wrapper existe (`src/yasma/align.py`, `shortstack_align`) pero está
  **comentado en `__init__.py`**; el comando que se registra es
  `nativealign.py`, un **bowtie1 nativo con pesado estilo ShortStack3**. Y no
  costaría el `-m 50`: `--max_multi` **vale 50 por defecto**. Tres pasadas por
  librería —`-m 1`, después `-m 50 -a --best --strata` con sorteo ponderado por
  cobertura única local, y los que se pasan de 50 al BAM como no mapeados con
  `XY:Z:H`—, `@RG` por corrida en la cabecera y por read, y
  `align/alignment.bam` ordenado. Detalle en `docs/yasma.md`.
  **Es el alineamiento del proyecto**, vía `scripts/align.sh`: reemplaza al
  bowtie de `orchestrate.sh`, que nunca llegó a este repo.
- **Un proyecto no se puede partir en varias llamadas a `yasma align`.**
  `nativealign.py` crea `unique_d` una sola vez (línea 439) y lo acumula sobre
  **todas** las librerías del proyecto antes de que la etapa `multi` lo use para
  pesar (línea 628): el peso de cada posición multimapeada sale de la cobertura
  única **agrupada de todo el proyecto**. Partirlo no cambia la contabilidad,
  cambia a qué locus se asigna cada read multimapeado. Así que la unidad
  reanudable es el **proyecto entero**, no la corrida — que es lo que decide el
  diseño de `20_alinear.ipynb`: de a un proyecto, del más chico al más grande, y
  el BAM a Drive apenas termina.
- **`yasma align` levanta un bowtie por librería que nadie lee, y eso es un
  OOM.** `bowtie_generator` arma `bowtie_call` en tres ramas —`unique`, `multi`
  y `over`— pero solo las dos primeras le agregan flags; después hace el `Popen`
  **sin mirar `mmap`** (líneas 324-336). O sea que en `over` también arranca un
  bowtie: sin `-v`, sin `-m`, sin `-S`, sobre la librería entera. Esa salida
  **no se usa** —la rama `over` lee `<RG>.max50.fq` del disco y `p.stdout` no se
  toca nunca— y como abajo dice `if mmap != 'over': p.wait()`, tampoco se
  espera: el proceso queda vivo, bloqueado escribiendo a un pipe que nadie lee,
  con el índice del genoma entero en RAM. **Uno por librería, todos a la vez.**
  Medido: `gadmo_duplicado` (12 librerías, genoma de 670 Mb) murió con `Killed`
  al **96.5%, en la etapa `over`**, que es justo cuando ya hay once huérfanos
  vivos. `sclsc_duplicado` había pasado porque son 2 librerías y 39 Mb.
  `galga_duplicado` son **95 librerías contra 1.05 Gb**: no hay máquina donde
  entre.
  Lo arregla `scripts/yasma_parche.py`, que saltea ese `Popen` y **nada más**:
  el proceso que deja de levantarse es exactamente el que no se lee, así que la
  salida no cambia. Es idempotente, deja respaldo, y **se niega si el fuente
  cambió** —si el ancla no aparece una sola vez, o si desaparece el
  `if mmap != 'over':` que es lo que vuelve seguro saltearlo— porque un parche
  que aplica a ciegas sobre otra versión es peor que no tenerlo. `align.sh
  correr` lo exige antes de alinear. Tiene banco (`tests/test_parche.py`, que
  corre el fixture de verdad y mira **qué procesos arranca**, no qué dice el
  código) y 5 mutaciones.
- **Lo que decide si un proyecto entra no es solo el disco: es la RAM, y sale
  del tamaño del genoma.** `nativealign.py` hace
  `unique_d[ref] = [0] * int(length)` por cada secuencia del `.fai`, o sea **un
  entero de Python por base del genoma**. Medido en CPython: **8 bytes por
  elemento** (los punteros de la lista), y 32 más por posición que pase de 256,
  que en sRNA son pocas. No depende de cuántos reads haya **ni de cuántas
  librerías**: depende solo del ensamblado, y se arma entero antes de alinear el
  primer read. Para `gadmo` son ~5.4 GB; para `galga`, ~8.4.
  §1 del notebook predecía **disco** nada más, así que un proyecto salía
  "entra" y moría por memoria a las seis horas — el mismo fallo que esa celda
  existe para evitar. Ahora `align.sh` lo dice en `genoma`, en `plan` (columna
  `RAM_GB`) y antes de cada proyecto en `correr`, con `10 B por base`: los 8
  medidos más ~1.5 estimados del índice de bowtie a `--offrate 3`, que **es
  estimación y va declarada como tal**.
  Y cuando igual se queda sin memoria, `correr` distingue el **código 137**
  (128+9, SIGKILL) y lo dice: sin eso, un problema de memoria se reportaba
  igual que un genoma que falta. El proceso no llega a imprimir traceback, así
  que el único rastro es la palabra `Killed` y el % donde se cortó.
- **El alineamiento en Colab lo limita el disco, no el tiempo.** El pico es
  `recortado + 2 × BAM`, porque `pysam.sort` escribe el BAM ordenado **antes**
  de borrar el sin ordenar. §1 del notebook lo mide contra el disco real de la
  VM antes de empezar y lo dice — se sabe en un segundo o a las seis horas. Los
  `.sra` no se copian a la VM: Drive está montado y `SRA_DEST` apunta al mount —
  la regla es no **escribir** grande al FUSE, leer está bien.

  **Los bytes por read eran una estimación y estaban 3× de más.** Decía ~45 B
  por alineamiento en BAM; `sclsc_duplicado` midió **14.1** (32 M reads, BAM de
  432 MB), y 21.8 en `.t.fq.gz` contra los 25 estimados. Todo lo que se derivaba
  de ese 45 estaba inflado: `galga_duplicado` pasó de un pico de ~174 GB a
  **~79**, `maldo_primario` de ~76 a ~34, y el total de BAMs del proyecto de los
  ~340 GB que se declaraban a **~100**. §1 usa 16 B/read y no 14 a propósito: se
  midió en una librería con **67% de reads sin alinear**, y un read sin alinear
  ocupa menos que uno colocado.

  **Y la VM de Colab Free no daba 78 GB sino 220**, así que la línea de "entran
  16 de los 18" era falsa por los dos lados. Con lo medido entran los 18 en esa
  VM, y 17 de 18 en una de 78 GB (`galga_duplicado`, con ~79 GB de pico, es el
  único que se queda afuera). Nada de esto se declara de nuevo a ojo: §1 lo
  recalcula contra `shutil.disk_usage('/content').free` en cada sesión, porque
  el disco que toque varía.
- **Lo que persiste entre sesiones de Colab es Drive, no `/content`, y por eso
  "qué falta" no se puede leer del disco local.** `align.sh estado` mira
  `PROY_DIR`, que en una VM nueva está vacío: dice que falta **todo**, incluido
  lo que se alineó la semana pasada. Costó 31 minutos re-alineando
  `sclsc_duplicado` entero, con la salida idéntica byte por byte, porque
  `PROYECTO` de §2 estaba cableado y nadie lo cambió. Ahora §1 mira
  `10_bam/<org>/<rol>.bam` en Drive —donde solo llega lo que pasó `verificar`
  en §5, así que un BAM ahí es trabajo terminado y verificado—, marca esos como
  `YA EN DRIVE` y deja `SIGUIENTE` con el primer pendiente **que además entra en
  el disco**; §2 hace `PROYECTO = SIGUIENTE`. Rehacer uno sigue siendo posible
  poniéndolo a mano. Tiene banco (4 escenarios en `tests/test_celda_alinear.py`)
  y dos mutaciones, una de ellas para que `SIGUIENTE` no pueda caer en un
  proyecto que no entra — mandar a gastar horas en algo que se queda sin disco a
  la mitad es peor que no proponer nada.
- **Paralelizar dentro de una máquina no sirve, y el techo entre máquinas es
  4×.** `yasma align` ya corre `bowtie -p <cores>`: en una VM de Colab Free con
  2 vCPU un proyecto solo ya las satura, así que lanzar dos no divide el tiempo,
  lo reparte — y encima **suma las dos `unique_d`**, que son 8 B por base del
  genoma **cada una**. Dos proyectos de `gadmo` juntos piden ~12.4 GB de los
  ~11.4 que hay: mueren los dos, y más tarde que uno solo. Y un proyecto tampoco
  se puede partir, por lo de `unique_d` acumulada. O sea que **la unidad
  paralelizable es el proyecto entero en otra máquina**.
  Entre máquinas son ~101 h de trabajo total, pero `galga_duplicado` solo son
  ~25 h: 1 máquina 101 h, 2 → 51, 3 → 34, 4 → 25, y **de ahí no baja con
  ninguna cantidad**. La quinta no hace nada. Eso se dice antes y no se descubre
  habiendo conseguido cinco.
  `scripts/reparto.py` hace el reparto con **LPT** —del más caro al más barato,
  cada uno a la máquina menos cargada que pueda hospedarlo— porque al revés
  todas terminan los baratos y una queda sola con `galga_duplicado` al final.
  **Dentro** de cada máquina el orden vuelve a ser de chico a grande: no cambia
  el makespan y deja algo terminado si la sesión se muere.
  Dos cosas del resumen que hacían falta y no estaban:
  el speedup se calcula **solo sobre lo repartido** —decía 4.5× mientras las
  34 h de `galga` no estaban asignadas a ninguna máquina, que es justo el número
  con el que se decide cuántas conseguir—, y lo que ninguna máquina puede correr
  sale nombrado con sus horas en vez de repartirse igual. La RAM que se le
  declara a una máquina es la **disponible** (`MemAvailable`), no la nominal:
  `galga` cae justo en esa diferencia —con 12 GB entra y con 11.4 no— y ese es
  el caso que decide si va a Colab o a la máquina local.
  La coordinación entre máquinas son claims en `Drive/90_claims/`, con TTL de
  12 h para que una sesión que se murió no bloquee el proyecto para siempre.
  **No es exclusión mutua de verdad** —el FUSE de Drive no da atomicidad—: relee
  lo que escribió y comprueba que siga siendo suyo, lo que atrapa el caso común
  y no el de dos escrituras simultáneas. El costo de una colisión son horas
  perdidas, no datos corruptos: el BAM se escribe local y se copia a Drive al
  final. El claim se suelta **cuando el BAM ya está en Drive**, no al terminar
  `align`, porque si no queda una ventana en la que nadie lo tiene y nadie lo
  hizo.
- **El cronograma salía de UN punto, y el genoma influye.** Los ~98 h se
  extrapolaban de `sclsc_duplicado`: 32 M reads contra un genoma de **39 Mb**.
  `galga` es **1.05 Gb**, 27× más grande. Y no es un detalle de borde: los
  cuatro proyectos más caros son **el 54% de los reads del set** y tres de ellos
  tienen genomas de 650 Mb o más, así que si el `s/M read` escala con el genoma
  el número se va — y conviene saberlo antes de empezar, no a mitad de camino.
  Por eso existe `data/calibracion.tsv`: §5 del notebook escribe una fila por
  proyecto alineado (reads, genoma, segundos de **reloj**) y §1 la lee. Con un
  punto usa un `s/M` plano; con dos o más sobre genomas distintos ajusta
  `s/M = a + b · genoma_Mb` y lo aplica **por proyecto**. Nadie vuelve a copiar
  un número a mano — que es exactamente cómo `B_BAM` llegó a estar en 45 sin
  haberse medido nunca cuando lo real son 14.1.
  Dos puntos no son un modelo, y la celda no finge que sí: imprime el estimado
  plano y el ajustado **lado a lado** y avisa si difieren más de un 25%, que es
  la señal de que hace falta un tercer punto y no de que el segundo sea verdad.
  **Cuál medir tampoco se elige a ojo**: §1b rankea los pendientes por
  *(Mb de genoma nuevo) ÷ (horas que cuesta)* y solo entre los que entran en
  esa VM. Con el estado de hoy propone `gadmo/duplicado` —37 M reads contra
  670 Mb, ~37 min— que además es el que ya había que rehacer por el kit mixto.
- **`galga` no entra en la RAM de Colab Free, y eso no lo decía nada.** Con
  1.05 Gb de ensamblado, `unique_d` pide ~9.8 GB de los ~11.4 disponibles: los
  dos proyectos de `galga` quedan fuera por **memoria**, no por disco —
  `galga_primario` entra holgado en disco (27 GB de pico contra 220 libres) y
  antes salía como "entra". Es marginal y depende del término estimado del
  índice de bowtie, así que `align.sh genoma` da el número exacto de `unique_d`
  antes de gastar la hora; pero el default es no mandar a alinear algo que va a
  morir con `Killed` a las horas.
- **`yasma align` tarda ~58 s por millón de reads recortados.** Medido dos
  veces sobre `sclsc_duplicado` (32.1 M reads, genoma de 39 Mb, Colab Free):
  26 min de bowtie más 4 de `pysam.sort`, o sea **31 min de reloj de punta a
  punta**. Extrapolado a los 6068 M reads que sobreviven al recorte en los 18
  proyectos son **~95-105 h**. Es un solo proyecto medido y el genoma influye:
  `galga` (1.1 Gb) y `maggi` (650 Mb) son más lentos por read que un hongo de
  39 Mb, así que §1 usa 60 y no 58. §4 se cronometra solo y §5 imprime los s/M
  read de cada proyecto para ir corrigiendo `S_POR_M` en vez de arrastrar este
  número.
  **El primer valor que se puso acá fueron 49 s y estaba mal por medir de menos:**
  salió de la línea `time elapsed` de YASMA más la de `Sorting`, y entre las dos
  hay una etapa —escribir la tabla de abundancia— que no aparece en ninguna. Por
  eso §4 ahora cronometra la celda entera en vez de sumar lo que el programa
  dice de sí mismo.
- **Colab ya no solo lee de GitHub: `scripts/colab_git.py` empuja.** Reemplazó
  los bloques de "copiá esta salida al repo" de §4 de `10_descarga_runs` y §6 de
  `descarga_genomas`. Tres reglas que tienen banco
  (`tests/test_colab_git.py`, contra un repo bare de verdad):
  **nunca `git add -A`** —una VM de Colab tiene Drive montado en
  `/content/drive` y un `add -A` es donde se cuela un `.sra`— y solo rutas bajo
  `data/`; **el token no aparece en ningún mensaje**, porque git mete la URL con
  el token adentro en sus errores; y **un push rechazado revienta** tras un
  reintento rebasando, porque uno que nadie mira deja el resultado en Drive y no
  en git. El destino sale del `origin` del clon y no de una URL cableada: con la
  URL cableada, cualquier clon empujaba al repo de verdad apenas hubiera un
  `GITHUB_TOKEN` en el entorno. Las celdas arrancan con `REVISAR_PRIMERO = True`.
- **El clon NO está en `colab_git.py`, y es a propósito.** Ese módulo vive
  adentro del repo, así que no se puede importar antes de clonarlo: el clon es
  el bootstrap y tiene que estar en la celda. Duplicarlo en el módulo dejaría
  dos implementaciones de lo mismo —la trampa de las tres rutas de los `.sra`,
  un nivel más arriba— y la del módulo no correría nunca. La celda está copiada
  en los cuatro notebooks que clonan y **tiene que ser idéntica**:
  `scripts/validate_notebooks.py` falla si una deriva.
- **El genoma tiene que estar ADENTRO del `-o` de `yasma align`, y por symlink
  al directorio.** `inputClass.check()` hace
  `value.relative_to(self.output_directory)` sin protegerlo: un genoma
  compartido fuera del proyecto no da un mensaje, tira `ValueError`. Copiar el
  FASTA a cada proyecto no hace falta —`validate_path` usa `.absolute()` y
  **no** `.resolve()`, así que un symlink queda lexicalmente adentro— pero el
  symlink tiene que ser **al directorio** del organismo y no al fichero: el
  índice se construye en `genome_file.with_suffix(".rev.1.ebwt")`, o sea **al
  lado del FASTA que le pasaste**. Con un symlink por fichero habría 18 índices
  en vez de 9 y `bowtie-build` correría dos veces por organismo — en un genoma
  de 1 Gb, horas. `align.sh` hace `ln -sfn <genomes>/<org> <proyecto>/genome`.
- **`pysam.FastaFile` no lee un FASTA comprimido con `gzip` plano.** Los 9
  ensamblados están en `70_genomas/<acc>.fna.gz` hechos con `gzip -c`, y
  `make_bam_header()` de `yasma align` los abre con `pysam.FastaFile`. Medido:
  `g.fna` OK, `g.fna.gz` → `OSError error when opening file`. (Con `bgzip`
  andaría, pero no es lo que tenemos.) Por eso `align.sh genoma` **verifica el
  `sha256` del `.gz` contra `data/genomas.sha256` y recién entonces
  descomprime**: el `.fna` y el índice `.ebwt` son derivados, no se respaldan y
  están en `.gitignore`, pero el ledger sigue siendo el ancla de qué se alineó.
- **Alinear contra el genoma equivocado no falla: sale con 0 y 0% alineado.**
  Medido contra el `yasma align` real, con un genoma del mismo tamaño y otro
  azar: código de salida 0, `alignment.bam` escrito, `estado` diciendo "2 de 2
  proyectos con BAM al día". Es el mismo modo de fallo que recortar con el
  adaptador equivocado, un paso más abajo, y lo único que lo delata es la
  fracción alineada. Por eso `./scripts/align.sh verificar` lee
  `align/library_stats.txt` —los conteos por read group que YASMA escribe— y
  grita `MUY BAJA` por debajo del 10% alineado. También reporta la corrida del
  manifiesto que **no tiene `@RG` en el BAM**: `tradeoff` agrega por read group,
  así que una librería que no llegó desaparece del análisis sin ruido. Y avisa
  —sin fallar— cuando más del 50% se pasa de `-m 50`, que es el fenómeno de los
  tRF medido en `danre`: hay que mirar la distribución de longitudes antes de
  tocar `-m`, no al revés.
- **La fracción alineada sola no distingue dos problemas opuestos, y por eso
  hay dos columnas.** `sclsc_duplicado` salió con **17.2% colocado y 67.1% sin
  ningún alineamiento válido**, y pasó como `ok` porque el único umbral miraba
  lo colocado y 17.2 > 10. Son diagnósticos distintos: lo que **no alinea en
  ninguna parte** (`XY:Z:N`) son reads que no son de este genoma —ensamblado
  equivocado, contaminación, o el huésped si el experimento es de infección—,
  mientras que lo que alinea **de más** (`XY:Z:H`, por encima de `-m 50`) es un
  genoma repetitivo, que es lo de los tRF. Con una sola columna los dos casos se
  leen igual: "poco alineado". `verificar` ahora emite `ALIN` y `SIN_AL` por
  separado y avisa —sin fallar, porque el BAM está bien escrito— cuando `SIN_AL`
  pasa del 50%. Lo que queda entre las dos (`XY:Z:Q`, por encima de
  `--max_random 3`) alineó pero no se colocó, así que `ALIN + SIN_AL` no suma
  100 y no tiene por qué.
- **El BAM vive en dos lugares a propósito, y son hard links.** `yasma align` lo
  deja en `<proyecto>/align/alignment.bam` y anota esa ruta absoluta en
  `inputs.json`, que es de donde `tradeoff` la lee; moverlo rompe la anotación.
  Pero `drive_push.sh bam <org>` sube `bams/<org>/`. `align.sh` enlaza
  (`ln`, no `cp`) a `bams/<org>/<rol>.bam`: mismo inodo, cero disco de más, las
  dos rutas válidas. Si el hard link no se puede —otro filesystem— copia y
  **avisa**, porque ~100 GB duplicados son una decisión y no un detalle.
- **`proyectos/<org>_<rol>/`, no `trim/<org>_<rol>/`.** El directorio guarda el
  proyecto YASMA entero —`trim/`, `align/`, `annotations/`— así que llamarlo
  `trim` pasó a mentir en cuanto `yasma align` escribió adentro. La ruta sale de
  `ruta_proyectos()` en `scripts/_drive_lib.sh`, por el mismo motivo que
  `ruta_local`: `trim.sh` escribe y `align.sh` lee, y si cada uno tuviera su
  idea, `align` diría "NADA RECORTADO" — que se lee como un problema del recorte
  y no de la ruta. `tests/test_rutas.sh` lo verifica preguntándoles a los dos.
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
- **El `read_count` del manifiesto está verificado contra el archivo, no solo
  declarado por la ENA.** `vdb-dump --info` informa `SEQ`, que es el conteo de
  reads del propio `.sra`; para `SRR317135` (14 311 812) y `SRR1066790`
  (11 760 772) coincide exacto con el manifiesto. Es una comprobación cruzada
  barata que `./scripts/fetch_runs.sh diag <RUN>` deja a mano si alguna corrida
  da un conteo sospechoso.
- **Dos corridas son SRA Lite y eso va en métodos.** `SRR317135` y
  `SRR1066790` (primario de `maggi`) solo existen en ese formato — `prefetch`
  dice explícitamente que prefiere el normalizado y que cae a lite *due to
  current file availability*, o sea que no hay versión full que pedir. Se
  aceptan a propósito, pero **SRA Lite guarda una sola calidad sintética para
  todas las bases**: en esas 2 corridas el filtro de calidad de `fastp` no
  descarta nada y en las otras 414 sí. La columna `formato` de
  `data/sra_md5.tsv` dice cuáles son; el archivo en Drive no, porque se guarda
  como `<RUN>.sra` igual que los demás. **La mitad del recorte está cerrada**:
  `yasma trim` llama a cutadapt sin `-q`, así que no hay filtro de calidad que
  distorsionar. **Y la del alineamiento también, si se usa `yasma align`**: sus
  dos pasadas de bowtie van con `-v 1`, que cuenta mismatches e **ignora las
  calidades**. Lo único que quedaría abierto es el bowtie de `orchestrate.sh`,
  que todavía no está en este repo: si corriera `-n`/`-e` (suma de calidades del
  seed) las calidades sintéticas sí cambiarían el alineamiento de esas 2
  corridas. Leer sus flags cuando llegue del `main` local.
- **Tres scripts creían que los `.sra` vivían en tres lugares distintos.**
  `drive_pull.sh` los dejaba en `<repo>/sra_cache/<org>/`, `fetch_runs.sh` los
  buscaba en **`/home/dev/sra_cache`** —un absoluto con un usuario que no existe
  en ninguna máquina de este proyecto— y `trim.sh` en
  `$HOME/tesis_data/80_sra/<org>/`. El síntoma habría sido `FALTA el .sra` en
  las 417 después de bajarlas bien. Ahora la ruta sale de `ruta_local()` en
  `scripts/_drive_lib.sh`, que ya era la fuente única del mapa de fases por
  exactamente el mismo motivo —"si push y pull se desincronizaran, subirías a
  una carpeta y bajarías de otra"—, solo que el razonamiento valía una vuelta
  más ancha. `tests/test_rutas.sh` lo verifica preguntándole a cada script en
  vez de leer su fuente.
- **`qc/`, `figures/` y `trim/` no estaban en `.gitignore`.** Son
  subdirectorios de data adentro del repo; uno que se escape hace que
  `git status` muestre cientos de GB. El banco de rutas exige que **todos** los
  del mapa de fases estén ignorados, así que agregar una fase nueva sin su
  entrada ahora falla.
- **Un OAuth que falla no vuelve a rclone: queda un remoto guardado sin
  token.** El error de Google —`Error 401: invalid_client` cuando no pudo
  resolver el `client_id`— muere en la pestaña del navegador, y `rclone config`
  sigue adelante y ofrece `Keep this remote? y`. Después **todos** los comandos
  fallan sin decir de dónde viene. Por eso `drive_check.sh` mira el `token` del
  dump: es el rastro que deja. Y mira el `client_id`, que es la causa casi
  siempre: tiene que terminar en `.apps.googleusercontent.com` y no traer
  espacios pegados del copiar. Los dos errores que se confunden con éste son de
  la consola de Google, no de rclone: `redirect_uri_mismatch` es haber creado el
  cliente como *Aplicación web* en vez de *de escritorio*, y `Acceso bloqueado`
  es la cuenta sin agregar a **Usuarios de prueba**. Todo en `docs/rclone.md`.
- **`read` descarta tabs a la izquierda aunque `IFS` sea solo tab.** Espacio,
  tab y salto de línea son *IFS whitespace* para bash pase lo que pase, así que
  `IFS=$'\t' read -r a b c` sobre `"\t\t24"` deja `a=24`, no `a=''`. El
  chequeo nuevo de `drive_check` nació con ese bug: un `client_id` vacío se leía
  como el tercer campo y se reportaba como "trae espacios pegados". Lo agarró su
  propio banco al primer escenario. Para campos que pueden venir vacíos:
  `mapfile -t` con una línea por campo, o que el productor emita una palabra de
  veredicto en vez del valor crudo.
- **Dos formas de configurar rclone mal que no dan error: el remoto lista
  vacío.** (a) elegir `scope = drive.file` en vez de `drive` — rclone solo ve
  los ficheros que él mismo creó, y las carpetas de `tesis/` se hicieron a mano
  en la web; (b) poner `root_folder_id` apuntando a `tesis/` **y** dejar
  `DRIVE_ROOT=tesis`, con lo que rclone busca `tesis/tesis/80_sra/`. En los dos
  casos `drive_pull.sh sra <org> --go` baja **0 ficheros y sale con código 0**.
  `./scripts/drive_check.sh` los distingue; la guía está en `docs/rclone.md`.
- **`${VAR:-x}` usa el default también cuando `VAR` está definida pero vacía**;
  `${VAR-x}` solo cuando no está definida. Los tres scripts de Drive usaban
  `:-` para `DRIVE_ROOT`, así que `DRIVE_ROOT=''` —la salida de escape
  documentada para un remoto que ya tiene `root_folder_id`— no existía. Lo
  encontró el banco de `drive_check`.
- **`purge` borraba la copia local sin comprobar que estuviera en Drive.** Su
  docstring decía "nunca toca Drive", que es verdad y no es el peligro: para
  `bam`, `yasma`, `features` y `modelos` la copia local es lo que se **produjo**
  acá, y si todavía no se subió, purgar la pierde. Los BAMs de un organismo son
  horas de alineamiento. Ahora corre `rclone check --one-way` contra el remoto y
  se niega si falta algo, salvo para `sra` y `genomas`, donde Drive es la fuente
  y la copia local es descartable.
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
| recorte y alineamiento (los que entran) | Colab | `notebooks/20_alinear.ipynb` |
| empujar a git desde Colab | Colab | `scripts/colab_git.py` + un PAT en Secrets |
| configurar rclone | máquina local | `docs/rclone.md` + `scripts/drive_check.sh` |
| traer `.sra` para alinear | máquina local | `scripts/drive_pull.sh sra <org> --go` |
| recorte | máquina local | `scripts/trim.sh plan/correr/verificar/rehacer` |
| alineamiento | máquina local | `scripts/align.sh genoma/plan/correr/verificar` |
| repartir entre máquinas | cualquiera | `scripts/reparto.py --maquinas ...` |
| BAMs a Drive | máquina local | `scripts/drive_push.sh` |

**Colab es el administrador de datos**: baja, valida y escribe a Drive sin pasar
por el disco local. Y desde `20_alinear.ipynb` también recorta y alinea — pero
solo los proyectos que entran en **disco y RAM** de esa VM, que §1 del notebook
mide en cada sesión. Con los 220 GB que dio la primera corrida real el disco
alcanza para los 18; lo que deja afuera a los dos `galga` es la **memoria**,
que pide ~9.8 GB de los ~11.4 disponibles. Esos van a la máquina local —que
trae los `.sra` de a un organismo con `drive_pull.sh` y los purga después— o a
Colab Pro high-RAM.

**Y se puede trabajar en varias máquinas a la vez.** `scripts/reparto.py` da la
tanda de cada una y §1c del notebook la toma; los claims en `Drive/90_claims/`
evitan que dos hagan el mismo proyecto. El techo es 4× y lo marca
`galga_duplicado`, que sola son ~25 h. Ver `docs/colab.md` y
`docs/plan_datos_colab.md`.

## Chequeos

Todo corre sin red y en segundos. Antes de cada push:

```bash
./tests/run_all.sh              # 21 bancos, 682 chequeos, binarios falsos en el PATH
./tests/mutar.py                # rompe el codigo y exige que algun banco grite
./scripts/check_docs.py         # lo que afirman los docs contra data/
./scripts/validate_notebooks.py # los .ipynb parsean y no hay duplicados
```

Los bancos encontraron **diez bugs** que ninguna lectura del código había visto,
y `check_docs.py` dos más.

**Un banco que pasa no prueba nada.** Prueba algo el día que se rompe lo que
cubre y el banco se queja, y la única forma de saberlo es romper el código a
propósito: eso es `tests/mutar.py`, 117 mutaciones que tienen que dar todas
`[OK]`. Un `[HUECO]` es un chequeo que falta; un `[VIEJA]` es una mutación cuyo
patrón ya no existe, que tampoco prueba nada. Así aparecieron los dos huecos que
ninguna otra cosa mostró — el veredicto de `perfil` que iba a la tabla sin estar
cubierto, y `estado` sin banco.

**Y una medición vale más que un umbral.** `align.sh verificar` daba `ok` a un
BAM con 67% de reads sin alinear, porque el único umbral miraba otra cosa. El
banco no lo iba a encontrar —el umbral hacía exactamente lo que decía hacer— y
la mutación tampoco. Lo encontró mirar el número al lado del veredicto, que es
la misma lección que la trampa del chequeo que mira el mensaje y no el dato,
un nivel más arriba: **el veredicto de una herramienta no reemplaza leer lo que
midió**, sobre todo la primera vez que se corre sobre datos reales.

`check_docs.py` pasa en verde: los 19 BioProjects de la spec están resueltos,
el manifiesto y el ledger reconcilian en las dos direcciones, y cada fichero que
los docs citan existe o está declarado como deuda.

## Entorno

micromamba, entorno `srna2` (se llama así porque el primer intento usó Python
3.11 y YASMA exige >= 3.12). Pin en `environment.yml` + `environment_pip.txt`.
YASMA se pinea por ref de git, no por número de versión — ver `docs/yasma.md`,
que documenta qué consume, qué decide y dónde no hay que confiarle.
`./check_env.sh` verifica que cada herramienta **arranque** de verdad, incluido
`import RNA` de ViennaRNA — el fallo clásico.
