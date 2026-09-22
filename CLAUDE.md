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
2. `yasma trim` (cutadapt, 15-50 nt) → `bowtie1` estilo ShortStack3 → `samtools`
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
    `PRE-TRIMMED`, así que YASMA la pasa de largo y **no le aplica ningún
    filtro**, ni el de longitud. Es la única de las 19 que entra al alineamiento
    sin pasar por la ventana 15-50.
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
  como `<RUN>.sra` igual que los demás. **Pendiente de confirmar contra
  `config.sh`**: si `bowtie` corre en modo `-v` (conteo de mismatches, el
  estilo ShortStack3) la calidad se ignora y el asunto se agota en `fastp`; si
  corriera `-n`/`-e` (suma de calidades del seed), las calidades sintéticas
  cambiarían el alineamiento de esas 2 corridas.
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
| configurar rclone | máquina local | `docs/rclone.md` + `scripts/drive_check.sh` |
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
./tests/run_all.sh              # 15 bancos, 366 chequeos, binarios falsos en el PATH
./tests/mutar.py                # rompe el codigo y exige que algun banco grite
./scripts/check_docs.py         # lo que afirman los docs contra data/
./scripts/validate_notebooks.py # los .ipynb parsean y no hay duplicados
```

Los bancos encontraron **diez bugs** que ninguna lectura del código había visto,
y `check_docs.py` dos más.

**Un banco que pasa no prueba nada.** Prueba algo el día que se rompe lo que
cubre y el banco se queja, y la única forma de saberlo es romper el código a
propósito: eso es `tests/mutar.py`, 46 mutaciones que tienen que dar todas
`[OK]`. Un `[HUECO]` es un chequeo que falta; un `[VIEJA]` es una mutación cuyo
patrón ya no existe, que tampoco prueba nada. Así aparecieron los dos huecos que
ninguna otra cosa mostró — el veredicto de `perfil` que iba a la tabla sin estar
cubierto, y `estado` sin banco.

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
