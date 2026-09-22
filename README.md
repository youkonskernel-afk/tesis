# tesis — descubrimiento de sRNAs no descritos con PU learning

Anotación de loci de sRNA en ~425 corridas sRNA-seq de 9 organismos —3 hongos,
3 plantas, 3 animales— y un clasificador **positive-unlabeled** para priorizar
loci que no corresponden a ningún sRNA descrito.

Cada organismo se descarga por duplicado: un BioProject **primario** para
descubrir y uno **independiente** para validar. Un candidato que aparece en
ambos vale mucho más que uno que solo existe en un experimento.

El planteamiento PU es deliberado: miRBase y Rfam dan positivos, pero lo **no
anotado no es negativo** — es justamente donde puede estar lo nuevo. Tratarlo
como negativo entrenaría al modelo en contra del objetivo.

## Dónde está cada cosa

| | |
| :-- | :-- |
| Código | este repo |
| Data | `Mi unidad/tesis/` en Drive — ver [`data/DRIVE.md`](data/DRIVE.md) |
| Descarga | Google Colab, que escribe directo a Drive — ver [`docs/colab.md`](docs/colab.md) |
| Cómputo largo | máquina local: alineamiento y YASMA |

La regla completa y las trampas del pipeline están en [`CLAUDE.md`](CLAUDE.md).

## Estructura

```
data/organismos.tsv       9 organismos × (primario + duplicado) — spec del dataset
data/genomas.tsv          ensamblado por organismo, con estado de verificación
data/srr_manifest_r1.tsv  169 corridas de la ronda anterior (referencia)
data/DRIVE.md             índice de Drive: qué hay, dónde, con qué ID
docs/positivos.md         cómo se arma el conjunto positivo y por qué
docs/yasma.md             qué consume YASMA, qué decide, y cómo pinearlo
docs/rclone.md            configurar rclone para Drive, paso a paso
scripts/drive_check.sh    verifica esa configuración
docs/colab.md             Colab como administrador de datos
docs/plan_datos_colab.md  el plan que implementa lo anterior
notebooks/00_setup.ipynb  monta Drive, clona, instala, verifica
notebooks/descarga_genomas.ipynb   verifica y baja ensamblados a Drive
notebooks/10_descarga_runs.ipynb   manifiesto + .sra a Drive, por tandas
notebooks/90_estado.ipynb          qué falta y cuánto ocupa
scripts/fetch_runs.sh     resuelve los 19 proyectos a corridas y las descarga
scripts/fetch_genomes.sh  resuelve y descarga los ensamblados
scripts/_drive_lib.sh     mapa de fases y ruta local — fuente única
scripts/trim.sh           recorte con yasma trim, por tandas, 18 proyectos
data/adaptadores.tsv      qué adaptador recortar en cada BioProject
scripts/drive_push.sh     sube a Drive vía rclone
scripts/drive_pull.sh     baja de Drive, y purga la copia local
scripts/loo_cv.py         validación dejando un organismo afuera
scripts/validate_notebooks.py   chequea los .ipynb del repo
scripts/check_docs.py     cruza lo que afirman los docs contra data/
tests/                    bancos de prueba, con binarios falsos en el PATH
tests/mutar.py            rompe el código y exige que algún banco grite
CLAUDE.md                 regla de ubicación + trampas conocidas
```

## Pipeline

```
prefetch → yasma trim (15-50 nt) → bowtie1 (-m 50) → samtools → yasma tradeoff
                                                              ↓
                                              features por locus → PU learning
```

**El upstream está acá; el alineamiento todavía no.** Lo que este repo corre hoy
es la mitad de arriba: resolver el manifiesto, bajar los `.sra` y los genomas, y
verificarlos. `config.sh`, `orchestrate.sh`, `verify.sh`, `check_env.sh` y
`environment.yml` están en el `main` local sin subir, así que los comandos que
los usan **no se pueden correr desde un clon de este repo**. Mientras no lleguen,
`scripts/check_docs.py` los lista como deuda declarada en vez de dejar que el
README los cite como si funcionaran.

## Reproducir desde cero

Lo que hoy se puede reproducir con lo que hay versionado:

```bash
./scripts/check_docs.py              # docs contra data/: ¿coincide todo?
./scripts/validate_notebooks.py      # los .ipynb parsean
tests/run_all.sh                     # los bancos, sin red

# desde Colab, que sí tiene salida a la ENA y a NCBI:
./scripts/fetch_runs.sh manifest     # los 19 proyectos -> corridas
./scripts/fetch_runs.sh perfil --proyectos   # ¿son sRNA-seq de verdad?
./scripts/fetch_runs.sh prefetch     # los .sra a Drive, por tandas
./scripts/fetch_genomes.sh fetch     # los ensamblados a Drive
./scripts/fetch_genomes.sh verificar # recalcula los sha256 contra el ledger
```

Y en la máquina local, una vez que rclone esté configurado (`docs/rclone.md`)
y los `.sra` bajados:

```bash
./scripts/drive_pull.sh sra galga --go   # traer un organismo
./scripts/trim.sh plan galga             # qué recortaría, y cuánto disco
./scripts/trim.sh correr galga           # recorta, en tandas
./scripts/trim.sh verificar galga        # ¿la retención da lo que perfil predijo?
```

El recorte produce **18 proyectos YASMA**, uno por organismo y rol
(`trim/galga_primario/`, `trim/galga_duplicado/`): el duplicado es la validación
independiente y no puede compartir anotación con el primario.

`verificar` no es opcional. cutadapt corre con `--trimmed-only`, así que recortar
con la secuencia equivocada **no da error**: deja el `.t.fq.gz` casi vacío y
`estado` sigue diciendo "0 faltan". Lo único que lo delata es la retención medida
contra la `retencion_est` de `data/adaptadores.tsv`.

El manifiesto y los ensamblados **no se editan a mano**: se regeneran. Sacar una
corrida va en `data/excluidas.tsv`, con el motivo medido; editar
`data/srr_manifest.tsv` a mano lo deshace la próxima regeneración, en silencio.

Falta, y va acá cuando llegue: `./check_env.sh`, `./orchestrate.sh` (las 4
fases, idempotente — relanzar tras un corte no rehace trabajo porque el
alineamiento salta las corridas cuyo BAM ya existe) y `./verify.sh`.
