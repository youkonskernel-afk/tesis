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
docs/colab.md             Colab como administrador de datos
docs/plan_datos_colab.md  el plan que implementa lo anterior
notebooks/00_setup.ipynb  monta Drive, clona, instala, verifica
notebooks/descarga_genomas.ipynb   verifica y baja ensamblados a Drive
notebooks/10_descarga_runs.ipynb   manifiesto + .sra a Drive, por tandas
notebooks/90_estado.ipynb          qué falta y cuánto ocupa
scripts/fetch_runs.sh     resuelve los 18 proyectos a corridas y las descarga
scripts/fetch_genomes.sh  resuelve y descarga los ensamblados
scripts/drive_push.sh     sube a Drive vía rclone
scripts/drive_pull.sh     baja de Drive, y purga la copia local
scripts/loo_cv.py         validación dejando un organismo afuera
scripts/validate_notebooks.py   chequea los .ipynb del repo
CLAUDE.md                 regla de ubicación + trampas conocidas
```

## Pipeline

```
prefetch → fastp (15-50 nt) → bowtie1 (-m 50) → samtools → YASMA v1.1.1
                                                              ↓
                                              features por locus → PU learning
```

Entrada: `./orchestrate.sh`. Estado: `./orchestrate.sh status`. Verificación:
`./verify.sh`. Es idempotente: relanzar tras un corte no rehace trabajo, porque
el paso de alineamiento salta las corridas cuyo BAM ya existe.

## Reproducir desde cero

Solo hacen falta tres ficheros, los tres versionados acá: `config.sh`,
`organismos.tsv` y `environment.yml`.

```bash
./check_env.sh      # línea base
# instalar micromamba + entorno srna2 (ver CLAUDE.md)
./check_env.sh      # verificar
./gen_manifest.sh   # regenerar el manifiesto desde la ENA Portal API
./orchestrate.sh    # lanza las 4 fases
```
