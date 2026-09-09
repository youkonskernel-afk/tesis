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
| Lecturas crudas | no se respaldan: `prefetch` las baja de SRA |

La regla completa y las trampas del pipeline están en [`CLAUDE.md`](CLAUDE.md).

## Estructura

```
data/organismos.tsv       9 organismos × (primario + duplicado) — spec del dataset
data/genomas.tsv          ensamblado por organismo, con estado de verificación
data/srr_manifest_r1.tsv  169 corridas de la ronda anterior (referencia)
data/DRIVE.md             índice de Drive: qué hay, dónde, con qué ID
docs/positivos.md         cómo se arma el conjunto positivo y por qué
docs/colab.md             bajar a Drive sin pasar por el disco local
notebooks/                cuaderno de Colab: NCBI -> Drive directo
scripts/fetch_runs.sh     resuelve los 18 proyectos a corridas y las descarga
scripts/fetch_genomes.sh  resuelve y descarga los ensamblados
scripts/drive_push.sh     sube resultados a Drive vía rclone
scripts/loo_cv.py         validación dejando un organismo afuera
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
