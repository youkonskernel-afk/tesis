# tesis — descubrimiento de sRNAs no descritos con PU learning

Anotación de loci de sRNA en 169 corridas sRNA-seq de 6 organismos (4 reinos), y
un clasificador **positive-unlabeled** para priorizar loci que no corresponden a
ningún sRNA descrito.

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
data/srr_manifest.tsv    169 corridas — fuente de verdad del dataset
data/DRIVE.md            índice de Drive: qué hay, dónde, con qué ID
scripts/drive_push.sh    sube resultados a Drive vía rclone
CLAUDE.md                regla de ubicación + trampas conocidas
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
`srr_manifest.tsv` y `environment.yml`.

```bash
./check_env.sh      # línea base
# instalar micromamba + entorno srna2 (ver CLAUDE.md)
./check_env.sh      # verificar
./gen_manifest.sh   # regenerar el manifiesto desde la ENA Portal API
./orchestrate.sh    # lanza las 4 fases
```
