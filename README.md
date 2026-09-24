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
notebooks/20_alinear.ipynb         recorte + alineamiento en Colab
notebooks/90_estado.ipynb          qué falta y cuánto ocupa
scripts/fetch_runs.sh     resuelve los 19 proyectos a corridas y las descarga
scripts/fetch_genomes.sh  resuelve y descarga los ensamblados
scripts/_drive_lib.sh     mapa de fases y ruta local — fuente única
scripts/trim.sh           recorte con yasma trim, por tandas, 18 proyectos
scripts/align.sh          alineamiento con yasma align (bowtie1 -v 1 -m 50)
data/adaptadores.tsv      qué adaptador recortar en cada BioProject
scripts/drive_push.sh     sube a Drive vía rclone
scripts/drive_pull.sh     baja de Drive, y purga la copia local
scripts/loo_cv.py         validación dejando un organismo afuera
scripts/validate_notebooks.py   chequea los .ipynb del repo
scripts/colab_git.py      commitea y empuja a GitHub desde Colab
scripts/check_docs.py     cruza lo que afirman los docs contra data/
tests/                    bancos de prueba, con binarios falsos en el PATH
tests/mutar.py            rompe el código y exige que algún banco grite
CLAUDE.md                 regla de ubicación + trampas conocidas
```

## Pipeline

```
prefetch → yasma trim (15-50 nt) → yasma align (bowtie1 -v 1 -m 50) → yasma tradeoff
                                                                            ↓
                                                    features por locus → PU learning
```

**El upstream, el recorte y el alineamiento están acá.** El alineamiento dejó de
depender del `orchestrate.sh` que nunca se subió: lo hace `yasma align`, que es
bowtie1 nativo con `--max_multi` en 50 por defecto — o sea el mismo `-m 50` que
el proyecto tenía medido. `config.sh`, `orchestrate.sh`, `verify.sh`,
`check_env.sh` y `environment.yml` siguen en el `main` local sin subir, y
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

./scripts/drive_pull.sh genomas galga --go
./scripts/align.sh genoma galga          # sha256 contra el ledger, y descomprime
./scripts/align.sh plan galga            # qué alinearía, contra qué ensamblado
./scripts/align.sh correr galga          # yasma align: bowtie1 -v 1 -m 50
./scripts/align.sh verificar galga       # ¿qué fracción alineó? ¿están los @RG?
./scripts/drive_push.sh bam galga --go   # los BAMs a Drive
./scripts/align.sh ledger                # data/alineamientos.tsv, para commitear
```

O en **Colab**, con `notebooks/20_alinear.ipynb`, para los 16 proyectos que
entran en el disco de una VM. §1 mide cuáles antes de empezar: el pico es
`recortado + 2 × BAM` y un proyecto no se puede partir sin cambiar el resultado.
Desde ahí Colab también **empuja a GitHub** (`scripts/colab_git.py`, con un PAT
en los Secrets), así que el manifiesto, los checksums y `alineamientos.tsv` ya no
se copian y pegan a mano. Ver `docs/colab.md`.

El recorte produce **18 proyectos YASMA**, uno por organismo y rol
(`proyectos/galga_primario/`, `proyectos/galga_duplicado/`): el duplicado es la validación
independiente y no puede compartir anotación con el primario.

Los dos `verificar` no son opcionales, y por el mismo motivo. cutadapt corre con
`--trimmed-only`, así que recortar con la secuencia equivocada **no da error**:
deja el `.t.fq.gz` casi vacío y `estado` sigue diciendo "0 faltan". Y alinear
contra el genoma equivocado tampoco: medido, `yasma align` sale con código 0 y
0% alineado. Lo único que delata cada caso es un número —la retención contra
`data/adaptadores.tsv`, la fracción alineada contra el sentido común— y eso es
lo que los dos `verificar` miran.

El BAM queda en `proyectos/<org>_<rol>/align/alignment.bam`, donde `tradeoff` lo
espera, y enlazado (hard link, no copia) a `bams/<org>/<rol>.bam`, que es de
donde `drive_push.sh` lo sube.

El manifiesto y los ensamblados **no se editan a mano**: se regeneran. Sacar una
corrida va en `data/excluidas.tsv`, con el motivo medido; editar
`data/srr_manifest.tsv` a mano lo deshace la próxima regeneración, en silencio.

Falta, y va acá cuando llegue: `./check_env.sh`, `./orchestrate.sh` (las 4
fases, idempotente — relanzar tras un corte no rehace trabajo porque el
alineamiento salta las corridas cuyo BAM ya existe) y `./verify.sh`.
