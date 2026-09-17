# Índice de Drive

Fuente de verdad de dónde vive la data. Si una ruta cambia, se actualiza **acá**
y en `scripts/drive_push.sh` (que direcciona por ruta, no por ID — los IDs de
esta tabla son para abrir la carpeta a mano y para las herramientas de Drive).

- **Cuenta**: `seb.ugazm@gmail.com` — al menos 1.6 TB de capacidad
- **Raíz**: `Mi unidad/tesis/`
- **ID raíz**: `1E_Q6XLg4_RD01TFgfHbGUSQxx2ExtuVP`
- https://drive.google.com/drive/folders/1E_Q6XLg4_RD01TFgfHbGUSQxx2ExtuVP

## Presupuesto de espacio

| | |
| :-- | --: |
| Capacidad de la cuenta | ≥ 1.6 TB |
| BAMs (~425 corridas) | ~340 GB |
| `.sra` crudos (**se respaldan**, `80_sra/`) | ~190 GB |
| Genomas + índices bowtie (9 organismos) | orden de 5 GB |
| YASMA, QC, features, modelos | orden de GB |

El espacio **no es la restricción**: BAMs + `.sra` + genomas ≈ 540 GB, un ~34%
de la cuenta. Lo que se excluye de abajo se excluye por criterio de reproducibilidad
—no se respalda lo que una base pública ya garantiza— no por falta de disco.

| Carpeta | ID | Contiene |
| :-- | :-- | :-- |
| `00_manifiestos/` | `1rZ6AQ7xu8D7F9WqR7qBfgzATWG_FrJWY` | Snapshots de `organismos.tsv`, `config.sh`, `environment.yml` por corrida |
| `10_bam/` | `1l_65g9VbWDqHik9UrBjgJCYNr7aurDB_` | `<org>/<RUN>.bam` + `.bai` — ~425 corridas, ~340 GB |
| `20_yasma/` | `14MGZizfNGREiX4CvHzZOPEnlwoTSu62s` | `<org>/annotations/<nombre>/loci.gff3`, counts |
| `30_qc/` | `11NXk3e0UFTw9RCCPMm-IpkqMsMVWnbsQ` | `<org>/` fastp `.json`/`.html`, flagstat, distribución de longitudes |
| `40_features/` | `1fHV2-QM6yvEW3yCCuyA8uWwvkEAeZS6x` | Matrices de features por locus para PU learning |
| `50_modelos/` | `1q59gopI-uKOGJhvQ_JbDJFx4qJL1wKO7` | Modelos entrenados, checkpoints, predicciones |
| `60_figuras/` | `13rOdpXqd_EbNi3RDamLLeIBUko-Zc7BC` | Figuras raster pesadas |
| `70_genomas/` | `1Ik2edbmuH6OLqjQDIOH6zLuUgYe_e4yy` | `<org>/<accession>.fna.gz` + `.sha256` — el FASTA exacto usado |
| `80_sra/` | `1RQYmPeaixM_-IXTzPeBCMg61xBkZlzL_` | `<org>/<RUN>.sra` — crudo validado, ~190 GB |

Dos corridas de `maggi` (`SRR317135`, `SRR1066790`) son **SRA Lite**: se
guardan como `<RUN>.sra` igual que las demás, así que el archivo no lo dice.
La columna `formato` de `data/sra_md5.tsv` es la que lo dice. Importa porque
las calidades de SRA Lite son sintéticas — ver `CLAUDE.md`.

Las subcarpetas por organismo (`rhirr`, `sclsc`, `cloro`, `phypa`, `prupe`,
`maldo`, `gadmo`, `galga`, `maggi`) no se crean a mano: `rclone copy` las crea
al subir.

## Qué NO está acá

- **Índices bowtie**. Se reconstruyen del FASTA en minutos.

Los **`.sra` pasaron a respaldarse** en `80_sra/`, revirtiendo la decisión
anterior. El motivo es el reparto con Colab: Colab los baja y una sesión que se
muere no puede perder la descarga, y el alineamiento local se desacopla de ella.
El caché local en `/home/dev/sra_cache` ya trae 57 corridas vigentes de R1
(`rhirr`, `sclsc` y `phypa` conservan su BioProject primario).

Los **genomas sí están** (`70_genomas/`), como excepción deliberada: un
ensamblado puede retirarse o reemplazarse, y sin el FASTA exacto el
alineamiento deja de ser reproducible. **Los 9 están fijados, bajados y
verificados**, con su `sha256` en `data/genomas.sha256`. Y la excepción se ganó
el lugar: **tres de los nueve tenían el ensamblado ya retirado por NCBI**
(`suppressed`), uno de ellos en uso desde R1. Ver `data/genomas.tsv` para el
porqué de cada elección, y `./scripts/fetch_genomes.sh verificar` para
comprobar que lo que hay en Drive es lo que dice el ledger.

Reconstruir todo desde cero necesita solo tres ficheros, y los tres están en
git: `config.sh`, `organismos.tsv` y `environment.yml`.

## Data heredada, fuera de este árbol

En la misma cuenta hay un árbol viejo del pipeline anterior (ShortStack:
`strucVis/`, `pudata/`, `puutils/`, y una carpeta por SRR bajo el id
`1-ocXZVs_z9Orfbk1IyIaheQfgX7HNW4v`). No se migró: `tesis/` arranca limpio.
Si algo de ahí resulta necesario, se copia explícitamente y se anota en esta tabla.
