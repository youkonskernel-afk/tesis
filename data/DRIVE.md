# Índice de Drive

Fuente de verdad de dónde vive la data. Si una ruta cambia, se actualiza **acá**
y en `scripts/drive_push.sh` (que direcciona por ruta, no por ID — los IDs de
esta tabla son para abrir la carpeta a mano y para las herramientas de Drive).

- **Cuenta**: `seb.ugazm@gmail.com`
- **Raíz**: `Mi unidad/tesis/`
- **ID raíz**: `1E_Q6XLg4_RD01TFgfHbGUSQxx2ExtuVP`
- https://drive.google.com/drive/folders/1E_Q6XLg4_RD01TFgfHbGUSQxx2ExtuVP

| Carpeta | ID | Contiene |
| :-- | :-- | :-- |
| `00_manifiestos/` | `1rZ6AQ7xu8D7F9WqR7qBfgzATWG_FrJWY` | Copia de respaldo de `srr_manifest.tsv`, `config.sh`, `environment.yml` |
| `10_bam/` | `1l_65g9VbWDqHik9UrBjgJCYNr7aurDB_` | `<org>/<RUN>.bam` + `.bai` — 169 corridas, ~100 GB |
| `20_yasma/` | `14MGZizfNGREiX4CvHzZOPEnlwoTSu62s` | `<org>/annotations/<nombre>/loci.gff3`, counts |
| `30_qc/` | `11NXk3e0UFTw9RCCPMm-IpkqMsMVWnbsQ` | `<org>/` fastp `.json`/`.html`, flagstat, distribución de longitudes |
| `40_features/` | `1fHV2-QM6yvEW3yCCuyA8uWwvkEAeZS6x` | Matrices de features por locus para PU learning |
| `50_modelos/` | `1q59gopI-uKOGJhvQ_JbDJFx4qJL1wKO7` | Modelos entrenados, checkpoints, predicciones |
| `60_figuras/` | `13rOdpXqd_EbNi3RDamLLeIBUko-Zc7BC` | Figuras raster pesadas |

Las subcarpetas por organismo (`rhirr`, `sclsc`, `arath`, `phypa`, `danre`,
`nemve`) no se crean a mano: `rclone copy` las crea al subir.

## Qué NO está acá

- **`.sra` crudos** (55 GB). Son públicos; `prefetch` los recupera desde SRA.
  El caché local vive en `/home/dev/sra_cache`.
- **Genomas e índices bowtie** (~2.5 GB). Se re-descargan desde Ensembl,
  EnsemblGenomes release-57 y NCBI; las URLs están en `config.sh`.

Reconstruir todo desde cero necesita solo tres ficheros, y los tres están en
git: `config.sh`, `srr_manifest.tsv` y `environment.yml`.

## Data heredada, fuera de este árbol

En la misma cuenta hay un árbol viejo del pipeline anterior (ShortStack:
`strucVis/`, `pudata/`, `puutils/`, y una carpeta por SRR bajo el id
`1-ocXZVs_z9Orfbk1IyIaheQfgX7HNW4v`). No se migró: `tesis/` arranca limpio.
Si algo de ahí resulta necesario, se copia explícitamente y se anota en esta tabla.
