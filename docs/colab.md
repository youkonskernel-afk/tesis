# Descargar a Drive sin pasar por el disco local

El problema: la sesión cloud de Claude tiene bloqueada por política la salida a
NCBI, Ensembl, la ENA, miRBase y MirGeneDB. Solo pasan repositorios de paquetes
y GitHub. Así que no puede bajar ni un genoma ni una corrida.

La salida es **Google Colab**: corre en infraestructura de Google, tiene salida
a NCBI, y monta Drive como sistema de archivos. El FASTA va de NCBI a
`Mi unidad/tesis/70_genomas/` sin tocar tu disco.

## Usarlo

1. Abrí `notebooks/descarga_genomas.ipynb` en https://colab.research.google.com
   (Archivo → Subir cuaderno, o desde GitHub si el repo es accesible).
2. Corré las celdas en orden. La celda 1 pide permiso para montar Drive —
   tiene que ser la cuenta `seb.ugazm@gmail.com`, la que tiene el árbol.
3. **Mirá la salida de la celda de verificación antes de seguir.** Contrasta
   cada candidato contra el ensamblado de referencia vigente que declara NCBI.
   Por defecto baja el vigente, no el candidato: si difieren, gana NCBI.
4. La última celda imprime las líneas para pegar en `data/genomas.tsv` y
   `data/genomas.sha256`, y commitearlas.

El notebook **se genera**, no se edita a mano:

    ./scripts/gen_colab_notebook.py          # regenerar desde data/genomas.tsv
    ./scripts/gen_colab_notebook.py --check  # falla si quedó desactualizado

Así la tabla de candidatos no se desincroniza de `data/genomas.tsv`.

## Por qué el paso de verificación vive acá

`scripts/fetch_genomes.sh resolve` hace lo mismo, pero necesita red hacia NCBI,
o sea tu máquina. En Colab hay red, así que verificar y descargar ocurren en el
mismo lugar y en el mismo orden. Es la única parte del proyecto donde eso se
puede hacer de corrido.

## Lo que Colab NO resuelve

- **Los BAMs.** Son del orden de 340 GB y las sesiones de Colab son efímeras y
  con límite de tiempo. Van por `rclone` desde tu máquina:
  `./scripts/drive_push.sh bam <org> --go`.
- **El alineamiento.** 30-40 h de bowtie no entran en una sesión de Colab.
- **Los `.sra`.** Se bajan con `prefetch` en la máquina donde se va a alinear;
  copiarlos a Drive para después volver a bajarlos es trabajo al pedo.

## Alternativas, por si Colab no sirve

- **Google Cloud Shell**: también hosteado por Google, con `rclone` instalable.
  Más shell y menos notebook, pero solo 5 GB de home persistente.
- **`rclone copyurl <url> gdrive-tesis:tesis/70_genomas/...`**: transmite de la
  URL al remoto sin escribir en disco local. El ancho de banda igual pasa por
  tu máquina, pero no gasta disco.
