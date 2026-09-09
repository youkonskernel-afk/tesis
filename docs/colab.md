# Colab como administrador de datos

La sesión cloud de Claude tiene bloqueada por política la salida a NCBI,
Ensembl, la ENA, miRBase y MirGeneDB: el gateway responde 403 al `CONNECT`.
Solo pasan repositorios de paquetes y GitHub. No puede bajar ni un genoma ni
una corrida.

**Colab es el puente.** Corre en infraestructura de Google, tiene salida a NCBI
y a la ENA, y monta Drive como sistema de archivos: el dato va de la fuente a
`Mi unidad/tesis/` sin tocar tu disco.

## Reparto

| tarea | dónde |
| :-- | :-- |
| setup de la sesión | Colab — `00_setup.ipynb` |
| ensamblados | Colab — `descarga_genomas.ipynb` |
| manifiesto y `.sra` | Colab — `10_descarga_runs.ipynb` |
| ver qué falta | Colab — `90_estado.ipynb` |
| traer `.sra` para alinear | local — `scripts/drive_pull.sh sra <org> --go` |
| alineamiento y YASMA | local — `orchestrate.sh` |
| BAMs a Drive | local — `scripts/drive_push.sh bam <org> --go` |

Con Colab Free el cómputo largo se queda local: son ~2 vCPU y sesiones de hasta
12 h con desconexión por inactividad, y el alineamiento son 30-40 h. El patrón
de cola reanudable de `10_descarga_runs` es el mismo que necesitaría un
`20_alinear.ipynb`, así que moverlo si pasás a Pro es incremental.

## Cómo se usa

1. Abrí los notebooks en https://colab.research.google.com (desde GitHub, o
   Archivo → Subir cuaderno).
2. **Siempre `00_setup.ipynb` primero.** Colab arranca de cero en cada sesión.
3. Montá con la cuenta `seb.ugazm@gmail.com`, que es la que tiene el árbol.
4. En `10_descarga_runs.ipynb`, ajustá `LIMITE` y corré la celda de descarga
   tantas veces como aguante la sesión.

## Cuatro reglas que no son opcionales

Esto es lo que separa un notebook que funciona de uno que pierde datos.

**Nunca escribir archivos grandes directo al FUSE de Drive.** `prefetch`
escribe en el disco efímero de la VM (`SRA_STAGING`), se valida, y recién ahí se
mueve a Drive (`SRA_DEST`). Escribir GB a través del mount es lento e inestable.

**Validar antes de mover.** `vdb-validate` después de cada `prefetch`. Un `.sra`
truncado **no falla ruidosamente**: alinea de menos, y el error aparece recién
en los resultados. El que no valida se descarta y se reintenta; nunca llega a
Drive.

**Morir es lo normal.** El estado del trabajo es qué archivos existen en Drive,
no un fichero de progreso que se pueda corromper a media escritura. Cualquier
notebook se puede re-ejecutar y retoma.

**Los checksums van a git.** `data/sra_md5.tsv` y `data/genomas.sha256`. El
checksum guardado solo al lado del dato en Drive no prueba nada: quien reemplace
el archivo reemplaza el checksum con él. El clon de Colab es efímero, así que
los notebooks imprimen el ledger para que lo copies y commitees.

## Los notebooks son delgados

Clonan el repo y **llaman a los scripts que ya existen**. El criterio de
selección de corridas (`library_source = TRANSCRIPTOMIC`, `SINGLE` para
RNA-Seq) vive en `scripts/fetch_runs.sh`, y la verificación de ensamblados en
`scripts/fetch_genomes.sh`. Duplicarlos en Python sería tener dos criterios de
selección: un riesgo de reproducibilidad, no una comodidad.

Por eso el repo es público — así el clon no necesita credenciales.

`scripts/validate_notebooks.py` chequea que los `.ipynb` del repo tengan JSON
válido y que cada celda de código parsee.

## Límites, dichos de frente

- **Colab Free no está pensado para trabajo desatendido largo.** El uso
  sostenido lleva a throttling. Bajar ~190 GB va a tomar varias sesiones y
  conviene espaciarlas.
- **Drive tiene un tope de ~750 GB/día de subida.** 190 GB entra, pero no de una.
- **El FUSE de Drive es lento con muchos archivos chicos.** Acá son ~425
  archivos grandes, que es el caso bueno; aun así, no generar temporales sobre
  el mount.
- **La versión de sra-tools en Colab puede no ser la 3.4.1 que pinea
  `environment.yml`.** No es grave porque Colab solo descarga y valida, no
  produce resultados que entren en la tesis.

## Alternativas

- **Google Cloud Shell**: también hosteado por Google, con `rclone` instalable.
  Más shell y menos notebook, pero solo 5 GB de home persistente.
- **`rclone copyurl <url> gdrive-tesis:tesis/...`**: transmite de la URL al
  remoto sin escribir en disco local. El ancho de banda pasa por tu máquina.
