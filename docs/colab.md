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
| recorte y alineamiento | Colab — `20_alinear.ipynb`, **si el proyecto entra** |
| traer `.sra` para alinear | local — `scripts/drive_pull.sh sra <org> --go` |
| recorte y alineamiento | local — `scripts/trim.sh` + `scripts/align.sh` |
| BAMs a Drive | local — `scripts/drive_push.sh bam <org> --go` |

El alineamiento está en los dos lados, y no es indecisión: hay proyectos que no
entran en una VM de Colab. Ver **Qué entra en Colab** más abajo.

## Cómo se usa

Links directos — Colab abre el notebook desde GitHub sin descargar nada:

| notebook | link |
| :-- | :-- |
| 00 setup | https://colab.research.google.com/github/youkonskernel-afk/tesis/blob/claude/github-google-drive-setup-cwapri/notebooks/00_setup.ipynb |
| genomas | https://colab.research.google.com/github/youkonskernel-afk/tesis/blob/claude/github-google-drive-setup-cwapri/notebooks/descarga_genomas.ipynb |
| 10 corridas | https://colab.research.google.com/github/youkonskernel-afk/tesis/blob/claude/github-google-drive-setup-cwapri/notebooks/10_descarga_runs.ipynb |
| 20 alinear | https://colab.research.google.com/github/youkonskernel-afk/tesis/blob/claude/github-google-drive-setup-cwapri/notebooks/20_alinear.ipynb |
| 90 estado | https://colab.research.google.com/github/youkonskernel-afk/tesis/blob/claude/github-google-drive-setup-cwapri/notebooks/90_estado.ipynb |

Ojo: la rama está en la URL. Cuando el default pase a `main` y esta rama se
mergee, hay que actualizar estos links.

1. Abrilos con los links de arriba, o desde https://colab.research.google.com
   (Archivo → Abrir cuaderno → GitHub).
2. **Siempre `00_setup.ipynb` primero.** Colab arranca de cero en cada sesión.
3. Montá con la cuenta `seb.ugazm@gmail.com`, que es la que tiene el árbol.
4. En `10_descarga_runs.ipynb`, ajustá `LIMITE` y corré la celda de descarga
   tantas veces como aguante la sesión.

## Qué entra en Colab, y qué no

**Un proyecto no se puede partir.** `yasma align` acumula la cobertura única de
*todas* las librerías del proyecto antes de pesar los reads multimapeados
(`unique_d` en `nativealign.py`), así que partirlo en varias llamadas cambia a
qué locus va cada uno. La unidad reanudable es el proyecto entero: si la sesión
se muere a la mitad de uno, ese se rehace; los que ya terminaron están en Drive.

El pico de disco es `recortado + 2 × BAM` — `pysam.sort` escribe el BAM ordenado
**antes** de borrar el sin ordenar, así que los dos conviven.

Los bytes por read están **medidos** en `sclsc_duplicado` (32.1 M reads): 21.8
en `.t.fq.gz` y 14.1 en BAM. Las estimaciones que había antes acá decían 25 y
**45**, o sea que todos los picos estaban inflados 3×. §1 usa 22 y 16 —el 16 en
vez del 14 porque se midió en una librería con 67% de reads sin alinear, y un
read sin alinear ocupa menos que uno colocado.

Con eso, los picos por proyecto quedan así:

| proyecto | runs | trim M | pico GB | horas |
| :-- | --: | --: | --: | --: |
| galga_duplicado | 95 | 1510 | **79** | ~23 |
| maldo_primario | 36 | 661 | 34 | ~10 |
| cloro_duplicado | 23 | 600 | 31 | ~9 |
| cloro_primario | 34 | 530 | 28 | ~8 |
| galga_primario | 27 | 505 | 26 | ~8 |
| phypa_primario | 30 | 444 | 23 | ~7 |
| … los otros 12 | | ≤433 | ≤23 | ≤7 |
| **total** | **417** | **6068** | | **~90** |

Cuánto disco da una VM **varía**: la primera corrida real dio **220 GB libres**
en Colab Free, no los ~78 que se daban por sentados. Con 220 entran los 18; con
78 entran 17 y se queda afuera `galga_duplicado`. Por eso el número no se
declara acá: **§1 de `20_alinear.ipynb` lo mide contra el disco de esa VM antes
de empezar** y ordena del más chico al más grande. Eso se sabe en un segundo o
a las seis horas.

Las horas salen de los **49 s por millón de reads** medidos en el mismo
proyecto (22:28 de bowtie + 3:35 de `pysam.sort`, genoma de 39 Mb, Colab Free).
§1 usa 55 de margen porque los genomas grandes —`galga` 1.1 Gb, `maggi`
650 Mb— son más lentos por read que un hongo. §4 se cronometra solo y §5
imprime los s/M read reales de cada proyecto.

Los `.sra` **no se copian** a la VM: Drive está montado, así que `SRA_DEST`
apunta al mount. La regla de abajo es no *escribir* archivos grandes al FUSE;
leerlos está bien, y ahorra ~190 GB de copia.

## De vuelta a GitHub

Colab ya no solo lee. `scripts/colab_git.py` commitea y empuja, y reemplazó los
bloques de "copiá esta salida al repo" de §4 de `10_descarga_runs` y §6 de
`descarga_genomas`.

Hace falta un **PAT con permiso de escritura** guardado como `GITHUB_TOKEN` en
los Secrets de Colab (la llave a la izquierda), habilitado para el notebook. El
token vive ahí: no va al repo ni a Drive.

Tres cosas que el módulo garantiza, y que tienen banco:

- **Nunca `git add -A`.** Solo las rutas que se le pasan, y solo bajo `data/`.
  Esta VM tiene Drive montado en `/content/drive`: un `add -A` es exactamente
  donde se cuela un `.sra`.
- **El token no aparece en ningún mensaje**, ni cuando git falla — git mete la
  URL, con el token adentro, en sus errores.
- **Un push rechazado falla fuerte.** Se reintenta una vez rebasando sobre lo
  que haya (dos sesiones de Colab sobre el mismo ledger divergen) y si vuelve a
  fallar, revienta. Un push rechazado que nadie mira deja el resultado en Drive
  y no en git.

Las celdas que empujan arrancan con `REVISAR_PRIMERO = True`: imprimen el diff y
**no** empujan. Es el mismo criterio que `--go` en `drive_push.sh`.

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
los notebooks imprimen lo que falta commitear para que lo copies.

**Pero el ledger de trabajo vive en Drive**, en `00_manifiestos/sra_md5.tsv`, no
adentro del clon. No es una excepción a lo anterior: git sigue siendo el
registro que vale. Son dos cosas distintas por dos motivos concretos.

El primero es que el clon **se resetea en cada corrida** (`git fetch` +
`git reset --hard`), así que nada escrito adentro sobrevive. El segundo es que
antes el ledger era un archivo **versionado** que el script modificaba, y eso
hacía que `git pull --ff-only` fallara para siempre después de la primera
descarga — sin hacer ruido, dejando al notebook corriendo con el código viejo.
Ese fallo mudo costó una ronda entera.

Que el ledger viva en Drive además lo salva de que se muera la sesión, que es
exactamente lo que pasó una vez y obligó a recalcular md5. La celda de setup lo
siembra desde la copia del repo cuando esta tiene más filas, así que el modo
`ledger` no recalcula las 416 corridas al pedo.

## Los notebooks son delgados

Clonan el repo y **llaman a los scripts que ya existen**. El criterio de
selección de corridas (`library_source = TRANSCRIPTOMIC`, `SINGLE` para
RNA-Seq) vive en `scripts/fetch_runs.sh`, y la verificación de ensamblados en
`scripts/fetch_genomes.sh`. Duplicarlos en Python sería tener dos criterios de
selección: un riesgo de reproducibilidad, no una comodidad.

### El clon funciona con el repo público o privado

La celda de clon intenta primero un `git clone` anónimo. Si el repo es público,
listo. Si falla, busca un token en los Secrets de Colab y reintenta. Así la
visibilidad del repo es una decisión tuya y no un bloqueo.

**Si querés dejar el repo privado**, cargá el token una vez:

1. En GitHub: Settings → Developer settings → Personal access tokens →
   Fine-grained tokens. Alcance mínimo: este repo, permiso *Contents: Read-only*.
2. En Colab: el ícono de llave en la barra izquierda → añadir secreto con
   nombre exacto `GITHUB_TOKEN`, y habilitar el acceso para el notebook.

Dos cuidados que la celda ya resuelve, y conviene no deshacer si alguien la
edita: **no se imprime el stderr de git en la rama con token** —git incluye la
URL en sus errores, y esa URL lleva el token— y **se resetea el remoto** a la
URL sin token después de clonar, para que no quede escrito en `.git/config`.

**Los notebooks clonan la rama por defecto del repo.** Hoy esa rama es
`claude/github-google-drive-setup-cwapri`, porque GitHub tomó como default la
primera que se empujó al repo vacío. Funciona, pero en cuanto exista `main` con
la historia mergeada hay que cambiar el default a `main` en Settings → General,
o los notebooks van a seguir clonando una rama que quedó atrás — y sin avisar,
que es lo peor de este tipo de fallas.

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
