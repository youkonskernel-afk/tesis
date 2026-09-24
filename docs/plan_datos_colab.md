# Plan: Colab como administrador de datos y puente de descarga

Estado: **implementado.** Queda como registro de por qué el diseño es así.

**El repo ya es público** (verificado por API: `"visibility":"public"`), así que
los notebooks clonan sin credenciales. Queda una sola cosa fuera de lo que Claude
puede hacer, en Settings de GitHub:

1. **Cambiar el default branch a `main`** una vez que exista con la historia
   mergeada. Hoy el default es `claude/github-google-drive-setup-cwapri` porque
   GitHub tomó la primera rama empujada al repo vacío, y los notebooks clonan la
   rama por defecto.

Y una que sí depende de vos pero no de Settings: para que Colab **empuje** hace
falta un PAT con permiso de escritura en los Secrets de la sesión
(`GITHUB_TOKEN`). Ver `docs/colab.md`.

## Contexto

La sesión cloud de Claude tiene bloqueada por política la salida a NCBI, Ensembl, la
ENA, miRBase y MirGeneDB: el gateway responde 403 al `CONNECT`. Todo lo que necesita
red quedó delegado a la máquina local, que pasa a ser cuello de botella para bajar
~190 GB de `.sra` y unos GB de genomas.

Colab resuelve eso: corre en infraestructura de Google, tiene salida a NCBI y a la
ENA, y monta Drive como sistema de archivos. El dato va de la fuente a
`Mi unidad/tesis/` **sin pasar por el disco local**. Ya hay prueba de concepto para
genomas (`notebooks/descarga_genomas.ipynb`, existe); esto lo generaliza a todo el
ingreso de datos.

### Decisiones tomadas

- **Alcance, con Colab Free**: Colab es el puente de descarga y el administrador de
  Drive. El alineamiento **también corre en Colab ahora** (`20_alinear.ipynb`),
  pero solo para los proyectos que entran en el disco de la VM: 16 de los 18.
  `maldo_primario` y `galga_duplicado` se quedan locales. No es por tiempo sino por disco — el pico
  es `recortado + 2 × BAM`, y un proyecto no se puede partir sin cambiar el
  resultado. Ver `docs/colab.md`.
- **Los `.sra` se guardan en Drive** (~190 GB). Rompe a propósito la regla de no
  respaldar lo público. Es lo que hace funcionar el reparto: una sesión de Colab que
  se muere no pierde la descarga, y el alineamiento local se desacopla de ella.
- **El repo pasa a público**, para que los notebooks clonen sin credenciales.
  Hecho.

## Arquitectura

| Pieza | Rol | Es |
| :-- | :-- | :-- |
| GitHub (público) | código, specs, manifiestos, checksums | fuente de verdad |
| Drive | genomas, `.sra`, BAMs, resultados | almacenamiento **y estado** |
| Colab | descarga, verificación, QC liviano | worker **efímero** |
| Máquina local | los proyectos grandes, YASMA | cómputo largo |

El estado del trabajo **es la existencia de los archivos en Drive**, no un archivo de
progreso aparte. Mismo criterio que ya usa el pipeline (03 saltea las corridas cuyo
BAM existe), y evita que un fichero de estado se corrompa cuando una sesión muere a
la mitad. Toda sesión de Colab es reanudable por construcción.

## Cambios en Drive

Carpeta nueva bajo `tesis/` (`1E_Q6XLg4_RD01TFgfHbGUSQxx2ExtuVP`):

```
80_sra/<org>/<RUN>.sra      ~190 GB — crudo descargado y validado
```

Presupuesto: ~340 GB de BAMs + ~190 GB de `.sra` + ~5 GB de genomas ≈ **540 GB de
1.6 TB (~34%)**. Hay que corregir `data/DRIVE.md`, que hoy lista los `.sra` como
"no está acá".

## El repo público permite borrar código

`scripts/gen_colab_notebook.py` embebe la tabla de candidatos dentro del `.ipynb`
porque el notebook no podía leer el repo. Con el repo público eso deja de hacer
falta: **los notebooks clonan y leen los TSV en vivo**, y el problema de
sincronización desaparece en vez de gestionarse. El generador y su `--check` se
retiran.

**Principio para todos los notebooks**: clonan el repo y **llaman a los scripts que
ya existen**, no reimplementan. El filtro de RNA (`library_source = TRANSCRIPTOMIC`,
`SINGLE` para RNA-Seq) vive en `scripts/fetch_runs.sh` y tiene que seguir viviendo en
un solo lugar — duplicarlo en Python sería crear dos criterios de selección, que es
un riesgo de reproducibilidad, no una comodidad.

## Notebooks

Preámbulo común: montar Drive, `git clone`, instalar herramientas, verificar.

| Notebook | Hace |
| :-- | :-- |
| `00_setup.ipynb` | Monta, clona, instala `sra-tools` + `curl`/`jq`, verifica que arranquen |
| `descarga_genomas.ipynb` | **Existe**; reescribir para leer `data/genomas.tsv` del clon |
| `10_descarga_runs.ipynb` | `fetch_runs.sh manifest`, después prefetch + validate + mover a `80_sra/` |
| `90_estado.ipynb` | Reconcilia Drive contra el manifiesto: qué falta, cuánto ocupa |

`20_alinear.ipynb` no entra ahora (Free). El patrón de cola reanudable de
`10_descarga_runs` es el mismo que necesitaría, así que agregarlo después es
incremental.

## Mecánica crítica

Esto separa un notebook que funciona de uno que pierde datos:

1. **Nunca escribir archivos grandes directo al FUSE de Drive.** `prefetch` al disco
   local de la VM, verificar, y recién ahí mover. Escribir GB a través del mount es
   lento e inestable.
2. **Validar antes de mover.** `vdb-validate` después de cada `prefetch`: las
   descargas truncadas de SRA son frecuentes y un `.sra` corto **no falla
   ruidosamente**, alinea de menos. El que no valida se descarta y se reintenta;
   nunca se sube.
3. **Guardia de disco.** Free da ~78 GB efímeros. Chequear espacio antes de cada
   descarga y borrar la copia local apenas se confirmó la de Drive.
4. **Morir es normal.** Cada notebook recorre la cola y sale limpio; re-ejecutar
   retoma. Ninguna celda asume que la anterior terminó.
5. **Orden por organismo.** Completar organismos enteros antes de empezar otros: un
   organismo completo se puede alinear, uno a medias no sirve para nada.
6. **Ledger de checksums en git**, `data/sra_md5.tsv`, mismo patrón que
   `data/genomas.sha256`. El checksum guardado solo al lado del dato en Drive no
   prueba nada.

## Archivos

**Nuevos**
- `notebooks/00_setup.ipynb`, `notebooks/10_descarga_runs.ipynb`, `notebooks/90_estado.ipynb`
- `scripts/drive_pull.sh` — espejo de `drive_push.sh`, que hoy solo sube. Trae los
  `.sra` de un organismo justo antes de alinearlo y permite borrarlos después, para
  que el disco local no tenga que aguantar los 190 GB de una vez.
- `data/sra_md5.tsv` — ledger, se llena al descargar

**Modificados**
- `notebooks/descarga_genomas.ipynb` — driver delgado que lee del clon
- `scripts/fetch_runs.sh` — `--limit N` (cuota por sesión) y raíz de Drive configurable
- `data/DRIVE.md` — carpeta `80_sra/`, presupuesto, sacar `.sra` de "no está acá"
- `CLAUDE.md` — regla de ubicación (los `.sra` ahora sí se respaldan y por qué) y
  tabla de la sección "Red"
- `docs/colab.md` — de "cómo bajar genomas" a la guía del rol completo de Colab
- `README.md`

**A retirar**
- `scripts/gen_colab_notebook.py` — su motivo desaparece con el repo público

## Verificación

1. **Notebooks bien formados**: chequeo que recorra `notebooks/*.ipynb`, valide el
   JSON y haga `ast.parse` de cada celda de código. Reemplaza lo que hoy hace el
   generador, y cubre todos los notebooks en vez de solo el generado.
2. **Filtro de RNA intacto**: repetir la prueba con `curl` simulado — de 8 corridas
   sintéticas tiene que retener 5 y descartar la GENOMIC, la RNA-Seq PAIRED y la de
   pocos reads.
3. **Idempotencia**: correr `10_descarga_runs` dos veces sobre el mismo organismo; la
   segunda no baja nada.
4. **Reconciliación**: `90_estado` da el mismo conteo por organismo que
   `data/srr_manifest.tsv`, y cero corridas en Drive que no estén en el manifiesto.
5. **`drive_pull.sh`**: dry-run sobre un organismo chico (`prupe`, 9+5 corridas)
   antes de correr con `--go`.
6. **Ida y vuelta real**: bajar `prupe` entero en Colab, traerlo con `drive_pull.sh`
   y verificar los md5 contra `data/sra_md5.tsv`.

## Riesgos

- **Colab Free no está pensado para trabajo desatendido largo.** El uso sostenido
  lleva a throttling. Bajar ~190 GB va a tomar varias sesiones y conviene
  espaciarlas. Si se vuelve intolerable: Pro, o volver a `prefetch` local.
- **Drive tiene un tope de ~750 GB/día de subida.** 190 GB entra, pero no de una.
- **El FUSE de Drive es lento con muchos archivos chicos.** Acá son ~425 archivos
  grandes, que es el caso bueno, pero conviene no generar temporales sobre el mount.
- **El repo pasa a público**: queda expuesto el diseño de la tesis antes de
  publicarla. Es reversible, pero no retroactivo — lo que se clonó, se clonó.
