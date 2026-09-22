# YASMA: qué consume, qué decide, y dónde no hay que confiarle

Medido contra el código, no leído de la documentación. Repo:
[`NateyJay/YASMA`](https://github.com/NateyJay/YASMA). Todo lo de acá está
verificado en el tag **`v1.1.1`** (commit `c15e19d`), que es el que pinea el
proyecto.

## Cómo pinearlo, porque el nombre engaña dos veces

**El branch por defecto no es la release.** El default del repo es
`library-scaling`, no `main`. Un `git clone` pelado o un
`pip install git+https://github.com/NateyJay/YASMA` te deja en 1.1.0, no en el
v1.1.1 que este proyecto usa.

**Y el tag no coincide con la versión declarada.** `pyproject.toml` dentro de
`v1.1.1` dice `version = "1.1.0"`. O sea que los metadatos del paquete
instalado reportan 1.1.0 incluso cuando instalaste el tag correcto: **la
versión no se puede verificar desde el paquete**, hay que pinear por tag o
commit.

En `environment_pip.txt` va la ref exacta, no un número:

```
git+https://github.com/NateyJay/YASMA@v1.1.1
```

Dependencias que no son obvias: Python **>= 3.12** (de ahí el nombre del
entorno `srna2`) y **ViennaRNA** — `yasma/__init__.py` importa `hairpin`, que
hace `import RNA` en el nivel superior, así que **sin ViennaRNA no arranca
ningún subcomando**, ni los que no tienen nada que ver con estructura. Es el
fallo que `check_env.sh` está hecho para atrapar.

## El punto de integración: YASMA anota nuestro BAM

`yasma tradeoff -a <BAM>` es el que escribe `loci.gff3`. Lee el BAM con `pysam`
y lo indexa si le falta el índice. **Nuestro alineamiento es el que se anota**,
así que la decisión de `bowtie -m 50` vale tal como está.

`yasma align` **sí es un camino viable, al contrario de lo que decía esta
página**. Ver la sección propia más abajo: en `v1.1.1` no envuelve a ShortStack.

### El BAM tiene que traer `@RG` o `tradeoff` revienta

`get_chromosomes()` hace `header['RG']` sin `.get()`. Medido:

| BAM | resultado |
| :-- | :-- |
| sin `@RG` | `KeyError: 'RG'` |
| con `@RG` | `librerias=['SRR1']` |

No degrada, no avisa: corta. Y no es cosmético — `tradeoff` agrega profundidad
por read group (`aggregate_by=['rg','chrom']`), así que sin `@RG` no puede
separar librerías aunque no se cayera. El paso de `samtools` de
`orchestrate.sh` **tiene que poner un `@RG` por corrida**; `yasma readgroups`
es la utilidad para inspeccionarlos.

### Tres trampas de rutas, las tres del mismo tipo

`yasma adapter` falló tres veces seguidas antes de correr, siempre por rutas:

1. `-o` tiene que **existir** de antes (`InputError: not found!`).
2. Las librerías tienen que estar **adentro** de `-o`, porque internamente hace
   `relative_to(output_directory)`.
3. `-o` tiene que ser **absoluto** si la librería lo es, o `relative_to`
   compara una ruta absoluta contra una relativa y tira `ValueError`.
4. Y hay que correrlo **parado dentro** de `-o`: guarda la ruta relativa pero
   después la abre desde el CWD.

O sea: `cd` al directorio del proyecto, librerías adentro, `-o .` absoluto.

## `yasma align` no es el wrapper de ShortStack — es bowtie1 nativo

Esta página decía que `yasma align` envolvía a `ShortStack --align_only --mmap u`
y que usarlo nos costaría el `-m 50`. **Las dos cosas son falsas en `v1.1.1`.**

El wrapper de ShortStack existe —`src/yasma/align.py`, función
`shortstack_align`— pero está **comentado en `__init__.py`**:

```python
# from .align import *
...
from .nativealign import *
```

El comando `yasma align` que se registra es `nativealign.py:align()`, un
alineador **bowtie1 nativo con el pesado estilo ShortStack3**. Los defaults:

| opción | default | qué es |
| :-- | :-- | :-- |
| `--max_multi` | **50** | el `-m` de bowtie — **es exactamente nuestro `-m 50`** |
| `--max_random` | 3 | sitios empatados por encima de los cuales el read queda sin mapear |
| `--unique_locality` | 50 | ventana en nt para pesar por cobertura única local |
| `--offrate` | 3 | de `bowtie-build`; bowtie usa 5 |
| `--min_length` / `--max_length` | 15 / 50 | filtro aplicado **en el alineamiento** |

### Tres etapas por librería, contra el mismo índice

1. **`unique`** — `bowtie -v 1 -p <cores> -S -m 1 --best --strata --offrate 3
   --max <RG>.max1.fq`. Lo que mapea a un solo sitio entra al BAM; el resto cae
   al fichero de `--max`.
2. **`multi`** — sobre ese fichero: `bowtie -v 1 -S -m <max_multi> -a --best
   --strata --offrate 3 --max <RG>.max50.fq`. Con todas las posiciones en mano,
   pesa cada una por la cobertura única en una ventana de `unique_locality/2` a
   cada lado y elige una por sorteo ponderado. Si todos los pesos empatan **y**
   hay más de `max_random` sitios, el read queda sin mapear. Tags: `XY:Z:U/R/P/Q`
   y `XZ:f:<prob>`.
3. **`over`** — los reads que pasaron `max_multi` **no se tiran**: entran al BAM
   como no mapeados con `XY:Z:H`. Esta etapa no corre bowtie, lee el fichero de
   `--max`.

El filtro de longitud y el de N se aplican acá (`XY:Z:F`), no antes.

### Lo que esto resuelve

- **`bowtie` corre en `-v 1` en las dos pasadas.** `-v` cuenta mismatches y
  **ignora las calidades**, así que la calidad sintética única de `SRR317135` y
  `SRR1066790` (las dos SRA Lite) no cambia el alineamiento. Era la mitad
  abierta de esa nota de métodos, y se cierra sin necesitar `config.sh`.
- **`-m 50` no se pierde: es el default.** El hallazgo de `danre` —que subirlo
  inunda la anotación de tRF-5— sigue valiendo y no hay que pelearlo.
- **Escribe `@RG` por librería**, tanto en la cabecera (`header['RG']`) como por
  read (`a.set_tag("RG", rg, "Z")`), que es justo lo que `tradeoff` exige sin
  `.get()`. El nombre del read group sale de `get_rg()`, que pela `.gz`, `.t`,
  `.fq`/`.fastq`: de `SRR123.t.fq.gz` sale **`SRR123`**, y de un PRE-TRIMMED
  `SRR123.fastq.gz` también. O sea **un `@RG` por corrida, ya resuelto**.
- Deja `align/alignment.bam` ordenado por coordenada, la tabla de profundidad
  que `tradeoff` consume, y `alignment_file` anotado en `inputs.json`.

### Dos cosas que hay que tener listas antes

- **El genoma no puede estar comprimido con `gzip`.** `make_bam_header()` hace
  `pysam.FastaFile(genome_file)`, y `bowtie-build` recibe el mismo fichero.
  Nuestros ensamblados están en `70_genomas/<acc>.fna.gz` hechos con `gzip -c`.
  Medido:

  | fichero | `pysam.FastaFile` |
  | :-- | :-- |
  | `g.fna` | OK |
  | `g.fna.gz` (gzip) | `OSError error when opening file` |

  Hay que descomprimirlo, o re-comprimirlo con `bgzip`. El `sha256` del ledger
  es el del `.gz`, así que la copia descomprimida es derivada y no se respalda.
- **El índice se construye solo** si falta `<genoma>.rev.1.ebwt`, con
  `bowtie-build --offrate 3`, **al lado del FASTA**. Coincide con la regla del
  proyecto de no respaldar índices.

### Lo que sigue sin resolver

`yasma align` alinea **un proyecto entero a un BAM único**, con un `@RG` por
librería. Eso encaja con `trim/<org>_<rol>/`: un BAM por organismo y rol, que es
la unidad que `tradeoff` anota. Lo que falta decidir es si reemplaza al bowtie de
`orchestrate.sh` —que todavía no está en este repo— o convive con él. No se puede
comparar hasta que ese script llegue del `main` local.

## `yasma adapter`: su criterio de "ya recortada" es más débil que el nuestro

`yasma adapter` marca una librería `PRE-TRIMMED`. El criterio, literal:

```python
pretrim = read_length_freq < 0.8 and best_perc < 0.10
```

- `read_length_freq`: fracción de reads que tienen el largo modal.
- `best_perc`: fracción de reads con el mejor adaptador conocido.

Y solo escribe `PRE-TRIMMED` si además `best == "None"`.

**Usa la dispersión del largo, no su magnitud.** Ahí está la diferencia con
`fetch_runs.sh perfil`, que usa el largo mediano. Corrí los dos sobre las
mismas 20 000 lecturas sintéticas de cuatro tipos:

| librería | largo modal | freq. modal | YASMA | `perfil` |
| :-- | :-- | :-- | :-- | :-- |
| recortada, largos variables (caso `cloro`) | 22 nt | 0.298 | **PRE-TRIMMED** | `YA RECORTADA` |
| recortada a **un solo largo** | 22 nt | 1.000 | `None` | `YA RECORTADA` |
| mRNA, 150 nt fijos (caso `SRR23277331`) | 150 nt | 1.000 | `None` | `NO PARECE` |
| sRNA-seq sin recortar | 150 nt | 1.000 | `TGGAATTC` (100%) | `PARECE sRNA-seq` |

Coinciden en las tres que importan hoy. **Divergen en la segunda fila**, y la
divergencia es exactamente el bug que `perfil` tuvo y se arregló: YASMA da
`None` tanto para una librería ya recortada a largo fijo como para mRNA. Son
dos situaciones opuestas —el read *es* el inserto contra el inserto es más
largo que el read— y el veredicto es el mismo, porque el largo del read no
entra en la decisión. YASMA lo imprime (`22` contra `150`) pero no lo usa.

**Ninguno de nuestros 19 proyectos está en ese estado**, así que YASMA acierta
en los 19. Pero por suerte, no por diseño: el único pre-recortado es
`cloro PRJEB43636`, con largos de 30-34 nt, y esa variación es justo lo que el
criterio de YASMA necesita. Si algún reemplazo de BioProject llegara recortado
a largo fijo, YASMA lo llamaría "sin adaptador" sin distinguirlo de mRNA.

**Conclusión operativa: el que decide es `perfil`, no `yasma adapter`.** El
recorte lo hace `yasma trim` (ver abajo), pero el adaptador se lo damos nosotros
desde `data/adaptadores.tsv`: `yasma adapter` no está en el camino crítico, y es
a propósito.

### Detalle menor, pero conviene saberlo

`best_perc` dio **100.1%** en la librería de sRNA sin recortar. Cuenta por
k-mer y puede pasarse de 100, así que `--min_adapter_content` compara contra un
número inflado. No cambia nada con un umbral de 0.1, pero no es una proporción.

## Cómo reproducir esta comparación

Sin red y sin datos reales: se generan las cuatro librerías sintéticas, se
corre `yasma adapter` sobre cada una y se compara contra `perfil`. El generador
usa los mismos tipos de librería que `tests/test_perfil.sh`.

```bash
python3.12 -m venv vy
./vy/bin/pip install ViennaRNA 'git+https://github.com/NateyJay/YASMA@v1.1.1'
# librerías adentro del outdir, -o absoluto, y correr parado adentro
mkdir -p p && cp lib.fastq p/
( cd p && yasma adapter -ul "$PWD/lib.fastq" -o "$PWD" --override -n 20000 )
```

## `yasma trim`: lo que le pasa a cutadapt, y lo que eso implica

`yasma trim` es un wrapper de cutadapt. La llamada, literal:

```
cutadapt -a <sec> --minimum-length 15 --maximum-length 50 -j <cores> -O 4 \
         --max-n 0 --trimmed-only -o <out> <in>
```

Los defaults de `--min_length`/`--max_length` son **15 y 50**, o sea la misma
ventana que el proyecto eligió a propósito. Coincidencia afortunada, no diseño
compartido: `scripts/trim.sh` los pasa explícitos igual, porque un default que
cambie en una versión nueva no avisa.

Cuatro consecuencias que hay que declarar en métodos:

- **`--trimmed-only` descarta los reads sin adaptador.** Es distinto de `fastp`,
  que los conservaría. Para sRNA-seq es lo correcto —un read sin adaptador tiene
  el inserto más largo que el read, o sea que no es un sRNA— pero significa que
  la retención esperada **es el `adapt_pct` de `data/adaptadores.tsv`**, no ~100%.
- **`--max-n 0`** tira cualquier read con una sola N.
- **No hay filtro de calidad**: cutadapt se llama sin `-q`. Eso **cierra la mitad
  abierta de la nota de SRA Lite**: en este paso la calidad sintética única de
  `SRR317135` y `SRR1066790` no distorsiona nada, porque no hay nada que
  distorsionar. Queda solo la pregunta del alineamiento (`-v` contra `-n`/`-e`).
- **Una librería `PRE-TRIMMED` no se filtra por longitud.** YASMA la pasa de
  largo tal cual, sin llamar a cutadapt. Para `cloro` da igual —sus reads de
  30-34 nt están todos dentro de 15-50— pero no es una regla general.

### Por qué el adaptador NO se lo dejamos a `yasma adapter`

`yasma trim` lee los adaptadores de `inputs.json`. Ante un adaptador `"None"`
**descarta la librería**: no la suma a `trimmed_libraries` y desaparece del
pipeline sin ningún error. Y `yasma adapter` devuelve `None` tanto para una
librería ya recortada a largo fijo como para mRNA — medido arriba. Una librería
mal clasificada se perdería en silencio.

Así que el adaptador se decide con `fetch_runs.sh perfil`, se versiona en
`data/adaptadores.tsv`, y `scripts/trim.sh` escribe `inputs.json` desde ahí. Un
proyecto que no esté en la tabla hace fallar el script, que es mejor que
recortar con una secuencia adivinada.

### El prefijo que detecta no es la secuencia que recorta

Los prefijos de `ADAPTADORES` sirven para **detectar**. Cotejados contra los 161
adaptadores de la tabla de YASMA, tres cosas:

| prefijo | qué es de verdad |
| :-- | :-- |
| `TGGAATTCTCGGG` | Illumina RNA 3p Adapter (RA3). Prefijo compartido por **50** entradas de la familia RPI: identifica la familia, no un adaptador |
| `GATCGTCGGACTG` | `RNA_Adapter_(RA5)` — un adaptador **5'**. Encontrarlo es dímero o quimera, **no** read-through 3': no sirve como `-a` |
| `CGCCTTGGCCGT` | **no aparece en ninguno de los 161.** Procedencia desconocida |

`ATCTCGTATGCCG` y `TCGTATGCCGTCTTCTGCTTG` sí son el adaptador small-RNA de 2011;
sus 110 coincidencias son constructos modernos que lo contienen aguas abajo.

Por eso `perfil` marca los dos primeros con `5p:` y `??:`, y `trim.sh` **se
niega a recortar** con ellos.

### Tres cosas del propio `trim` que obligan a envolverlo

Las tres medidas contra el binario de `v1.1.1`, no leídas:

**1. `trimmed_libraries` se pisa, no se acumula.** `trim()` hace
`ic.inputs['trimmed_libraries'] = []` al entrar y `= <lo de esta llamada>` al
salir. Comprobado con dos librerías en tandas separadas:

| | `trimmed_libraries` tras la llamada | qué hay en `trim/` |
| :-- | :-- | :-- |
| tanda 1 (RUNA) | `['trim/RUNA.t.fq.gz']` | `RUNA.t.fq.gz` |
| tanda 2 (RUNB) | `['trim/RUNB.t.fq.gz']` | `RUNA.t.fq.gz`, `RUNB.t.fq.gz` |

El fichero de la tanda 1 **sigue en disco y desaparece del registro**. Como
`trim.sh` lee de ahí para saber qué está recortado, sin corregirlo cada tanda
desmentiría a la anterior: `estado` reportaría la corrida como faltante y
`correr` la volvería a volcar y recortar — con 417 corridas, para siempre. De
ahí el ledger `recortadas.tsv`, que es acumulativo, y la re-escritura de
`inputs.json` con la lista completa después de cada tanda (que es lo que
`tradeoff` y compañía van a leer).

Pasarle **todas** las librerías en cada llamada tampoco sirve: no saltea las que
ya tienen salida, las re-recorta desde cero.

**2. `trim/log.txt` se trunca en cada llamada.** `Logger.__init__` abre el
fichero con `"w"`. Las estadísticas de cutadapt de las tandas anteriores se
pierden, así que cada tanda se loguea aparte y de esos logs salen los conteos
de `trim.sh verificar`.

**3. `--cleanup` no se puede usar.** Itera `ic.inputs['srrs']`, que en un
`inputs.json` que no escribió `yasma download` queda en `None` → `TypeError`. Y
además vacía `untrimmed_libraries`, que para una librería `PRE-TRIMMED` es
justamente donde está la salida. El borrado del fastq sin recortar lo hace
`trim.sh`.

### Dos nombres de salida que no se pueden adivinar

Probado contra el YASMA real, no leído:

- El fichero es **`<RUN>.t.fq.gz`**, no `<RUN>.tfq.gz`: YASMA hace
  `'.t' + library_format` y ese formato ya viene con punto (`.fq`).
- Una librería `PRE-TRIMMED` **no produce fichero nuevo**.

Con el nombre adivinado nada contaba como recortado: ni la idempotencia ni el
`estado`. `trim.sh` lee `trimmed_libraries` de `inputs.json`, que es el registro
que YASMA deja de lo que produjo.

Y un detalle que cuesta una corrida entera: **`yasma trim` no tiene
`--override`** (sí lo tiene `yasma adapter`). Pasárselo por analogía aborta.
