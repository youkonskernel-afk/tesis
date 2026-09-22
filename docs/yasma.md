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

`yasma align` **no está en nuestro camino y no debería estarlo**: envuelve a
`ShortStack` (`--align_only --mmap u`) y pide `-tl/--trimmed_libraries`. Usarlo
reemplazaría nuestro bowtie y con él el `-m 50`, que está medido y justificado.

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

**Conclusión operativa: el que decide es `perfil`, no `yasma adapter`.** Nuestro
pipeline recorta con `fastp`, no con YASMA, así que `yasma adapter` no está en
el camino crítico — y está bien que no esté.

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
