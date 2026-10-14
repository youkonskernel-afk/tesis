# La matriz de features por locus

Una fila por locus y una columna por feature. Es la entrada del PU learning
(`scripts/loo_cv.py`). La arma `scripts/matriz.py` y se corre en Colab con
`notebooks/30_matriz.ipynb`. El banco es `tests/test_matriz.py`.

## De dónde salen los loci

Los loci los anota **`yasma tradeoff` sobre nuestro BAM**, el de `yasma align`
con `-m 50`. El alineamiento define los loci, igual que en ShortStack.

Lo que sigue se verificó contra el código de YASMA v1.1.1 (commit `c15e19d`), no
contra su documentación.

### Dónde escribe y qué escribe

`tradeoff` escribe en `<-o>/tradeoff/`, o en `tradeoff_<nombre>/` si se le pasa
`-n`.

| fichero | qué es |
| :-- | :-- |
| `loci.txt` | Una fila por locus, con las columnas de `assessClass.header` en `generics.py`: `Locus`, `Name`, `Length`, `Reads`, `RPM`, `UniqueReads`, `FracTop`, `Strand`, `MajorRNA`, `MajorRNAReads`, `Complexity`, `Gap`, `skew`, `size_1n`…`size_3n_depth`, `sizecall`. `matriz.py` exige esa cabecera exacta y corta si cambia. |
| `loci.gff3` | Los mismos loci en GFF. |
| `reads.txt` | Los reads más abundantes de cada locus. |
| `stats.txt` | Totales de la anotación. |
| `params.json`, `log.txt` | Los parámetros y el registro de la corrida. |

La anotación **no necesita el genoma**: lee los largos de los cromosomas de la
cabecera del BAM. Por eso se puede anotar directo desde los BAMs de Drive.

### `--max_length` no filtra reads

- El filtro por largo de `samtools_view` está **comentado** en v1.1.1.
- La llamada de tamaño (`sizecall`) usa **15-30 fijo**: está escrito en
  `sizeClass.update()`.
- O sea que todos los reads del BAM, de 15 a 50 nt, cuentan para cobertura,
  profundidad y RNA mayoritario. Solo la clase de tamaño mira 15-30.
- Igual se pasa `--min_length 15 --max_length 50`. El motivo es que `tradeoff` escribe
  ese valor en `inputs.json` y `params.json`, y ahí tiene que quedar la ventana del
  proyecto y no el 30 por defecto, que se leería como un filtro que no existe.

### Parámetros

Todos van explícitos (`TRADEOFF` en `matriz.py`), aunque coincidan con los
defaults de v1.1.1:

| parámetro | valor |
| :-- | :-- |
| `--genome_scaling_factor` | 0.4 |
| `--merge_dist` | 500 |
| `--merge_strand_similarity` | 0.5 |
| `--max_skew` | 0.9 |
| `--skew_ignore_range` | 20-24 |
| `--min_complexity` | 10 |
| `--min_abundance` | 50 |
| `--min_abundance_density` | 100 |
| `--min_size` | 19 |

Dos de estos filtros hay que nombrarlos en métodos:
- **`skew` mayor que 0.9 descarta el locus**: más del 90% de los reads son una sola
  secuencia. **Salvo** que su clase de tamaño caiga en 20-24, que es justo lo que
  protege a los miRNAs, muy sesgados por definición.
- **Una clase de tamaño menor a 19 nt descarta el locus**: se toma como
  degradación.

### Qué reads cuentan para un locus

YASMA cuenta todo read **mapeado** que solape el intervalo
(`bamf.fetch(chrom, start, stop)`, coordenadas base 0). `matriz.py` usa la misma
regla y guarda los dos conteos:
- `reads_yasma`: el de `loci.txt`;
- `reads_bam`: el propio.

Si difieren, no se está mirando lo mismo que miró YASMA (por ejemplo, se anotó
otro BAM). En ese caso `features` sale con código 1 y el notebook no sube la
matriz.

La secuencia de cada read es `get_forward_sequence()`: la molécula tal como se
secuenció, no la hebra del genoma. El RNA mayoritario viene así, y así se
compara contra los maduros.

## Las features

Son 31, en el orden de `FEATURES`. Todas son números. Se estandarizan **dentro de
cada organismo** en `loo_cv.py`, para que la profundidad de cada proyecto no entre
como señal.

| feature | de dónde | por qué |
| :-- | :-- | :-- |
| `log_largo` | `Length` | Los loci de miRNA son cortos (~22 nt de maduro en ~60-150 nt). Los clusters de siRNA son largos. |
| `log_reads`, `log_rpm` | `Reads`, `RPM` | Abundancia, en log10(x+1). |
| `log_unicos` | `UniqueReads` | Diversidad de secuencias. Un miRNA tiene pocas, un siRNA muchas. |
| `complejidad` | `Complexity` | Lecturas únicas por nt del locus. |
| `skew` | `skew` | Fracción del read más abundante. Es alta en miRNAs, que se procesan con precisión. |
| `sesgo_hebra` | `FracTop` | \|2·FracTop − 1\|: 1 si todo está en una hebra, 0 si está repartido. Los miRNAs son de una hebra; los siRNAs de dsRNA, de las dos. |
| `frac_15_17` … `frac_31_50` | reads | Distribución de largos: 15-17, cada uno de 18 a 25, 26-30 y 31-50. Separa miRNA (21-22), siRNA vegetal (24), piRNA (26-31) y tRF (30-40). |
| `frac_5A/C/G/U` | reads | Nucleótido 5'. Lo eligen las Argonautas: U para miRNA/AGO1, A para AGO2/AGO4 en plantas. |
| `precision_5p` | reads | En la hebra dominante, qué fracción de reads empieza en la posición 5' más común. Dicer corta preciso y la degradación no. En la hebra − el 5' de un read es su último nt en el genoma. |
| `frac_multimapeo` | marcas `XY` | Fracción de reads colocados por el sorteo ponderado (`P`) o al azar (`R`), sobre los colocados (`U` + `P` + `R`). Alta en loci de repetidos y transposones. |
| `frac_librerias` | `@RG` | En qué fracción de las librerías del proyecto aparece el locus. Lo real se repite entre réplicas. |
| `largo_mayoritario`, `gc_mayoritario`, `mayoritario_5U` | `MajorRNA` | El producto dominante del locus. |
| `mfe`, `mfe_por_nt` | genoma + ViennaRNA | Energía del pliegue de una ventana de ±80 nt alrededor del RNA mayoritario, en su hebra. Un precursor de miRNA pliega en horquilla estable. |
| `frac_pareado_mayoritario` | genoma + ViennaRNA | Qué fracción del RNA mayoritario queda apareada en ese pliegue. Un miRNA maduro vive en el tallo. |

### Por qué ±80 nt

Un pre-miRNA animal mide ~60-70 nt y el maduro puede estar en cualquiera de los
dos brazos. Con 80 nt de cada lado del mayoritario, la horquilla entra entera sin
importar el brazo. En plantas el precursor puede ser más largo, pero la horquilla
local alrededor del maduro sigue entrando.

Junto a un extremo del cromosoma la ventana se recorta. Si el RNA mayoritario no
aparece entre los reads del locus, el pliegue da NaN y `exportar` se niega a usar
esa fila.

### Lo que no es feature y no puede serlo

Si el locus reaparece en el duplicado. El duplicado es la validación
independiente, y usarlo como señal la anula.

## Etiquetado

Por secuencia, según `docs/positivos.md`: el RNA mayoritario contra los maduros
de la base.

- **Semilla:** las primeras **8 nt exactas**, sin desplazar el 5' (un corrimiento
  de una base ya es otro miRNA).
- **Cuerpo:** **a lo sumo 1 diferencia** en el resto de la parte común.
- **Largo:** una diferencia de **hasta 3 nt**, por los isomiRs del 3'.

Los tres umbrales son parámetros (`--semilla`, `--max-mm`, `--holgura`). La matriz
guarda la mejor coincidencia aunque quede fuera del umbral:
- `y_ref`: el maduro;
- `y_mm`: las diferencias;
- `y_dlen`: la diferencia de largo.

Así, medir la sensibilidad del resultado al umbral, que es lo que pide
`docs/positivos.md`, es releer la tabla y no volver a comparar.

`y = 0` significa **sin etiquetar**, no negativo.

Los positivos de entrenamiento salen de **MirGeneDB** para `gadmo`, `galga` y
`maggi`; `maggi` figura allí como *Crassostrea gigas*. Los códigos de especie y
las URLs de §3 del notebook están a confirmar en la primera corrida: esta sesión
de desarrollo no llega a mirgenedb.org. Lo bajado se valida antes de usarlo,
porque una URL que cambió devuelve una página HTML con código 200, y etiquetar
contra eso daría todo sin etiquetar sin ningún error.

## Exportar

`matriz.py exportar` produce `org, locus_id, y` y las 31 features, que es el
formato de `loo_cv.py`.

- Saca **solo los primarios de los organismos de entrenamiento**: `set_modelo` en
  `data/organismos.tsv`.
- El duplicado y los organismos de aplicación quedan afuera, y el resumen dice
  cuántos y por qué.
- Corta si:
  - una matriz no está etiquetada;
  - una feature viene vacía o NaN;
  - un `locus_id` se repite.
- No exporta nada de `y_ref`, `y_mm` ni `y_dlen`: le darían la etiqueta al modelo.

## Dónde queda cada cosa

| qué | dónde |
| :-- | :-- |
| Anotación | `Drive/20_yasma/<org>/<rol>/tradeoff/` |
| Matriz por proyecto | `Drive/40_features/<org>/<rol>.matriz.tsv` |
| Maduros de MirGeneDB | `Drive/40_features/maduros/<org>.fa` |
| Lo que entra a `loo_cv.py` | `Drive/40_features/features_entrenamiento.tsv` |
| Registro en git: qué matriz se armó | `data/matrices.tsv` (loci, positivos, features, **sha256**) |
| Registro en git: de dónde salieron los positivos | `data/maduros.tsv` (fuente, URL, sha256) |

## Lo que queda para después

- **`yasma hairpin`.** Evalúa los loci con criterios de miRNA tipo ShortStack:
  estructura, precisión, estrella. Daría features más finas, pero son 2700 líneas
  que no se verificaron. Va como agregado posterior, no en la primera matriz.
- **El brazo estrella (miRNA\*).** Si aparece el read complementario en el otro
  brazo de la horquilla, es la evidencia más fuerte de procesamiento por Dicer.
  Requiere ubicar el dúplex en el pliegue.
- **Las matrices de aplicación.** Se arman igual, pero se etiquetan contra
  miRBase solo para evaluar cuántos miRNAs ya descritos recupera el modelo, nunca
  para entrenar.
