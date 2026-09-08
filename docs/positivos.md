# Construcción del conjunto positivo

Decisión: **los positivos se asignan por secuencia, no por coordenada.**

Un locus de YASMA se marca positivo si el sRNA que produce coincide con un sRNA
ya descrito. No si su intervalo genómico se solapa con una anotación de
referencia.

## Por qué

Etiquetar por coordenada ata el proyecto a que nuestro ensamblado sea el mismo
que usa la base de referencia. Si no coinciden —y para `galga` no pudimos
siquiera averiguar cuál usa MirGeneDB— los miRNAs conocidos no caen sobre
nuestros loci, quedan como *unlabeled*, y el clasificador aprende que un miRNA
real es un candidato novedoso. Es el modo de falla exacto que PU learning
existe para evitar, y no produce ningún error visible: solo un modelo peor.

Etiquetar por secuencia elimina de un golpe tres problemas:

- la dependencia del ensamblado y el liftover entre versiones,
- la diferencia de nombres de cromosoma entre Ensembl (`1`) y NCBI (`NC_006088.5`),
- y la pregunta abierta de MirGeneDB para `galga`.

Además gana algo: un miRNA conservado que en nuestro organismo cae en un locus
que la anotación de referencia no tiene, igual se reconoce como positivo.

**Consecuencia sobre los genomas**: el ensamblado ya no tiene que coincidir con
el de ninguna base externa, así que se elige por sus propios méritos —
contigüidad, completitud, que sea la referencia vigente.

## Qué se compara contra qué

Esta es la parte que hay que hacer bien. **No** se compara la secuencia
genómica completa del locus contra el miRNA maduro: un locus puede medir
cientos de nt y el maduro mide ~22, así que cualquier medida de identidad
global sobre el locus entero no significa nada.

Lo que se compara es el **RNA mayoritario del locus** (el que YASMA reporta
como producto dominante) contra las **secuencias maduras** de la base.
Secundariamente, el precursor del locus contra el hairpin de referencia.

## Criterio de coincidencia

Un umbral de identidad global es la elección obvia y es la equivocada, porque
trata todas las posiciones como equivalentes. En miRNAs no lo son:

- El **extremo 5'** define la semilla y determina los blancos. Está bajo fuerte
  presión selectiva y es preciso: un corrimiento de una base ya es otro miRNA.
- El **extremo 3'** varía de forma rutinaria por isomiRs, adiciones no
  templadas y recorte. Exigir identidad ahí genera falsos negativos.

Criterio propuesto, a validar:

1. **5' exacto** en las primeras ~8 nt (semilla más posición 1).
2. **≤1-2 mismatches** en el cuerpo.
3. **Holgura en 3'**: ±2-3 nt de longitud, sin penalizar.

## El parámetro que hay que justificar

Cualquier umbral tiene un costo asimétrico y hay que declararlo en métodos:

- **Demasiado laxo**: parálogos y miembros de la misma familia entran como
  positivos. Las familias de miRNA comparten semilla por definición, así que un
  criterio flojo colapsa la familia entera en un positivo y el modelo aprende
  la familia, no el gen.
- **Demasiado estricto**: variantes reales quedan como *unlabeled*. Eso
  reintroduce exactamente el sesgo que motivó todo esto, solo que por otra
  puerta.

La forma honesta de fijarlo es medir la sensibilidad del resultado al umbral y
reportarla, en vez de elegir un número y defenderlo.

## Fuentes por organismo

| org | fuente principal de positivos | nota |
| :-- | :-- | :-- |
| gadmo | MirGeneDB | curado, declarado |
| galga | MirGeneDB | curado, declarado |
| maggi | MirGeneDB | curado, declarado; buscar como *Crassostrea gigas* |
| rhirr, sclsc, cloro | miRBase / Rfam | hongos: cobertura pobre en miRBase |
| phypa, prupe, maldo | miRBase / Rfam / PmiREN | plantas: revisar solapamiento entre bases |

Para los seis organismos sin MirGeneDB el conjunto positivo es más ruidoso.
miRBase tiene una tasa conocida de falsos positivos, que es precisamente lo que
MirGeneDB se propuso corregir (Fromm et al., 2019). Un positivo falso en PU
learning es peor que un unlabeled: contamina la clase que se asume limpia.

## Alcance: entrenamiento vs aplicación

**Decidido**: el modelo se entrena **solo** en los tres organismos con positivos
curados por MirGeneDB (`gadmo`, `galga`, `maggi`). Los otros seis son
**conjunto de aplicación**: se predice sobre ellos, no se entrena.

La razón es que en PU learning un falso positivo es peor que un *unlabeled*.
El método asume que la clase positiva está limpia y que el ruido vive en la no
etiquetada; meter entradas dudosas de miRBase en los positivos invierte ese
supuesto y no hay forma de recuperarse después.

### El costo, que hay que declarar

Los tres organismos de entrenamiento son **los tres animales**. Los seis de
aplicación son tres hongos y tres plantas. O sea que esto no es solo un cambio
de especie: es **transferencia entre reinos**, y es la parte más frágil del
diseño.

Los sRNA de plantas y de animales difieren en cosas que probablemente sean
justo las features del modelo:

- Los precursores de plantas son más largos y mucho más heterogéneos en
  longitud y estructura de hairpin que los de animales.
- Las plantas tienen clases de tamaño propias — 21 nt y 24 nt — mientras el
  pico animal es ~22 nt. Un modelo que aprendió "22 nt es señal de miRNA"
  puede descartar la clase de 24 nt de plantas por construcción.
- En hongos el panorama es todavía más distinto: varios linajes tienen vías de
  RNAi divergentes o reducidas, y lo que se describe suele ser milRNA y siRNA
  Dicer-dependiente antes que miRNA canónico. No por nada el proyecto primario
  de `cloro` son mutantes Dicer-like.

Si el modelo aprende la firma de un miRNA animal, aplicarlo a plantas y hongos
puede dar pocos candidatos no porque no los haya, sino porque busca la forma
equivocada.

### Dos cosas que hay que hacer para que el diseño se sostenga

**1. miRBase como evaluación, nunca como entrenamiento.** La decisión de no
entrenar con miRBase no obliga a ignorarlo. En los seis organismos de
aplicación, los miRNAs de miRBase sirven como *control de recuperación*:
¿cuántos de los ya descritos vuelve a encontrar el modelo? Es una medida
directa de si la transferencia entre reinos funciona, cuesta nada, y no
contamina nada porque no toca el entrenamiento. Si la recuperación es baja en
plantas y hongos, eso **es** un resultado, no un fracaso — pero hay que poder
medirlo.

**2. Validación cruzada dejando un organismo afuera, entre los tres curados.**
Antes de cruzar a otro reino, medir si el modelo transfiere entre taxones.
`maggi` es un molusco y `gadmo`/`galga` son vertebrados: entrenar en dos y
predecir en el tercero da una estimación honesta de la capacidad de
transferencia. Si el modelo no puede ir de pez a molusco, no va a ir de pez a
musgo, y conviene saberlo antes y no después de correr los nueve organismos.

Implementado en `scripts/loo_cv.py`:

```
./scripts/loo_cv.py features.tsv        # org, locus_id, y, <features>
./scripts/loo_cv.py --self-test         # verifica la lógica sin datos
```

### Cómo se mide, sin negativos

No hay negativos conocidos, así que no se reporta precisión ni AUC contra
"negativos". Las métricas son las que se sostienen con solo positivos y no
etiquetados:

- **enriquecimiento** = recall en el top 5% dividido por 5%. **1.0 es azar.**
  Es la métrica que responde la pregunta: si da ~1 en el organismo dejado
  afuera, el modelo no transfiere.
- **Lee & Liu** = recall² / P(predicho positivo). Criterio estándar de PU que
  no necesita negativos.
- **percentil mediano** de los positivos conocidos. 0.5 es azar.

### Un resultado del propio código, que hay que declarar en métodos

`scripts/loo_cv.py` implementa tres métodos, y el self-test verifica algo que
conviene saber antes de escribir nada:

**Elkan-Noto da exactamente el mismo ranking que tratar los no etiquetados como
negativos.** Su estimación es P(y=1|x) = g(x)/c, y como c es una constante,
dividir por ella no reordena nada. Elkan-Noto sirve para *calibrar* la
probabilidad y elegir umbral — no para mejorar el orden.

O sea que decir "usamos PU learning" no cambia por sí solo qué candidatos
salen priorizados. Para que el orden cambie hace falta un método que toque la
pérdida (nnPU) o el muestreo (**bagging PU**, Mordelet & Vert, que es el que
está implementado y sí reordena). Es una distinción que un revisor puede
preguntar y conviene tenerla resuelta.

### Normalización por organismo

Las features se estandarizan **dentro de cada organismo** antes de entrenar.
Sin eso, la diferencia de profundidad entre organismos —`galga` declara 32 M
spots y el primario de `maggi` ~353 M— entra como señal y el modelo aprende a
reconocer el organismo en vez del sRNA. El self-test verifica que un
corrimiento de escala artificial no rompe el resultado.

### El prior de clase

Casi todos los métodos de PU learning (Elkan-Noto, nnPU y parientes) necesitan
π, la proporción de positivos entre los no etiquetados. π estimado en animales
**no** es válido en plantas ni hongos: la proporción de loci que son miRNA real
sobre el total de loci anotados depende del organismo y del reino. Hay que
estimar π por organismo en el conjunto de aplicación, o reportar los resultados
como ranking en vez de como probabilidad calibrada.

### Alternativa, si la transferencia falla

Si la validación del punto 2 muestra que el modelo no transfiere ni entre
animales, el diseño de entrenar-en-tres-y-aplicar-a-seis no se sostiene y hay
que replantearlo: por ejemplo entrenar un modelo por reino, aceptando positivos
de miRBase filtrados en plantas y hongos, y declarando el ruido. Peor, pero
honesto. Conviene tener esto resuelto antes de escribir métodos.

## Referencia

Fromm, B., Domanska, D., Høye, E., et al. (2019). MirGeneDB 2.0: the metazoan
microRNA complement. *Nucleic Acids Research*, 48(D1), D132–D141.
https://doi.org/10.1093/nar/gkz885 (erratum: https://doi.org/10.1093/nar/gkz1016)
