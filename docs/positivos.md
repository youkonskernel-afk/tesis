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

**Pendiente**: decidir si los organismos sin MirGeneDB se etiquetan con miRBase
tal cual, con miRBase filtrado por algún criterio de calidad, o si el modelo se
entrena solo donde el positivo es curado y los demás se usan como conjunto de
aplicación. No es una decisión menor: cambia qué se puede afirmar del resultado.

## Referencia

Fromm, B., Domanska, D., Høye, E., et al. (2019). MirGeneDB 2.0: the metazoan
microRNA complement. *Nucleic Acids Research*, 48(D1), D132–D141.
https://doi.org/10.1093/nar/gkz885 (erratum: https://doi.org/10.1093/nar/gkz1016)
