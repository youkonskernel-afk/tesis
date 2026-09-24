# Guion de presentación — el upstream

Todo el contenido de la presentación en texto corrido, para leer, recortar o
pasar a otro formato. Está ordenado por las **10 láminas** del deck, y cada
bloque trae lo que se muestra, lo que se dice, y las preguntas probables.

Las cifras salen de `data/` y de mediciones registradas en `CLAUDE.md` y
`docs/yasma.md`. Donde hay una hipótesis sin comprobar, está dicho.

Duración estimada leyendo el guion completo: **35-45 min**. Para 20 min,
recortar los bloques marcados *(ampliación)*.

---

## 1 · Portada

**Se muestra:** título, y cuatro cifras — 9 organismos, 3 reinos, 417 corridas,
7.90 G spots.

**Se dice:**

Esta es la presentación del *upstream*: la parte bioinformática que produce los
loci candidatos sobre los que después corre el clasificador. No es el aporte de
la tesis, es lo que lo hace posible, y es donde está casi todo el trabajo hecho
hasta ahora.

Voy a recorrer cinco cosas: qué se busca y por qué el método tiene que ser
positive-unlabeled; cómo se eligieron los nueve organismos y sus datos; cómo se
eligieron y verificaron los genomas; qué herramientas corren y por qué YASMA es
la central; y qué problemas aparecieron.

Hay un hilo que atraviesa todo y conviene tenerlo presente desde ahora: **en
este pipeline casi ningún error se manifiesta como un error**. Por eso la mitad
de lo construido son verificaciones y no funcionalidad.

---

## 2 · El foco

**Se muestra:** qué se busca; por qué PU y no clasificación binaria; el supuesto
de la clase positiva limpia; el etiquetado por secuencia.

### Qué se busca

Descubrir **sRNAs que todavía no están descritos**, priorizando candidatos con
un clasificador entrenado sobre los que sí lo están.

Cuando se secuencia RNA pequeño de un organismo y se alinea contra su genoma,
aparecen miles de regiones con lecturas acumuladas. A esas regiones las llamamos
**loci**. Una fracción chica corresponde a sRNAs ya descritos en las bases
públicas. El resto es una mezcla: sRNAs reales todavía no caracterizados,
fragmentos de degradación, productos de otras rutas.

La pregunta de la tesis es **cuáles de esos loci no anotados valen una
validación experimental**. Y la respuesta tiene que ser un **ranking**, no un sí
o un no, porque lo que sigue es trabajo de laboratorio con presupuesto finito.

El aporte no es el pipeline: es el planteo.

### Por qué positive-unlabeled

Tenemos un conjunto **positivo** —los sRNAs curados en las bases públicas— y un
conjunto **sin etiqueta** —todos los demás loci. Lo que **no** tenemos es un
conjunto **negativo**: nadie publica "esta región, mirada a fondo, no es un
sRNA".

La objeción natural es: ¿por qué no entrenar un clasificador binario normal, con
los anotados como positivos y el resto como negativos? **Porque el resto
contiene exactamente lo que buscamos.** Un modelo binario aprendería a llamar
negativo a un sRNA no descrito, que es la definición de nuestro objetivo. El
error no sería un bug: sería el comportamiento correcto de un método mal
elegido.

PU learning trata al conjunto sin etiqueta como **una mezcla** de positivos y
negativos en proporción desconocida. Un locus sin etiqueta con probabilidad alta
no es un error del modelo: **es el resultado**.

### El supuesto que hay que proteger

En PU learning **un falso positivo es peor que un no etiquetado**. El método
asume que lo marcado como positivo lo es; una entrada dudosa invierte ese
supuesto sin vuelta atrás.

Consecuencia práctica: no todas las bases públicas de sRNA tienen el mismo nivel
de curación. Hay bases que agregan lo que se publica, y bases que curan
manualmente cada entrada verificando que cumpla los criterios de biogénesis.
**Sólo las segundas sirven para entrenar bajo PU.** Las primeras son excelentes
para **evaluar**, y ahí se usan.

Tres cosas que van en métodos:

- **Elkan-Noto no cambia el ranking** respecto de tratar lo no etiquetado como
  negativo: divide por una constante. Sirve para calibrar, no para reordenar.
  Verificado. Como lo que entregamos es un ranking, no aporta.
- Para que el orden cambie hace falta **bagging PU** o **nnPU**.
- El **prior de clase** se estima **por organismo**: el de animales no vale en
  plantas ni en hongos. Si no se puede estimar bien, se reporta ranking y no
  probabilidad calibrada.

### El etiquetado va por secuencia, no por coordenada

Lo intuitivo sería marcar positivo el locus cuyo **intervalo se solapa** con una
anotación publicada. Está mal.

Lo que se hace: marcar positivo el locus cuyo **RNA mayoritario coincide** con
un sRNA descrito.

El motivo es un fallo silencioso. Por coordenada, el etiquetado depende del
ensamblado. Si alineamos contra un genoma y la base publicó sobre otro, los
sRNAs conocidos **no caerían sobre nuestros loci**: quedarían sin etiqueta, y el
modelo aprendería que un sRNA real es un candidato novedoso. **Es el modo de
falla exacto que PU learning evita** — y sin ningún error visible.

Dos cosas que se siguen y es fácil hacer mal:

- **No** comparar la secuencia genómica completa del locus contra el maduro: el
  locus mide cientos de nucleótidos y el maduro unos 22. La identidad global no
  significa nada.
- **No** usar un umbral de identidad plano. En los miRNAs la región que
  determina a qué transcripto se unen son los nucleótidos 2 a 8 desde el extremo
  5'; el extremo 3' varía de rutina entre isoformas. Un umbral plano trataría
  las dos regiones como equivalentes. El criterio va **anclado en 5'** con
  holgura en 3'.

**Beneficio secundario, y es grande:** al no depender de coordenadas, la
elección de ensamblado deja de estar atada a qué genoma usó cada base pública.
Podemos elegir por contigüidad y completitud, que es lo que conviene para
alinear.

---

## 3 · Selección de organismos

**Se muestra:** el diseño 3×3, la tabla de los nueve, primario/duplicado,
entrenamiento vs aplicación, el costo y las salvaguardas.

### El diseño

Nueve organismos, **tres por reino**: hongos, plantas y animales.

- **Por qué tres reinos:** las rutas de biogénesis de sRNA difieren mucho entre
  reinos. Es lo que hace la pregunta interesante **y** lo que la hace difícil.
- **Por qué tres por reino:** con uno solo, cualquier resultado puede ser una
  particularidad de esa especie. Con tres hay variación interna contra la que
  contrastar.
- **Por qué no más:** el cuello de botella no es el cómputo, es la
  verificación. Cada BioProject hay que perfilarlo para saber si de verdad es
  sRNA-seq y con qué adaptador se recorta; cada genoma hay que resolverlo y
  verificarlo. Nueve es lo que entra con verificación seria de cada paso.

Y dentro de cada organismo, **dos BioProjects independientes**. 9 × 2 = **18
proyectos**.

### Los nueve

| org | especie | reino | primario | duplicado |
| :-- | :-- | :-- | :-- | :-- |
| rhirr | *Rhizophagus irregularis* | Fungi · Glomeromycota | PRJEB29180 | PRJNA722321 |
| sclsc | *Sclerotinia sclerotiorum* | Fungi · Ascomycota | PRJNA477286 | PRJNA1135930 |
| cloro | *Clonostachys rosea* | Fungi · Ascomycota | PRJEB43636 | PRJEB51338 |
| phypa | *Physcomitrium patens* | Plantae · briofita | PRJNA222997 | PRJNA277372 |
| prupe | *Prunus persica* | Plantae · rosácea | PRJNA929031 | PRJNA780811 |
| maldo | *Malus domestica* | Plantae · rosácea | PRJNA681626 | PRJNA784097 |
| gadmo | *Gadus morhua* | Animalia · pez | PRJNA284846 | PRJNA328800 |
| galga | *Gallus gallus* | Animalia · ave | PRJEB12164 | PRJNA694114 |
| maggi | *Magallana gigas* | Animalia · molusco | PRJNA154615 + PRJNA232734 | PRJNA1254880 |

Tres cosas para señalar sin leer la tabla entera:

1. **Distancia filogenética deliberada dentro de cada reino.** En animales hay un
   pez, un ave y un molusco: si el modelo no transfiere de pez a molusco, no va a
   transferir de pez a musgo, y eso se puede saber temprano.
2. **Las dos rosáceas son un control interno casi ideal.** Durazno y manzana: si
   el modelo funciona en una y no en la otra, el problema no es la distancia
   evolutiva.
3. **19 accessions en 18 slots.** El primario de la ostra son dos BioProjects
   combinados porque ninguno solo daba profundidad suficiente. Eso trajo un
   problema propio: son de 2011 y 2014, o sea de eras de kit distintas.

### Primario y duplicado

- **Primario · descubrimiento.** Sobre él se anotan loci y se corre el modelo.
- **Duplicado · validación.** Un experimento independiente, de otro laboratorio
  y otra condición, para preguntar si un candidato **reaparece**.

Un locus que el modelo prioriza en el primario **y** que reaparece en el
duplicado es mucho más defendible que uno que sólo existe en un experimento. Es
la diferencia entre un candidato y un artefacto de una librería.

Por eso **el duplicado no entra al entrenamiento**. Si entrara, la validación
dejaría de ser independiente y el criterio de reaparición no significaría nada.

La separación está **impuesta por la estructura de directorios**, no por
acordarse: primario y duplicado son proyectos distintos aguas abajo, y la unidad
de anotación es el directorio. Si compartieran directorio, la independencia
quedaría en "acordarse de filtrar por una columna" — y se perdería el día que
alguien no se acuerde. Son **18 directorios de trabajo, no 9**.

### Entrenamiento vs aplicación

- **Entrenamiento · 3:** `gadmo`, `galga`, `maggi`. Los tres con positivos
  curados de alta calidad.
- **Aplicación · 6:** `rhirr`, `sclsc`, `cloro`, `phypa`, `prupe`, `maldo`. Se
  predice sobre ellos, no se entrena.

Está explícito en una columna del fichero de especificación versionado. No es
una convención que haya que recordar al analizar.

### El costo del diseño, declarado

**Los tres de entrenamiento son los tres animales**, y los seis de aplicación
son hongos y plantas. Es transferencia entre reinos y es la parte más frágil del
diseño. Lo digo yo antes de que lo pregunten.

- En **plantas** los precursores son más largos y heterogéneos, y las clases de
  tamaño son **21 y 24 nt** contra el pico animal de ~22.
- En **hongos** predominan milRNA y siRNA dependiente de Dicer sobre el miRNA
  canónico.
- **El riesgo concreto:** un modelo que aprendió la forma animal puede no
  encontrar nada en plantas por buscar la forma equivocada — y eso se confunde
  fácil con "no hay nada que encontrar".

El origen del problema no es una mala elección: es que las bases curadas a mano
cubren sobre todo animales. No hay un equivalente de la misma calidad para
plantas o para hongos. **La elección real era entre transferencia entre reinos
con positivos limpios, o entrenamiento dentro del reino con positivos sucios.**
Se eligió lo primero.

### Las dos salvaguardas

1. **Las anotaciones públicas como evaluación, nunca como entrenamiento.** En
   los seis de aplicación, **cuántos sRNAs ya descritos recupera el modelo** es
   la medida directa de si la transferencia funciona. No contamina nada: esos
   datos nunca entran al ajuste. Si recupera muchos, transfiere. Si no recupera
   ninguno, no transfiere — y eso es un resultado interpretable, no una lista
   vacía ambigua.
2. **Validación cruzada dejando un organismo afuera, entre los tres curados.**
   La ostra es molusco; los otros dos son vertebrados. Si el modelo no
   transfiere de pez a molusco, no va a transferir de pez a musgo, y conviene
   saberlo **antes** de correr los nueve organismos. Ya está implementada, con
   una prueba que verifica su lógica sin necesitar datos reales.

La segunda convierte un posible fracaso en un resultado temprano y barato, en
vez de en cien horas de cómputo y una lista vacía.

---

## 4 · Selección y filtrado de los datos

**Se muestra:** de 19 accessions a 417 corridas; el filtro duro; el perfilado;
los adaptadores; la retención; las exclusiones; la escala.

### De la intención a los datos

Cuatro pasos, reproducibles:

1. **La especificación.** Un fichero versionado con los 19 BioProjects, su rol y
   para qué sirve cada organismo. Lo escribe una persona.
2. **La consulta.** Se resuelve contra el repositorio europeo de nucleótidos,
   que devuelve las corridas de cada proyecto con sus metadatos.
3. **El filtro.** Se descarta lo que no puede ser sRNA-seq.
4. **El manifiesto.** 417 corridas con identificador, proyecto, rol y conteos.
   Versionado en git.

**El manifiesto es la verdad; la especificación es de cuando se eligió el
proyecto.** Difieren en 7 de los 19 proyectos, y eso está bien.

*Lección que costó cara:* la especificación tiene una columna con el número
aproximado de corridas. Se tomó ese número estimado y se escribió como un hecho
de métodos — se declaró tres veces que la validación de un organismo era "más
débil" porque tenía una sola corrida, y el manifiesto real devolvió dos. **Lo
que va a métodos sale del dato resuelto, no de la nota de cuando se eligió.**

### El filtro duro

- **Sólo datos de RNA:** se exige `library_source = TRANSCRIPTOMIC`. Hay
  BioProjects que agrupan todo el trabajo de un laboratorio sobre un organismo,
  incluyendo secuenciación de genoma. Sin este filtro se bajan gigabytes de
  lecturas de DNA para después descubrir que no tienen adaptador de sRNA.
- **Sólo lecturas simples:** para las etiquetadas `RNA-Seq` se exige `SINGLE`.
  Un proyecto de RNA-Seq con lecturas apareadas es, casi siempre, RNA mensajero.
  Los sRNA se secuencian en lectura simple porque el inserto es más corto que la
  lectura.

Este filtro es **estructural** —opera sobre metadatos, es barato— y **no
alcanza**.

### El perfilado: la etiqueta engaña, el adaptador no

Una corrida etiquetada **miRNA-Seq**, `TRANSCRIPTOMIC` y `SINGLE`, pasó todos
los filtros formales y **era RNA mensajero**. Se midió: **0 de 40 000 lecturas
con adaptador 3'**.

Por qué cero es concluyente: si el inserto fuera corto, la secuenciación leería
el inserto entero y seguiría leyendo el adaptador que viene después. Eso se ve
como el adaptador apareciendo dentro de la lectura. Probamos contra una lista
que **incluye el adaptador universal de Illumina**, así que sea cual sea el kit,
algo tendría que aparecer. Cero en cuarenta mil lecturas significa que ningún
inserto terminó antes de los 150 nucleótidos.

Y el largo promedio no distingue: un sRNA de 22 nt corrido en 2×150 da un
promedio de 273 nt, **igual que un mRNA**.

Lo único que separa los dos casos es **dónde empieza el adaptador**, y eso hay
que medirlo abriendo el fichero. Por eso existe un comando de perfilado que baja
una muestra de lecturas y mide dónde aparece cada adaptador conocido. **Se corre
antes de adoptar un proyecto, no después de bajarlo entero.** Perfilar cuesta
minutos; bajar un proyecto para descubrir que no sirve cuesta días.

### El resultado del perfilado

**19 proyectos perfilados. 18 dieron "parece sRNA-seq"; 1 ya viene recortado de
origen.**

Pero cuatro estaban marcados como sospechosos al principio, y **los cuatro eran
defectos de la herramienta, no de los datos**:

- La lista de adaptadores del perfilador era **sólo moderna**. Dos proyectos de
  2011 y 2014 daban cero por ciento. Al agregar las secuencias de la época, uno
  pasó de **0% a 92%** y otro de **0% a 97%**. Eran sRNA-seq perfectamente
  normales — y uno es el primario de un organismo de entrenamiento.
- El proyecto pre-recortado: lecturas de 30-34 nt donde **la lectura es el
  inserto**. La herramienta lo llamaba "no parece sRNA-seq". **Casi se tiran 34
  corridas buenas.** Cero por ciento de adaptador con lecturas cortas significa
  una cosa; cero por ciento con lecturas largas significa la opuesta.

### Con qué adaptador cada uno

| familia | proyectos | nota |
| :-- | --: | :-- |
| RA3 | 14 | el más común en kits de sRNA |
| universal de Illumina | 2 | se usa el prefijo común de 21 nt |
| small-RNA de 2011 | 2 | sin esto daban 0% |
| ninguno | 1 | ya viene recortado |

**Sobre el prefijo de 21 nt:** TruSeq y el kit de sRNA de NEBNext comparten los
primeros 21 nucleótidos y divergen después. Si uno le pasa la versión larga del
kit equivocado, el recorte igual funciona — pero porque el recortador rechaza la
coincidencia larga por tener demasiadas diferencias (6 en 27 nt = 22%, sobre el
10% por defecto) y recién entonces acepta la corta. Eso depende del umbral de
error, que es un parámetro. **Con el prefijo común no hay ninguna diferencia
posible y no depende de nada.**

Dos reglas operativas:

- La tabla **la emite la herramienta**, no se transcribe. Pasar 19 filas de una
  tabla formateada a un fichero es una tarea donde un error no produce ningún
  síntoma inmediato.
- Un proyecto que no esté en la tabla **hace fallar el recorte**, que es mejor
  que recortar con una secuencia adivinada.

### La retención esperada no es el porcentaje de adaptador

El recortador descarta **dos cosas**: las lecturas sin adaptador, y las que
quedan fuera de la ventana 15-50 nt. La retención es la **intersección**, y la
brecha puede ser enorme.

| proyecto | con adaptador | retiene | por qué la brecha |
| :-- | --: | --: | :-- |
| galga PRJEB12164 | 95% | **52%** | el 20% de sus insertos mide 6-7 nt |
| gadmo PRJNA328800 | 84% | **51%** | el 28% tiene inserto de 10 nt |
| cloro PRJEB51338 | 62% | **57%** | poco adaptador de entrada |
| sclsc PRJNA1135930 | 98% | 97% | sin brecha: insertos de 22 nt |

Esos insertos de 6-10 nt son **dímeros de adaptador**: el adaptador 5' pegado al
3' sin inserto en el medio, un artefacto conocido de la preparación.

**Valor práctico:** sin este número, una retención del 51% se lee como "algo
falló". Con él, es lo esperado — y una del 2% es una alarma real.

### Cómo se saca una corrida

Editar el manifiesto a mano no sirve: la próxima vez que se regenera, la corrida
vuelve a aparecer, en silencio. Por eso hay un **fichero de exclusiones
versionado** que el generador aplica y reporta.

Una exclusión sólo se agrega **con evidencia medida**. La columna `motivo` dice
**qué se midió**, no qué se sospechó. Hoy tiene **una sola fila**.

Dos casos del proceso funcionando:

- **Un candidato a duplicado descartado:** 6 corridas, 758 M spots, etiquetado
  `miRNA-Seq`. El perfilado dio **1% de adaptador** con lecturas de 100 nt e
  insertos de 71-87 nt. Es mRNA. Los 126 M spots por corrida ya lo hacían
  sospechoso, y perfilarlo antes de adoptarlo evitó cambiar un duplicado que no
  servía por otro que tampoco.
- **Un duplicado reemplazado:** el original de *S. sclerotiorum* era RNA-Seq
  apareado y caía **entero** en el filtro, dejando al organismo sin validación.
  Se buscó un reemplazo con el mismo filtro, se perfiló, y recién entonces se
  adoptó. Con eso **los nueve organismos tienen primario y duplicado**.

### La escala resultante

**417 corridas · 7.90 G spots · 193 en primarios y 224 en duplicados · ~190 GB
de ficheros crudos respaldados.**

| organismo | corridas | spots |
| :-- | --: | --: |
| galga | 122 | 2 577 M |
| cloro | 57 | 1 582 M |
| maldo | 65 | 1 206 M |
| maggi | 57 | 809 M |
| phypa | 40 | 730 M |
| rhirr | 22 | 357 M |
| prupe | 15 | 287 M |
| sclsc | 20 | 216 M |
| gadmo | 19 | 133 M |

El set está **muy desbalanceado**, y es consecuencia de elegir proyectos
públicos reales en vez de diseñar un experimento. El ave tiene 122 corridas y el
pez 19.

Dos implicaciones. Para el cómputo: el proyecto más grande solo son ~23 horas de
alineamiento, y hay que procesarlo entero de una vez por una limitación del
alineador. Para el modelo: **la profundidad por organismo varía en un factor de
veinte**, y eso hay que tenerlo en cuenta al comparar cuántos loci se anotan en
cada uno — más profundidad detecta loci de menor expresión.

Contra la ronda anterior del proyecto —169 corridas, 2.18 G spots— es **3.4×
más**.

---

## 5 · Selección y verificación de genomas

**Se muestra:** el criterio; los tres retirados; la cepa de los datos; las
trampas del catálogo.

### El criterio

**Se respaldan, aunque sean públicos.** Rompe a propósito la regla general del
proyecto de no guardar lo que se puede volver a bajar, porque "se puede volver a
bajar" resultó ser **falso**: los ensamblados se retiran, se renumeran, se
reemplazan. Guardar el fichero exacto con su huella digital es la única forma de
poder decir, dentro de dos años, contra qué se alineó.

**Y se verifican recalculando.** El verificador recalcula la huella contra un
registro versionado, en vez de confiar en el tamaño del fichero.

El criterio de elección es **contigüidad y completitud**. Gracias al etiquetado
por secuencia, no hace falta que coincida con el ensamblado que usó ninguna base
pública.

| organismo | ensamblado | accession |
| :-- | :-- | :-- |
| rhirr | ASM2621079v1 | GCF_026210795.1 |
| sclsc | ASM14694v2 | GCF_000146945.2 |
| phypa | Phypa V5 | GCF_000002425.5 |
| cloro | C_rosea_IK726 | GCA_902827195.2 |
| prupe | Prunus_persica_NCBIv2 | GCF_000346465.2 |
| maldo | GDT2T_hap1 | GCF_042453785.1 |
| gadmo | gadMor3.0 | GCF_902167405.1 |
| galga | bGalGal1.mat.broiler.GRCg7b | GCF_016699485.2 |
| maggi | xbMagGiga1.1 | GCF_963853765.1 |

### Tres de los nueve estaban retirados

**Tres** ensamblados que el proyecto tenía elegidos estaban marcados como
retirados por el repositorio. **Dos son de los más citados de su especie.** Y
uno **estaba en uso** desde la ronda anterior.

**Descargarlos habría funcionado sin error visible.** El fichero se baja igual;
lo que cambió es que el repositorio lo retiró del catálogo vigente. Se retiran
por contaminación detectada después, por un ensamblado mejor del mismo material,
o por problemas de calidad encontrados tras publicarlo.

Qué se hizo: el resolvedor ahora **grita** cuando un accession está retirado, y
el descargador no baja nada que no esté verificado.

Dos detalles que importan:

- **El del molusco importa más porque es un organismo de entrenamiento.** El
  ensamblado viejo retuvo **haplotigos** —la misma región dos veces—. Si una
  región está dos veces en el genoma, toda lectura que venga de ahí mapea a dos
  sitios en vez de a uno, y eso cambia qué lecturas se colocan bajo el límite de
  multimapeo. El efecto contaminaría el **modelo**, no sólo la aplicación.
- **Los tres heredados de la ronda anterior no tenían accession registrado.** Era
  el hueco más serio de reproducibilidad: había resultados producidos contra un
  genoma del que no quedaba constancia de cuál era. Se resolvieron desde el
  nombre del ensamblado.

### Un organismo usa la cepa de los datos, no la referencia

El experimento primario de *C. rosea* son mutantes de genes de la maquinaria de
RNA de interferencia sobre una cepa concreta. Se eligió **el ensamblado de esa
cepa**, no la referencia de la especie.

Por qué: si alineáramos contra la referencia, las diferencias entre cepas
aparecerían mezcladas con el efecto de las mutaciones. La coincidencia de cepa
no se recupera de ninguna otra forma. El resolvedor marca ese ensamblado como
**diferente** de la referencia, y está bien que lo diga.

**Y hay una anomalía sin explicar:** ese ensamblado mide **70.7 Mb** contra una
mediana de **55.2 Mb** en las otras 20 cepas, y ~58 Mb que declara su propia
publicación.

En vez de ignorarlo o de descartar el ensamblado sin fundamento, quedó anotado
con una **comprobación concreta**: si es contenido duplicado, después de anotar
los loci van a aparecer **pares de loci distintos con el mismo RNA
mayoritario**. Es una predicción falsable y barata de verificar.

### Dos trampas del catálogo, y son la misma idea

**Un listado paginado que se trunca miente en silencio.** La consulta pedía 20
ensamblados y devolvió exactamente 20 —el límite— e imprimía el parcial como si
fuera todo. La especie tiene **78**, y el de la cepa correcta era **uno de los
58 que no se veían**. Se estuvo a un paso de elegir otra cepa por un defecto de
la herramienta, no por los datos. Ahora pagina hasta agotar y avisa si el total
declarado no coincide con lo que trajo.

**Un campo que no se lee se ve igual que un campo vacío.** La cepa se leía de
**un solo lugar** del registro. Cuando estaba en otro, la columna salía vacía y
se leía como "este ensamblado no declara aislado". Ahora cae en cascada por los
tres lugares posibles y, cuando de verdad no hay, imprime `cepa=?`.

Las dos son **una limitación de la consulta presentada como un hecho sobre el
mundo**.

Y una tercera, de nomenclatura: **el género del molusco se renombró**. Buena
parte de las bases todavía usa el nombre viejo. Buscar el genoma por el nombre
nuevo no lo encuentra.

---

## 6 · Las herramientas, y por qué YASMA

**Se muestra:** el stack completo; las tres razones; el contrafáctico.

### El stack

| herramienta | para qué | por qué esa |
| :-- | :-- | :-- |
| **YASMA v1.1.1** | recorte, alineamiento y **anotación de loci** | la central |
| bowtie 1 | el alineador que YASMA invoca | end-to-end; el estándar para 15-50 nt |
| cutadapt | el recorte real, dentro de `yasma trim` | lo elige YASMA, no nosotros |
| ViennaRNA | dependencia de importación de YASMA | sin ella no arranca ningún subcomando |
| SRA Toolkit | descarga y validación de las 417 corridas | es la vía oficial del repositorio |
| pysam · samtools | lectura y ordenamiento del resultado | dentro de YASMA |
| rclone | respaldo del dato pesado | con credencial propia, no la compartida |
| micromamba | el entorno, fijado por versión | YASMA exige Python ≥ 3.12 |

Ocho piezas, y tres de ellas las usa YASMA por dentro. Todo lo demás son scripts
propios del proyecto.

Dos elecciones que conviene poder defender:

- **bowtie 1 y no bowtie 2 ni BWA.** Para lecturas de 15 a 50 nucleótidos,
  bowtie 1 alineando de extremo a extremo es el estándar del campo y es lo que
  usan las herramientas de referencia de RNA pequeño. Los alineadores modernos
  están diseñados para lecturas largas y su estrategia de semillas no tiene
  sentido a este tamaño.
- **rclone con credencial propia.** El cliente por defecto de la herramienta
  está compartido por todos sus usuarios y la cuota del servicio se aplica por
  cliente, así que está permanentemente saturado. Con 190 GB que mover, la
  diferencia es entre horas y días.

### Por qué YASMA

**1 · Cubre tres de las cuatro etapas.** Recorte, alineamiento y anotación de
loci con una sola herramienta y un solo formato intermedio. No hay pegamento
entre etapas, y cada conversión de formato entre herramientas es un lugar donde
se pierde información sin avisar.

**2 · Anota loci de novo, que es la parte difícil.** Alinear lecturas cortas es
un problema resuelto; hay varias herramientas buenas. **Decidir que una
acumulación de lecturas constituye un locus, dónde empieza, dónde termina, y si
es uno o dos loci solapados, es un problema abierto** con criterios que dependen
del tipo de RNA. YASMA implementa ese paso, y es el paso que produce las filas
de nuestra matriz.

**3 · Está hecha para sRNA.** No es un pipeline genérico adaptado. Sus
decisiones por defecto —ventana de tamaño, límite de multimapeo, pesado de
multimapeadores— ya son las del campo.

**Y la razón práctica: es mucho más simple de usar.** Se instala con una línea,
se corre con tres comandos, y **el resultado del alineamiento es exactamente el
que consume su propio anotador**. Sin conversiones, sin reordenar, sin adaptar
cabeceras.

### Qué costaría no usarla

| etapa | con YASMA | sin YASMA |
| :-- | :-- | :-- |
| recorte | un comando | llamar al recortador y decidir todos los parámetros |
| alineamiento | un comando, con pesado de multimapeadores incluido | **implementar el pesado por cobertura única local** |
| grupos de lectura | los escribe solo, uno por corrida | un paso extra con otra herramienta, y acordarse |
| anotación de loci | un comando | **el proyecto entero** |

El pesado de multimapeadores **no es opcional** en sRNA: los sRNA que buscamos
en hongos derivan de retrotransposones, o sea que son multicopia por definición.
Reimplementarlo bien sería un trabajo en sí mismo.

### La tensión honesta

**Es mucho más simple de usar que armar el pipeline a mano, y al mismo tiempo la
mayor parte del trabajo de esta etapa fue caracterizarla.** Las dos cosas son
ciertas y no se contradicen.

Es simple de usar porque la interfaz es de tres comandos y hace lo correcto por
defecto. Hubo que caracterizarla porque es software académico reciente, su
documentación está incompleta, y varias de sus conductas importantes —cómo lleva
el registro de lo hecho, qué pasa cuando no puede detectar un adaptador, si un
proyecto se puede procesar por partes— **sólo se descubren midiendo**.

¿Eso la descalifica? No. Cualquier herramienta de este campo tiene conductas no
documentadas. **La diferencia es haberlas medido y dejarlas escritas.**

---

## 7 · YASMA por dentro

**Se muestra:** fijar la versión; el recorte; por qué no usamos su detector de
adaptador; el alineamiento; el límite que no se puede partir; la anotación.

### Fijar la versión: el nombre engaña dos veces

- **El branch por defecto no es la release.** Un clon pelado deja una versión
  anterior, no la que el proyecto usa.
- **Y el tag declara otro número.** El fichero de metadatos **dentro** del tag
  correcto declara la versión anterior.

Consecuencia: **la versión no se puede verificar desde el paquete instalado**.
Si escribiéramos en la tesis el número que reporta el instalador, estaríamos
declarando la versión equivocada de buena fe. Se fija por **referencia de git**:

```
git+https://github.com/NateyJay/YASMA@v1.1.1
```

Y hay un verificador de entorno que comprueba que cada herramienta **arranque**
de verdad, no que esté instalada. YASMA importa en su inicialización un módulo
de estructura secundaria de RNA que a su vez importa una biblioteca externa. Si
esa biblioteca falta, **ningún subcomando arranca**, ni los que no tienen nada
que ver con estructura — y el error que tira no la menciona de forma obvia.

### El recorte

Es un envoltorio del recortador. La llamada, literal:

```
cutadapt -a <sec> --minimum-length 15 --maximum-length 50 \
         -j <cores> -O 4 --max-n 0 --trimmed-only
```

Cuatro consecuencias que van en métodos:

- **`--trimmed-only` descarta las lecturas sin adaptador.** La retención
  esperada **no** es ~100%: es el porcentaje de adaptador medido. Para sRNA-seq
  es lo correcto —una lectura sin adaptador tiene el inserto más largo que la
  lectura, o sea que no es un sRNA— pero hay que declararlo.
- **`--max-n 0`** tira cualquier lectura con una sola base indeterminada.
- **No hay filtro de calidad**: se llama sin `-q`. Eso cierra la mitad del
  asunto de las dos corridas en formato reducido, que guardan una única calidad
  sintética para todas las bases: **no hay calidad que distorsionar**.
- **Una librería pre-recortada no pasa por la ventana 15-50.** Se saltea entera.
  Para el proyecto que tenemos da igual —sus lecturas de 30-34 nt están dentro—
  pero no es una regla general. Donde sí se filtran es en el alineamiento.

*Sobre la ventana 15-50:* es más ancha que la típica de miRNA porque queremos
conservar fragmentos de RNA de transferencia, que miden 30-40 nt, y piRNAs, que
miden 24-32. Los valores por defecto coinciden con la ventana que el proyecto
eligió por ese motivo — **coincidencia afortunada, no diseño compartido**. Se
pasan explícitos igual: un default que cambia de versión no avisa.

### Tres cosas del recorte que rompen el trabajo por tandas

Las tres se encontraron **corriendo el binario y mirando qué dejaba**, no
leyendo documentación:

1. **Pisa su propio registro.** Lo vacía al entrar y lo reescribe sólo con lo de
   esta llamada al salir. **La tanda 2 borra el registro de la tanda 1.** El
   síntoma: el pipeline reporta la corrida como faltante, la vuelve a volcar y
   recortar, **para siempre**. El fichero sigue en disco y desapareció del
   registro.
2. **Trunca su propio log.** Se abre en modo escritura: las estadísticas del
   recortador de las tandas previas se pierden.
3. **La limpieza automática no se puede usar.** Recorre un campo que en nuestro
   caso está vacío y revienta. Y borraría la salida de una librería
   pre-recortada.

La solución: un **registro propio acumulativo**, y reescribir el de la
herramienta con la lista completa después de cada tanda.

**Y las tandas hacen falta:** volcar las 417 corridas a texto sin comprimir de
una vez son **~1.1 TB**, y un solo organismo son 376 GB. No entra en ningún
disco razonable junto con los datos crudos y los resultados. Se procesa por
tandas acotadas por un presupuesto de disco, borrando el texto sin recortar
apenas la tanda termina — que es exactamente lo que el primer bug rompía.

### Por qué no usamos su detector de adaptador

| librería de prueba | largo modal | YASMA dice | el nuestro dice |
| :-- | --: | :-- | :-- |
| recortada, largos variables | 22 nt | PRE-TRIMMED | YA RECORTADA |
| recortada a **un solo largo** | 22 nt | **None** | YA RECORTADA |
| mRNA (no es sRNA-seq) | 150 nt | **None** | NO PARECE |
| sRNA-seq sin recortar | 150 nt | el adaptador | PARECE sRNA-seq |

Su criterio usa la **dispersión** del largo de lectura, no su **magnitud**. Por
eso da el mismo veredicto —`None`— para una librería ya recortada y para mRNA:
**situaciones opuestas**.

Acierta en nuestros 19 proyectos, pero **por suerte, no por diseño**: el único
pre-recortado tiene largos variables, que es justo lo que su criterio necesita.
Si un reemplazo llegara recortado a largo fijo, lo llamaría "sin adaptador".

**Y ahí está el peligro real:** el recorte, ante un adaptador nulo, **descarta la
librería**. No la suma al registro y desaparece del pipeline **sin ningún
error**. Una librería mal clasificada se perdería en silencio.

Por eso el adaptador se mide una vez, se versiona en git, y se le pasa.

### El alineamiento: no envuelve a ShortStack

**Lo que creíamos:** que el comando de alineamiento llamaba a ShortStack con una
opción que nos costaría el límite de multimapeo que habíamos justificado con
datos. Íbamos a **escribir nuestro propio alineamiento** para no perderlo.

**Lo que dice el código:** el envoltorio de ShortStack existe pero está
**comentado**. El comando que se registra es **bowtie 1 nativo** con el pesado
estilo ShortStack3. Y el límite que queríamos **es el valor por defecto**.

*Sobre el límite de 50 sitios:* se midió en un organismo que después salió del
set. Descartaba el 76% de las lecturas, pero **el 98.8% de lo que se recuperaría
al subirlo eran fragmentos de RNA de transferencia de 34 nt exactos, en cientos
de copias génicas**. Subirlo inunda la anotación con una sola especie
repetitiva.

### Tres pasadas por librería

1. **únicas** — lo que mapea a un solo sitio entra al resultado; el resto se
   aparta.
2. **múltiples** — con todas las posiciones en mano, pesa cada una por la
   **cobertura única local** y elige por **sorteo ponderado**.
3. **excedidas** — las que pasaron el límite **no se tiran**: entran al
   resultado como no mapeadas, **etiquetadas**.

| etiqueta | qué significa |
| :-- | :-- |
| U | única — mapeó a un solo sitio |
| P | multimapeada, **colocada con peso** |
| Q | alineó, pero por encima del máximo de sorteo: **no se colocó** |
| H | por encima del límite de 50 sitios |
| N | **sin ningún alineamiento válido** |
| F | descartada por el filtro de tamaño o de bases indeterminadas |

**La idea del pesado por cobertura única local** es la contribución central de
ShortStack. Una lectura que mapea a diez sitios no se reparte entre los diez
—inflaría la cobertura— ni se tira —perdería señal real—. Se coloca en **uno
solo**, sorteado con probabilidad proporcional a cuánta cobertura de lecturas
**únicas** hay alrededor de cada sitio candidato. La intuición: si un locus ya
tiene lecturas que sólo pueden venir de ahí, es más probable que las ambiguas
también vengan de ahí.

Que las excedidas no se tiren es una decisión buena: **el conteo cierra siempre
y todo se puede auditar**.

### Un proyecto no se puede partir

El diccionario de cobertura única se crea **una sola vez** y se acumula sobre
**todas** las librerías del proyecto **antes** de que la etapa de multimapeo lo
use para pesar. El peso de cada posición sale de la cobertura única **agrupada
de todo el proyecto**.

**Partirlo no cambia la contabilidad: cambia a qué locus se asigna cada lectura
multimapeada.** Es un resultado distinto, no un resultado parcial.

Consecuencia: la unidad reanudable es el **proyecto entero**, no la corrida. Una
sesión de cómputo que se muere pierde un proyecto completo. Por eso se procesa
de a un proyecto, **del más chico al más grande**, y el resultado se respalda
apenas termina — así una sesión interrumpida dejó proyectos terminados en vez de
nada.

La tentación natural, con un proyecto de 95 corridas y sesiones que se cortan,
es procesar de a pedazos. Eso daría un resultado que parece correcto, con los
mismos conteos totales, pero con las multimapeadas repartidas distinto. **Y las
multimapeadas no son un detalle acá.**

### Cuatro trampas de rutas

El genoma tiene que estar **adentro** del directorio de salida: uno compartido
afuera **no da un mensaje**, revienta desde dentro de la librería, porque la
validación hace una operación de ruta relativa sin protegerla.

La solución ingenua —copiar el genoma a cada proyecto— funciona y es mala:
**duplica el índice**. Como cada organismo tiene dos proyectos, construiría el
índice dos veces; en un genoma de 1 Gb son horas por nada.

La correcta es un **enlace simbólico al directorio** del organismo. La
validación interna usa la ruta absoluta **sin resolver los enlaces**, así que el
enlace queda lexicalmente dentro aunque apunte afuera. Y tiene que ser **al
directorio, no al fichero**: el índice se construye **al lado del genoma que le
pasaste**.

**Y el genoma no puede ir comprimido:** se abre con una librería que no lee un
comprimido común. Medido: sin comprimir funciona, comprimido tira un error de
apertura. Se **verifica la huella del comprimido** contra el registro versionado
y **recién entonces** se descomprime. El descomprimido y el índice son derivados
y no se respaldan.

### La anotación

El resultado tiene que traer un **grupo de lectura por corrida**. Sin ellos, el
anotador **corta** con un error de clave faltante: la consulta a la cabecera se
hace sin protección.

No es cosmético: **agrega profundidad por grupo**. Sin grupos no podría separar
librerías aunque no se cayera — y separar librerías es lo que permite decir si
un locus aparece en un experimento o en varios.

Lo bueno: el alineamiento ya los escribe, uno por corrida, en la cabecera y en
cada lectura.

**El caso peligroso no es que falten todos** —eso se cae ruidosamente— **sino
que falte uno**. Ahí no se cae, anota normalmente, y el resultado es un análisis
al que le falta una librería entera sin que nada lo diga. Por eso el verificador
cruza la lista de corridas del manifiesto contra los grupos presentes.

---

## 8 · Los problemas

**Se muestra:** el fallo silencioso; el catálogo; los verificadores; el testing
por mutación.

### El hilo conductor

**Casi ningún error de este pipeline se manifiesta como un error.**

| si te equivocás en... | lo que pasa | lo único que lo delata |
| :-- | :-- | :-- |
| la secuencia del adaptador | código 0, salida **casi vacía** | la retención medida |
| el genoma | código 0, **0% alineado**, resultado escrito | la fracción alineada |
| el adaptador de **una** corrida | la librería **se descarta** sin sumarse al registro | la corrida que falta |
| el grupo de lectura de una | anota igual, **sin esa librería** | el manifiesto contra el resultado |
| el ensamblado retirado | baja y alinea **normalmente** | el catálogo, consultado a tiempo |

Ninguno tira una excepción. Ninguno escribe una advertencia. **Los cinco
producen un resultado plausible obtenido de la forma equivocada** — que es el
peor tipo de resultado que un pipeline puede dar.

Un pipeline que se cae es un pipeline fácil: leés el error y lo arreglás. Un
pipeline que produce un número plausible por el camino equivocado **te hace
escribir una tesis sobre ruido**.

El del genoma se comprobó a propósito: se corrió el alineador real contra un
genoma del mismo tamaño pero de secuencia aleatoria. Salió con **código cero**,
escribió el resultado, y el estado del pipeline decía "dos de dos proyectos al
día".

### El catálogo de problemas

| problema | cómo se vio | qué costó |
| :-- | :-- | :-- |
| Una corrida etiquetada mal era mRNA | perfilado: 0 de 40 000 con adaptador | se excluyó con evidencia |
| 4 proyectos marcados como sospechosos | eran defectos de la herramienta | casi 34 corridas buenas |
| 3 ensamblados retirados | consulta al catálogo antes de bajar | reproducibilidad, si no |
| Un listado que se truncaba en silencio | 20 de 78 devueltos, sin avisar | casi la cepa equivocada |
| El registro del recorte se pisaba | medido corriendo dos tandas | bucle infinito |
| Un proyecto mezcla **dos kits** | 6 de 12 corridas vacías | el verificador lo paró |
| El 67% de una librería no alinea | dos columnas en vez de una | **pendiente de explicar** |
| Una estimación de disco, 3× de más | el primer proyecto real | 2 proyectos mal descartados |

**Ninguno se manifestó como un error.** Todos se encontraron midiendo, y cada
uno dejó un chequeo automático detrás.

### Los verificadores

- **Tras el recorte:** compara la retención **medida** contra la esperada. Grita
  `VACÍA` por debajo del 5% y `DESVIADA` más allá de 15 puntos.
- **Tras el alineamiento:** lee los conteos por grupo. Grita `MUY BAJA` por
  debajo del 10% colocado, reporta la corrida que no llegó, y avisa si mucho se
  pasa del límite de multimapeo.
- **Sobre los documentos:** cruza lo que afirman los documentos contra los
  ficheros de datos, y **falla** si un documento dice algo que los datos no
  dicen. Encontró dos errores que ninguna lectura había visto.

Y antes de respaldar nada: **verificar antes de subir**. Un resultado alineado
contra el genoma equivocado se ve igual que uno bueno una vez guardado.

**584 chequeos automáticos · 19 bancos de prueba · 89 mutaciones.**

Ninguno de estos chequeos se agregó por precaución abstracta. **Cada uno existe
porque una cosa concreta falló en silencio.**

### Un banco de pruebas que pasa no prueba nada

Un banco verde puede estar verde porque el código está bien, **o porque el banco
no mira lo que importa**. Desde afuera no se distinguen.

La prueba: se **rompe el código a propósito, 89 veces**, y se exige que algún
banco se ponga rojo en cada una.

Así aparecieron **dos huecos** que ninguna otra cosa mostró. El caso concreto es
didáctico: el verificador imprimía un **mensaje** para la persona y por separado
calculaba una **variable** que decidía el resultado. El banco miraba el mensaje.
Se podía romper la decisión —reintroduciendo el bug que casi tira 34 corridas
buenas— y los chequeos seguían **todos en verde**.

**Regla que quedó: al agregar un chequeo, afirmar sobre el valor que el programa
usa, no sobre el texto que imprime al lado.**

Es el mismo principio que aplicamos a las herramientas externas, aplicado al
propio código: **no confiar en lo que un programa dice de sí mismo, medir lo que
hace**.

---

## 9 · Resultados

**Se muestra:** el estado; el primer alineamiento; el 67% sin alinear; el kit
mixto.

### Lo que está cerrado

- **417 corridas** descargadas, validadas y con huella versionada en git.
  Reconciliadas contra el manifiesto **en las dos direcciones**.
- **9 ensamblados** fijados con huella verificada. Tres estaban retirados y hubo
  que reemplazarlos.
- **19 proyectos** perfilados: los 19 son sRNA-seq, y sabemos con qué adaptador
  recortar cada uno.

**Empezando:** el alineamiento, con **1 proyecto** hecho como calibración.

**Por delante:** los otros 17, la anotación de loci, la extracción de
características y el modelo. **~100 h** de alineamiento estimadas, con un número
medido.

La extracción y el entrenamiento son **órdenes de magnitud más baratos** que el
alineamiento: el cuello de botella es esta etapa.

### El primer alineamiento

Se corrió el proyecto más chico como calibración: 32 millones de lecturas,
genoma de 39 Mb.

**La estimación de tiempo, bien:** 26 min de alineador más 4 de ordenamiento.
De ahí salen las ~100 h del total.

**La de espacio, mal por 3×:**

| | decía | midió |
| :-- | --: | --: |
| bytes por alineamiento | 45 | **14.1** |
| pico de disco del proyecto más grande | ~174 GB | **~79 GB** |
| total de resultados del proyecto | ~340 GB | **~100 GB** |
| disco libre de la máquina de cómputo | ~78 GB | **220 GB** |

Estaba escrito en cuatro documentos como si fuera un hecho. La consecuencia
práctica: **dos proyectos estaban declarados como "no entran en la máquina de
cómputo" y entran perfectamente**.

*Detalle honesto:* el 14.1 se midió en una librería donde el 67% de las lecturas
no alinea, y una lectura sin alinear ocupa menos que una colocada. Por eso el
pipeline usa **16 y no 14**. Es un margen declarado, no un número redondeado a
ojo.

**Un número inventado no falla: manda el trabajo al lugar equivocado.**

### El 67% que no alinea

| etiqueta | qué es | % |
| :-- | :-- | --: |
| U | única | 4.1 |
| P | multimapeada, colocada con peso | 12.9 |
| Q | alineó pero no se colocó | 15.6 |
| N | **sin ningún alineamiento válido** | **67.1** |

**Y pasó el umbral automático:** el único criterio miraba la fracción colocada,
y 17.2% está por encima del 10% que se exigía. **Pasar un umbral no es lo mismo
que estar bien.**

- **No es el genoma:** contra un ensamblado ajeno da ~0%, no 33%.
- **No es el recorte:** 99% de retención y **cero** lecturas descartadas por
  tamaño.
- **La hipótesis:** que el experimento sea sobre **tejido vegetal infectado**. Es
  un hongo necrótrofo, su sRNA se estudia sobre todo por RNAi entre reinos, y una
  librería de tejido infectado sería mayoritariamente del huésped.

**El arreglo es conceptual, no cosmético.** Había una sola columna —"fracción
alineada"— que colapsaba dos problemas **opuestos**:

- Lo que **no alinea en ninguna parte** son lecturas que no son de este genoma:
  ensamblado equivocado, contaminación, o el huésped.
- Lo que alinea **de más**, por encima del límite de multimapeo, es un genoma
  repetitivo.

Con una sola columna los dos casos se leen igual: "poco alineado". **Ahora hay
dos.**

Por qué no se afirma ya que es sobre tejido infectado: **porque no está
medido**. La comprobación es barata —alinear el experimento de cultivo puro del
mismo hongo y comparar— y hasta que no esté hecha, la hipótesis se nombra como
hipótesis.

*Dato que sí es esperable:* el multimapeado **triplica** al único. Los sRNA
efectores descritos en este hongo derivan de retrotransposones, que son
multicopia por definición. **Eso es señal, no ruido.**

### Un proyecto público puede mezclar kits

- **6 de 12 corridas** retuvieron **0.6% – 2.0%** donde se esperaba 51%.
- **Las otras 6** retuvieron **38% – 58%**. Normales.
- **Todas** del mismo BioProject, descargadas juntas, etiquetadas igual.

El perfilado medía **una** corrida por proyecto y extendía esa conclusión a las
demás. Le tocó una de las buenas, así que **las otras once nunca se midieron**.

El supuesto estaba escrito en un comentario del código: *"se busca por proyecto
porque el adaptador es del kit, no de la corrida"*. **El kit sí es del proyecto.
Lo que no es cierto es que un proyecto público use un solo kit.**

**Lo que funcionó:** el verificador las marcó vacías y el pipeline **se
detuvo**. Nada llegó al alineamiento. Si no existiera, el resultado habría sido
un alineamiento con seis librerías casi vacías, la verificación diría que está
bien —porque la fracción alineada de lo poco que quedó es normal— y la anotación
de ese organismo se haría **sobre la mitad de la profundidad real** sin que nada
lo indicara.

**Lo que se agregó:** excepción **por corrida** en la tabla de adaptadores, un
modo que mide **todas** las corridas una por una, y una forma de **rehacer** sólo
las malas.

Y una cosa más: el paso siguiente eran **seis identificadores transcritos a
mano** desde una tabla. Eso también se automatizó, porque transcribir a mano es
exactamente donde se cuela el próximo error silencioso.

---

## 10 · Qué puede salir, y cierre

**Se muestra:** los desenlaces del upstream; los tres desenlaces del modelo; el
cierre.

### Lo abierto del upstream

| pregunta abierta | si sale que sí | si sale que no |
| :-- | :-- | :-- |
| ¿El 67% sin alinear es el experimento? | genoma y recorte confirmados; ese experimento de validación tiene **un tercio** de la profundidad nominal | hay que revisar el ensamblado o algo aguas arriba |
| ¿Las 6 corridas son de otro kit? | se recortan con su secuencia y el proyecto se recupera entero | si vienen ya recortadas, se marcan y se saltean |
| ¿El ensamblado anómalo tiene contenido duplicado? | aparecerían **loci duplicados con el mismo RNA mayoritario** | la anomalía queda sin explicar pero sin efecto |
| ¿Cuántos loci por organismo? | fija la escala del problema de clasificación | — |

Los cuatro se contestan **con los datos que ya tenemos**. Dos de ellos, con una
corrida de alineamiento de un par de horas.

**No son riesgos vagos: son preguntas acotadas con plan de resolución.**

### Los tres desenlaces del modelo

**1 · La transferencia funciona.** El modelo entrenado en los tres animales
curados recupera sRNAs ya descritos en plantas y hongos, y prioriza candidatos
nuevos que **reaparecen** en el experimento independiente. Es el resultado
buscado.

**2 · La transferencia falla.** No recupera casi nada fuera de animales. **Es un
resultado interpretable**, no una lista vacía ambigua: sabríamos que falló
porque la firma animal no transfiere, y no porque no haya nada que encontrar.

**3 · Funciona a medias.** Recupera en plantas pero no en hongos, o al revés.
Acota **dónde** está el límite de la transferencia entre reinos, que es una
contribución por sí misma.

**Lo que hace que los tres sean interpretables son las dos salvaguardas:** las
anotaciones públicas como evaluación y nunca como entrenamiento, y la validación
cruzada entre los tres curados **antes** de correr los nueve.

**Un resultado negativo también es un resultado.**

### Cierre

Cuatro cosas para llevarse:

1. **El diseño sale del método.** Nueve organismos, dos experimentos cada uno, y
   la partición entrenamiento/aplicación salen del supuesto de PU, no de lo que
   había disponible.
2. **El upstream está instrumentado.** No sólo corre: **mide qué hace** y se
   detiene cuando no coincide con lo esperado.
3. **YASMA está caracterizada.** Qué consume, qué decide y **dónde no hay que
   confiarle**, medido contra el código de la versión fijada.
4. **Los hallazgos son del tipo que no avisa.** Los modos de fallo encontrados
   salen con **código cero**. Ninguno se habría notado leyendo la salida.

Todo sale de una sola observación: **en este pipeline los errores no se
manifiestan como errores**. De ahí salen los verificadores, de ahí sale medir
las herramientas en vez de creerles, y de ahí sale romper el propio código a
propósito para saber si los chequeos sirven.

**La regla que ordena todo el trabajo: el veredicto de una herramienta no
reemplaza leer lo que midió** — sobre todo la primera vez que corre sobre datos
reales.

Y el aporte de la tesis es el planteo positive-unlabeled sobre loci anotados de
novo en nueve organismos de tres reinos, con un protocolo de etiquetado que no
depende del ensamblado y una validación por experimento independiente. **El
upstream es el medio.**

---

## Preguntas probables, con respuesta

**¿Esto no es simplemente correr herramientas?**
En parte sí, y el trabajo está en elegirlas con criterio, **medir qué hacen de
verdad**, y construir los chequeos que hacen falta porque ninguna avisa cuando
se equivoca. Las ocho entradas del catálogo de problemas son eso.

**¿Por qué no diez o veinte organismos?**
Porque el cuello de botella no es el cómputo, es la verificación. Cada
BioProject hay que perfilarlo y cada genoma resolverlo. Nueve es lo que entra
con verificación seria.

**¿Y si la transferencia entre reinos falla?**
Se reporta como resultado negativo interpretable gracias a las dos
salvaguardas, no como lista vacía ambigua. Y la validación dejando un organismo
afuera da la señal **antes** de gastar cien horas de cómputo.

**¿Por qué no hay candidatos validados experimentalmente todavía?**
Porque el upstream recién se está cerrando, y **la priorización es justamente lo
que define qué valdría la pena validar**.

**¿Cómo se elegirían los candidatos a validar?**
Por ranking, con el criterio de **reaparición en el experimento independiente**.

**¿Por qué bowtie 1 y no un alineador moderno?**
Para lecturas de 15-50 nt, bowtie 1 end-to-end es el estándar del campo y es lo
que usan las herramientas de referencia de sRNA. Los modernos están diseñados
para lecturas largas y su estrategia de semillas no tiene sentido a este tamaño.

**¿El límite de multimapeo no descarta demasiado?**
Descarta el 76% de las lecturas en el organismo donde se midió, y el 98.8% de lo
que se recuperaría al subirlo son fragmentos de RNA de transferencia de un largo
exacto en cientos de copias génicas. **Subirlo inunda la anotación con una sola
especie repetitiva.**

**¿YASMA es confiable si tiene tantas conductas no documentadas?**
Cualquier herramienta de este campo tiene conductas no documentadas. **La
diferencia es haberlas medido y dejarlas escritas**, que es lo que permite
declararlas en métodos en vez de descubrirlas en la defensa.
