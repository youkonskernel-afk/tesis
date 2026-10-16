#!/usr/bin/env bash
# Recorte de adaptador con YASMA (que envuelve a cutadapt).
#
#   ./scripts/trim.sh plan      [<org>[/<rol>]]   que se haria, sin hacerlo
#   ./scripts/trim.sh correr    [<org>[/<rol>]]   recorta, en tandas
#   ./scripts/trim.sh estado    [<org>[/<rol>]]   que esta recortado y que falta
#   ./scripts/trim.sh verificar [<org>[/<rol>]]   la retencion medida vs la esperada
#   ./scripts/trim.sh rehacer <org>/<rol> RUN...  saca esas del registro y las re-recorta
#   ./scripts/trim.sh muestra <org>/<rol> [N] [RUN...]
#                                    N lecturas de cada corrida, recortadas igual que
#                                    yasma trim, en $PROY_DIR/_muestras/<org>_<rol>/
#
# El filtro es un organismo (`galga`) o un organismo y un rol (`galga/primario`).
#
# Variables: SRA_DEST, PROY_DIR, ADAPTADORES_TSV, EXTREMOS_TSV, MANIFEST, CORES,
# PRESUPUESTO_GB, SOLAPAR.
#
# LA UNIDAD DE TRABAJO ES <org>_<rol>, NO <org>
# Un proyecto YASMA por organismo Y rol: trim/galga_primario/, trim/galga_duplicado/.
# El primario es el set de descubrimiento y el duplicado la validacion
# independiente, y esa separacion hay que respetarla o la validacion deja de ser
# independiente. Mezclarlos en un solo directorio deja la separacion en "hay que
# acordarse de filtrar por la columna rol del manifiesto" en vez de en el disco.
# Aguas abajo (`yasma tradeoff`) el directorio de proyecto ES la unidad, asi que
# separarlos aca es lo que hace que anotar el primario y el duplicado por
# separado salga solo.
#
# `maggi_primario` junta PRJNA154615 y PRJNA232734, que son de 2011 y 2014 y
# llevan adaptadores distintos. No es un problema: el `adapters` de inputs.json
# es un dict POR FICHERO, no uno por proyecto.
#
# POR QUE NO SE USA `yasma adapter`
# `yasma trim` lee los adaptadores de inputs.json, que normalmente escribe
# `yasma adapter`. Dos razones para no dejarselo a el:
#   1. Ante un adaptador "None", `yasma trim` DESCARTA la libreria — no la suma
#      a trimmed_libraries y desaparece del pipeline sin ningun error. Y
#      `yasma adapter` da "None" tanto para una libreria ya recortada a largo
#      fijo como para mRNA, que son situaciones opuestas (medido: docs/yasma.md).
#   2. El adaptador es del kit, o sea del proyecto. Redetectarlo corrida por
#      corrida es tirar una decision ya tomada y medida.
# Asi que el adaptador sale de data/adaptadores.tsv y este script escribe
# inputs.json. Un proyecto que no este en esa tabla hace fallar el script.
#
# LO QUE YASMA LE PASA A CUTADAPT, Y SUS CONSECUENCIAS
#   cutadapt -a <sec> --minimum-length 15 --maximum-length 50 -j <cores> -O 4
#            --max-n 0 --trimmed-only
#   - 15-50 son los defaults de yasma trim y coinciden con la ventana del
#     proyecto. No hay que pasarlos, pero se pasan explicitos igual: un default
#     que cambia en una version nueva no avisa.
#   - `--trimmed-only` DESCARTA los reads sin adaptador. Es distinto de fastp,
#     que los conservaria. Para sRNA-seq es lo correcto —un read sin adaptador
#     tiene el inserto mas largo que el read, o sea que no es un sRNA— pero hay
#     que declararlo, porque la retencion va a ser el adapt_pct de la tabla.
#   - `--max-n 0` tira cualquier read con una N.
#   - NO hay filtro de calidad: cutadapt se llama sin -q. De ahi que las dos
#     corridas SRA Lite, con su calidad sintetica unica, no distorsionen nada en
#     este paso.
#
# TRES COSAS DE `yasma trim` v1.1.1 MEDIDAS CONTRA EL BINARIO, NO LEIDAS
#   1. PISA `trimmed_libraries` EN VEZ DE ACUMULARLO. Hace
#      `ic.inputs['trimmed_libraries'] = []` al entrar y `= <lo de esta llamada>`
#      al salir. O sea que la tanda 2 BORRA el registro de la tanda 1: el
#      .t.fq.gz sigue en disco pero desaparece de inputs.json, `estado` reporta
#      la corrida como faltante y `correr` la vuelve a volcar y recortar. Por eso
#      existe el ledger `recortadas.tsv`, que es acumulativo, y por eso despues
#      de cada tanda se vuelve a escribir inputs.json con la lista completa —
#      que es lo que los comandos YASMA de aguas abajo van a leer.
#   2. `trim/log.txt` SE TRUNCA EN CADA LLAMADA: `Logger` abre el fichero con
#      "w". Las estadisticas de cutadapt de las tandas anteriores se pierden, asi
#      que cada tanda se loguea aparte en `logs/tanda_NN.log` y de ahi salen los
#      conteos de `verificar`.
#   3. `--cleanup` NO SE PUEDE USAR: itera `ic.inputs['srrs']`, que en nuestro
#      inputs.json es None -> TypeError. Ademas vacia `untrimmed_libraries`, que
#      para una libreria PRE-TRIMMED es justamente donde esta la salida. El
#      borrado del fastq sin recortar lo hace este script.
#
# KITS 4N: LO UNICO QUE NO RECORTA `yasma trim`
# rhirr_duplicado alineo 0.17% y cortando 4 nt de cada punta, 82.07%: adaptadores
# con bases aleatorias en las dos uniones (NEXTflex y los "HD"), que quedan
# pegadas al inserto. `yasma trim` le pasa a cutadapt `-a` y la ventana, nada
# mas, asi que las corridas con fila en data/extremos.tsv (la escribe
# `align.sh prueba`, ver scripts/extremos.py) no van por yasma: van por la MISMA
# llamada a cutadapt, en dos pasadas:
#   cutadapt -u <5'> -a <sec> -j <cores> -O 4 --max-n 0 --trimmed-only
#     | cutadapt -u -<3'> --minimum-length 15 --maximum-length 50 -j <cores>
# Dos porque cutadapt aplica -u ANTES de buscar el adaptador, y las bases del 3'
# estan entre el inserto y el adaptador. La ventana va al final, sobre el inserto
# ya cortado. Lo cortado queda en la columna `extremos` del ledger, y lo
# recortado con otro corte que el que pide la tabla HOY se rehace solo: la fila
# llega, por construccion, despues de recortar y ver que no alinea.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${MANIFEST:-$ROOT/data/srr_manifest.tsv}"
ADAPTADORES_TSV="${ADAPTADORES_TSV:-$ROOT/data/adaptadores.tsv}"
EXTREMOS_TSV="${EXTREMOS_TSV:-$ROOT/data/extremos.tsv}"
EXTREMOS_PY="$ROOT/scripts/extremos.py"
# shellcheck source=_drive_lib.sh
. "$ROOT/scripts/_drive_lib.sh"
# El mismo directorio donde drive_pull.sh deja los .sra, no uno propio.
SRA_DEST="${SRA_DEST:-$(ruta_local sra)}"
# Las lecturas recortadas NO se respaldan: se re-generan de forma determinista
# desde los .sra y data/adaptadores.tsv, que si estan respaldados. Por eso no
# hay fase 'trim' en el mapa de Drive. El directorio es el PROYECTO YASMA
# completo —trim/, align/, annotations/— y la ruta sale de _drive_lib.sh, que es
# de donde la lee align.sh.
PROY_DIR="${PROY_DIR:-$(ruta_proyectos)}"
CORES="${CORES:-$(nproc 2>/dev/null || echo 1)}"
# Cuanto fastq sin recortar se permite tener en disco a la vez. Las 417 corridas
# son ~1.1 TB descomprimidas: volcarlas todas antes de recortar no entra en
# ningun disco de este proyecto, y `galga` solo ya son 376 GB. Con el solapado
# activo el pico es ~2x este numero, porque la tanda que se vuelca convive con
# la que se esta recortando.
PRESUPUESTO_GB="${PRESUPUESTO_GB:-40}"
# Volcar la tanda siguiente mientras se recorta la actual. fasterq-dump es de
# disco y cutadapt de CPU, asi que se tapan bastante bien. SOLAPAR=0 lo apaga.
SOLAPAR="${SOLAPAR:-1}"
# Sin -t los temporales de fasterq-dump van al CWD; llegaron a 109 GB.
TMP_FASTERQ="${TMP_FASTERQ:-$PROY_DIR/.tmp}"

die() { echo "ERROR: $*" >&2; exit 1; }

# --- lectura de las specs -----------------------------------------------------

sin_comentarios() {
  grep -vE '^[[:space:]]*#' "$1" | grep -vE '^[[:space:]]*$'
}

# org+bioproject[+run] -> una columna de la tabla.
#
# Se buscaba SOLO por proyecto, con el razonamiento de que "el adaptador es del
# kit, no de la corrida". El kit si es del proyecto; lo que no es cierto es que
# un BioProject use un solo kit. Medido en gadmo PRJNA328800: 6 de 12 corridas
# retuvieron 0.6-2.0% contra el 51% esperado, y las otras 6 dieron 38-58%. Una
# fila por proyecto no puede expresar eso, y el modo de falla es el silencioso
# de siempre — cutadapt con --trimmed-only descarta lo que no matchea y sale
# con 0.
#
# Por eso la tabla lleva columna `run`: una fila con el RUN exacto le gana a la
# fila `-` del proyecto. El fallback a `-` es lo que mantiene las 19 filas
# existentes valiendo sin tocarlas.
columna_de() {
  local col="$1" org="$2" proy="$3" run="${4:-}"
  sin_comentarios "$ADAPTADORES_TSV" | tail -n +2 \
    | awk -F'\t' -v c="$col" -v o="$org" -v p="$proy" -v r="$run" '
        $1 != o || $2 != p { next }
        $3 == r && r != "" { print $c; hecho = 1; exit }  # la de esta corrida
        $3 == "-"          { gen = $c }                   # la del proyecto, de reserva
        # OJO: en awk `exit` TAMBIEN corre el END, asi que sin el flag la fila
        # del proyecto se imprimiria pegada a la de la corrida y el llamador
        # recibiria dos lineas donde espera una.
        END { if (!hecho && gen != "") print gen }'
}

secuencia_de()          { columna_de 5 "$@"; }
familia_de()            { columna_de 4 "$@"; }
retencion_esperada_de() { columna_de 8 "$@"; }

# run -> "<5'> <3'>" de cada corrida del proyecto, de data/extremos.tsv. La
# busqueda (la fila de la corrida le gana a la del proyecto) es la de
# scripts/extremos.py, la misma que usan align.sh y de_quien.py: una sola.
declare -A EXT=()
cargar_extremos() {
  local salida r a b
  salida=$("$EXTREMOS_PY" --tabla "$EXTREMOS_TSV" --manifest "$MANIFEST" requeridos "$1" "$2") \
    || die "no pude leer $EXTREMOS_TSV (ver arriba)"
  EXT=()
  while IFS=$'\t' read -r r a b; do
    [[ -n "$r" ]] && EXT[$r]="$a $b"
  done <<<"$salida"
  return 0
}

corte_de()  { echo "${EXT[$1]:-0 0}"; }
etiqueta()  { if [[ "$1 $2" == "0 0" ]]; then echo "-"; else echo "$1/$2"; fi; }

# Lo recortado con otro corte que el que pide la tabla: run <TAB> tiene <TAB> pide.
desactualizadas_de() {
  "$EXTREMOS_PY" --tabla "$EXTREMOS_TSV" desactualizadas "$1" "$2"
}

# El .sra en el destino. Las tres disposiciones que deja prefetch, y -s para que
# un fichero truncado no cuente como bajado.
ruta_sra() {
  local org="$1" run="$2" c
  for c in "$SRA_DEST/$org/$run.sra" "$SRA_DEST/$run/$run.sra" "$SRA_DEST/$run.sra"; do
    [[ -s "$c" ]] && { echo "$c"; return 0; }
  done
  return 1
}

# --- el manifiesto ------------------------------------------------------------

# Los proyectos (org + rol) que toca el filtro, en el orden del manifiesto.
proyectos() {
  local filtro="${1:-}" org rol=""
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: fetch_runs.sh manifest"
  if [[ "$filtro" == */* ]]; then org="${filtro%%/*}"; rol="${filtro##*/}"; else org="$filtro"; fi
  awk -F'\t' -v o="$org" -v r="$rol" '
    NR>1 && (o=="" || $1==o) && (r=="" || $4==r) {
      k = $1 "\t" $4
      if (!(k in visto)) { visto[k]=1; print k }
    }' "$MANIFEST"
}

# run, bioproject, read_count, base_count de un proyecto.
corridas() {
  awk -F'\t' -v o="$1" -v r="$2" 'NR>1 && $1==o && $4==r {print $2"\t"$3"\t"$6"\t"$7}' "$MANIFEST"
}

# --- el registro de lo recortado ----------------------------------------------
#
# `recortadas.tsv` es acumulativo y lo escribe este script; inputs.json lo pisa
# YASMA en cada llamada (ver cabecera). Se leen los dos y se unen, para que un
# directorio recortado por una version anterior del script —que solo tenia
# inputs.json— siga contando.
#
# El nombre del fichero NO se reconstruye: sale de `trimmed_libraries`, que es el
# registro que YASMA deja de lo que produjo. Adivinarlo fallo de dos formas, las
# dos silenciosas: es <RUN>.t.fq.gz y no <RUN>.tfq.gz ('.t' + un formato que ya
# trae punto), y una libreria PRE-TRIMMED no escribe nada en trim/ — anota la
# ruta del fichero original, o sea que la salida ES la entrada.
LEDGER=recortadas.tsv

recortadas_de() {
  local dir="$1"
  [[ -d "$dir" ]] || return 0
  LEDGER="$LEDGER" python3 - "$dir" <<'PY'
import json, os, pathlib, sys
d = pathlib.Path(sys.argv[1])
rutas = []
led = d / os.environ['LEDGER']
if led.is_file():
    for ln in led.read_text().splitlines()[1:]:
        c = ln.split('\t')
        if len(c) > 2:
            rutas.append(c[2])
f = d / 'inputs.json'
if f.is_file():
    try:
        rutas += json.loads(f.read_text()).get('trimmed_libraries') or []
    except (json.JSONDecodeError, OSError):
        pass
visto = set()
for r in rutas:
    q = pathlib.Path(r)
    if not q.is_absolute():
        q = d / q
    # Solo cuenta si el fichero esta y no esta vacio: una corrida cortada a
    # mitad deja la entrada en el registro y el fichero truncado.
    if q.is_file() and q.stat().st_size > 0:
        run = q.name.split('.')[0]
        if run not in visto:
            visto.add(run)
            print(run)
PY
}

# --- validacion ---------------------------------------------------------------

# Que cada proyecto del manifiesto tenga fila en la tabla, y que su secuencia
# sea usable. Corre ANTES de tocar nada: descubrir a la corrida 300 que falta un
# adaptador es tarde.
validar() {
  local filtro="${1:-}"
  local org rol run proy sec fam _rc _bc
  local -a sin_fila=() no_recortable=() pre_cortada=()
  local vistos=""

  while IFS=$'\t' read -r org rol; do
    cargar_extremos "$org" "$rol"
    # Se recorria deduplicando por proyecto, que era correcto cuando la tabla
    # tenia una sola fila por BioProject. Con filas por corrida eso deja de
    # cubrir: una corrida con su propia fila mala no se miraria nunca. Ahora se
    # miran TODAS las corridas y lo que se deduplica es el mensaje, asi que un
    # proyecto entero sin fila sigue siendo una linea y una excepcion por
    # corrida nombra la corrida.
    local origen etiqueta
    while IFS=$'\t' read -r run proy _rc _bc; do
      sec=$(secuencia_de "$org" "$proy" "$run")
      fam=$(familia_de "$org" "$proy" "$run")
      # La columna `run` de la fila que aplico ES la clave: "-" si vino la del
      # proyecto, el RUN si gano una excepcion.
      origen=$(columna_de 3 "$org" "$proy" "$run")
      if [[ "$origen" == "-" || -z "$origen" ]]; then
        etiqueta="$org $proy"
      else
        etiqueta="$org $proy $run"
      fi
      # Una PRE-TRIMMED no pasa por cutadapt: no hay donde cortarle las puntas.
      if [[ "$sec" == "PRE-TRIMMED" && "$(corte_de "$run")" != "0 0" ]]; then
        pre_cortada+=("$org $proy $run")
      fi
      [[ " $vistos " == *" $etiqueta|"* ]] && continue
      vistos="$vistos $etiqueta|"
      if [[ -z "$sec" ]]; then
        sin_fila+=("$etiqueta")
        continue
      fi
      # Un 5p o un sin_identificar no puede ir como -a de cutadapt: recortariamos
      # con la secuencia equivocada y el resultado no falla ruidosamente.
      case "$fam" in
        5p:*|'??:'*) no_recortable+=("$etiqueta ($fam)") ;;
      esac
    done < <(corridas "$org" "$rol")
  done < <(proyectos "$filtro")

  if [[ ${#sin_fila[@]} -gt 0 ]]; then
    echo "Proyectos sin fila en $(basename "$ADAPTADORES_TSV"):" >&2
    printf '  %s\n' "${sin_fila[@]}" >&2
    echo "  Medilos con: ./scripts/fetch_runs.sh perfil --proyectos" >&2
    echo "  y agregalos a la tabla. No recorto con un adaptador adivinado." >&2
    return 1
  fi
  if [[ ${#no_recortable[@]} -gt 0 ]]; then
    echo "Proyectos cuyo adaptador mayoritario no sirve para recortar:" >&2
    printf '  %s\n' "${no_recortable[@]}" >&2
    echo "  Un 5p es dimero o quimera, no read-through; un sin_identificar no" >&2
    echo "  esta en los 161 adaptadores de YASMA. Los dos necesitan mirarse." >&2
    return 1
  fi
  if [[ ${#pre_cortada[@]} -gt 0 ]]; then
    echo "Corridas PRE-TRIMMED con corte de extremos en $(basename "$EXTREMOS_TSV"):" >&2
    printf '  %s\n' "${pre_cortada[@]}" >&2
    echo "  Una PRE-TRIMMED no pasa por cutadapt, y cortarle las puntas no está" >&2
    echo "  implementado: hay que mirarla antes de inventar cómo." >&2
    return 1
  fi
  return 0
}

# --- modos --------------------------------------------------------------------

cmd_plan() {
  local filtro="${1:-}"
  echo "manifiesto : $MANIFEST"
  echo "adaptadores: $ADAPTADORES_TSV"
  echo "origen     : $SRA_DEST"
  echo "destino    : $PROY_DIR"
  echo "presupuesto: $PRESUPUESTO_GB GB de fastq sin recortar por tanda"
  echo
  validar "$filtro" || return 1

  local org rol run proy rc bc dir sec ya est viejas desact
  local total=0 listo=0 hay=0 sin_sra=0
  printf '%-18s %-14s %-12s %-22s %-5s %8s  %s\n' PROYECTO BIOPROJECT CORRIDA SECUENCIA CORTE FASTQ_GB ESTADO
  while IFS=$'\t' read -r org rol; do
    dir="$PROY_DIR/${org}_${rol}"
    ya=" $(recortadas_de "$dir" | tr '\n' ' ')"
    cargar_extremos "$org" "$rol"
    desact=$(desactualizadas_de "$dir" "$org") || die "no pude leer $EXTREMOS_TSV"
    viejas=" $(cut -f1 <<<"$desact" | tr '\n' ' ')"
    while IFS=$'\t' read -r run proy rc bc; do
      sec=$(secuencia_de "$org" "$proy" "$run")
      if [[ "$ya" == *" $run "* && "$viejas" == *" $run "* ]]; then
        est="recortada con otro corte: se rehace"; hay=$((hay+1))
      elif [[ "$ya" == *" $run "* ]]; then
        est="ya recortada"; listo=$((listo+1))
      elif ruta_sra "$org" "$run" >/dev/null; then
        est="por recortar"; hay=$((hay+1))
      else
        est="FALTA el .sra"; sin_sra=$((sin_sra+1))
      fi
      printf '%-18s %-14s %-12s %-22s %-5s %8.1f  %s\n' \
        "${org}_${rol}" "$proy" "$run" "${sec:0:22}" "$(etiqueta $(corte_de "$run"))" \
        "$(awk -v b="$bc" -v r="$rc" 'BEGIN{printf "%.1f", (b*2+r*35)/1e9}')" "$est"
      total=$((total+1))
    done < <(corridas "$org" "$rol")
  done < <(proyectos "$filtro")
  [[ $total -gt 0 ]] || die "el filtro '${filtro:-(todo)}' no encontró ningún proyecto"
  echo
  echo "$total corridas: $listo ya recortadas, $hay por recortar, $sin_sra sin .sra"
  [[ $sin_sra -eq 0 ]] || echo "Traé los .sra con: ./scripts/drive_pull.sh sra <org> --go" >&2
}

cmd_estado() {
  local filtro="${1:-}"
  local org rol run proy rc bc dir ya n_ok n_no t_ok=0 t_no=0 desact viejas
  printf '%-18s %10s %8s\n' PROYECTO RECORTADAS FALTAN
  while IFS=$'\t' read -r org rol; do
    dir="$PROY_DIR/${org}_${rol}"
    ya=" $(recortadas_de "$dir" | tr '\n' ' ')"
    # Recortada con otro corte de extremos que el de la tabla = falta.
    desact=$(desactualizadas_de "$dir" "$org") || die "no pude leer $EXTREMOS_TSV"
    viejas=" $(cut -f1 <<<"$desact" | tr '\n' ' ')"
    n_ok=0; n_no=0
    while IFS=$'\t' read -r run proy rc bc; do
      if [[ "$ya" == *" $run "* && "$viejas" != *" $run "* ]]; then
        n_ok=$((n_ok+1))
      else
        n_no=$((n_no+1))
      fi
    done < <(corridas "$org" "$rol")
    printf '%-18s %10d %8d\n' "${org}_${rol}" "$n_ok" "$n_no"
    t_ok=$((t_ok+n_ok)); t_no=$((t_no+n_no))
  done < <(proyectos "$filtro")
  echo
  echo "total: $t_ok recortadas, $t_no faltan"
}

# inputs.json para yasma trim. Las librerias tienen que estar ADENTRO del
# output_directory y las rutas se guardan relativas a el: YASMA hace
# relative_to(output_directory) al leer y open() desde el CWD al usar, asi que
# tambien hay que correrlo parado adentro. Ver docs/yasma.md.
escribir_inputs() {
  local dir="$1"; shift
  python3 - "$dir" "$@" <<'PY'
import json, pathlib, sys
dir_proy = pathlib.Path(sys.argv[1])
# los pares vienen como <nombre_fichero> <secuencia> ...
pares = sys.argv[2:]
libs, adapters = [], {}
for i in range(0, len(pares), 2):
    nombre, sec = pares[i], pares[i + 1]
    libs.append(f"untrimmed/{nombre}")
    adapters[nombre] = sec

f = dir_proy / 'inputs.json'
datos = json.loads(f.read_text()) if f.is_file() else {}
datos.update({
    'project_name': dir_proy.name,
    'untrimmed_libraries': libs,
    'adapters': adapters,
    # Explicitos aunque coincidan con los defaults de yasma trim: un default que
    # cambie en una version nueva no avisa, y la ventana 15-50 esta elegida para
    # no truncar tRFs (30-40 nt) ni dejar los piRNAs (24-32) sin margen.
    'min_length': 15,
    'max_length': 50,
})
f.write_text(json.dumps(datos, indent=1) + '\n')
print(f"   inputs.json: {len(libs)} librerías, {len(set(adapters.values()))} adaptador(es)")
PY
}

# Lo que YASMA acaba de producir -> al ledger acumulativo, y el ledger completo
# de vuelta a inputs.json. Los conteos de cutadapt salen del log de ESTA tanda,
# porque trim/log.txt lo trunca la llamada siguiente.
cosechar() {
  local dir="$1" log="$2" lista="$3"
  LEDGER="$LEDGER" python3 - "$dir" "$log" "$lista" <<'PY'
import json, os, pathlib, re, sys, time
d, log = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
# run -> bioproject de ESTA tanda. Una tanda puede cruzar dos BioProjects
# (maggi_primario junta PRJNA154615 y PRJNA232734, con 62% y 87% de retencion
# esperada), asi que tomar el del primero desviaria justo lo que `verificar` mide.
de_proyecto = {}
for ln in pathlib.Path(sys.argv[3]).read_text().splitlines():
    c = ln.split('\t')
    if len(c) >= 2:
        de_proyecto[c[0]] = c[1]

# --- conteos de cutadapt, por libreria ---------------------------------------
# El bloque va precedido por 'trimming: <ruta>'. Los numeros de cutadapt vienen
# con separador de miles, asi que hay que sacarle las comas antes de int().
# Una corrida 4N trae un segundo bloque, 'extremos: <salida> cortados 5/3', con
# el reporte de la segunda pasada: de ahi sale lo que SOBREVIVIO, y lo que
# entro sigue siendo lo de la primera.
stats, cur, fase = {}, None, None
num = lambda s: int(s.replace(',', ''))
for ln in log.read_text(errors='replace').splitlines() if log.is_file() else []:
    m = re.match(r'\s*trimming:\s*(\S.*)$', ln)
    if m:
        cur = pathlib.Path(m.group(1).strip()).name.split('.')[0]
        stats.setdefault(cur, {})
        fase = 'adaptador'
        continue
    m = re.match(r'\s*extremos:\s*(\S+)\s+cortados\s+(\d+)/(\d+)', ln)
    if m:
        cur = pathlib.Path(m.group(1)).name.split('.')[0]
        stats.setdefault(cur, {})['extremos'] = f"{m.group(2)}/{m.group(3)}"
        fase = 'extremos'
        continue
    if cur is None:
        continue
    m = re.match(r'\s*Total reads processed:\s*([\d,]+)', ln)
    if m:
        if fase != 'extremos':
            stats[cur]['in'] = num(m.group(1))
        continue
    m = re.match(r'\s*Reads written \(passing filters\):\s*([\d,]+)', ln)
    if m:
        stats[cur]['out'] = num(m.group(1))

# --- lo que YASMA dice que produjo -------------------------------------------
producidas = json.loads((d / 'inputs.json').read_text()).get('trimmed_libraries') or []

led = d / os.environ['LEDGER']
# `extremos` es lo que se le corto a cada punta ('-' nada). Un ledger de antes
# de esa columna no corto nada, asi que sus filas se completan con '-'.
cab = 'run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc\textremos\n'
previas, orden = {}, []
if led.is_file():
    for ln in led.read_text().splitlines()[1:]:
        c = ln.split('\t')
        if len(c) >= 7:
            c = c[:8] + ['-'] * (8 - len(c[:8]))
            if c[0] not in previas:
                orden.append(c[0])
            previas[c[0]] = c

ahora = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
nuevas = 0
for lib in producidas:
    q = pathlib.Path(lib)
    rel = str(q) if not q.is_absolute() else (
        str(q.relative_to(d)) if q.is_relative_to(d) else str(q))
    run = q.name.split('.')[0]
    s = stats.get(run, {})
    # Una PRE-TRIMMED no pasa por cutadapt: YASMA la anota tal cual y no hay
    # conteos que leer. Retiene el 100% por definicion, y ademas NO pasa por el
    # filtro de longitud: es la unica que entra al alineamiento sin la ventana.
    if 'in' in s and 'out' in s and s['in']:
        pct = '%.1f' % (100.0 * s['out'] / s['in'])
        fila = [run, de_proyecto.get(run, '-'), rel, str(s['in']), str(s['out']), pct, ahora,
                s.get('extremos', '-')]
    else:
        fila = [run, de_proyecto.get(run, '-'), rel, '-', '-', 'PRE-TRIMMED', ahora, '-']
    if run not in previas:
        orden.append(run)
        nuevas += 1
    previas[run] = fila

with open(led, 'w') as f:
    f.write(cab)
    for run in orden:
        f.write('\t'.join(previas[run]) + '\n')

# --- y el ledger completo de vuelta a inputs.json ----------------------------
# yasma trim lo dejo con SOLO lo de esta tanda. Los comandos de aguas abajo leen
# de aca, asi que tiene que estar la lista entera.
datos = json.loads((d / 'inputs.json').read_text())
completas = []
for run in orden:
    rel = previas[run][2]
    q = pathlib.Path(rel)
    if not q.is_absolute():
        q = d / q
    if q.is_file() and q.stat().st_size > 0:
        completas.append(rel)
datos['trimmed_libraries'] = completas
(d / 'inputs.json').write_text(json.dumps(datos, indent=1) + '\n')
print(f"   recortadas: +{nuevas} (acumulado {len(completas)})")
PY
}

# Una corrida de kit 4N: la llamada a cutadapt de yasma trim, partida en dos
# para cortar las puntas (ver la cabecera). Deja la salida donde la dejaria
# yasma (trim/<RUN>.t.fq.gz), los dos reportes en el log de la tanda con el
# formato que lee cosechar(), y la anota en trimmed_libraries.
#
# El .t.fq.gz aparece de una vez, con mv, y recien despues se anota: una sesion
# que muere a la mitad no deja un fichero a medio cortar contando como hecho.
recortar_extremos() {
  local dir="$1" log="$2" fq_nombre="$3" sec="$4" c5="$5" c3="$6"
  local run="${fq_nombre%%.*}"
  local sal="trim/$run.t.fq.gz" tmp="$dir/trim/.$run.extremos.fq.gz"
  local r1="$dir/logs/.$run.pasada1" r2="$dir/logs/.$run.pasada2"
  local -a u5=() u3=()
  if [[ "$c5" -gt 0 ]]; then u5=(-u "$c5"); fi
  if [[ "$c3" -gt 0 ]]; then u3=(-u "-$c3"); fi
  echo "   $run: cortando $c5 nt del 5' y $c3 del 3' (data/extremos.tsv)"
  if ! { cutadapt ${u5[@]+"${u5[@]}"} -a "$sec" -j "$CORES" -O 4 --max-n 0 --trimmed-only \
           "$dir/untrimmed/$fq_nombre" 2> "$r1" </dev/null \
         | cutadapt ${u3[@]+"${u3[@]}"} --minimum-length 15 --maximum-length 50 -j "$CORES" \
           -o "$tmp" - > "$r2" 2>&1; }; then
    tail -20 "$r1" "$r2" >&2 || true
    rm -f "$tmp"
    die "cutadapt falló cortando los extremos de $run (logs: $r1, $r2)"
  fi
  mv -f "$tmp" "$dir/$sal"
  { echo "trimming: untrimmed/$fq_nombre"; cat "$r1"
    echo "extremos: $sal cortados $c5/$c3"; cat "$r2"; } >> "$log"
  rm -f "$r1" "$r2"
  python3 - "$dir" "$sal" <<'PY'
import json, pathlib, sys
f = pathlib.Path(sys.argv[1]) / 'inputs.json'
d = json.loads(f.read_text()) if f.is_file() else {}
libs = d.get('trimmed_libraries') or []
if sys.argv[2] not in libs:
    libs.append(sys.argv[2])
d['trimmed_libraries'] = libs
f.write_text(json.dumps(d, indent=1) + '\n')
PY
}

# trimmed_libraries = [] : lo que cosechar() lee como "producido" tiene que ser
# lo de ESTA tanda, que es lo que yasma trim deja. Cuando la tanda es toda 4N
# yasma no corre, y sin esto cosechar veria la lista completa de las tandas
# anteriores, sin conteos, y las reescribiria como PRE-TRIMMED.
sin_registradas() {
  python3 - "$1" <<'PY'
import json, pathlib, sys
f = pathlib.Path(sys.argv[1]) / 'inputs.json'
d = json.loads(f.read_text()) if f.is_file() else {}
d['trimmed_libraries'] = []
f.write_text(json.dumps(d, indent=1) + '\n')
PY
}

# Vuelca a fastq las corridas de un fichero de tanda. Cada linea:
#   run <TAB> bioproject <TAB> secuencia <TAB> nombre_fichero
volcar() {
  local org="$1" dir="$2" lista="$3"
  local run proy sec fq_nombre src destino
  while IFS=$'\t' read -r run proy sec fq_nombre; do
    destino="$dir/untrimmed/$fq_nombre"
    [[ -s "$destino" ]] && continue
    src=$(ruta_sra "$org" "$run") || { echo "   FALTA el .sra: $run" >&2; continue; }
    echo "   fastq: $run"
    # El fastq va SIN sufijo _1: check_paired_end de YASMA descarta el _2 de
    # cualquier par, y un nombre que termina en _1 lo hace tratar la corrida
    # como paired aunque sea SINGLE.
    # El `</dev/null` no es decorativo: el `while ... done < "$lista"` de abajo
    # le deja como stdin el fichero de tanda, y un comando que leyera de ahi se
    # comeria las corridas siguientes sin que nada lo dijera.
    fasterq-dump --concatenate-reads -e "$CORES" -t "$TMP_FASTERQ" \
      -O "$dir/untrimmed" -o "${fq_nombre%.gz}" "$src" >/dev/null </dev/null
    # Una PRE-TRIMMED se guarda comprimida porque NO se borra: YASMA la anota
    # como su propia salida, asi que ese fastq se queda en disco para siempre.
    # `cloro PRJEB43636` son 53 GB sin comprimir contra ~14 con gzip, y YASMA lee
    # .gz sin problema (get_library_format lo contempla).
    if [[ "$fq_nombre" == *.gz ]]; then
      if command -v pigz >/dev/null; then
        pigz -p "$CORES" -f "$dir/untrimmed/${fq_nombre%.gz}"
      else
        gzip -f "$dir/untrimmed/${fq_nombre%.gz}"
      fi
    fi
  done < "$lista"
}

cmd_correr() {
  local filtro="${1:-}"
  command -v yasma >/dev/null || die "no está yasma en el PATH (ver docs/yasma.md)"
  command -v cutadapt >/dev/null || die "no está cutadapt: yasma trim lo envuelve"
  command -v fasterq-dump >/dev/null || die "no está fasterq-dump"
  validar "$filtro" || return 1

  mkdir -p "$TMP_FASTERQ"
  local org rol dir ya run proy rc bc sec fq_nombre
  local n_proy=0
  while IFS=$'\t' read -r org rol; do
    n_proy=$((n_proy+1))
    dir="$PROY_DIR/${org}_${rol}"
    mkdir -p "$dir/untrimmed" "$dir/trim" "$dir/logs"
    echo "== ${org}_${rol}  ($dir)"

    ya=" $(recortadas_de "$dir" | tr '\n' ' ')"
    cargar_extremos "$org" "$rol"
    # Recortadas con otro corte de extremos que el que pide la tabla hoy: se
    # rehacen. Es el caso normal, no una rareza: la fila de un kit 4N llega
    # despues de recortar, alinear y ver que no alinea.
    local desact viejas r_ a_ q_
    desact=$(desactualizadas_de "$dir" "$org") || die "no pude leer $EXTREMOS_TSV"
    viejas=" $(cut -f1 <<<"$desact" | tr '\n' ' ')"
    if [[ -n "$desact" ]]; then
      echo "   recortadas con otro corte de extremos que el de $(basename "$EXTREMOS_TSV"), se rehacen:"
      while IFS=$'\t' read -r r_ a_ q_; do
        echo "     $r_: tiene $a_, la tabla pide $q_"
      done <<<"$desact"
    fi

    # --- armar las tandas, acotadas por el presupuesto de disco ---------------
    local pend="$dir/logs/.pendientes" nt=0
    : > "$pend"
    while IFS=$'\t' read -r run proy rc bc; do
      [[ "$ya" == *" $run "* && "$viejas" != *" $run "* ]] && continue   # idempotente
      ruta_sra "$org" "$run" >/dev/null || { echo "   FALTA el .sra: $run" >&2; continue; }
      sec=$(secuencia_de "$org" "$proy" "$run")
      # La PRE-TRIMMED es la unica que se guarda comprimida: ver volcar().
      if [[ "$sec" == "PRE-TRIMMED" ]]; then fq_nombre="$run.fastq.gz"; else fq_nombre="$run.fastq"; fi
      printf '%s\t%s\t%s\t%s\t%s\n' "$run" "$proy" "$sec" "$fq_nombre" \
        "$(awk -v b="$bc" -v r="$rc" 'BEGIN{printf "%d", b*2+r*35}')" >> "$pend"
    done < <(corridas "$org" "$rol")

    if [[ ! -s "$pend" ]]; then
      echo "   nada que recortar"
      rm -f "$pend"
      continue
    fi

    rm -f "$dir/logs"/tanda_*.lst
    nt=$(awk -F'\t' -v dir="$dir/logs" -v pres="$PRESUPUESTO_GB" '
      BEGIN { lim = pres * 1e9; n = 1; acum = 0 }
      {
        # Una corrida que sola pasa el presupuesto va igual, sola: no se puede
        # partir un .sra a la mitad.
        if (acum > 0 && acum + $5 > lim) { n++; acum = 0 }
        printf "%s\t%s\t%s\t%s\n", $1, $2, $3, $4 > sprintf("%s/tanda_%02d.lst", dir, n)
        acum += $5
      }
      END { print n }' "$pend")
    rm -f "$pend"
    echo "   $(cat "$dir/logs"/tanda_*.lst | wc -l) corridas por recortar en $nt tanda(s) de <= $PRESUPUESTO_GB GB"

    # --- volcar / recortar, solapando una tanda con la siguiente --------------
    local i lista lista_sig log pid=""
    for ((i = 1; i <= nt; i++)); do
      lista=$(printf '%s/logs/tanda_%02d.lst' "$dir" "$i")
      if [[ -n "$pid" ]]; then wait "$pid" || die "falló el volcado de la tanda $i en ${org}_${rol}"; pid=""; fi
      [[ -s "$lista" ]] || continue
      volcar "$org" "$dir" "$lista"

      if [[ "$SOLAPAR" == "1" && $i -lt $nt ]]; then
        lista_sig=$(printf '%s/logs/tanda_%02d.lst' "$dir" "$((i+1))")
        volcar "$org" "$dir" "$lista_sig" & pid=$!
      fi

      # Las corridas sin corte de extremos van a yasma trim; las 4N, a la misma
      # llamada a cutadapt partida en dos (recortar_extremos).
      local -a pares=() cortes=()
      local j
      while IFS=$'\t' read -r run proy sec fq_nombre; do
        [[ -s "$dir/untrimmed/$fq_nombre" ]] || continue
        if [[ "$(corte_de "$run")" == "0 0" ]]; then
          pares+=("$fq_nombre" "$sec")
        else
          cortes+=("$fq_nombre" "$sec" "$(corte_de "$run")")
        fi
      done < "$lista"
      [[ ${#pares[@]} -gt 0 || ${#cortes[@]} -gt 0 ]] || continue

      echo "   -- tanda $i/$nt: $(( ${#pares[@]} / 2 + ${#cortes[@]} / 3 )) librería(s)"
      log=$(printf '%s/logs/tanda_%02d.log' "$dir" "$i")
      if [[ ${#pares[@]} -gt 0 ]]; then
        escribir_inputs "$dir" "${pares[@]}"
        # Parado adentro del directorio del proyecto, con -o absoluto.
        # Sin --override: `yasma trim` no tiene esa opcion (si la tiene `yasma
        # adapter`, y asumirlo por analogia hacia fallar la corrida entera).
        # Sin --cleanup: iteraria ic.inputs['srrs'], que es None -> TypeError.
        ( cd "$dir" && yasma trim -o "$dir" --cores "$CORES" \
            --min_length 15 --max_length 50 </dev/null ) > "$log" 2>&1 \
          || { tail -20 "$log" >&2; die "yasma trim falló en ${org}_${rol}, tanda $i (log: $log)"; }
      else
        : > "$log"
        sin_registradas "$dir"
      fi
      for ((j = 0; j < ${#cortes[@]}; j += 3)); do
        # shellcheck disable=SC2086  # "5 3" -> dos argumentos, a proposito
        recortar_extremos "$dir" "$log" "${cortes[j]}" "${cortes[j+1]}" ${cortes[j+2]}
      done
      cosechar "$dir" "$log" "$lista"

      # El fastq sin recortar ya no hace falta: se re-genera del .sra. La
      # PRE-TRIMMED no se toca, porque ES la salida que YASMA anoto.
      while IFS=$'\t' read -r run proy sec fq_nombre; do
        [[ "$sec" == "PRE-TRIMMED" ]] && continue
        rm -f "$dir/untrimmed/$fq_nombre"
      done < "$lista"
    done
    if [[ -n "$pid" ]]; then wait "$pid" || true; fi
  done < <(proyectos "$filtro")
  [[ $n_proy -gt 0 ]] || die "el filtro '${filtro:-(todo)}' no encontró ningún proyecto"
  echo
  cmd_estado "$filtro"
}

# La retencion medida contra la que `perfil` predijo. Es el chequeo que atrapa el
# unico modo de fallo que el recorte tiene y que NO hace ruido: recortar con la
# secuencia equivocada. cutadapt corre con --trimmed-only, asi que un adaptador
# que no corresponde no da error — deja un .t.fq.gz casi vacio y el pipeline
# sigue. `maggi` primario es el caso de manual: dos BioProjects de 2011 y 2014
# con kits distintos, y una sola fila de adaptador para los dos habria vaciado el
# que no corresponde.
cmd_verificar() {
  local filtro="${1:-}" org rol dir esp fallas=0 filas=0 n_desact=0
  printf '%-18s %-14s %-12s %10s %10s %11s %7s  %s\n' \
    PROYECTO BIOPROJECT CORRIDA READS_IN READS_OUT MEDIDA ESPERA VEREDICTO
  while IFS=$'\t' read -r org rol; do
    dir="$PROY_DIR/${org}_${rol}"
    [[ -f "$dir/$LEDGER" ]] || { echo "   sin recortar: ${org}_${rol}" >&2; continue; }
    local run proy fich rin rout pct fecha ext desact r_ a_ q_
    local -A pide=()
    desact=$(desactualizadas_de "$dir" "$org") || die "no pude leer $EXTREMOS_TSV"
    while IFS=$'\t' read -r r_ a_ q_; do [[ -n "$r_" ]] && pide[$r_]="$q_"; done <<<"$desact"
    while IFS=$'\t' read -r run proy fich rin rout pct fecha ext; do
      [[ "$run" == "run" ]] && continue
      filas=$((filas+1))
      ext="${ext:--}"
      esp=$(retencion_esperada_de "$org" "$proy" "$run")
      local ver
      if [[ -n "${pide[$run]:-}" ]]; then
        # Primero: un recorte que no corta lo que pide la tabla alinea como
        # rhirr_duplicado (0.17%), y eso no lo dice la retencion.
        ver="EXTREMOS DESACTUALIZADOS: cortado ${ext}, la tabla pide ${pide[$run]}"
        n_desact=$((n_desact+1))
      elif [[ "$pct" == "PRE-TRIMMED" ]]; then
        ver="pre-trimmed (sin recorte ni filtro de largo)"
      elif [[ "$ext" != "-" ]]; then
        # `retencion_est` la midio perfil sin cortar las puntas: su inserto
        # incluye las 8 bases y la ventana 15-50 se le aplicaba corrida. Contra
        # eso un DESVIADA no significa nada; VACIA si.
        ver=$(awk -v m="$pct" -v x="$ext" 'BEGIN{
          if (m < 5) print "VACIA — casi seguro el adaptador equivocado"
          else       print "ok, extremos " x " cortados (retencion_est es sin cortar: solo se mira VACIA)" }')
      elif [[ -z "$esp" || "$esp" == "-" ]]; then
        ver="sin retencion_est en la tabla"
      else
        # 15 puntos de holgura: `retencion_est` se midio sobre una corrida por
        # proyecto, y dentro de un mismo proyecto las corridas varian.
        ver=$(awk -v m="$pct" -v e="$esp" 'BEGIN{
          d = m - e; if (d < 0) d = -d
          if (m < 5)       print "VACIA — casi seguro el adaptador equivocado"
          else if (d > 15) print "DESVIADA " (m>e?"+":"-") int(d) " pts"
          else             print "ok" }')
      fi
      [[ "$ver" == ok* || "$ver" == pre-trimmed* ]] || fallas=$((fallas+1))
      # El % se pega solo a lo que es un numero: una PRE-TRIMMED lleva la
      # palabra en la columna, y "PRE-TRIMMED%" se lee como un dato roto.
      local m e
      [[ "$pct" =~ ^[0-9.]+$ ]] && m="$pct%" || m="$pct"
      [[ "$esp" =~ ^[0-9.]+$ ]] && e="$esp%" || e="${esp:--}"
      printf '%-18s %-14s %-12s %10s %10s %11s %7s  %s\n' \
        "${org}_${rol}" "$proy" "$run" "$rin" "$rout" "$m" "$e" "$ver"
    done < "$dir/$LEDGER"
  done < <(proyectos "$filtro")
  echo
  if [[ $fallas -eq 0 ]]; then
    echo "$filas corridas verificadas, ninguna desviada"
    return 0
  fi
  echo "$filas corridas verificadas, $fallas fuera de lo esperado" >&2
  if [[ $n_desact -gt 0 ]]; then
    echo "EXTREMOS DESACTUALIZADOS: recortadas antes de su fila en data/extremos.tsv." >&2
    echo "  ./scripts/trim.sh correr <org>/<rol> las rehace con el corte." >&2
  fi
  echo "Una retención muy por debajo de retencion_est es, casi siempre, la" >&2
  echo "secuencia equivocada en data/adaptadores.tsv: --trimmed-only descarta" >&2
  echo "todo lo que no matchea y no da error. Re-medí con:" >&2
  echo "  ./scripts/fetch_runs.sh perfil --proyectos" >&2
  return 1
}

# Saca corridas del registro para que `correr` las vuelva a recortar.
#
# Hace falta porque el recorte es idempotente por diseno —lo que esta en
# recortadas.tsv no se vuelve a tocar— y eso es exactamente lo que estorba
# cuando lo que hay que rehacer es un recorte MALO. El caso: gadmo/duplicado
# recorto 12 corridas, 6 con la secuencia equivocada; esas 6 figuran como
# hechas, con su .t.fq.gz casi vacio en disco, y `correr` las saltea para
# siempre.
#
# No toca los .sra ni el untrimmed/: solo el registro y las salidas que se van
# a rehacer. Y pide las corridas explicitas — rehacer un proyecto entero por
# error son horas.
cmd_rehacer() {
  local filtro="${1:-}"; shift || true
  [[ "$filtro" == */* ]] || die "uso: $0 rehacer <org>/<rol> RUN [RUN...]"
  [[ $# -ge 1 ]] || die "decime QUE corridas rehacer: $0 rehacer $filtro RUN [RUN...]"
  local org="${filtro%%/*}" rol="${filtro##*/}"
  local dir="$PROY_DIR/${org}_${rol}"
  [[ -d "$dir" ]] || die "no existe $dir"

  LEDGER="$LEDGER" python3 - "$dir" "$@" <<'PY'
import json, os, pathlib, sys
d = pathlib.Path(sys.argv[1])
runs = set(sys.argv[2:])
led = d / os.environ['LEDGER']

fuera, quedan = [], []
if led.is_file():
    lineas = led.read_text().splitlines()
    cab, filas = lineas[0], lineas[1:]
    for ln in filas:
        c = ln.split('\t')
        (fuera if c and c[0] in runs else quedan).append(ln)
    led.write_text('\n'.join([cab] + quedan) + '\n')

# El .t.fq.gz se borra para que no quede un fichero casi vacio pareciendo
# salida buena. La ruta sale del ledger, no se reconstruye: los nombres de
# yasma trim no se pueden adivinar, y una PRE-TRIMMED apunta a su ENTRADA
# —borrarla seria tirar el fastq original.
borrados = 0
for ln in fuera:
    c = ln.split('\t')
    if len(c) < 3 or c[5:6] == ['PRE-TRIMMED']:
        continue
    q = pathlib.Path(c[2])
    if not q.is_absolute():
        q = d / q
    # Solo lo que esta en trim/, que es lo que YASMA produjo. Una PRE-TRIMMED
    # apunta a untrimmed/, o sea a su propia ENTRADA: borrarla seria tirar el
    # fastq original y dejar la corrida sin forma de rehacerse.
    if q.is_file() and q.parent.name == 'trim':
        q.unlink(); borrados += 1

# Y inputs.json, que es de donde leen los comandos de aguas abajo.
f = d / 'inputs.json'
if f.is_file():
    datos = json.loads(f.read_text())
    prev = datos.get('trimmed_libraries') or []
    datos['trimmed_libraries'] = [r for r in prev
                                  if not any(run in r for run in runs)]
    f.write_text(json.dumps(datos, indent=1) + '\n')

faltaron = runs - {ln.split('\t')[0] for ln in fuera}
if faltaron:
    print(f"   no estaban en el registro: {' '.join(sorted(faltaron))}")
print(f"   sacadas del registro: {len(fuera)}, .t.fq.gz borrados: {borrados}")
print(f"   ahora `correr` las vuelve a recortar con lo que diga la tabla")
PY
}

# --- una muestra, para diagnosticar sin recortar el proyecto entero ----------
#
# Existe para la prueba 4N de rhirr_duplicado: hasta ahora solo se podia correr
# despues de recortar y alinear el proyecto entero (~2.5 h) para ver que no
# alineaba. Con N lecturas de cada corrida la prueba son minutos.
#
# El recorte es la MISMA llamada a cutadapt que hace yasma trim v1.1.1
# (trim.py, linea 153), con la secuencia que sale de data/adaptadores.tsv por
# corrida —la misma busqueda que `correr`, incluida la fila por corrida que le
# gana a la del proyecto—. Una PRE-TRIMMED no se recorta, igual que en YASMA.
# Deja un recortadas.tsv con las columnas del de verdad, asi que de_quien.py la
# lee sin saber que es una muestra.
#
# Todas las corridas y no una: un BioProject puede mezclar kits, y ya mezclo.
#
# NO corta extremos aunque data/extremos.tsv lo pida, y es a proposito: esta
# muestra es el instrumento que MIDE ese corte (align.sh prueba -> de_quien.py
# --registrar). Cortada, la medicion saldria relativa a lo ya cortado, y lo que
# se registrara seria la diferencia y no el corte.
cmd_muestra() {
  local filtro="${1:-}" n="${2:-50000}"
  [[ "$filtro" == */* ]] || die "muestra pide <org>/<rol>"
  [[ "$n" =~ ^[0-9]+$ && "$n" -gt 0 ]] || die "N tiene que ser un entero positivo (dio: $n)"
  shift 2 2>/dev/null || shift $#
  local org="${filtro%%/*}" rol="${filtro##*/}"
  command -v fastq-dump >/dev/null || die "no está fastq-dump (sra-tools)"
  command -v cutadapt >/dev/null || die "no está cutadapt"
  local dir="$PROY_DIR/_muestras/${org}_${rol}"
  rm -rf "$dir"; mkdir -p "$dir/trim" "$dir/untrimmed" "$dir/logs"
  local led="$dir/$LEDGER" run proy _rc _bc sec src raw fich rin rout hechas=0
  local -a pedidas=("$@")
  printf 'run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc\n' > "$led"
  while IFS=$'\t' read -r run proy _rc _bc; do
    if [[ ${#pedidas[@]} -gt 0 && " ${pedidas[*]} " != *" $run "* ]]; then continue; fi
    sec=$(secuencia_de "$org" "$proy" "$run")
    if [[ -z "$sec" || "$sec" == "-" ]]; then
      echo "   $run: sin adaptador recortable en $ADAPTADORES_TSV, la salteo" >&2; continue
    fi
    src=$(ruta_sra "$org" "$run") || { echo "   FALTA el .sra: $run" >&2; continue; }
    raw="$dir/untrimmed/$run.fq"
    # -X N: solo las primeras N; -Z a stdout. fasterq-dump no tiene como cortar.
    fastq-dump -X "$n" -Z "$src" > "$raw" 2> "$dir/logs/$run.fastq-dump.log" </dev/null \
      || { echo "   $run: fastq-dump falló (log: $dir/logs/$run.fastq-dump.log)" >&2; continue; }
    rin=$(( $(wc -l < "$raw") / 4 ))
    if [[ "$sec" == "PRE-TRIMMED" ]]; then
      gzip -f "$raw"; fich="untrimmed/$run.fq.gz"; rout=$rin
    else
      fich="trim/$run.t.fq.gz"
      cutadapt -a "$sec" --minimum-length 15 --maximum-length 50 -j "$CORES" -O 4 \
        --max-n 0 --trimmed-only -o "$dir/$fich" "$raw" > "$dir/logs/$run.cutadapt.log" 2>&1 </dev/null \
        || { echo "   $run: cutadapt falló (log: $dir/logs/$run.cutadapt.log)" >&2; continue; }
      rm -f "$raw"
      rout=$(( $(gzip -dc "$dir/$fich" | wc -l) / 4 ))
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$run" "$proy" "$fich" "$rin" "$rout" \
      "$(awk -v a="$rout" -v b="$rin" 'BEGIN { printf "%.1f", (b ? 100*a/b : 0) }')" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$led"
    printf '   %-14s %8d -> %8d   %s\n' "$run" "$rin" "$rout" "${sec:0:21}"
    hechas=$((hechas + 1))
  done < <(corridas "$org" "$rol")
  [[ $hechas -gt 0 ]] || die "no salió ninguna corrida de $filtro (ver arriba)"
  echo "muestra de $hechas corrida(s) en $dir"
}

[[ $# -ge 1 ]] || { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
case "$1" in
  plan)      shift; cmd_plan      "${1:-}" ;;
  estado)    shift; cmd_estado    "${1:-}" ;;
  correr)    shift; cmd_correr    "${1:-}" ;;
  verificar) shift; cmd_verificar "${1:-}" ;;
  rehacer)   shift; cmd_rehacer   "$@" ;;
  muestra)   shift; cmd_muestra   "$@" ;;
  -h|--help) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//' ;;
  *) die "modo desconocido: $1 (plan|correr|estado|verificar|rehacer|muestra)" ;;
esac
