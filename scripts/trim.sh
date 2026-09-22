#!/usr/bin/env bash
# Recorte de adaptador con YASMA (que envuelve a cutadapt).
#
#   ./scripts/trim.sh plan      [<org>[/<rol>]]   que se haria, sin hacerlo
#   ./scripts/trim.sh correr    [<org>[/<rol>]]   recorta, en tandas
#   ./scripts/trim.sh estado    [<org>[/<rol>]]   que esta recortado y que falta
#   ./scripts/trim.sh verificar [<org>[/<rol>]]   la retencion medida vs la esperada
#
# El filtro es un organismo (`galga`) o un organismo y un rol (`galga/primario`).
#
# Variables: SRA_DEST, TRIM_DIR, ADAPTADORES_TSV, MANIFEST, CORES,
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
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${MANIFEST:-$ROOT/data/srr_manifest.tsv}"
ADAPTADORES_TSV="${ADAPTADORES_TSV:-$ROOT/data/adaptadores.tsv}"
# shellcheck source=_drive_lib.sh
. "$ROOT/scripts/_drive_lib.sh"
# El mismo directorio donde drive_pull.sh deja los .sra, no uno propio.
SRA_DEST="${SRA_DEST:-$(ruta_local sra)}"
# Las lecturas recortadas NO se respaldan: se re-generan de forma determinista
# desde los .sra y data/adaptadores.tsv, que si estan respaldados. Por eso no
# hay fase 'trim' en el mapa de Drive.
TRIM_DIR="${TRIM_DIR:-$ROOT/trim}"
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
TMP_FASTERQ="${TMP_FASTERQ:-$TRIM_DIR/.tmp}"

die() { echo "ERROR: $*" >&2; exit 1; }

# --- lectura de las specs -----------------------------------------------------

sin_comentarios() {
  grep -vE '^[[:space:]]*#' "$1" | grep -vE '^[[:space:]]*$'
}

# org+bioproject -> secuencia. Se busca por proyecto porque el adaptador es del
# kit, no de la corrida.
secuencia_de() {
  local org="$1" proy="$2"
  sin_comentarios "$ADAPTADORES_TSV" | tail -n +2 \
    | awk -F'\t' -v o="$org" -v p="$proy" '$1==o && $2==p {print $4; exit}'
}

familia_de() {
  local org="$1" proy="$2"
  sin_comentarios "$ADAPTADORES_TSV" | tail -n +2 \
    | awk -F'\t' -v o="$org" -v p="$proy" '$1==o && $2==p {print $3; exit}'
}

retencion_esperada_de() {
  local org="$1" proy="$2"
  sin_comentarios "$ADAPTADORES_TSV" | tail -n +2 \
    | awk -F'\t' -v o="$org" -v p="$proy" '$1==o && $2==p {print $7; exit}'
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
  local -a sin_fila=() no_recortable=()
  local vistos=""

  while IFS=$'\t' read -r org rol; do
    while IFS=$'\t' read -r run proy _rc _bc; do
      [[ " $vistos " == *" $org/$proy "* ]] && continue
      vistos="$vistos $org/$proy"
      sec=$(secuencia_de "$org" "$proy")
      fam=$(familia_de "$org" "$proy")
      if [[ -z "$sec" ]]; then
        sin_fila+=("$org $proy")
        continue
      fi
      # Un 5p o un sin_identificar no puede ir como -a de cutadapt: recortariamos
      # con la secuencia equivocada y el resultado no falla ruidosamente.
      case "$fam" in
        5p:*|'??:'*) no_recortable+=("$org $proy ($fam)") ;;
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
  return 0
}

# --- modos --------------------------------------------------------------------

cmd_plan() {
  local filtro="${1:-}"
  echo "manifiesto : $MANIFEST"
  echo "adaptadores: $ADAPTADORES_TSV"
  echo "origen     : $SRA_DEST"
  echo "destino    : $TRIM_DIR"
  echo "presupuesto: $PRESUPUESTO_GB GB de fastq sin recortar por tanda"
  echo
  validar "$filtro" || return 1

  local org rol run proy rc bc dir sec ya est
  local total=0 listo=0 hay=0 sin_sra=0
  printf '%-18s %-14s %-12s %-22s %8s  %s\n' PROYECTO BIOPROJECT CORRIDA SECUENCIA FASTQ_GB ESTADO
  while IFS=$'\t' read -r org rol; do
    dir="$TRIM_DIR/${org}_${rol}"
    ya=" $(recortadas_de "$dir" | tr '\n' ' ')"
    while IFS=$'\t' read -r run proy rc bc; do
      sec=$(secuencia_de "$org" "$proy")
      if [[ "$ya" == *" $run "* ]]; then
        est="ya recortada"; listo=$((listo+1))
      elif ruta_sra "$org" "$run" >/dev/null; then
        est="por recortar"; hay=$((hay+1))
      else
        est="FALTA el .sra"; sin_sra=$((sin_sra+1))
      fi
      printf '%-18s %-14s %-12s %-22s %8.1f  %s\n' \
        "${org}_${rol}" "$proy" "$run" "${sec:0:22}" \
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
  local org rol run proy rc bc dir ya n_ok n_no t_ok=0 t_no=0
  printf '%-18s %10s %8s\n' PROYECTO RECORTADAS FALTAN
  while IFS=$'\t' read -r org rol; do
    dir="$TRIM_DIR/${org}_${rol}"
    ya=" $(recortadas_de "$dir" | tr '\n' ' ')"
    n_ok=0; n_no=0
    while IFS=$'\t' read -r run proy rc bc; do
      if [[ "$ya" == *" $run "* ]]; then n_ok=$((n_ok+1)); else n_no=$((n_no+1)); fi
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
stats, cur = {}, None
num = lambda s: int(s.replace(',', ''))
for ln in log.read_text(errors='replace').splitlines() if log.is_file() else []:
    m = re.match(r'\s*trimming:\s*(\S.*)$', ln)
    if m:
        cur = pathlib.Path(m.group(1).strip()).name.split('.')[0]
        stats.setdefault(cur, {})
        continue
    if cur is None:
        continue
    m = re.match(r'\s*Total reads processed:\s*([\d,]+)', ln)
    if m:
        stats[cur]['in'] = num(m.group(1)); continue
    m = re.match(r'\s*Reads written \(passing filters\):\s*([\d,]+)', ln)
    if m:
        stats[cur]['out'] = num(m.group(1))

# --- lo que YASMA dice que produjo -------------------------------------------
producidas = json.loads((d / 'inputs.json').read_text()).get('trimmed_libraries') or []

led = d / os.environ['LEDGER']
cab = 'run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc\n'
previas, orden = {}, []
if led.is_file():
    for ln in led.read_text().splitlines()[1:]:
        c = ln.split('\t')
        if len(c) >= 7:
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
        fila = [run, de_proyecto.get(run, '-'), rel, str(s['in']), str(s['out']), pct, ahora]
    else:
        fila = [run, de_proyecto.get(run, '-'), rel, '-', '-', 'PRE-TRIMMED', ahora]
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
    dir="$TRIM_DIR/${org}_${rol}"
    mkdir -p "$dir/untrimmed" "$dir/trim" "$dir/logs"
    echo "== ${org}_${rol}  ($dir)"

    ya=" $(recortadas_de "$dir" | tr '\n' ' ')"

    # --- armar las tandas, acotadas por el presupuesto de disco ---------------
    local pend="$dir/logs/.pendientes" nt=0
    : > "$pend"
    while IFS=$'\t' read -r run proy rc bc; do
      [[ "$ya" == *" $run "* ]] && continue          # idempotente
      ruta_sra "$org" "$run" >/dev/null || { echo "   FALTA el .sra: $run" >&2; continue; }
      sec=$(secuencia_de "$org" "$proy")
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

      local -a pares=()
      while IFS=$'\t' read -r run proy sec fq_nombre; do
        [[ -s "$dir/untrimmed/$fq_nombre" ]] || continue
        pares+=("$fq_nombre" "$sec")
      done < "$lista"
      [[ ${#pares[@]} -gt 0 ]] || continue

      echo "   -- tanda $i/$nt: $(( ${#pares[@]} / 2 )) librería(s)"
      escribir_inputs "$dir" "${pares[@]}"
      log=$(printf '%s/logs/tanda_%02d.log' "$dir" "$i")
      # Parado adentro del directorio del proyecto, con -o absoluto.
      # Sin --override: `yasma trim` no tiene esa opcion (si la tiene `yasma
      # adapter`, y asumirlo por analogia hacia fallar la corrida entera).
      # Sin --cleanup: iteraria ic.inputs['srrs'], que es None -> TypeError.
      ( cd "$dir" && yasma trim -o "$dir" --cores "$CORES" \
          --min_length 15 --max_length 50 </dev/null ) > "$log" 2>&1 \
        || { tail -20 "$log" >&2; die "yasma trim falló en ${org}_${rol}, tanda $i (log: $log)"; }
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
  local filtro="${1:-}" org rol dir esp fallas=0 filas=0
  printf '%-18s %-14s %-12s %10s %10s %11s %7s  %s\n' \
    PROYECTO BIOPROJECT CORRIDA READS_IN READS_OUT MEDIDA ESPERA VEREDICTO
  while IFS=$'\t' read -r org rol; do
    dir="$TRIM_DIR/${org}_${rol}"
    [[ -f "$dir/$LEDGER" ]] || { echo "   sin recortar: ${org}_${rol}" >&2; continue; }
    local run proy fich rin rout pct fecha
    while IFS=$'\t' read -r run proy fich rin rout pct fecha; do
      [[ "$run" == "run" ]] && continue
      filas=$((filas+1))
      esp=$(retencion_esperada_de "$org" "$proy")
      local ver
      if [[ "$pct" == "PRE-TRIMMED" ]]; then
        ver="pre-trimmed (sin recorte ni filtro de largo)"
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
  echo "Una retención muy por debajo de retencion_est es, casi siempre, la" >&2
  echo "secuencia equivocada en data/adaptadores.tsv: --trimmed-only descarta" >&2
  echo "todo lo que no matchea y no da error. Re-medí con:" >&2
  echo "  ./scripts/fetch_runs.sh perfil --proyectos" >&2
  return 1
}

[[ $# -ge 1 ]] || { sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
case "$1" in
  plan)      shift; cmd_plan      "${1:-}" ;;
  estado)    shift; cmd_estado    "${1:-}" ;;
  correr)    shift; cmd_correr    "${1:-}" ;;
  verificar) shift; cmd_verificar "${1:-}" ;;
  -h|--help) sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//' ;;
  *) die "modo desconocido: $1 (plan|correr|estado|verificar)" ;;
esac
