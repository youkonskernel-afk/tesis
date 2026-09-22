#!/usr/bin/env bash
# Recorte de adaptador con YASMA (que envuelve a cutadapt).
#
#   ./scripts/trim.sh plan [<org>]      que se haria, sin hacerlo
#   ./scripts/trim.sh correr [<org>]    recorta
#   ./scripts/trim.sh estado [<org>]    que esta recortado y que falta
#
# Variables: SRA_DEST (donde estan los .sra), TRIM_DIR (donde deja la salida),
# ADAPTADORES_TSV, MANIFEST, CORES.
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
#   cutadapt -a <sec> --minimum-length 15 --maximum-length 50 -O 4 --max-n 0
#            --trimmed-only
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
CORES="${CORES:-1}"
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

# El .sra en el destino. Las tres disposiciones que deja prefetch, y -s para que
# un fichero truncado no cuente como bajado.
ruta_sra() {
  local org="$1" run="$2" c
  for c in "$SRA_DEST/$org/$run.sra" "$SRA_DEST/$run/$run.sra" "$SRA_DEST/$run.sra"; do
    [[ -s "$c" ]] && { echo "$c"; return 0; }
  done
  return 1
}

# Que corridas de un organismo estan YA recortadas, una por linea.
#
# No se reconstruye el nombre del fichero: se lee `trimmed_libraries` de
# inputs.json, que es el registro que YASMA deja de lo que produjo. Adivinarlo
# fallo de dos formas distintas, las dos silenciosas:
#   - el nombre es <RUN>.t.fq.gz, no <RUN>.tfq.gz: yasma hace '.t' + el formato
#     detectado, y ese formato ya viene con punto ('.fq').
#   - para una libreria PRE-TRIMMED yasma NO escribe nada en trim/: anota la
#     ruta del fichero original. La salida ES la entrada.
# Con el nombre mal, nada contaba como recortado: ni idempotencia ni estado.
recortadas_de() {
  local dir="$TRIM_DIR/$1"
  [[ -f "$dir/inputs.json" ]] || return 0
  python3 - "$dir" <<'PY'
import json, pathlib, sys
d = pathlib.Path(sys.argv[1])
try:
    libs = json.loads((d / 'inputs.json').read_text()).get('trimmed_libraries') or []
except (json.JSONDecodeError, OSError):
    sys.exit(0)
for lib in libs:
    q = pathlib.Path(lib)
    if not q.is_absolute():
        q = d / q
    # Solo cuenta si el fichero esta y no esta vacio: una corrida cortada a
    # mitad deja la entrada en el json y el fichero truncado.
    if q.is_file() and q.stat().st_size > 0:
        # <RUN>.t.fq.gz -> RUN ; <RUN>.fastq -> RUN
        print(q.name.split('.')[0])
PY
}

# Las corridas de un organismo (o todas), con su proyecto.
corridas() {
  local filtro="${1:-}"
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: fetch_runs.sh manifest"
  awk -F'\t' -v o="$filtro" 'NR>1 && (o=="" || $1==o) {print $1"\t"$2"\t"$3}' "$MANIFEST"
}

# --- validacion ---------------------------------------------------------------

# Que cada proyecto del manifiesto tenga fila en la tabla, y que su secuencia
# sea usable. Corre ANTES de tocar nada: descubrir a la corrida 300 que falta un
# adaptador es tarde.
validar() {
  local filtro="${1:-}"
  local org run proy sec fam
  local -a sin_fila=() no_recortable=()
  local vistos=""

  while IFS=$'\t' read -r org run proy; do
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
  done < <(corridas "$filtro")

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
  echo
  validar "$filtro" || return 1

  local org run proy sec hay=0 falta=0 sin_sra=0 listo=0 ult_org="" ya=""
  printf '%-7s %-14s %-12s %-22s %s\n' ORG PROYECTO CORRIDA SECUENCIA ESTADO
  while IFS=$'\t' read -r org run proy; do
    # inputs.json se lee una vez por organismo, no una por corrida.
    [[ "$org" != "$ult_org" ]] && { ya=" $(recortadas_de "$org" | tr '\n' ' ')"; ult_org="$org"; }
    sec=$(secuencia_de "$org" "$proy")
    local est
    if [[ "$ya" == *" $run "* ]]; then
      est="ya recortada"; listo=$((listo+1))
    elif ruta_sra "$org" "$run" >/dev/null; then
      est="por recortar"; hay=$((hay+1))
    else
      est="FALTA el .sra"; sin_sra=$((sin_sra+1))
    fi
    printf '%-7s %-14s %-12s %-22s %s\n' "$org" "$proy" "$run" "${sec:0:22}" "$est"
    falta=$((falta+1))
  done < <(corridas "$filtro")
  echo
  echo "$falta corridas: $listo ya recortadas, $hay por recortar, $sin_sra sin .sra"
  [[ $sin_sra -eq 0 ]] || echo "Traé los .sra con: ./scripts/drive_pull.sh sra <org> --go" >&2
}

cmd_estado() {
  local filtro="${1:-}"
  local org run proy ult_org="" ya=""
  declare -A ok no
  while IFS=$'\t' read -r org run proy; do
    [[ "$org" != "$ult_org" ]] && { ya=" $(recortadas_de "$org" | tr '\n' ' ')"; ult_org="$org"; }
    if [[ "$ya" == *" $run "* ]]; then
      ok[$org]=$(( ${ok[$org]:-0} + 1 ))
    else
      no[$org]=$(( ${no[$org]:-0} + 1 ))
    fi
  done < <(corridas "$filtro")
  printf '%-8s %8s %8s\n' ORG RECORTADAS FALTAN
  local o t_ok=0 t_no=0
  for o in $(corridas "$filtro" | cut -f1 | sort -u); do
    printf '%-8s %8d %8d\n' "$o" "${ok[$o]:-0}" "${no[$o]:-0}"
    t_ok=$(( t_ok + ${ok[$o]:-0} )); t_no=$(( t_no + ${no[$o]:-0} ))
  done
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
print(f"  inputs.json: {len(libs)} librerias, {len(set(adapters.values()))} adaptador(es)")
PY
}

cmd_correr() {
  local filtro="${1:-}"
  command -v yasma >/dev/null || die "no está yasma en el PATH (ver docs/yasma.md)"
  command -v cutadapt >/dev/null || die "no está cutadapt: yasma trim lo envuelve"
  command -v fasterq-dump >/dev/null || die "no está fasterq-dump"
  validar "$filtro" || return 1

  mkdir -p "$TMP_FASTERQ"
  local org
  for org in $(corridas "$filtro" | cut -f1 | sort -u); do
    local dir="$TRIM_DIR/$org"
    mkdir -p "$dir/untrimmed" "$dir/trim"
    echo "== $org  ($dir)"

    local run proy sec src
    local -a pares=()
    local ya=" $(recortadas_de "$org" | tr '\n' ' ')"
    while IFS=$'\t' read -r _o run proy; do
      [[ "$ya" == *" $run "* ]] && continue   # idempotente
      src=$(ruta_sra "$org" "$run") || { echo "   FALTA el .sra: $run" >&2; continue; }

      # El fastq va SIN sufijo _1: check_paired_end de YASMA descarta el _2 de
      # cualquier par, y un nombre que termina en _1 lo hace tratar la corrida
      # como paired aunque sea SINGLE.
      local fq="$dir/untrimmed/$run.fastq"
      if [[ ! -s "$fq" ]]; then
        echo "   fastq: $run"
        fasterq-dump --concatenate-reads -t "$TMP_FASTERQ" -O "$dir/untrimmed" \
          -o "$run.fastq" "$src" >/dev/null
      fi
      sec=$(secuencia_de "$org" "$proy")
      pares+=("$run.fastq" "$sec")
    done < <(corridas "$org")

    if [[ ${#pares[@]} -eq 0 ]]; then
      echo "   nada que recortar"
      continue
    fi
    escribir_inputs "$dir" "${pares[@]}"

    # Parado adentro del directorio del proyecto, con -o absoluto.
    # Sin --override: `yasma trim` no tiene esa opcion (si la tiene `yasma
    # adapter`, y asumirlo por analogia hacia fallar la corrida entera).
    ( cd "$dir" && yasma trim -o "$dir" --cores "$CORES" \
        --min_length 15 --max_length 50 </dev/null ) \
      || die "yasma trim falló en $org"
  done
  echo
  cmd_estado "$filtro"
}

[[ $# -ge 1 ]] || { sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
case "$1" in
  plan)   shift; cmd_plan   "${1:-}" ;;
  estado) shift; cmd_estado "${1:-}" ;;
  correr) shift; cmd_correr "${1:-}" ;;
  -h|--help) sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//' ;;
  *) die "modo desconocido: $1 (plan|correr|estado)" ;;
esac
