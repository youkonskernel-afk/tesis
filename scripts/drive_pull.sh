#!/usr/bin/env bash
# Trae data de Mi unidad/tesis/ a la máquina local. Espejo de drive_push.sh.
#
#   ./scripts/drive_pull.sh <fase> [ORG]        # dry-run, no baja nada
#   ./scripts/drive_pull.sh <fase> [ORG] --go   # baja de verdad
#   ./scripts/drive_pull.sh purge <fase> ORG --go   # borra la copia LOCAL
#
#   fases: bam | yasma | qc | features | modelos | figuras | genomas | sra
#   ORG:   rhirr | sclsc | cloro | phypa | prupe | maldo | gadmo | galga | maggi
#
# Para qué: los .sra viven en Drive (80_sra/, ~190 GB). Bajarlos todos de una
# no entra cómodo en disco, así que se trae un organismo, se alinea, y se
# purga antes de pasar al siguiente.
#
#   ./scripts/drive_pull.sh sra prupe --go
#   ./03_stream_align.sh prupe
#   ./scripts/drive_pull.sh purge sra prupe --go
#
# 'purge' borra SOLO la copia local. Nunca toca Drive.

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$DIR/_drive_lib.sh"

REMOTE="${DRIVE_REMOTE:-$DRIVE_REMOTE_DEFAULT}"
# OJO con `-` y no `:-`: DRIVE_ROOT='' es una configuracion valida
# (remoto con root_folder_id ya apuntando a tesis/), y `:-` la pisaria.
DRIVE_ROOT="${DRIVE_ROOT-tesis}"
LOCAL_ROOT="${LOCAL_ROOT:-$(cd "$DIR/.." && pwd)}"

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }

[[ $# -ge 1 ]] || usage
case "$1" in -h|--help) usage 0 ;; esac

PURGE=0
if [[ "$1" == "purge" ]]; then PURGE=1; shift; fi
[[ $# -ge 1 ]] || usage
parse_args_drive "$@" || usage

SRC_PATH="$REMOTE:${DRIVE_ROOT:+$DRIVE_ROOT/}$DST"
DST_PATH="$LOCAL_ROOT/$SRC"

if [[ $PURGE -eq 1 ]]; then
  # Exigir organismo: 'purge sra' sin ORG borraría el caché entero.
  [[ -n "$ORG" ]] || { echo "purge necesita un ORG explícito" >&2; exit 2; }
  [[ -d "$DST_PATH" ]] || { echo "nada que purgar: $DST_PATH no existe"; exit 0; }
  local_n=$(find "$DST_PATH" -type f | wc -l)
  local_sz=$(du -sh "$DST_PATH" 2>/dev/null | cut -f1)

  # "Drive intacto" no es lo mismo que "seguro". Para sra y genomas, Drive es
  # la fuente y la copia local es descartable. Para el resto, la copia local es
  # lo que se PRODUJO acá: si todavía no se subió, purgar la borra y no hay de
  # dónde recuperarla. Los BAMs de un organismo son horas de alineamiento.
  #
  # `rclone check` compara por checksum contra el remoto y sale != 0 si falta
  # algo o difiere. Es lo mismo que haría un push, pero sin escribir.
  case "$FASE" in
    sra|genomas) ;;   # Drive es la fuente; nada que confirmar
    *)
      echo ">> confirmando que $SRC_PATH tiene lo que hay acá, antes de borrar"
      if ! rclone check "$DST_PATH" "$SRC_PATH" --checksum --one-way; then
        echo >&2
        echo "NO purgo: Drive no tiene todo lo que hay en $DST_PATH." >&2
        echo "  Es data producida acá; borrarla ahora la pierde." >&2
        echo "  Subila primero: ./scripts/drive_push.sh $FASE $ORG --go" >&2
        exit 1
      fi
      echo "   está todo en Drive."
      ;;
  esac

  if [[ $GO -eq 1 ]]; then
    echo ">> borrando LOCAL $DST_PATH ($local_n archivos, $local_sz)"
    rm -rf "$DST_PATH"
    echo "   listo. Drive intacto."
  else
    echo ">> DRY-RUN borraría LOCAL $DST_PATH ($local_n archivos, $local_sz)"
    echo "   agregá --go para borrar. Drive no se toca."
  fi
  exit 0
fi

ARGS=(copy "$SRC_PATH" "$DST_PATH" --checksum --progress --transfers 4)

if [[ $GO -eq 1 ]]; then
  mkdir -p "$DST_PATH"   # solo al bajar de verdad: un dry-run no deja rastro
  echo ">> bajando $SRC_PATH -> $DST_PATH"
  rclone "${ARGS[@]}"
else
  echo ">> DRY-RUN $SRC_PATH -> $DST_PATH   (agregá --go para bajar)"
  rclone "${ARGS[@]}" --dry-run
fi
