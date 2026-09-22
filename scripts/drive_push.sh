#!/usr/bin/env bash
# Sube resultados a Mi unidad/tesis/ en Drive. Ver data/DRIVE.md.
#
#   ./scripts/drive_push.sh <fase> [ORG]        # dry-run, no sube nada
#   ./scripts/drive_push.sh <fase> [ORG] --go   # sube de verdad
#
#   fases: bam | yasma | qc | features | modelos | figuras | genomas | sra
#   ORG:   rhirr | sclsc | cloro | phypa | prupe | maldo | gadmo | galga | maggi
#
# Por qué rclone y no el conector de Drive: son cientos de GB. El conector
# mueve manifiestos, no alineamientos.
#
# Requiere un remoto rclone de tipo drive apuntando a seb.ugazm@gmail.com:
#   rclone config       # nombre del remoto: gdrive-tesis
# El rclone.conf NO va al repo (está en .gitignore).

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$DIR/_drive_lib.sh"

REMOTE="${DRIVE_REMOTE:-$DRIVE_REMOTE_DEFAULT}"
# OJO con `-` y no `:-`: DRIVE_ROOT='' es una configuracion valida
# (remoto con root_folder_id ya apuntando a tesis/), y `:-` la pisaria.
DRIVE_ROOT="${DRIVE_ROOT-tesis}"
LOCAL_ROOT="${LOCAL_ROOT:-$(cd "$DIR/.." && pwd)}"

usage() { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }

[[ $# -ge 1 ]] || usage
case "$1" in -h|--help) usage 0 ;; esac
parse_args_drive "$@" || usage

SRC_PATH="$LOCAL_ROOT/$SRC"
DST_PATH="$REMOTE:${DRIVE_ROOT:+$DRIVE_ROOT/}$DST"
[[ -d "$SRC_PATH" ]] || { echo "no existe: $SRC_PATH" >&2; exit 1; }

# --checksum en vez de --size-only: un BAM truncado por un corte de red tiene
# tamaño distinto, pero uno re-generado con otro orden de sort puede tener el
# mismo tamaño y contenido distinto.
ARGS=(copy "$SRC_PATH" "$DST_PATH" --checksum --progress --transfers 4 --drive-chunk-size 64M)

if [[ $GO -eq 1 ]]; then
  echo ">> subiendo $SRC_PATH -> $DST_PATH"
  rclone "${ARGS[@]}"
else
  echo ">> DRY-RUN $SRC_PATH -> $DST_PATH   (agregá --go para subir)"
  rclone "${ARGS[@]}" --dry-run
fi
