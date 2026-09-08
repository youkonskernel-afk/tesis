#!/usr/bin/env bash
# Sube resultados a Mi unidad/tesis/ en Drive. Ver data/DRIVE.md.
#
#   ./scripts/drive_push.sh <fase> [ORG]        # dry-run, no sube nada
#   ./scripts/drive_push.sh <fase> [ORG] --go   # sube de verdad
#
#   fases: bam | yasma | qc | features | modelos | figuras
#   ORG:   rhirr | sclsc | arath | phypa | danre | nemve   (opcional)
#
# Por qué rclone y no el conector de Drive: son ~100 GB de BAMs. El conector
# mueve manifiestos, no alineamientos.
#
# Requiere un remoto rclone de tipo drive apuntando a seb.ugazm@gmail.com:
#   rclone config       # nombre del remoto: gdrive-tesis
# El rclone.conf NO va al repo (está en .gitignore).

set -euo pipefail

REMOTE="${DRIVE_REMOTE:-gdrive-tesis}"
DRIVE_ROOT="${DRIVE_ROOT:-tesis}"
LOCAL_ROOT="${LOCAL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

usage() { sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }

[[ $# -ge 1 ]] || usage
case "$1" in -h|--help) usage 0 ;; esac
FASE="$1"; shift

ORG=""
GO=0
for arg in "$@"; do
  case "$arg" in
    --go) GO=1 ;;
    -h|--help) usage 0 ;;
    *) ORG="$arg" ;;
  esac
done

# fase -> (subdir local, subdir en Drive)
case "$FASE" in
  bam)      SRC="bams";        DST="10_bam" ;;
  yasma)    SRC="yasma_out";   DST="20_yasma" ;;
  qc)       SRC="qc";          DST="30_qc" ;;
  features) SRC="features";    DST="40_features" ;;
  modelos)  SRC="models";      DST="50_modelos" ;;
  figuras)  SRC="figures";     DST="60_figuras" ;;
  *) echo "fase desconocida: $FASE" >&2; usage ;;
esac

if [[ -n "$ORG" ]]; then
  case "$ORG" in
    rhirr|sclsc|arath|phypa|danre|nemve) ;;
    *) echo "organismo desconocido: $ORG" >&2; exit 2 ;;
  esac
  SRC="$SRC/$ORG"; DST="$DST/$ORG"
fi

SRC_PATH="$LOCAL_ROOT/$SRC"
DST_PATH="$REMOTE:$DRIVE_ROOT/$DST"

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
