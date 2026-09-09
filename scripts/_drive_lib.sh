# Compartido por drive_push.sh y drive_pull.sh. No se ejecuta solo.
#
# El mapa de fases y la lista de organismos viven acá y en un solo lugar: si
# push y pull se desincronizaran, subirías a una carpeta y bajarías de otra.

DRIVE_REMOTE_DEFAULT="gdrive-tesis"
ORGS_VALIDOS="rhirr sclsc cloro phypa prupe maldo gadmo galga maggi"

# fase -> "<subdir local>|<subdir en Drive>"
fase_a_rutas() {
  case "$1" in
    bam)      echo "bams|10_bam" ;;
    yasma)    echo "yasma_out|20_yasma" ;;
    qc)       echo "qc|30_qc" ;;
    features) echo "features|40_features" ;;
    modelos)  echo "models|50_modelos" ;;
    figuras)  echo "figures|60_figuras" ;;
    genomas)  echo "genomes|70_genomas" ;;
    sra)      echo "sra_cache|80_sra" ;;
    *) return 1 ;;
  esac
}

org_valido() {
  local o
  for o in $ORGS_VALIDOS; do [[ "$o" == "$1" ]] && return 0; done
  return 1
}

# Parsea "<fase> [ORG] [--go]" y deja FASE, ORG, GO, SRC, DST.
parse_args_drive() {
  [[ $# -ge 1 ]] || return 2
  FASE="$1"; shift
  ORG=""; GO=0
  local arg
  for arg in "$@"; do
    case "$arg" in
      --go) GO=1 ;;
      *) ORG="$arg" ;;
    esac
  done

  local rutas
  rutas=$(fase_a_rutas "$FASE") || { echo "fase desconocida: $FASE" >&2; return 2; }
  SRC="${rutas%%|*}"; DST="${rutas##*|}"

  if [[ -n "$ORG" ]]; then
    org_valido "$ORG" || { echo "organismo desconocido: $ORG" >&2; return 2; }
    SRC="$SRC/$ORG"; DST="$DST/$ORG"
  fi
  return 0
}
