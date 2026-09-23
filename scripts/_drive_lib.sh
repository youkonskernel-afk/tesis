# shellcheck shell=bash
# Compartido por drive_push.sh, drive_pull.sh, fetch_runs.sh, trim.sh y align.sh.
# No se ejecuta solo.
#
# El mapa de fases, la lista de organismos y LA RUTA LOCAL viven acá y en un
# solo lugar: si push y pull se desincronizaran, subirías a una carpeta y
# bajarías de otra.
#
# Eso mismo paso con los .sra, una vuelta mas ancha: drive_pull los dejaba en
# <repo>/sra_cache/<org>/, fetch_runs los buscaba en /home/dev/sra_cache —una
# ruta absoluta con un usuario que no existe en ninguna maquina de este
# proyecto— y trim.sh en $HOME/tesis_data/80_sra/<org>/. Tres lugares para la
# misma data, y el sintoma habria sido "FALTA el .sra" en las 417 despues de
# bajarlas bien.

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

# Raiz de la data local. La usan los cuatro scripts; se cambia con LOCAL_ROOT.
# Por defecto es la raiz del repo, y .gitignore cubre todos los subdirectorios
# del mapa de arriba.
ruta_local() {
  local fase="$1" org="${2:-}" rutas raiz
  rutas=$(fase_a_rutas "$fase") || return 1
  raiz="${LOCAL_ROOT:-$(cd "${BASH_SOURCE[0]%/*}/.." && pwd)}"
  printf '%s/%s%s\n' "$raiz" "${rutas%%|*}" "${org:+/$org}"
}

# Raiz de los PROYECTOS YASMA: un directorio por organismo y rol, con inputs.json,
# trim/, align/ y annotations/ adentro. No es una fase de Drive —lo de adentro se
# re-genera de los .sra y del genoma, que si estan respaldados— pero vive aca por
# el mismo motivo que ruta_local: trim.sh escribe y align.sh lee, y si cada uno
# tuviera su idea de donde estan, align no encontraria nada recortado.
#
# Se llamaba trim/ cuando solo guardaba el recorte. Con `yasma align` escribiendo
# align/alignment.bam adentro, ese nombre pasaba a mentir.
ruta_proyectos() {
  local raiz="${LOCAL_ROOT:-$(cd "${BASH_SOURCE[0]%/*}/.." && pwd)}"
  printf '%s/proyectos%s\n' "$raiz" "${1:+/$1}"
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
