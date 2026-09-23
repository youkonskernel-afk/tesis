#!/usr/bin/env bash
# Que los cinco scripts crean que la data local vive en el MISMO lugar.
#
# Por que existe: habia tres respuestas distintas para donde estan los .sra.
#   drive_pull.sh los dejaba en  <repo>/sra_cache/<org>/
#   fetch_runs.sh los buscaba en /home/dev/sra_cache
#   trim.sh       los buscaba en $HOME/tesis_data/80_sra/<org>/
# Ninguna falla ruidosamente: el sintoma es "FALTA el .sra" en las 417 despues
# de haberlas bajado bien, o un `estado` que dice que no hay nada.
#
# `/home/dev/` ademas es un absoluto con un usuario que no existe en ninguna
# maquina de este proyecto, asi que no habia forma de que funcionara.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
igual(){ [[ "$2" == "$3" ]] && ok "$1" || mal "$1
         uno dice: $2
         otro    : $3"; }

# Lo que cada script resuelve, preguntandoselo a el y no leyendo su fuente.
ruta_de_pull() {   # a donde BAJA drive_pull.sh
  LOCAL_ROOT="$1" bash -c '
    . "'"$RAIZ"'/scripts/_drive_lib.sh"
    parse_args_drive "$@" >/dev/null
    echo "${LOCAL_ROOT}/${SRC}"' _ "$2" "${3:-}"
}
ruta_de_script() { # donde BUSCA un script, via su variable interna
  LOCAL_ROOT="$1" bash -c '
    set -a
    src=$(sed -n "/^'"$3"'=/p" "'"$RAIZ"'/scripts/'"$2"'")
    ROOT="'"$RAIZ"'"
    . "'"$RAIZ"'/scripts/_drive_lib.sh"
    eval "$src"
    eval echo "\$'"$3"'"'
}

L="$TMP/data"

echo "== 1. los .sra: donde los deja pull y donde los busca cada uno"
P=$(ruta_de_pull "$L" sra)
F=$(ruta_de_script "$L" fetch_runs.sh CACHE)
T=$(ruta_de_script "$L" trim.sh SRA_DEST)
echo "   drive_pull deja en : $P"
echo "   fetch_runs busca en: $F"
echo "   trim.sh    busca en: $T"
igual "fetch_runs coincide con drive_pull" "$F" "$P"
igual "trim.sh coincide con drive_pull"    "$T" "$P"

echo "== 2. ninguna ruta absoluta con un usuario cableado"
# Solo en codigo: los comentarios NOMBRAN la ruta vieja a proposito, para que
# se entienda por que existe esta regla.
_cableadas=0
for f in "$RAIZ"/scripts/*.sh; do
  hit=$(grep -vE '^\s*#' "$f" | grep -oE '/home/(dev|ubuntu)/[a-zA-Z_/]*' | head -1)
  [[ -n "$hit" ]] && { mal "$(basename "$f") cablea $hit"; _cableadas=1; }
done
[[ $_cableadas -eq 0 ]] && ok "ningun /home/<usuario> en el codigo de scripts/"

echo "== 3. los genomas: pull y fetch_genomes"
PG=$(ruta_de_pull "$L" genomas)
FG=$(ruta_de_script "$L" fetch_genomes.sh DEST)
igual "fetch_genomes coincide con drive_pull" "$FG" "$PG"

echo "== 4. LOCAL_ROOT manda sobre todos por igual"
P2=$(ruta_de_pull "$TMP/otra" sra)
F2=$(ruta_de_script "$TMP/otra" fetch_runs.sh CACHE)
igual "pull respeta LOCAL_ROOT"       "$P2" "$TMP/otra/sra_cache"
igual "y fetch_runs tambien"          "$F2" "$P2"

echo "== 5. cada subdir del mapa esta en .gitignore"
# Son directorios de data adentro del repo: si uno se escapa, `git status`
# muestra cientos de GB y algo termina commiteado.
GI="$RAIZ/.gitignore"
for fase in bam yasma qc features modelos figuras genomas sra; do
  sub=$(bash -c '. "'"$RAIZ"'/scripts/_drive_lib.sh"; fase_a_rutas "$1" | cut -d"|" -f1' _ "$fase")
  if grep -qxF "$sub/" "$GI"; then ok "$sub/ ignorado"
  else mal "$sub/ NO esta en .gitignore"; fi
done
# y el de los proyectos YASMA, que no es fase de Drive pero tambien es data
TD=$(ruta_de_script "$L" trim.sh PROY_DIR)
sub_proy=$(basename "$TD")
grep -qxF "$sub_proy/" "$GI" && ok "$sub_proy/ ignorado" \
  || mal "$sub_proy/ NO esta en .gitignore (PROY_DIR=$TD)"

echo "== 6. trim.sh y align.sh escriben y leen el MISMO proyecto"
# trim.sh deja <org>_<rol>/trim/*.t.fq.gz y align.sh los busca ahi. Si cada uno
# tuviera su idea de donde esta el directorio, align no encontraria nada
# recortado y lo diria como "NADA RECORTADO", que se lee como un problema del
# recorte y no de la ruta.
AD=$(ruta_de_script "$L" align.sh PROY_DIR)
igual "align.sh coincide con trim.sh"  "$AD" "$TD"

echo "== 7. align.sh toma el genoma y el BAM del mapa de fases"
AG=$(ruta_de_script "$L" align.sh GENOMES_DIR)
igual "el genoma, de donde lo deja fetch_genomes" "$AG" "$(ruta_de_pull "$L" genomas)"
AB=$(ruta_de_script "$L" align.sh BAM_DIR)
igual "y el BAM, de donde lo sube drive_push"     "$AB" "$(ruta_de_pull "$L" bam)"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
