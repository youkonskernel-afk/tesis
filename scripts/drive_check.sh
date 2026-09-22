#!/usr/bin/env bash
# Comprueba que rclone este configurado como los scripts de este repo esperan.
#
#   ./scripts/drive_check.sh
#
# Existe porque "segui los pasos" no es lo mismo que "funciona", y los dos modos
# de fallar mas comunes NO dan un error claro:
#
#   1. Elegir el scope `drive.file` en vez de `drive`. rclone solo ve entonces
#      los ficheros que el mismo creo, asi que las carpetas de tesis/ —que se
#      crearon a mano en la web— son INVISIBLES. El sintoma es un remoto que
#      lista vacio, no un error de permisos.
#   2. Poner root_folder_id apuntando a tesis/ Y dejar DRIVE_ROOT=tesis. rclone
#      busca entonces tesis/tesis/80_sra, que no existe. Tambien lista vacio.
#
# En los dos casos `drive_pull.sh sra <org> --go` baja 0 ficheros y sale con
# codigo 0. Este script los distingue.
set -uo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=_drive_lib.sh
. "$DIR/_drive_lib.sh"

REMOTE="${DRIVE_REMOTE:-$DRIVE_REMOTE_DEFAULT}"
# OJO con `-` y no `:-`: DRIVE_ROOT='' es una configuracion valida
# (remoto con root_folder_id ya apuntando a tesis/), y `:-` la pisaria.
DRIVE_ROOT="${DRIVE_ROOT-tesis}"
CUENTA_ESPERADA="${DRIVE_CUENTA:-seb.ugazm@gmail.com}"

FALLAS=0
ok()  { printf '  [ok]  %s\n' "$1"; }
mal() { printf '  [MAL] %s\n' "$1"; FALLAS=$((FALLAS+1)); }
nota(){ printf '        %s\n' "$1"; }

echo "remoto esperado : $REMOTE"
echo "raiz en Drive   : $REMOTE:$DRIVE_ROOT"
echo

echo "== 1. rclone instalado"
if command -v rclone >/dev/null; then
  ok "$(rclone version 2>/dev/null | head -1)"
else
  mal "no esta en el PATH"
  nota "Debian/Ubuntu: sudo -v && curl https://rclone.org/install.sh | sudo bash"
  nota "macOS:         brew install rclone"
  echo; echo "$FALLAS problema(s). Ver docs/rclone.md"; exit 1
fi

echo "== 2. el remoto existe y es de tipo drive"
if rclone listremotes 2>/dev/null | grep -qx "$REMOTE:"; then
  tipo=$(rclone config dump 2>/dev/null \
         | python3 -c "import json,sys; print(json.load(sys.stdin).get('$REMOTE',{}).get('type','?'))")
  if [[ "$tipo" == "drive" ]]; then
    ok "existe y es type=drive"
  else
    mal "existe pero type=$tipo, no drive"
  fi
else
  mal "no hay un remoto llamado '$REMOTE'"
  nota "los que hay: $(rclone listremotes 2>/dev/null | tr '\n' ' ')"
  nota "creá uno con: rclone config   (ver docs/rclone.md)"
  echo; echo "$FALLAS problema(s). Ver docs/rclone.md"; exit 1
fi

echo "== 3. el scope permite ver lo que ya existe"
# drive.file solo ve lo que rclone creo. Las carpetas de tesis/ se crearon a
# mano, asi que con ese scope el remoto lista vacio y nada lo dice.
scope=$(rclone config dump 2>/dev/null \
        | python3 -c "import json,sys; print(json.load(sys.stdin).get('$REMOTE',{}).get('scope','(sin declarar)'))")
case "$scope" in
  drive|'(sin declarar)') ok "scope=$scope" ;;
  drive.file)
    mal "scope=drive.file: rclone solo ve lo que el mismo creo"
    nota "las carpetas de tesis/ se crearon a mano, asi que son invisibles"
    nota "arreglo: rclone config -> editar $REMOTE -> scope 1 (drive)" ;;
  *) mal "scope=$scope, se espera drive" ;;
esac

echo "== 4. root_folder_id no choca con DRIVE_ROOT"
rfid=$(rclone config dump 2>/dev/null \
       | python3 -c "import json,sys; print(json.load(sys.stdin).get('$REMOTE',{}).get('root_folder_id',''))")
if [[ -z "$rfid" ]]; then
  ok "sin root_folder_id (los scripts direccionan por ruta)"
elif [[ -n "$DRIVE_ROOT" ]]; then
  mal "root_folder_id=$rfid Y DRIVE_ROOT=$DRIVE_ROOT"
  nota "rclone buscaria $DRIVE_ROOT/$DRIVE_ROOT/..., que no existe"
  nota "arreglo: sacá el root_folder_id, o exportá DRIVE_ROOT=''"
else
  ok "root_folder_id=$rfid con DRIVE_ROOT vacio"
fi

echo "== 5. la cuenta es la correcta"
cuenta=$(rclone about "$REMOTE:" --json 2>/dev/null >/dev/null && \
         rclone backend drives "$REMOTE:" 2>/dev/null | head -1)
email=$(rclone config userinfo "$REMOTE:" 2>/dev/null \
        | sed -n 's/.*"emailAddress": *"\([^"]*\)".*/\1/p' | head -1)
if [[ -z "$email" ]]; then
  nota "no pude leer el email (rclone config userinfo no respondio); sigo"
elif [[ "$email" == "$CUENTA_ESPERADA" ]]; then
  ok "$email"
else
  mal "autenticado como $email, se espera $CUENTA_ESPERADA"
fi

echo "== 6. la raiz $DRIVE_ROOT/ se ve"
if ! salida=$(rclone lsd "$REMOTE:${DRIVE_ROOT:-}" 2>&1); then
  mal "no pude listar $REMOTE:$DRIVE_ROOT"
  nota "$(head -2 <<<"$salida" | tr '\n' ' ')"
else
  n=$(grep -c . <<<"$salida")
  if [[ "$n" -eq 0 ]]; then
    mal "$REMOTE:$DRIVE_ROOT lista VACIO"
    nota "casi siempre es el scope (chequeo 3) o el root_folder_id (chequeo 4)"
    nota "una raiz vacia hace que drive_pull baje 0 ficheros y salga con 0"
  else
    ok "$n subcarpetas"
  fi
fi

echo "== 7. estan las carpetas que el mapa de fases espera"
for fase in bam yasma qc features modelos figuras genomas sra; do
  sub=$(fase_a_rutas "$fase" | cut -d'|' -f2)
  if rclone lsd "$REMOTE:${DRIVE_ROOT:+$DRIVE_ROOT/}$sub" >/dev/null 2>&1; then
    ok "$sub/"
  else
    mal "$sub/ no se ve (fase '$fase')"
  fi
done

echo "== 8. espacio"
if libre=$(rclone about "$REMOTE:" 2>/dev/null | sed -n 's/^Free: *//p'); then
  [[ -n "$libre" ]] && nota "libre: $libre  (BAMs + .sra + genomas son ~540 GB)"
fi

echo
if [[ $FALLAS -eq 0 ]]; then
  echo "TODO OK — probá un dry-run real:"
  echo "  ./scripts/drive_pull.sh sra prupe      # sin --go, no baja nada"
  exit 0
fi
echo "$FALLAS problema(s). Ver docs/rclone.md"
exit 1
