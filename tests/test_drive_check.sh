#!/usr/bin/env bash
# drive_check.sh contra un rclone falso, en los escenarios que importan.
#
# Los dos modos de fallo que este script existe para distinguir NO dan un error:
#   - scope=drive.file  -> rclone solo ve lo que el mismo creo, y las carpetas
#     de tesis/ se crearon a mano. El remoto lista VACIO.
#   - root_folder_id apuntando a tesis/ con DRIVE_ROOT=tesis -> busca
#     tesis/tesis/... Tambien lista vacio.
# En los dos casos `drive_pull.sh sra <org> --go` baja 0 ficheros y sale con 0.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
tiene(){ grep -qF -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
notiene(){ grep -qF -- "$2" <<<"$3" && mal "$1 — no debia: $2" || ok "$1"; }

mkdir -p "$TMP/bin"
# rclone falso. Lo que responde se controla con variables:
#   FAKE_SCOPE, FAKE_RFID, FAKE_EMAIL, FAKE_VACIO (1 = la raiz lista vacia)
cat > "$TMP/bin/rclone" <<'RC'
#!/usr/bin/env bash
case "$1 ${2:-}" in
  "version"*) echo "rclone v1.65.0"; exit 0 ;;
  "listremotes"*) echo "gdrive-tesis:"; exit 0 ;;
  "config dump")
    printf '{"gdrive-tesis":{"type":"drive","scope":"%s","root_folder_id":"%s"}}\n' \
      "${FAKE_SCOPE:-drive}" "${FAKE_RFID:-}"
    exit 0 ;;
  "config userinfo")
    printf '{ "emailAddress": "%s" }\n' "${FAKE_EMAIL:-seb.ugazm@gmail.com}"; exit 0 ;;
  "about"*) echo "Free: 1.100T"; exit 0 ;;
  "lsd"*)
    [[ "${FAKE_VACIO:-0}" == "1" ]] && exit 0
    # solo conoce las carpetas del mapa
    case "$*" in
      *"tesis/10_bam"*|*"tesis/20_yasma"*|*"tesis/30_qc"*|*"tesis/40_features"*|\
      *"tesis/50_modelos"*|*"tesis/60_figuras"*|*"tesis/70_genomas"*|*"tesis/80_sra"*)
        echo "          -1 2026-01-01 00:00:00        -1 prupe"; exit 0 ;;
      *"tesis"*) printf '%s\n' "  -1 x 10_bam" "  -1 x 20_yasma" "  -1 x 30_qc" \
          "  -1 x 40_features" "  -1 x 50_modelos" "  -1 x 60_figuras" \
          "  -1 x 70_genomas" "  -1 x 80_sra"; exit 0 ;;
      *) echo "directory not found" >&2; exit 3 ;;
    esac ;;
  "backend drives") exit 0 ;;
esac
exit 0
RC
chmod +x "$TMP/bin/rclone"
corre() { env PATH="$TMP/bin:$PATH" "$@" bash "$RAIZ/scripts/drive_check.sh" 2>&1; }

echo "== 1. todo bien"
S=$(corre); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
tiene "ve el remoto"        "existe y es type=drive"  "$S"
tiene "ve las 8 carpetas"   "80_sra/"                 "$S"
tiene "y cierra bien"       "TODO OK"                 "$S"

echo "== 2. scope=drive.file: el fallo que no da error"
S=$(corre FAKE_SCOPE=drive.file); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo nombra"       "scope=drive.file"                      "$S"
tiene "y explica por que" "solo ve lo que el mismo creo"        "$S"
tiene "y como arreglarlo" "scope 1 (drive)"                     "$S"

echo "== 3. root_folder_id chocando con DRIVE_ROOT"
S=$(corre FAKE_RFID=1E_Q6XLg4); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "explica el doble prefijo" "tesis/tesis/"   "$S"

echo "== 3b. pero con DRIVE_ROOT vacio, root_folder_id esta bien"
S=$(corre FAKE_RFID=1E_Q6XLg4 DRIVE_ROOT=)
tiene "lo acepta" "root_folder_id=1E_Q6XLg4 con DRIVE_ROOT vacio" "$S"

echo "== 4. la cuenta equivocada"
S=$(corre FAKE_EMAIL=otro@gmail.com); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "dice cual es"  "autenticado como otro@gmail.com"  "$S"

echo "== 5. la raiz lista vacia"
S=$(corre FAKE_VACIO=1); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo dice"            "lista VACIO"                       "$S"
tiene "y manda a los dos culpables" "scope (chequeo 3)"        "$S"
tiene "y avisa del sintoma" "baje 0 ficheros y salga con 0"    "$S"

echo "== 6. sin el remoto, no sigue adelante"
cat > "$TMP/bin/rclone" <<'RC2'
#!/usr/bin/env bash
case "$1" in version) echo "rclone v1.65.0" ;; listremotes) echo "otro:" ;; esac
exit 0
RC2
chmod +x "$TMP/bin/rclone"
S=$(corre); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "dice cual falta"  "no hay un remoto llamado"  "$S"
tiene "y lista los que hay" "otro:"                  "$S"
notiene "y no sigue chequeando" "80_sra/"            "$S"

echo "== 7. sin rclone, dice como instalarlo"
rm "$TMP/bin/rclone"
S=$(corre); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo dice"            "no esta en el PATH"  "$S"
tiene "y como instalarlo"  "rclone.org/install.sh" "$S"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
