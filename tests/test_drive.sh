#!/usr/bin/env bash
# drive_push.sh y drive_pull.sh con un rclone falso.
#
# No tenian banco, y son los dos scripts que mueven cientos de GB. El que mas
# importa es `purge`: borra la copia LOCAL, y para las fases que se PRODUCEN
# aca —bam, yasma, features, modelos— esa copia puede ser la unica que existe.
# "Drive intacto" no es lo mismo que "seguro": los BAMs de un organismo son
# horas de alineamiento.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
tiene(){ grep -qF -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
notiene(){ grep -qF -- "$2" <<<"$3" && mal "$1 — no debia: $2" || ok "$1"; }

mkdir -p "$TMP/bin"
# rclone falso: registra como lo llamaron, y `check` falla o no segun RCLONE_CHECK
cat > "$TMP/bin/rclone" <<'RC'
#!/usr/bin/env bash
echo "rclone $*" >> "$LOG_RC"
if [[ "${1:-}" == "check" ]]; then
  exit "${RCLONE_CHECK:-0}"
fi
exit 0
RC
chmod +x "$TMP/bin/rclone"
export PATH="$TMP/bin:$PATH" LOG_RC="$TMP/rc.log"

L="$TMP/data"
corre() { : > "$LOG_RC"; LOCAL_ROOT="$L" bash "$RAIZ/scripts/$1" "${@:2}" 2>&1; }

echo "== 1. push: dry-run por defecto, sube solo con --go"
mkdir -p "$L/bams/gadmo"; echo x > "$L/bams/gadmo/a.bam"
S=$(corre drive_push.sh bam gadmo)
tiene "dice que es dry-run" "DRY-RUN"          "$S"
tiene "y se lo pasa a rclone" "--dry-run"      "$(cat "$LOG_RC")"
S=$(corre drive_push.sh bam gadmo --go)
tiene "con --go sube"      "subiendo"          "$S"
notiene "y sin --dry-run"  "--dry-run"         "$(cat "$LOG_RC")"

echo "== 2. push: --checksum, no --size-only"
# Un BAM re-generado con otro orden de sort puede pesar igual y ser distinto.
tiene "usa --checksum"     "--checksum"        "$(cat "$LOG_RC")"
notiene "y no --size-only" "--size-only"       "$(cat "$LOG_RC")"

echo "== 3. push: si no existe el origen, no llama a rclone"
S=$(corre drive_push.sh features gadmo); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo dice"            "no existe"         "$S"
[[ ! -s "$LOG_RC" ]] && ok "no llamo a rclone" || mal "no llamo a rclone"

echo "== 4. pull: baja de Drive al mismo sitio donde los otros scripts buscan"
S=$(corre drive_pull.sh sra prupe --go)
tiene "destino local" "$L/sra_cache/prupe"     "$(cat "$LOG_RC")"
tiene "origen remoto" "tesis/80_sra/prupe"     "$(cat "$LOG_RC")"

echo "== 5. purge exige ORG: 'purge sra' borraria el cache entero"
S=$(corre drive_pull.sh purge sra --go); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo dice"  "necesita un ORG"             "$S"

echo "== 6. purge de una fase PRODUCIDA aca: confirma contra Drive antes"
# Es el hallazgo que motivo este banco: borraba sin comprobar nada.
mkdir -p "$L/bams/gadmo"; echo x > "$L/bams/gadmo/a.bam"
RCLONE_CHECK=1
S=$(RCLONE_CHECK=1 LOCAL_ROOT="$L" bash "$RAIZ/scripts/drive_pull.sh" purge bam gadmo --go 2>&1); RC=$?
[[ $RC -ne 0 ]] && ok "si Drive no lo tiene, exit != 0" || mal "si Drive no lo tiene, exit != 0 (rc=$RC)"
tiene "explica por que"       "borrarla ahora la pierde" "$S"
tiene "y como arreglarlo"     "drive_push.sh bam gadmo"  "$S"
[[ -f "$L/bams/gadmo/a.bam" ]] && ok "NO borro nada" || mal "NO borro nada"

echo "== 7. y si Drive SI lo tiene, purga"
S=$(RCLONE_CHECK=0 LOCAL_ROOT="$L" bash "$RAIZ/scripts/drive_pull.sh" purge bam gadmo --go 2>&1)
tiene "lo confirma"  "está todo en Drive"      "$S"
[[ ! -d "$L/bams/gadmo" ]] && ok "y borro la copia local" || mal "y borro la copia local"

echo "== 8. purge de sra no confirma: Drive es la fuente, es descartable"
mkdir -p "$L/sra_cache/prupe"; echo x > "$L/sra_cache/prupe/a.sra"
S=$(RCLONE_CHECK=1 LOCAL_ROOT="$L" bash "$RAIZ/scripts/drive_pull.sh" purge sra prupe --go 2>&1)
notiene "no confirma"  "confirmando"           "$S"
[[ ! -d "$L/sra_cache/prupe" ]] && ok "y purga igual" || mal "y purga igual"

echo "== 9. purge sin --go no borra"
mkdir -p "$L/sra_cache/gadmo"; echo x > "$L/sra_cache/gadmo/a.sra"
S=$(corre drive_pull.sh purge sra gadmo)
tiene "dice que es dry-run"  "DRY-RUN"         "$S"
[[ -f "$L/sra_cache/gadmo/a.sra" ]] && ok "y no borro" || mal "y no borro"

echo "== 10. una fase o un organismo que no existen fallan"
tiene "fase mala"  "fase desconocida"  "$(corre drive_push.sh nosequé gadmo 2>&1 || true)"
tiene "org malo"   "organismo desconocido" "$(corre drive_push.sh bam nosequé 2>&1 || true)"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
