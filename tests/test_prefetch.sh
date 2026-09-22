#!/usr/bin/env bash
# Banco de pruebas de cmd_prefetch con un prefetch falso en el PATH.
# Misma tecnica que ya se uso para el filtro de RNA y el flujo de clon.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$RAIZ"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0

ok()   { printf '  ok   %s\n' "$1"; }
mal()  { printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
chk()  { if [[ "$2" == "$3" ]]; then ok "$1"; else mal "$1 (esperaba '$3', obtuve '$2')"; fi; }
tiene(){ if grep -qF -- "$2" <<<"$3"; then ok "$1"; else mal "$1 — no aparece: $2"; fi; }
# columna del ledger para una corrida: fila(caso) campo
campo(){ awk -F'	' -v r="$2" -v n="$3" '$2==r {print $n}' "$TMP/$1/data/ledger.tsv" 2>/dev/null; }
notiene(){ if grep -qF -- "$2" <<<"$3"; then mal "$1 — no debia aparecer: $2"; else ok "$1"; fi; }

# --- escenario: monta un PATH con un prefetch falso y corre fetch_runs.sh ---
correr() {
  local nombre="$1" cuerpo="$2" exitcode="${3:-0}"
  ESC="$TMP/$nombre"
  mkdir -p "$ESC/bin" "$ESC/dest" "$ESC/staging" "$ESC/data"

  printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n' > "$ESC/data/man.tsv"
  printf 'maggi\tSRR317135\tPRJNA154615\tprimario\tentrenamiento\t14311812\t701278788\t49\tRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n' >> "$ESC/data/man.tsv"

  cat > "$ESC/bin/prefetch" <<PREFETCH
#!/usr/bin/env bash
# args: --output-directory DIR --max-size u RUN
OUT="\$2"; RUN="\${@: -1}"
$cuerpo
exit $exitcode
PREFETCH
  cat > "$ESC/bin/vdb-validate" <<'VAL'
#!/usr/bin/env bash
exit 0
VAL
  chmod +x "$ESC/bin/prefetch" "$ESC/bin/vdb-validate"

  PATH="$ESC/bin:$PATH" \
  SRA_DEST="$ESC/dest" SRA_STAGING="$ESC/staging" \
  SRA_LEDGER="$ESC/data/ledger.tsv" MANIFEST="$ESC/data/man.tsv" \
  MIN_LIBRE_GB=0 \
    "$REPO/scripts/fetch_runs.sh" prefetch 2>&1
}

echo "== 1. prefetch sale 0 y no escribe nada"
SAL=$(correr caso1 'mkdir -p "$OUT/$RUN"' 0)
tiene    "reporta el fallo mudo"        "FALLO: prefetch salió 0 pero no dejó .sra" "$SAL"
tiene    "muestra el log de prefetch"   "prefetch no imprimió nada"                 "$SAL"
tiene    "lista el staging"             "quedó en el staging:"                      "$SAL"
tiene    "el staging estaba vacio"      "(nada)"                                    "$SAL"
tiene    "apunta al diagnostico"        "diag SRR317135"                            "$SAL"
tiene    "cuenta el fallo"              "bajadas=0 fallos=1"                        "$SAL"
if [[ -d "$TMP/caso1/staging/SRR317135" ]]; then mal "staging limpiado"; else ok "staging limpiado"; fi

echo "== 2. prefetch deja un .sralite"
SAL=$(correr caso2 'mkdir -p "$OUT/$RUN"; head -c 4096 /dev/zero > "$OUT/$RUN/$RUN.sralite"; echo "downloaded"' 0)
tiene    "se acepta"                    "bajadas=1 fallos=0"                        "$SAL"
tiene    "avisa del formato lite"       "vino en formato lite"                      "$SAL"
if [[ -s "$TMP/caso2/dest/maggi/SRR317135.sra" ]]; then ok "movido a DEST como .sra"; else mal "movido a DEST como .sra"; fi
if grep -q "SRR317135" "$TMP/caso2/data/ledger.tsv" 2>/dev/null; then ok "md5 en el ledger"; else mal "md5 en el ledger"; fi
chk "el ledger lo marca sralite"   "$(campo caso2 SRR317135 4)" "sralite"
chk "cabecera con formato"         "$(head -1 "$TMP/caso2/data/ledger.tsv")" "$(printf 'org	run	md5	formato	fecha_utc')"

echo "== 3. prefetch deja los originales del envio"
SAL=$(correr caso3 'mkdir -p "$OUT/$RUN"; head -c 4096 /dev/zero > "$OUT/$RUN/$RUN.fastq.gz"; echo "downloaded original format"' 0)
tiene    "se rechaza"                   "FALLO: prefetch salió 0 pero no dejó .sra" "$SAL"
tiene    "muestra el fastq encontrado"  "SRR317135.fastq.gz"                        "$SAL"
tiene    "muestra el log"               "downloaded original format"                "$SAL"
if [[ -e "$TMP/caso3/dest/maggi/SRR317135.sra" ]]; then mal "no se movio nada a DEST"; else ok "no se movio nada a DEST"; fi

echo "== 4. camino normal: .sra (no regresion)"
SAL=$(correr caso4 'mkdir -p "$OUT/$RUN"; head -c 4096 /dev/zero > "$OUT/$RUN/$RUN.sra"' 0)
tiene    "se baja"                      "bajadas=1 fallos=0"                        "$SAL"
notiene  "sin aviso de lite"            "formato lite"                              "$SAL"
if [[ -s "$TMP/caso4/dest/maggi/SRR317135.sra" ]]; then ok "movido a DEST"; else mal "movido a DEST"; fi
chk "el ledger lo marca sra"       "$(campo caso4 SRR317135 4)" "sra"

echo "== 5. prefetch sale != 0"
SAL=$(correr caso5 'echo "err: connection failed" >&2' 1)
tiene    "reporta exit != 0"            "FALLO prefetch (exit != 0)"                "$SAL"
tiene    "muestra el stderr capturado"  "connection failed"                         "$SAL"
tiene    "cuenta el fallo"              "bajadas=0 fallos=1"                        "$SAL"

echo "== 6. idempotencia: si ya esta en DEST, no se rebaja"
SAL=$(correr caso4 'exit 9' 0)
tiene    "nada pendiente"               "no hay nada pendiente"                     "$SAL"

echo "== 7. idempotencia del ledger: registrar dos veces no duplica"
ESC="$TMP/caso7"; mkdir -p "$ESC/data"
(
  set -euo pipefail
  LEDGER="$ESC/data/ledger.tsv"
  source <(sed -n '/^registrar_md5()/,/^}/p' "$REPO/scripts/fetch_runs.sh")
  registrar_md5 maggi SRR317135 aaa sralite
  registrar_md5 maggi SRR317135 bbb sralite
  registrar_md5 maggi SRR1066790 ccc sralite
)
chk "3 lineas (cabecera + 2)"      "$(wc -l < "$ESC/data/ledger.tsv")" "3"
chk "se quedo con el md5 nuevo"    "$(campo caso7 SRR317135 3)" "bbb"
chk "formato preservado"           "$(campo caso7 SRR317135 4)" "sralite"
chk "default sra sin 4o argumento" "$(LEDGER="$ESC/d2.tsv"; source <(sed -n '/^registrar_md5()/,/^}/p' "$REPO/scripts/fetch_runs.sh"); registrar_md5 maggi SRRX zzz; awk -F'	' 'NR==2{print $4}' "$ESC/d2.tsv")" "sra"

echo "== 8. el ledger se poda contra el manifiesto"
ESC="$TMP/caso8"; mkdir -p "$ESC/d/xx"
printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n' > "$ESC/man.tsv"
printf 'xx\tSRR_OK\tP1\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n' >> "$ESC/man.tsv"
printf 'org\trun\tmd5\tformato\tfecha_utc\nxx\tSRR_OK\taaa\tsra\t2026-09-18\nxx\tEXCLUIDA\tbbb\tsra\t2026-09-18\n' > "$ESC/led.tsv"
echo x > "$ESC/d/xx/SRR_OK.sra"
S8=$(SRA_DEST="$ESC/d" SRA_LEDGER="$ESC/led.tsv" MANIFEST="$ESC/man.tsv" "$REPO/scripts/fetch_runs.sh" ledger 2>&1)
tiene  "avisa que la saca"         "ya no estan en el manifiesto"  "$S8"
tiene  "la nombra"                 "EXCLUIDA"                      "$S8"
grep -q "EXCLUIDA" "$ESC/led.tsv" && mal "la saca del ledger" || ok "la saca del ledger"
grep -q "SRR_OK"   "$ESC/led.tsv" && ok "no toca la que si esta" || mal "no toca la que si esta"
chk "quedan 1 fila + cabecera"     "$(wc -l < "$ESC/led.tsv")" "2"

echo
if [[ $FALLAS -eq 0 ]]; then echo "TODO OK"; else echo "$FALLAS fallas"; exit 1; fi
