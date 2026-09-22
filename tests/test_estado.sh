#!/usr/bin/env bash
# 'estado' contra un DEST falso: cuenta lo bajado por organismo y lo que falta.
#
# Por que este banco existe: 'estado' es el modo que contesta "¿ya esta todo?",
# y ya tuvo un bug silencioso — el awk del resumen usaba -F'\t' cuando 'filas'
# emite \x1f, asi que la linea de las que faltaban no se imprimia nunca. Un modo
# de reporte que reporta de menos no falla: solo dice que esta todo bien.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R="$RAIZ/scripts/fetch_runs.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
tiene(){ grep -qF -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
fila(){ grep -E "^$1 " <<<"$2" | tr -s ' '; }

MAN="$TMP/man.tsv"
{ printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n'
  printf 'aa\tSRR_A1\tPRJ_A\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRR_A2\tPRJ_A\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRR_A3\tPRJ_A\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'bb\tSRR_B1\tPRJ_B\tduplicado\tentr\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'cc\tSRR_C1\tPRJ_C\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
} > "$MAN"

DEST="$TMP/dest"; mkdir -p "$DEST"
# Las tres disposiciones que ruta_sra acepta, una cada una: asi el banco falla
# si alguien saca alguna del bucle.
mkdir -p "$DEST/aa";     head -c 2048 /dev/zero > "$DEST/aa/SRR_A1.sra"
mkdir -p "$DEST/SRR_A2"; head -c 2048 /dev/zero > "$DEST/SRR_A2/SRR_A2.sra"
head -c 2048 /dev/zero > "$DEST/SRR_A3.sra"
# bb: existe pero VACIO. ruta_sra pide -s, o un fichero truncado por una sesion
# de Colab que se murio a mitad contaria como bajado.
mkdir -p "$DEST/bb";     : > "$DEST/bb/SRR_B1.sra"
# cc: no bajada.

S=$(MANIFEST="$MAN" SRA_DEST="$DEST" bash "$R" estado 2>&1)
echo "$S" | sed -n '/^ORG/,$p'

echo "== 1. cuenta por organismo"
[[ "$(fila aa "$S")" == "aa 3 0" ]] && ok "aa: 3 bajadas, 0 faltan (las 3 disposiciones)" \
  || mal "aa: 3 bajadas (dijo '$(fila aa "$S")')"
[[ "$(fila bb "$S")" == "bb 0 1" ]] && ok "bb: el .sra vacio NO cuenta como bajado" \
  || mal "bb: el .sra vacio NO cuenta (dijo '$(fila bb "$S")')"
[[ "$(fila cc "$S")" == "cc 0 1" ]] && ok "cc: 0 bajadas, 1 falta" \
  || mal "cc: 0 bajadas, 1 falta (dijo '$(fila cc "$S")')"

echo "== 2. el total no se contradice con la tabla"
tiene "total correcto" "total: 3 bajadas, 2 faltan" "$S"

echo "== 3. dice contra que manifiesto y que destino corrio"
tiene "nombra el manifiesto" "$MAN"  "$S"
tiene "nombra el destino"    "$DEST" "$S"

echo "== 4. un .sralite en DEST no cuenta, y es a proposito"
# Al moverlo a Drive el fichero se guarda como <RUN>.sra venga normalizado o
# lite (la columna 'formato' del ledger es lo que los distingue), asi que un
# .sralite en DEST seria un fichero que el pipeline no espera. Si algun dia
# ruta_sra empieza a aceptarlo, este chequeo lo avisa en vez de que el cambio
# pase inadvertido.
D2="$TMP/dest2"; mkdir -p "$D2/cc"; head -c 2048 /dev/zero > "$D2/cc/SRR_C1.sralite"
S4=$(MANIFEST="$MAN" SRA_DEST="$D2" bash "$R" estado 2>&1)
[[ "$(fila cc "$S4")" == "cc 0 1" ]] && ok "sigue contando como que falta" \
  || mal "sigue contando como que falta (dijo '$(fila cc "$S4")')"

echo "== 5. sin manifiesto no inventa un estado vacio"
S5=$(MANIFEST="$TMP/noexiste.tsv" SRA_DEST="$DEST" bash "$R" estado 2>&1); RC5=$?
[[ $RC5 -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC5)"
tiene "y dice como arreglarlo" "manifest" "$S5"

echo "== 6. manifiesto entero bajado"
D3="$TMP/dest3"; mkdir -p "$D3"
for r in SRR_A1 SRR_A2 SRR_A3 SRR_B1 SRR_C1; do head -c 2048 /dev/zero > "$D3/$r.sra"; done
S6=$(MANIFEST="$MAN" SRA_DEST="$D3" bash "$R" estado 2>&1)
tiene "total sin faltantes" "total: 5 bajadas, 0 faltan" "$S6"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
