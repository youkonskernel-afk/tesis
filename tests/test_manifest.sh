#!/usr/bin/env bash
# Regresion del filtro del manifiesto tras sacarlo a pasa_filtro().
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R="$RAIZ/scripts/fetch_runs.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }

SPEC="$TMP/organismos.tsv"
{ printf 'org\tespecie\treino\tclado\trol\tbioproject\truns\tspots_M\tstrategy\tassembly\tnota\tset_modelo\n'
  printf 'test\tEspecie test\tFungi\tAsco\tprimario\tPRJNA000001\t8\t100\tmiRNA-Seq\tA\tn\taplicacion\n'
} > "$SPEC"

mkdir -p "$TMP/bin"
cat > "$TMP/bin/curl" <<'CURL'
#!/usr/bin/env bash
printf 'run_accession\tread_count\tbase_count\tlibrary_strategy\tlibrary_layout\tlibrary_source\n'
printf 'SRR1\t9000000\t200000000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'   # entra
printf 'SRR2\t8000000\t180000000\tncRNA-Seq\tPAIRED\tTRANSCRIPTOMIC\n'   # entra (ncRNA-Seq acepta PAIRED)
printf 'SRR3\t7000000\t160000000\tRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'     # entra
printf 'SRR4\t6000000\t140000000\tmiRNA-Seq\tPAIRED\tTRANSCRIPTOMIC\n'   # entra
printf 'SRR5\t5000000\t120000000\tmiRNA-Seq\tSINGLE\tGENOMIC\n'          # FUERA: GENOMIC
printf 'SRR6\t60000000\t9000000000\tRNA-Seq\tPAIRED\tTRANSCRIPTOMIC\n'   # FUERA: RNA-Seq PAIRED
printf 'SRR7\t1000\t20000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'          # FUERA: pocas reads
printf 'SRR8\t9000000\t200000000\tWGS\tSINGLE\tTRANSCRIPTOMIC\n'         # FUERA: estrategia
exit 0
CURL
chmod +x "$TMP/bin/curl"

# EXCLUIDAS explicito: sin esto el test lee data/excluidas.tsv del repo, que es
# una dependencia oculta que lo hace fallar cuando alguien agregue una fila.
: > "$TMP/sin_excluidas.tsv"
PATH="$TMP/bin:$PATH" ORGANISMOS="$SPEC" EXCLUIDAS="$TMP/sin_excluidas.tsv" \
  MANIFEST="$TMP/man.tsv" "$R" manifest >"$TMP/out" 2>&1
S=$(cat "$TMP/out"); echo "$S" | tail -6

RET=$(( $(wc -l < "$TMP/man.tsv") - 1 ))
[[ $RET -eq 4 ]] && ok "retiene 4 de 8" || mal "retiene 4 de 8 (retuvo $RET)"
for r in SRR1 SRR2 SRR3 SRR4; do
  grep -q "	$r	" "$TMP/man.tsv" && ok "$r entra" || mal "$r entra"
done
for r in SRR5 SRR6 SRR7 SRR8; do
  grep -q "	$r	" "$TMP/man.tsv" && mal "$r NO debe entrar" || ok "$r descartado"
done
grep -q "corridas vistas=8  descartadas=4  excluidas=0  retenidas=4" <<<"$S" \
  && ok "el resumen cuadra" || mal "el resumen cuadra"

echo "== una corrida excluida no vuelve al regenerar"
{ printf 'run\tmotivo\tfecha_utc\n'; printf 'SRR1\tprueba\t2026-09-18\n'; } > "$TMP/excl.tsv"
PATH="$TMP/bin:$PATH" ORGANISMOS="$SPEC" EXCLUIDAS="$TMP/excl.tsv" \
  MANIFEST="$TMP/man2.tsv" "$R" manifest >"$TMP/out2" 2>&1
S2=$(cat "$TMP/out2")
grep -q "	SRR1	" "$TMP/man2.tsv" && mal "SRR1 excluida no entra" || ok "SRR1 excluida no entra"
grep -q "EXCLUIDA SRR1" <<<"$S2" && ok "lo dice en la salida" || mal "lo dice en la salida"
grep -q "corridas vistas=8  descartadas=4  excluidas=1  retenidas=3" <<<"$S2" \
  && ok "las categorias suman 8" || mal "las categorias suman 8"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
