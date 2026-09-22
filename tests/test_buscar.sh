#!/usr/bin/env bash
# 'buscar' con un curl falso, y la prueba de que usa el mismo filtro que
# 'manifest' — que es lo que evita proponer un proyecto que despues no entra.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R="$RAIZ/scripts/fetch_runs.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
tiene(){ grep -qE "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
notiene(){ grep -qE "$2" <<<"$3" && mal "$1 — no debia: $2" || ok "$1"; }

SPEC="$TMP/organismos.tsv"
{ printf 'org\tespecie\treino\tclado\trol\tbioproject\truns\tspots_M\tstrategy\tassembly\tnota\tset_modelo\n'
  printf 'sclsc\tSclerotinia sclerotiorum\tFungi\tAsco\tprimario\tPRJNA477286\t18\t184\tmiRNA-Seq\tA\tn\taplicacion\n'
} > "$SPEC"

mkdir -p "$TMP/bin"
cat > "$TMP/bin/curl" <<'CURL'
#!/usr/bin/env bash
printf 'run_accession\tstudy_accession\tread_count\tbase_count\tlibrary_strategy\tlibrary_layout\tlibrary_source\n'
# PRJNA477286: el primario que YA esta en la spec (18 -> aca 2, alcanza)
printf 'SRR1\tPRJNA477286\t6726018\t150000000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
printf 'SRR2\tPRJNA477286\t6772389\t150000000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
# PRJNA985401: el duplicado que no sirve — RNA-Seq PAIRED, tiene que caer entero
printf 'SRR3\tPRJNA985401\t60000000\t9000000000\tRNA-Seq\tPAIRED\tTRANSCRIPTOMIC\n'
printf 'SRR4\tPRJNA985401\t61000000\t9000000000\tRNA-Seq\tPAIRED\tTRANSCRIPTOMIC\n'
# PRJNA999001: candidato bueno, miRNA-Seq SINGLE, NO esta en la spec
printf 'SRR5\tPRJNA999001\t9000000\t200000000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
printf 'SRR6\tPRJNA999001\t8000000\t180000000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
printf 'SRR7\tPRJNA999001\t7000000\t160000000\tncRNA-Seq\tPAIRED\tTRANSCRIPTOMIC\n'
# PRJNA999002: RNA-Seq pero SINGLE -> entra (el filtro solo exige SINGLE en RNA-Seq)
printf 'SRR8\tPRJNA999002\t5000000\t120000000\tRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
# PRJNA999003: todas por debajo del minimo de reads -> cae entero
printf 'SRR9\tPRJNA999003\t1000\t20000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
# PRJNA999004: GENOMIC -> el filtro duro lo saca
printf 'SRR10\tPRJNA999004\t50000000\t7000000000\tmiRNA-Seq\tSINGLE\tGENOMIC\n'
exit 0
CURL
chmod +x "$TMP/bin/curl"

S=$(PATH="$TMP/bin:$PATH" "$R" buscar 'Sclerotinia sclerotiorum' 2>&1)
echo "$S"
echo "-------------------------------------------------------"

echo "== el filtro es el mismo que el del manifiesto"
tiene   "cuenta crudas y filtradas"  "10 corridas TRANSCRIPTOMIC .* 6 pasan"  "$S"
notiene "descarta RNA-Seq PAIRED"    "PRJNA985401"                            "$S"
notiene "descarta GENOMIC"           "PRJNA999004"                            "$S"
notiene "descarta pocas reads"       "PRJNA999003"                            "$S"
tiene   "acepta RNA-Seq SINGLE"      "PRJNA999002"                            "$S"
tiene   "acepta ncRNA-Seq PAIRED"    "PRJNA999001 +3"                         "$S"

echo "== marca lo que ya esta en la spec"
tiene   "el primario aparece marcado" "PRJNA477286.*sclsc primario"           "$S"
tiene   "el candidato nuevo con -"    "PRJNA999001.*-$"                       "$S"

echo "== ordenado por cantidad de corridas"
[[ "$(grep -oE 'PRJNA[0-9]+' <<<"$S" | head -1)" == "PRJNA999001" ]] \
  && ok "el de mas corridas primero" || mal "orden por corridas"

echo "== el criterio vive en un solo lugar (chequeo estructural)"
if grep -q 'grep -qw -- "\$strat"' "$R"; then
  mal "cmd_manifest ya no tiene el filtro inline"
else ok "cmd_manifest ya no tiene el filtro inline"; fi
# La invariante no es "exactamente 2 llamadas" sino "ningun sitio filtra por su
# cuenta". Fijarlo en 2 hacia fallar el test al agregar un tercer consumidor
# legitimo, que es justo lo contrario de lo que se quiere proteger.
N=$(grep -c 'pasa_filtro "' "$R")
[[ $N -ge 2 ]] && ok "todos los consumidores llaman a pasa_filtro ($N)" \
               || mal "pasa_filtro se usa en $N lugares, esperaba 2 o mas"

echo "== errores"
tiene   "exige la especie"           "uso:"  "$("$R" buscar 2>&1 || true)"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
