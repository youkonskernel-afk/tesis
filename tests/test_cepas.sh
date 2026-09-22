#!/usr/bin/env bash
# Paginacion de listar_cepas, --taxon y --grep, con un curl falso.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
G="$RAIZ/scripts/fetch_genomes.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
tiene(){ grep -qF -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
notiene(){ grep -qF -- "$2" <<<"$3" && mal "$1 — no debia: $2" || ok "$1"; }

SPEC="$TMP/genomas.tsv"
{ printf 'org\tespecie\tfuente\tassembly\taccession\testado\tconfianza\tnota\n'
  printf 'cloro\tClonostachys rosea\t?\t?\t\tcandidato\tnula\tn\n'
  printf 'loop\tLoop especie\t?\t?\t\tcandidato\tnula\tn\n'
  printf 'vacio\tVacio especie\t?\t?\t\tcandidato\tnula\tn\n'
} > "$SPEC"

mkdir -p "$TMP/bin"
cat > "$TMP/bin/curl" <<'CURL'
#!/usr/bin/env bash
dst=""; url=""
while [[ $# -gt 0 ]]; do
  case "$1" in -o) dst="$2"; shift 2 ;; --max-time) shift 2 ;; -*) shift ;; *) url="$1"; shift ;; esac
done
echo "$url" >> "$URLLOG"
rep(){ printf '{"accession":"%s","assembly_info":{"assembly_name":"%s","assembly_level":"Scaffold","assembly_status":"current"},"assembly_stats":{"scaffold_n50":%s,"total_sequence_length":"58000000"},"organism":{"organism_name":"X","infraspecific_names":{"strain":"%s"}}}' "$1" "$2" "$3" "$4"; }
env(){ printf '{"total_count":%s,"reports":[%s]%s}\n' "$1" "$2" "${3:+,\"next_page_token\":\"$3\"}" > "$dst"; }

case "$url" in
  # 3 paginas. El mas contiguo (N50 9 Mb) esta en la pagina 3, a proposito.
  *Clonostachys*page_token=p3*) env 9 "$(rep GCA_3.1 ASM_C3 9000000 IK726)" "" ;;
  *Clonostachys*page_token=p2*) env 9 "$(rep GCA_2.1 ASM_C2 2000000 NF-06),$(rep GCA_2.2 ASM_C2b 1000000 CR-7)" p3 ;;
  *Clonostachys*)               env 9 "$(rep GCA_1.1 ASM_C1 3000000 MLY32),$(rep GCA_1.2 ASM_C1b 1500000 CanS41)" p2 ;;
  # Bionectria: el sinonimo SI tiene IK726
  *Bionectria*) env 1 "$(rep GCA_BIO.1 ASM_IK726_v1 4000000 IK726)" "" ;;
  # devuelve siempre el mismo token -> tiene que cortar por el tope
  *Loop*) env 99 "$(rep GCA_L.1 ASM_L 1000000 L1)" siempre ;;
  *Vacio*) env 0 "" "" ;;
  *) env 0 "" "" ;;
esac
exit 0
CURL
chmod +x "$TMP/bin/curl"
export PATH="$TMP/bin:$PATH" GENOMES_SPEC="$SPEC" URLLOG="$TMP/urls.txt"
: > "$URLLOG"

echo "== 1. REGRESION: pagina las 3 paginas, no se queda con la primera"
S=$(bash "$G" cepas cloro 2>&1); echo "$S"
tiene "lista los 5 de las 3 paginas" "ensamblados disponibles (5"  "$S"
tiene "trae el de la pagina 3"       "cepa=IK726"                  "$S"
tiene "trae el de la pagina 2"       "cepa=NF-06"                  "$S"
tiene "trae el de la pagina 1"       "cepa=MLY32"                  "$S"
[[ $(grep -c page_token "$URLLOG") -eq 2 ]] && ok "pidio 3 paginas" || mal "pidio 3 paginas ($(wc -l < "$URLLOG") urls)"
tiene "page_size ya no es 20"        "page_size=100"               "$(cat "$URLLOG")"

echo "== 2. el orden es global, no por pagina"
[[ "$(grep -o 'cepa=[A-Za-z0-9-]*' <<<"$S" | head -1)" == "cepa=IK726" ]] \
  && ok "el mas contiguo (pagina 3) sale primero" || mal "orden global"

echo "== 3. un listado incompleto se anuncia"
: > "$URLLOG"; S3=$(bash "$G" cepas cloro 2>&1)
tiene "avisa que la API dice 9"      "la API dice 9 en total — listado incompleto" "$S3"

echo "== 4. --taxon consulta el sinonimo"
: > "$URLLOG"; S4=$(bash "$G" cepas cloro --taxon 'Bionectria ochroleuca' 2>&1); echo "$S4"
tiene "marca que forzaste el taxon"  "(taxon forzado)"             "$S4"
tiene "encuentra IK726"              "ASM_IK726_v1"                "$S4"
tiene "pidio Bionectria"             "Bionectria%20ochroleuca"     "$(cat "$URLLOG")"
notiene "no pidio Clonostachys"      "Clonostachys"                "$(cat "$URLLOG")"

echo "== 5. --grep filtra sin ocultar el total"
S5=$(bash "$G" cepas cloro --grep IK726 2>&1); echo "$S5"
tiene "dice cuantos hay y cuantos matchean" "5 ensamblados; 1 coinciden con 'IK726'" "$S5"
tiene "muestra el que matchea"       "cepa=IK726"                  "$S5"
notiene "no muestra los otros"       "cepa=MLY32"                  "$S5"
S5b=$(bash "$G" cepas cloro --grep NOEXISTE 2>&1)
tiene "un grep sin matches lo dice"  "NINGUNO coincide con 'NOEXISTE'" "$S5b"

echo "== 6. no se cuelga si la API pagina para siempre"
S6=$(timeout 60 bash "$G" cepas loop 2>&1); RC=$?
[[ $RC -eq 0 ]] && ok "termina (no timeout)" || mal "termina (rc=$RC)"
tiene "avisa del corte"              "corte tras 20 paginas"       "$S6"

echo "== 7. especie sin ensamblados"
S7=$(bash "$G" cepas vacio 2>&1)
tiene "lo dice"                      "no tiene ningun ensamblado"  "$S7"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
