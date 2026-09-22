#!/usr/bin/env bash
# Banco de pruebas de fetch_genomes.sh contra una spec sintetica y un curl falso
# que imita la API de NCBI Datasets v2alpha.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
G="$RAIZ/scripts/fetch_genomes.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok()  { printf '  ok   %s\n' "$1"; }
mal() { printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
tiene()   { if grep -qF -- "$2" <<<"$3"; then ok "$1"; else mal "$1 — falta: $2"; fi; }
notiene() { if grep -qF -- "$2" <<<"$3"; then mal "$1 — no debia: $2"; else ok "$1"; fi; }

# ---- spec sintetica: un organismo por desenlace ----
SPEC="$TMP/genomas.tsv"
{
  printf '# comentario que tiene que sobrevivir\n'
  printf 'org\tespecie\tfuente\tassembly\taccession\testado\tconfianza\tnota\n'
  printf 'here\tHeredada especie\tX\tASM1\t\theredado\talta\tn\n'
  printf 'coinc\tCoincide especie\tNCBI\tASM_OK\tGCF_111.1\tcandidato\talta\tn\n'
  printf 'difie\tDifiere especie\tNCBI\tASM_VIEJO\tGCF_222.1\tcandidato\talta\tn\n'
  printf 'supre\tSuprimida especie\tNCBI\tASM_SUP\tGCF_333.1\tcandidato\tmedia\tn\n'
  printf 'cepa\tCepa especie\t?\t?\t\tcandidato\tnula\tn\n'
  printf 'nored\tNored especie\tNCBI\tASM_N\tGCF_444.1\tcandidato\talta\tn\n'
  printf 'listo\tLista especie\tNCBI\tASM_L\tGCF_555.1\tverificado\talta\tn\n'
} > "$SPEC"

mkdir -p "$TMP/bin"
cat > "$TMP/bin/curl" <<'CURL'
#!/usr/bin/env bash
dst=""; url=""
while [[ $# -gt 0 ]]; do
  case "$1" in -o) dst="$2"; shift 2 ;; --max-time) shift 2 ;; -*) shift ;; *) url="$1"; shift ;; esac
done
rep() {
  local cepa=""
  [[ -n "${6:-}" ]] && cepa=',"infraspecific_names":{"strain":"'"$6"'"}'
  printf '{"accession":"%s","assembly_info":{"assembly_name":"%s","assembly_level":"%s","assembly_status":"%s"},"assembly_stats":{"scaffold_n50":%s,"total_sequence_length":"%s"},"organism":{"organism_name":"X"%s}}' \
    "$1" "$2" "$3" "${7:-current}" "$4" "$5" "$cepa"
}
env() { printf '{"reports":[%s]}\n' "$1" > "$dst"; }
case "$url" in
  *Nored*|*GCF_444.1*)  echo "curl: (28) Operation timed out" >&2; exit 28 ;;
  *accession/GCF_111.1/*|*taxon/Coincide*) env "$(rep GCF_111.1 ASM_OK Chromosome 27400000 227400000)" ;;
  *accession/GCF_222.1/*) env "$(rep GCF_222.1 ASM_VIEJO Scaffold 1000000 700000000)" ;;
  *taxon/Difiere*)        env "$(rep GCF_999.9 ASM_NUEVO Chromosome 37700000 650800000)" ;;
  *accession/GCF_333.1/*) env "$(rep GCF_333.1 ASM_SUP Chromosome 37600000 703000000 "" suppressed)" ;;
  *taxon/Suprimida*)      env "$(rep GCF_888.8 ASM_VIGENTE Chromosome 37700000 650800000)" ;;
  # el caso cloro: la especie SI tiene referencia, pero es de otra cepa
  *taxon/Cepa*reference_only*) env "$(rep GCA_000111.1 ASM_NF06 "Complete Genome" 7900000 56900000 NF-06)" ;;
  *taxon/Cepa*)                env "$(rep GCA_000222.1 ASM_IK726 Scaffold 1200000 58300000 IK726),$(rep GCA_000111.1 ASM_NF06 "Complete Genome" 7900000 56900000 NF-06)" ;;
  # descarga
  *download*) printf 'ZIPFALSO' > "$dst" ;;
  *) env "" ;;
esac
exit 0
CURL
cat > "$TMP/bin/unzip" <<'UNZIP'
#!/usr/bin/env bash
printf '>chr1\nACGTACGTAC\n'
UNZIP
chmod +x "$TMP/bin/curl" "$TMP/bin/unzip"
export PATH="$TMP/bin:$PATH" GENOMES_SPEC="$SPEC" GENOMES_LEDGER="$TMP/led.tsv"

SAL=$(bash "$G" resolve 2>&1)
pro() { sed -n "/^== $1 /,/^\$/p" <<<"$SAL"; }

echo "== veredictos (no regresion)"
tiene   "coinc COINCIDE"        ">>> COINCIDE"                   "$(pro coinc)"
tiene   "difie DIFIERE"         ">>> DIFIERE — la referencia vigente es GCF_999.9" "$(pro difie)"
tiene   "nored error de red"    "ERROR DE RED"                   "$(pro nored)"
tiene   "nored sin veredicto"   ">>> SIN RESPUESTA UTIL"         "$(pro nored)"
tiene   "here se salta"         "here: heredado, se salta"       "$SAL"

echo "== NUEVO: un ensamblado retirado se ve"
tiene   "supre marca RETIRADO"  "!!! RETIRADO por NCBI (suppressed) — no usar" "$(pro supre)"
tiene   "supre igual da DIFIERE" ">>> DIFIERE"                   "$(pro supre)"
notiene "coinc sin RETIRADO"    "RETIRADO"                       "$(pro coinc)"

echo "== NUEVO: sin candidato pero CON referencia, igual lista las cepas"
tiene   "avisa del riesgo de cepa" "la referencia puede ser de otra cepa" "$(pro cepa)"
tiene   "lista los ensamblados"    "ensamblados disponibles (2"          "$(pro cepa)"
tiene   "muestra IK726"            "cepa=IK726"                          "$(pro cepa)"
tiene   "muestra NF-06"            "cepa=NF-06"                          "$(pro cepa)"
tiene   "veredicto mira la cepa"   ">>> SIN CANDIDATO — elegi de la lista de arriba, mirando la cepa" "$(pro cepa)"

echo "== NUEVO: subcomando cepas"
C=$(bash "$G" cepas cepa 2>&1)
tiene   "cepas suelto anda"        "cepa=IK726"                          "$C"
C2=$(bash "$G" cepas noexiste 2>&1); RC=$?
tiene   "rechaza org desconocido"  "organismo desconocido: noexiste"     "$C2"
[[ $RC -ne 0 ]] && ok "sale con error" || mal "sale con error"
C3=$(bash "$G" cepas 2>&1 || true)
tiene   "exige un ORG"             "uso:"                                "$C3"

echo "== el gate de fetch sigue firme"
F=$(GENOMES_DIR="$TMP/g" bash "$G" fetch 2>&1)
tiene   "baja el verificado"       "listo: bajando GCF_555.1"            "$F"
tiene   "saltea el candidato"      "SALTO cepa: 'candidato'"             "$F"
tiene   "saltea el heredado"       "SALTO here: heredado"                "$F"
tiene   "cuenta bien"              "bajados=1 saltados=6"                "$F"
[[ -s "$TMP/g/listo/GCF_555.1.fna.gz" ]] && ok "el FASTA quedo escrito" || mal "el FASTA quedo escrito"
tiene   "sha256 al ledger"         "GCF_555.1"                           "$(cat "$TMP/led.tsv" 2>/dev/null)"

echo
if [[ $FALLAS -eq 0 ]]; then echo "TODO OK"; else echo "$FALLAS fallas"; exit 1; fi
