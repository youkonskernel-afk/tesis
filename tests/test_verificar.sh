#!/usr/bin/env bash
# Pruebas de 'verificar' con archivos .gz REALES, no falsos: asi gzip -t y
# sha256sum se ejercitan de verdad.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
G="$RAIZ/scripts/fetch_genomes.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
tiene(){ grep -qF "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
fila(){ grep -E "^$1 " <<<"$2"; }

SPEC="$TMP/genomas.tsv"; LED="$TMP/genomas.sha256"; D="$TMP/genomes"
{ printf 'org\tespecie\tfuente\tassembly\taccession\testado\tconfianza\tnota\n'
  printf 'heredo\tHeredada sp\tX\tASM_H\t\theredado\talta\tn\n'
  printf 'bueno\tBuena sp\tNCBI\tASM_B\tGCF_1.1\tverificado\talta\tn\n'
  printf 'nobaj\tSinbajar sp\tNCBI\tASM_N\tGCF_2.1\tverificado\talta\tn\n'
  printf 'shamal\tShamala sp\tNCBI\tASM_S\tGCF_3.1\tverificado\talta\tn\n'
  printf 'trunc\tTruncada sp\tNCBI\tASM_T\tGCF_4.1\tverificado\talta\tn\n'
  printf 'nofa\tNofasta sp\tNCBI\tASM_F\tGCF_5.1\tverificado\talta\tn\n'
  printf 'nosid\tNosidecar sp\tNCBI\tASM_C\tGCF_6.1\tverificado\talta\tn\n'
} > "$SPEC"

fasta() { printf '>chr1 fixture\n%s\n' "$(head -c 2000 /dev/zero | tr '\0' 'ACGT')"; }

# bueno: gz valido, sha en el ledger, sidecar presente
mkdir -p "$D/bueno"; fasta | gzip -c > "$D/bueno/GCF_1.1.fna.gz"
SHA_B=$(sha256sum "$D/bueno/GCF_1.1.fna.gz" | cut -d' ' -f1)
echo "$SHA_B  GCF_1.1.fna.gz" > "$D/bueno/GCF_1.1.fna.gz.sha256"
# shamal: el ledger dice otro sha
mkdir -p "$D/shamal"; fasta | gzip -c > "$D/shamal/GCF_3.1.fna.gz"
echo "deadbeef  GCF_3.1.fna.gz" > "$D/shamal/GCF_3.1.fna.gz.sha256"
# trunc: gz cortado, y SIN entrada en el ledger
mkdir -p "$D/trunc"; fasta | gzip -c > "$D/trunc/full.gz"
head -c $(( $(stat -c%s "$D/trunc/full.gz") / 2 )) "$D/trunc/full.gz" > "$D/trunc/GCF_4.1.fna.gz"
rm "$D/trunc/full.gz"; touch "$D/trunc/GCF_4.1.fna.gz.sha256"; echo x > "$D/trunc/GCF_4.1.fna.gz.sha256"
# nofa: gz valido pero adentro no hay FASTA
mkdir -p "$D/nofa"; printf 'esto no es fasta\n' | gzip -c > "$D/nofa/GCF_5.1.fna.gz"
echo x > "$D/nofa/GCF_5.1.fna.gz.sha256"
# nosid: gz valido, sha en el ledger, pero falta el sidecar
mkdir -p "$D/nosid"; fasta | gzip -c > "$D/nosid/GCF_6.1.fna.gz"
SHA_C=$(sha256sum "$D/nosid/GCF_6.1.fna.gz" | cut -d' ' -f1)

{ printf 'org\taccession\tassembly\tsha256\tfecha_utc\n'
  printf 'bueno\tGCF_1.1\tASM_B\t%s\t2026-09-17\n' "$SHA_B"
  printf 'shamal\tGCF_3.1\tASM_S\tdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeef\t2026-09-17\n'
  printf 'nofa\tGCF_5.1\tASM_F\t%s\t2026-09-17\n' "$(sha256sum "$D/nofa/GCF_5.1.fna.gz" | cut -d' ' -f1)"
  printf 'nosid\tGCF_6.1\tASM_C\t%s\t2026-09-17\n' "$SHA_C"
} > "$LED"

export GENOMES_SPEC="$SPEC" GENOMES_LEDGER="$LED" GENOMES_DIR="$D"
S=$(bash "$G" verificar 2>&1); RC=$?
echo "$S"
echo "-------------------------------------------------------"

tiene "bueno pasa"                  "OK"                                  "$(fila bueno "$S")"
tiene "el heredado sin accession"   "SIN RESPALDO"                        "$(fila heredo "$S")"
tiene "el no bajado"                "FALTA — corre:"                      "$(fila nobaj "$S")"
tiene "sha que no coincide"         "sha256 NO coincide con el ledger"    "$(fila shamal "$S")"
tiene "gz truncado (sin ledger)"    "gzip corrupto o truncado"            "$(fila trunc "$S")"
tiene "y dice que no esta en ledger" "sin entrada en genomas.sha256"      "$(fila trunc "$S")"
tiene "gz valido que no es FASTA"   "no arranca con '>' (no es FASTA)"    "$(fila nofa "$S")"
tiene "falta el sidecar"            "falta el .sha256 al lado"            "$(fila nosid "$S")"
tiene "cuenta bien"                 "ok=1  con problemas=4  sin bajar=1  sin respaldo=1" "$S"
tiene "explica los sin respaldo"    "nunca se registro su accession"      "$S"
[[ $RC -ne 0 ]] && ok "exit != 0 con problemas" || mal "exit != 0 con problemas (rc=$RC)"

echo "== solo un organismo"
S1=$(bash "$G" verificar bueno 2>&1); RC1=$?
tiene "chequea bueno"   "OK"    "$(fila bueno "$S1")"
[[ $(grep -c "^shamal " <<<"$S1") -eq 0 ]] && ok "no toca los otros" || mal "no toca los otros"
[[ $RC1 -eq 0 ]] && ok "exit 0 cuando todo ok" || mal "exit 0 cuando todo ok (rc=$RC1)"

echo "== --rapido no lee el archivo, asi que no ve el sha malo"
S2=$(bash "$G" verificar shamal --rapido 2>&1)
tiene "avisa que no verifica checksums" "modo rapido"                     "$S2"
[[ $(grep -c "sha256 NO coincide" <<<"$S2") -eq 0 ]] && ok "no compara sha" || mal "no compara sha"
tiene "pero igual detecta lo barato"    "falta el .sha256 al lado"        "$(fila nosid "$(bash "$G" verificar nosid --rapido 2>&1)")"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
