#!/usr/bin/env bash
# align.sh con yasma/bowtie falsos en el PATH.
#
# Lo que este banco protege, por orden de gravedad:
#  1. Que no se alinee contra un genoma cuyo sha256 NO coincide con
#     data/genomas.sha256. Ese ledger es el unico registro de contra que se
#     alineo; si se alinea contra otra cosa, no hay forma de saberlo despues.
#  2. Que `verificar` atrape una fraccion alineada muy baja. Alinear contra el
#     genoma equivocado NO falla: bowtie alinea mal y sale con 0. Medido contra
#     el yasma real — con un genoma de otro azar, 0.0% alineado y exit 0.
#  3. Que un @RG que falta en el BAM se reporte: una libreria que se perdio en
#     el camino no la dice nadie, y `yasma tradeoff` agrega por read group.
#  4. Que el genoma que se le pasa a yasma este ADENTRO del -o. `ic.check()`
#     hace relative_to(output_directory) sin protegerlo: un genoma compartido
#     fuera del proyecto tira ValueError, no un mensaje.
#  5. Que ese "adentro" sea un symlink AL DIRECTORIO del organismo, para que el
#     indice .ebwt sea uno por organismo y no uno por proyecto: con un genoma de
#     1 Gb, bowtie-build de mas son horas.
#  6. Que no se alinee un proyecto a medio recortar: sale un BAM al que le
#     faltan corridas y nada lo dice.
#  7. Que los parametros de bowtie vayan explicitos. -m 50 esta medido en danre.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
A="$RAIZ/scripts/align.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
# El `--` es imprescindible: sin el, un patron que empieza con '-' (como
# "-m 50") lo toma grep como opcion, aborta, y en notiene() el exit nonzero cae
# en la rama de exito -> el chequeo pasa SIEMPRE.
tiene(){ grep -qF -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
notiene(){ grep -qF -- "$2" <<<"$3" && mal "$1 — no deberia estar: $2" || ok "$1"; }

MAN="$TMP/man.tsv"
{ printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n'
  printf 'aa\tSRR_P1\tPRJ_A\tprimario\tapl\t1000\t50000\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRR_P2\tPRJ_A\tprimario\tapl\t1000\t50000\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRR_D1\tPRJ_B\tduplicado\tapl\t1000\t50000\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
} > "$MAN"

GEN="$TMP/genomes"; mkdir -p "$GEN/aa"
printf '>chr1\nACGTACGTACGTACGTACGT\n' | gzip -c > "$GEN/aa/GCF_TEST.1.fna.gz"
SHA=$(sha256sum "$GEN/aa/GCF_TEST.1.fna.gz" | cut -d' ' -f1)
LED="$TMP/genomas.sha256"
ledger_ok(){ printf 'org\taccession\tassembly\tsha256\tfecha_utc\naa\tGCF_TEST.1\tTestAsm\t%s\t2026-09-23\n' "$SHA" > "$LED"; }
ledger_ok

# Proyectos ya recortados, como los deja trim.sh.
PROY="$TMP/proy"
sembrar_trim() {
  rm -rf "$PROY"
  local p r
  for p in aa_primario:SRR_P1,SRR_P2 aa_duplicado:SRR_D1; do
    local dir="$PROY/${p%%:*}"
    mkdir -p "$dir/trim"
    printf 'run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc\n' > "$dir/recortadas.tsv"
    for r in ${p##*:}; do :; done
    local IFS=,
    for r in ${p##*:}; do
      printf 'x' | gzip -c > "$dir/trim/$r.t.fq.gz"
      printf '%s\tPRJ\ttrim/%s.t.fq.gz\t1000\t1000\t100.0\t2026-09-23T00:00:00Z\n' "$r" "$r" >> "$dir/recortadas.tsv"
    done
    printf '{"project_name":"%s","trimmed_libraries":[]}\n' "${p%%:*}" > "$dir/inputs.json"
  done
}
sembrar_trim

mkdir -p "$TMP/bin"
cat > "$TMP/bin/bowtie" <<'BW'
#!/usr/bin/env bash
echo "bowtie version 1.3.1"
BW
cp "$TMP/bin/bowtie" "$TMP/bin/bowtie-build"
cat > "$TMP/bin/yasma" <<'YA'
#!/usr/bin/env bash
# Imita a `yasma align` de nativealign.py v1.1.1 en lo que importa:
#  - exige el genoma ADENTRO del output_directory (ic.check hace relative_to)
#  - escribe align/alignment.bam, su .bai y align/library_stats.txt con los
#    conteos por read group: umap mmap_wg mmap_nw xmap_nw xmap_ma xmap_nv xmap_fr
#  - el read group sale de get_rg(): pela .gz/.t/.fq -> SRR.t.fq.gz da SRR
echo "CWD=$PWD" >> "$LOG_YA"
echo "ARGS=$*" >> "$LOG_YA"
OUT=""; GEN=""
while [[ $# -gt 0 ]]; do
  case "$1" in -o) OUT="$2"; shift 2;; -g) GEN="$2"; shift 2;; *) shift;; esac
done
case "$GEN" in
  "$OUT"/*) echo "GENOMA_ADENTRO=$GEN" >> "$LOG_YA" ;;
  *) echo "ValueError: '$GEN' is not in the subpath of '$OUT'" >&2; exit 1 ;;
esac
[[ -e "$GEN" ]] || { echo "InputError: genome not found: $GEN" >&2; exit 1; }
mkdir -p "$OUT/align"
printf 'BAMFALSO' > "$OUT/align/alignment.bam"
printf 'BAIFALSO' > "$OUT/align/alignment.bam.bai"
{ printf 'project\tlibrary\tumap\tmmap_wg\tmmap_nw\txmap_nw\txmap_ma\txmap_nv\txmap_fr\n'
  for f in "$OUT"/trim/*.t.fq.gz; do
    [[ -e "$f" ]] || continue
    rg=$(basename "$f"); rg=${rg%%.*}
    printf '%s\t%s\t%s\n' "$(basename "$OUT")" "$rg" "${FAKE_COUNTS:-1000	0	0	0	0	0	0}"
  done
} > "$OUT/align/library_stats.txt"
YA
chmod +x "$TMP/bin"/*
export PATH="$TMP/bin:$PATH" LOG_YA="$TMP/ya.log"
: > "$LOG_YA"

# Un paquete yasma de mentira, para que yasma_parche.py tenga que mirar. El
# binario falso del PATH no alcanza: el parche es sobre el FUENTE instalado, y
# `correr` se niega a alinear sin el. Ver el escenario 0.
PYLIB="$TMP/pylib"; mkdir -p "$PYLIB/yasma"
: > "$PYLIB/yasma/__init__.py"
nativealign_sin_parche() {
  { printf '\t\tprint(f"stage: {mmap}", file=errf)\n\n'
    printf '\t\tif ".gz" in lib.suffixes:\n\t\t\tp = 1\n\n'
    printf "\t\tif mmap != 'over':\n\t\t\tp.wait()\n"
  } > "$PYLIB/yasma/nativealign.py"
}
nativealign_sin_parche
export PYTHONPATH="$PYLIB"
"$RAIZ/scripts/yasma_parche.py" >/dev/null 2>&1

corre() {
  MANIFEST="$MAN" GENOMAS_LEDGER="$LED" GENOMES_DIR="$GEN" PROY_DIR="$PROY" \
  BAM_DIR="$TMP/bams" CORES=2 bash "$A" "$@" 2>&1
}

echo "== 1. genoma: el sha256 manda"
S=$(corre genoma); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0 con el ledger al día" || mal "exit 0 (rc=$RC)"
tiene "dice que coincide"       "sha256 ok"          "$S"
tiene "y lo descomprime"        "descomprimiendo"    "$S"
[[ -s "$GEN/aa/GCF_TEST.1.fna" ]] && ok "deja el .fna (pysam no lee gzip plano)" \
  || mal "deja el .fna"

echo "== 1b. un sha256 que no coincide no se alinea"
printf 'org\taccession\tassembly\tsha256\tfecha_utc\naa\tGCF_TEST.1\tTestAsm\tdeadbeef\t2026-09-23\n' > "$LED"
S=$(corre genoma); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo dice"                 "no coincide con el ledger"  "$S"
tiene "y muestra los dos"       "deadbeef"                   "$S"
S=$(corre correr); RC=$?
[[ $RC -ne 0 ]] && ok "y correr tampoco arranca" || mal "y correr tampoco arranca (rc=$RC)"
notiene "no llamó a yasma"      "ARGS="                      "$(cat "$LOG_YA")"
ledger_ok

echo "== 2. plan: un proyecto por organismo y rol, con su accession"
S=$(corre plan); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
tiene "el primario"             "aa_primario"     "$S"
tiene "y el duplicado, aparte"  "aa_duplicado"    "$S"
tiene "nombra el genoma"        "GCF_TEST.1"      "$S"
tiene "los parámetros de bowtie" "-m 50"          "$S"
tiene "y dónde va el BAM"       "drive_push.sh bam" "$S"

echo "== 3. correr: cómo se llama a yasma"
rm -f "$GEN/aa/GCF_TEST.1.fna"   # para que el aviso de descompresión se dispare
S=$(corre correr); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
LY=$(cat "$LOG_YA")
tiene "parado en el dir del proyecto" "CWD=$PROY/aa_primario"       "$LY"
tiene "con -o absoluto"               "-o $PROY/aa_primario"        "$LY"
tiene "-m 50 explícito"               "--max_multi 50"              "$LY"
tiene "max_random explícito"          "--max_random 3"              "$LY"
tiene "unique_locality explícito"     "--unique_locality 50"        "$LY"
tiene "offrate explícito"             "--offrate 3"                 "$LY"
tiene "la ventana 15-50"              "--min_length 15 --max_length 50" "$LY"

echo "== 3b. el genoma va ADENTRO del -o, por symlink al directorio"
# ic.check() hace relative_to(output_directory) sin protegerlo: un genoma de
# afuera tira ValueError. Y tiene que ser un symlink AL DIRECTORIO, o
# bowtie-build correria una vez por proyecto en vez de una por organismo.
tiene "yasma lo vio adentro" "GENOMA_ADENTRO=$PROY/aa_primario/genome/" "$LY"
[[ -L "$PROY/aa_primario/genome" ]] && ok "es un symlink" || mal "es un symlink"
[[ "$(readlink "$PROY/aa_primario/genome")" == "$GEN/aa" ]] \
  && ok "y apunta al DIRECTORIO del organismo (índice compartido)" \
  || mal "y apunta al DIRECTORIO del organismo (apunta a '$(readlink "$PROY/aa_primario/genome")')"

echo "== 3c. la ruta del genoma llega limpia, sin los avisos pegados"
# preparar_genoma imprime "sha256 ok" y "descomprimiendo"; el llamador la captura
# con fna=$(...). Si esos avisos salieran por stdout se le pegarian a la ruta y
# el -g llegaria con el texto adentro. Se afirma sobre el valor que el programa
# USA —lo que yasma recibio— y no sobre lo que se imprime al lado.
GVIO=$(grep -F -- 'GENOMA_ADENTRO=' "$LOG_YA" | head -1)
[[ "$GVIO" == "GENOMA_ADENTRO=$PROY/aa_primario/genome/GCF_TEST.1.fna" ]] \
  && ok "el -g es exactamente la ruta" \
  || mal "el -g es exactamente la ruta (recibió '${GVIO#GENOMA_ADENTRO=}')"

echo "== 4. el BAM queda enlazado donde drive_push.sh lo busca"
[[ -s "$TMP/bams/aa/primario.bam" ]] && ok "bams/aa/primario.bam" || mal "bams/aa/primario.bam"
[[ -s "$TMP/bams/aa/duplicado.bam" ]] && ok "bams/aa/duplicado.bam" || mal "bams/aa/duplicado.bam"
[[ -s "$TMP/bams/aa/primario.bam.bai" ]] && ok "y su .bai" || mal "y su .bai"
i1=$(stat -c %i "$PROY/aa_primario/align/alignment.bam")
i2=$(stat -c %i "$TMP/bams/aa/primario.bam")
[[ "$i1" == "$i2" ]] && ok "hard link: no duplica el disco" || mal "hard link: no duplica el disco"

echo "== 5. queda registrado contra QUÉ se alineó"
# inputs.json guarda la RUTA del genoma, que dentro de dos años no dice nada.
L=$(cat "$PROY/aa_primario/alineado.tsv")
tiene "el accession"  "GCF_TEST.1"  "$L"
tiene "y el sha256"   "$SHA"        "$L"
tiene "y el -m usado" "$(printf '\t50\t3\t50\t3\t')" "$L"

echo "== 6. idempotencia, y un BAM más viejo que el recorte no cuenta"
: > "$LOG_YA"
S=$(corre correr)
tiene "no rehace"            "ya alineado"   "$S"
notiene "ni llama a yasma"   "ARGS="         "$(cat "$LOG_YA")"
sleep 1; touch "$PROY/aa_primario/trim/SRR_P1.t.fq.gz"
S=$(corre estado)
tiene "un recorte nuevo lo marca" "DESACTUALIZADO" "$S"

echo "== 7. no se alinea un proyecto a medio recortar"
sembrar_trim
rm -f "$PROY/aa_primario/trim/SRR_P2.t.fq.gz"
sed -i '/SRR_P2/d' "$PROY/aa_primario/recortadas.tsv"
: > "$LOG_YA"
S=$(corre correr); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "dice cuántas faltan"  "recorte incompleto: 1 de 2"  "$S"
tiene "y manda a trim.sh"    "trim.sh correr aa/primario"  "$S"
notiene "sin llamar a yasma" "ARGS="                       "$(cat "$LOG_YA")"

echo "== 8. verificar: todo bien"
sembrar_trim; corre correr >/dev/null 2>&1
S=$(corre verificar); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
tiene "reporta la fracción alineada" "100.0%"  "$S"
tiene "las 3 librerías"              "SRR_D1"  "$S"

echo "== 9. verificar atrapa el genoma equivocado (yasma sale con 0)"
# Medido contra el yasma real: con un genoma de otro azar da 0.0% alineado y
# exit 0. Nada mas lo delata.
sembrar_trim
FAKE_COUNTS=$(printf '0\t0\t0\t0\t0\t1000\t0') corre correr >/dev/null 2>&1
S=$(corre verificar); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo nombra"        "MUY BAJA"                  "$S"
tiene "y sugiere qué es" "el genoma correcto"        "$S"
tiene "y a dónde ir"     "fetch_genomes.sh verificar" "$S"

echo "== 9b. mucho por encima de -m 50 avisa, pero no es una falla"
# El hallazgo de danre: los tRF multimapean y -m 50 descarta el 76%. Es
# esperable, no un error — pero hay que verlo antes de tocar -m.
sembrar_trim
FAKE_COUNTS=$(printf '300\t0\t0\t0\t700\t0\t0') corre correr >/dev/null 2>&1
S=$(corre verificar); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0: no es falla" || mal "exit 0: no es falla (rc=$RC)"
tiene "pero lo dice"     "se pasa de -m 50"   "$S"

echo "== 9c. mucho sin alineamiento posible avisa, pero no es una falla"
# Medido en sclsc_duplicado contra el yasma real: 17.2% colocado y 67.1% sin
# ningun alineamiento valido. Pasaba como "ok" porque el unico umbral miraba la
# fraccion colocada, y 17.2 > 10. Son diagnosticos distintos: lo que no alinea
# en ninguna parte no es un genoma repetitivo, son reads que no son de este
# genoma. No es falla —el BAM esta bien escrito— pero hay que mirarlo.
sembrar_trim
FAKE_COUNTS=$(printf '41\t129\t2\t156\t1\t671\t0') corre correr >/dev/null 2>&1
S=$(corre verificar); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0: no es falla" || mal "exit 0: no es falla (rc=$RC)"
tiene "pero lo dice"        "no alinea en ninguna parte"  "$S"
tiene "con la fracción"     "67%"                         "$S"
tiene "y la columna SIN_AL" "SIN_AL"                      "$S"
# Y no se confunde con el caso de los tRF, que es el otro umbral.
notiene "no es el de -m"    "se pasa de -m"               "$S"

echo "== 9d. colocado bajo pero todo alinea en alguna parte: es el caso de -m"
# 30% colocado, 0% sin alineamiento, 60% por encima de -m. Que ALIN sea bajo no
# alcanza para decir cual de los dos problemas es.
sembrar_trim
FAKE_COUNTS=$(printf '300\t0\t0\t100\t600\t0\t0') corre correr >/dev/null 2>&1
S=$(corre verificar); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0: no es falla" || mal "exit 0: no es falla (rc=$RC)"
tiene "es el de -m"      "se pasa de -m 50"            "$S"
notiene "no el de SIN_AL" "no alinea en ninguna parte" "$S"

echo "== 10. una corrida sin @RG en el BAM se reporta"
# yasma tradeoff agrega por read group: una libreria que no llego al BAM
# desaparece del analisis y nadie lo dice.
sembrar_trim; corre correr >/dev/null 2>&1
sed -i '/SRR_P2/d' "$PROY/aa_primario/align/library_stats.txt"
S=$(corre verificar); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "nombra la corrida" "SRR_P2"                 "$S"
tiene "y qué le pasa"     "no tiene @RG en el BAM" "$S"

echo "== 10b. sin el parche de yasma no se alinea"
# La etapa `over` de nativealign.py levanta un bowtie por libreria, no lee su
# salida y no lo espera: quedan todos vivos con el indice en RAM. Eso mato a
# gadmo_duplicado (12 librerias, 670 Mb) con un `Killed` al 96.5%, y el mensaje
# de entonces decia solo "yasma align fallo". Sin el parche no se arranca.
sembrar_trim; : > "$LOG_YA"
nativealign_sin_parche
S=$(corre correr); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "dice que falta el parche"  "falta el parche de yasma"   "$S"
tiene "y cómo aplicarlo"          "./scripts/yasma_parche.py"  "$S"
notiene "sin llamar a yasma"      "CWD="                       "$(cat "$LOG_YA")"
"$RAIZ/scripts/yasma_parche.py" >/dev/null 2>&1

echo "== 10c. la RAM se dice ANTES de alinear, no después del Killed"
# `unique_d` es un entero de Python por base del genoma (8 B medidos) y se arma
# entero al principio: el numero se sabe en cuanto existe el FASTA, sin leer un
# solo read. El fallo que esto evita no da un mensaje, da `Killed` a las horas.
S=$(corre genoma)
tiene "genoma reporta la RAM"     "RAM de yasma align"         "$S"
sembrar_trim
S=$(corre correr)
tiene "correr también"            "RAM: pide"                  "$S"
S=$(corre plan)
tiene "y plan trae la columna"    "RAM_GB"                     "$S"

echo "== 11. errores"
tiene "modo desconocido"  "modo desconocido"  "$(corre nosequé 2>&1 || true)"
S=$(corre plan noexiste 2>&1 || true)
tiene "un filtro vacío falla" "no encontró ningún proyecto" "$S"
S=$(MANIFEST="$TMP/nada.tsv" GENOMAS_LEDGER="$LED" GENOMES_DIR="$GEN" PROY_DIR="$PROY" \
    bash "$A" plan 2>&1 || true)
tiene "sin manifiesto no sigue" "no existe" "$S"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
