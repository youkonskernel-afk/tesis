#!/usr/bin/env bash
# trim.sh con yasma/cutadapt/fasterq-dump falsos en el PATH.
#
# Lo que este banco protege, por orden de gravedad:
#  1. Que un proyecto sin adaptador en la tabla HAGA FALLAR el recorte, en vez
#     de recortar con una secuencia adivinada.
#  2. Que un adaptador 5p o sin identificar tambien lo haga fallar: no sirven
#     como -a de cutadapt.
#  3. Que la tanda N no borre el registro de la N-1. `yasma trim` v1.1.1 PISA
#     trimmed_libraries con lo de su propia llamada, asi que sin el ledger
#     acumulativo cada tanda desmiente a la anterior y todo se re-recorta.
#  4. Que el primario y el duplicado vayan a proyectos YASMA distintos: el
#     duplicado es la validacion independiente y no puede compartir anotacion.
#  5. Que dos BioProjects dentro del mismo rol lleven CADA UNO su adaptador
#     (maggi primario son 2011 y 2014, con kits distintos).
#  6. Que `verificar` atrape una retencion vacia, que es como se ve recortar con
#     la secuencia equivocada: cutadapt corre con --trimmed-only y no da error.
#  7. Que inputs.json quede como YASMA lo espera —rutas relativas al
#     output_directory, adapters como dict— porque si queda mal `yasma trim`
#     descarta librerias sin error.
#  8. Que el fastq no se llame <RUN>_1.fastq: check_paired_end de YASMA lo
#     trataria como par y descartaria el _2 de cualquier cosa.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T="$RAIZ/scripts/trim.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
# El `--` es imprescindible: sin el, un patron que empieza con '-' (como
# "-o /ruta") lo toma grep como opcion, aborta, y en notiene() el exit nonzero
# cae en la rama de exito -> el chequeo pasa SIEMPRE. Estaba en los 14 helpers.
tiene(){ grep -qF -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
notiene(){ grep -qF -- "$2" <<<"$3" && mal "$1 — no deberia estar: $2" || ok "$1"; }

RA3=TGGAATTCTCGGGTGCCAAGG
UNIV=AGATCGGAAGAGCACACGTCT

MAN="$TMP/man.tsv"
{ printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n'
  # aa/primario cruza DOS BioProjects, como maggi. base_count ~1 GB por corrida.
  printf 'aa\tSRR_A1\tPRJ_A\tprimario\tapl\t1000\t1000000000\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRR_A2\tPRJ_A\tprimario\tapl\t1000\t1000000000\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRR_A3\tPRJ_A2\tprimario\tapl\t1000\t1000000000\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRR_D1\tPRJ_D\tduplicado\tapl\t1000\t1000000000\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'bb\tSRR_B1\tPRJ_B\tprimario\tapl\t1000\t32000000\t35\tncRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
} > "$MAN"

TABLA="$TMP/adaptadores.tsv"
tabla_completa() {
  { printf '# comentario que se tiene que saltear\n'
    printf 'org\tbioproject\trun\tfamilia\tsecuencia\tadapt_pct\tinserto_modal\tretencion_est\tveredicto\tfecha_utc\n'
    printf 'aa\tPRJ_A\t-\tRA3\t%s\t98\t22\t95\tPARECE sRNA-seq\t2026-09-22\n'  "$RA3"
    printf 'aa\tPRJ_A2\t-\tIllumina_universal\t%s\t90\t22\t88\tPARECE sRNA-seq\t2026-09-22\n' "$UNIV"
    printf 'aa\tPRJ_D\t-\tRA3\t%s\t96\t22\t93\tPARECE sRNA-seq\t2026-09-22\n'  "$RA3"
    printf 'bb\tPRJ_B\t-\t-\tPRE-TRIMMED\t0\t-\t100\tYA RECORTADA\t2026-09-22\n'
  } > "$TABLA"
}

DEST="$TMP/sra"; mkdir -p "$DEST/aa" "$DEST/bb"
for r in SRR_A1 SRR_A2 SRR_A3 SRR_D1; do head -c 2048 /dev/zero > "$DEST/aa/$r.sra"; done
head -c 2048 /dev/zero > "$DEST/bb/SRR_B1.sra"

mkdir -p "$TMP/bin"
cat > "$TMP/bin/fasterq-dump" <<'FD'
#!/usr/bin/env bash
# Deja un fastq con el nombre que le piden, y anota como fue llamado.
O="."; N="out.fastq"
while [[ $# -gt 0 ]]; do
  case "$1" in
    -O) O="$2"; shift 2 ;;
    -o) N="$2"; shift 2 ;;
    -t) echo "TMP=$2" >> "$LOG_FQ"; shift 2 ;;
    -e) echo "HILOS=$2" >> "$LOG_FQ"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$O"; printf '@r1\nACGTACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIIIIIII\n' > "$O/$N"
echo "SALIDA=$O/$N" >> "$LOG_FQ"
FD
cat > "$TMP/bin/yasma" <<'YA'
#!/usr/bin/env bash
# Imita a yasma trim v1.1.1, verificado contra el binario real:
#  - el nombre de salida es <RUN>.t.fq.gz ('.t' + el formato, que trae punto)
#  - una libreria PRE-TRIMMED NO produce fichero nuevo: se anota la ruta original
#  - PISA trimmed_libraries con lo de ESTA llamada (no acumula). Es el defecto
#    contra el que existe el ledger, asi que el falso tiene que reproducirlo o
#    el banco no protege nada.
#  - imprime el bloque de cutadapt por libreria, de donde salen los conteos
echo "CWD=$PWD" >> "$LOG_YA"
echo "ARGS=$*" >> "$LOG_YA"
python3 - "$PWD" "${SEC_BUENA:-}" 2>>"$LOG_YA" <<'PY'
import json, pathlib, sys
d = pathlib.Path(sys.argv[1]); sec_buena = sys.argv[2]
j = json.loads((d / 'inputs.json').read_text())
log = []
log.append("LIBS=" + ",".join(j['untrimmed_libraries']))
log.append("ADAPTERS=" + json.dumps(j['adapters'], sort_keys=True))
log.append("VENTANA=%s-%s" % (j.get('min_length'), j.get('max_length')))
(d / 'trim').mkdir(exist_ok=True)
trimmed = []
for lib in j['untrimmed_libraries']:
    nombre = pathlib.Path(lib).name
    sec = j['adapters'][nombre]
    base = nombre.split('.')[0]
    log.append(f"TRIM {base} adapter={sec}")
    print(f"trimming: {lib}")
    if sec == 'PRE-TRIMMED':
        print(f"  {lib} -> PRE-TRIMMED")
        trimmed.append(lib)
        continue
    # --trimmed-only: si la secuencia no es la de la libreria, no sobrevive nada
    salen = 1000 if (not sec_buena or sec == sec_buena) else 0
    print(f"Total reads processed:                   1,000")
    print(f"Reads written (passing filters):         {salen:,}")
    out = d / 'trim' / (base + '.t.fq.gz')
    out.write_bytes(b'x' * 64)
    trimmed.append(f"trim/{out.name}")
j['trimmed_libraries'] = trimmed          # PISA, no acumula
(d / 'inputs.json').write_text(json.dumps(j, indent=1))
sys.stderr.write("\n".join(log) + "\n")
PY
YA
cat > "$TMP/bin/cutadapt" <<'CA'
#!/usr/bin/env bash
echo "cutadapt 4.0"
CA
chmod +x "$TMP/bin"/*
export PATH="$TMP/bin:$PATH"
export LOG_FQ="$TMP/fq.log" LOG_YA="$TMP/ya.log"
: > "$LOG_FQ"; : > "$LOG_YA"

# PRESUPUESTO_GB=2 con corridas de ~2 GB estimados -> una tanda por corrida, que
# es lo que ejercita el ledger. SOLAPAR=0 para que el log de fasterq sea
# determinista; el solapado tiene su propio chequeo.
corre() {
  MANIFEST="$MAN" ADAPTADORES_TSV="$TABLA" SRA_DEST="$DEST" PROY_DIR="$TMP/trim" \
  PRESUPUESTO_GB="${PRES:-2}" SOLAPAR="${SOL:-0}" CORES=2 bash "$T" "$@" 2>&1
}

echo "== 1. un proyecto sin fila en la tabla hace fallar, no adivina"
{ printf 'org\tbioproject\trun\tfamilia\tsecuencia\tadapt_pct\tinserto_modal\tretencion_est\tveredicto\tfecha_utc\n'
  printf 'aa\tPRJ_A\t-\tRA3\t%s\t98\t22\t95\tPARECE sRNA-seq\t2026-09-22\n' "$RA3"
} > "$TABLA"
S=$(corre plan); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "nombra el proyecto que falta" "bb PRJ_B"                       "$S"
tiene "y dice como medirlo"          "perfil --proyectos"             "$S"
notiene "no lista corridas"          "SRR_A1"                         "$S"

echo "== 2. un adaptador 5p tampoco: no sirve como -a"
tabla_completa
sed -i "s/^bb\tPRJ_B\t-\t-\tPRE-TRIMMED/bb\tPRJ_B\t-\t5p:RA5\tGATCGTCGGACTGTAGAACTCTGAAC/" "$TABLA"
S=$(corre plan); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "explica por que"              "dimero o quimera"               "$S"

echo "== 2b. y un sin_identificar igual"
sed -i "s/^bb\tPRJ_B\t-\t5p:RA5\t[A-Z]*/bb\tPRJ_B\t-\t??:sin_identificar\tCGCCTTGGCCGT/" "$TABLA"
S=$(corre plan); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo nombra"                    "sin_identificar"                "$S"

echo "== 2c. una corrida puede tener su propia fila, y le gana a la del proyecto"
# gadmo PRJNA328800: 6 de 12 corridas retuvieron 0.6-2.0% contra el 51%
# esperado. Un BioProject PUEDE mezclar kits, y `perfil --proyectos` mide una
# sola corrida, asi que las otras 11 nunca se miraron.
tabla_completa
printf 'aa\tPRJ_A\tSRR_A2\tIllumina_universal\t%s\t91\t22\t77\tPARECE sRNA-seq\t2026-09-24\n' \
  "$UNIV" >> "$TABLA"
S=$(corre plan); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
tiene "SRR_A2 lleva la suya"    "SRR_A2       $UNIV"  "$S"
tiene "y SRR_A1 sigue con RA3"  "SRR_A1       $RA3"   "$S"

echo "== 2d. y si la fila de la corrida no sirve, se nombra la CORRIDA"
# Con la dedup por proyecto que habia antes, esta fila no se miraba nunca:
# la primera corrida del proyecto ya habia marcado PRJ_A como visto.
tabla_completa
printf 'aa\tPRJ_A\tSRR_A2\t5p:RA5\tGATCGTCGGACTGTAGAACTCTGAAC\t91\t22\t77\tPARECE sRNA-seq\t2026-09-24\n' \
  >> "$TABLA"
S=$(corre plan); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "nombra la corrida"  "aa PRJ_A SRR_A2"  "$S"
tiene "y por que"          "dimero o quimera" "$S"

echo "== 3. plan: proyectos por organismo Y rol"
tabla_completa
S=$(corre plan); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
tiene "el primario de aa"            "aa_primario"                    "$S"
tiene "y su duplicado, aparte"       "aa_duplicado"                   "$S"
tiene "lista las 5 corridas"         "SRR_B1"                         "$S"
tiene "muestra la secuencia"         "$RA3"                           "$S"
tiene "el PRE-TRIMMED aparece"       "PRE-TRIMMED"                    "$S"
tiene "estima el fastq en disco"     "FASTQ_GB"                       "$S"
tiene "cuenta bien"                  "5 corridas: 0 ya recortadas, 5 por recortar, 0 sin .sra" "$S"

echo "== 3b. el filtro acepta <org>/<rol>"
S=$(corre plan aa/duplicado)
tiene "trae el duplicado"            "aa_duplicado"                   "$S"
notiene "y NO el primario"           "aa_primario"                    "$S"
S=$(corre plan noexiste); RC=$?
[[ $RC -ne 0 ]] && ok "un filtro que no existe falla" || mal "un filtro que no existe falla (rc=$RC)"

echo "== 4. correr: como se llama a yasma y que ve en inputs.json"
S=$(corre correr); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
LY=$(cat "$LOG_YA")
tiene "corre parado en el dir del proyecto" "CWD=$TMP/trim/aa_primario"    "$LY"
tiene "y le pasa -o absoluto"               "-o $TMP/trim/aa_primario"     "$LY"
tiene "la ventana va explicita"             "--min_length 15 --max_length 50" "$LY"
tiene "inputs.json trae la ventana"         "VENTANA=15-50"                "$LY"
tiene "rutas relativas al outdir"           "LIBS=untrimmed/SRR_A1.fastq"  "$LY"
notiene "y NO absolutas"                    "LIBS=$TMP"                    "$LY"
tiene "adapters es un dict por fichero"     "\"SRR_A1.fastq\": \"$RA3\""   "$LY"
tiene "el PRE-TRIMMED llega tal cual"       '"SRR_B1.fastq.gz": "PRE-TRIMMED"' "$LY"
notiene "y sin --cleanup (iteraria srrs=None)" "--cleanup"                 "$LY"

echo "== 4b. dos BioProjects en el mismo rol, cada uno con SU adaptador"
tiene "PRJ_A va con RA3"        "TRIM SRR_A1 adapter=$RA3"   "$LY"
tiene "PRJ_A2 con el universal" "TRIM SRR_A3 adapter=$UNIV"  "$LY"

echo "== 5. el fastq no lleva sufijo _1 (check_paired_end lo tomaria por par)"
LF=$(cat "$LOG_FQ")
tiene "se pide como <RUN>.fastq"    "SALIDA=$TMP/trim/aa_primario/untrimmed/SRR_A1.fastq" "$LF"
notiene "nada de _1.fastq"          "_1.fastq"                     "$LF"

echo "== 6. fasterq-dump con -t, o los temporales van al CWD (109 GB una vez)"
tiene "le pasa -t"                  "TMP=$TMP/trim/.tmp"           "$LF"
tiene "y -e para no ir de a un hilo" "HILOS=2"                     "$LF"

echo "== 6b. lo recortado se lee de trimmed_libraries, no de un nombre adivinado"
tiene "el nombre es .t.fq.gz"       "SRR_A1.t.fq.gz"  "$(ls "$TMP/trim/aa_primario/trim")"
notiene "y no .tfq.gz"              "SRR_A1.tfq.gz"   "$(ls "$TMP/trim/aa_primario/trim")"
[[ ! -e "$TMP/trim/bb_primario/trim/SRR_B1.t.fq.gz" ]] \
  && ok "la PRE-TRIMMED no produce fichero nuevo" \
  || mal "la PRE-TRIMMED no produce fichero nuevo"
S6=$(corre estado)
grep -E '^bb_primario ' <<<"$S6" | tr -s ' ' | grep -q '^bb_primario 1 0' \
  && ok "y aun asi cuenta como recortada" \
  || mal "y aun asi cuenta como recortada (dijo '$(grep -E '^bb_primario ' <<<"$S6" | tr -s ' ')')"

echo "== 7. la tanda N no borra el registro de la N-1"
# yasma trim PISA trimmed_libraries con lo de su llamada. Con 3 corridas de ~2 GB
# y un presupuesto de 2 GB son 3 tandas: si el script le creyera a inputs.json tal
# como yasma lo deja, solo contaria la ultima.
tiene "3 tandas para aa_primario"  "en 3 tanda(s)"  "$S"
J=$(python3 -c "import json;print(' '.join(json.load(open('$TMP/trim/aa_primario/inputs.json'))['trimmed_libraries']))")
for r in SRR_A1 SRR_A2 SRR_A3; do
  tiene "inputs.json conserva $r"  "$r.t.fq.gz"  "$J"
done
L=$(cat "$TMP/trim/aa_primario/recortadas.tsv")
tiene "el ledger tiene cabecera"   "retencion_pct"  "$L"
[[ $(grep -c . <<<"$L") -eq 4 ]] && ok "y 3 filas" || mal "y 3 filas (tiene $(( $(grep -c . <<<"$L") - 1 )))"

echo "== 7b. el bioproject del ledger es por corrida, no el primero de la tanda"
tiene "SRR_A1 -> PRJ_A"   "$(printf 'SRR_A1\tPRJ_A\t')"   "$L"
tiene "SRR_A3 -> PRJ_A2"  "$(printf 'SRR_A3\tPRJ_A2\t')"  "$L"

echo "== 8. el fastq sin recortar se borra; el de una PRE-TRIMMED NO"
[[ -z "$(ls -A "$TMP/trim/aa_primario/untrimmed")" ]] \
  && ok "untrimmed/ de aa_primario quedo vacio" \
  || mal "untrimmed/ de aa_primario quedo vacio (hay: $(ls "$TMP/trim/aa_primario/untrimmed"))"
# Para una PRE-TRIMMED la salida ES la entrada: borrarla es borrar el resultado.
[[ -s "$TMP/trim/bb_primario/untrimmed/SRR_B1.fastq.gz" ]] \
  && ok "la PRE-TRIMMED sigue ahi, y comprimida" \
  || mal "la PRE-TRIMMED sigue ahi, y comprimida"

echo "== 9. es idempotente: relanzar no rehace lo hecho"
: > "$LOG_YA"; : > "$LOG_FQ"
S=$(corre correr)
tiene "lo dice"                     "nada que recortar"            "$S"
notiene "no vuelve a llamar fasterq" "SALIDA="                     "$(cat "$LOG_FQ")"

echo "== 10. estado cuenta por proyecto"
S=$(corre estado)
grep -E '^aa_primario ' <<<"$S" | tr -s ' ' | grep -q '^aa_primario 3 0' \
  && ok "aa_primario: 3 recortadas" || mal "aa_primario: 3 recortadas (dijo '$(grep -E '^aa_primario ' <<<"$S" | tr -s ' ')')"
tiene "total"                       "total: 5 recortadas, 0 faltan" "$S"

echo "== 11. verificar: la retencion medida contra la esperada"
S=$(corre verificar); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0 con todo bien" || mal "exit 0 con todo bien (rc=$RC)"
tiene "reporta la medida"           "100.0%"                        "$S"
tiene "y la esperada de la tabla"   "95%"                           "$S"
tiene "la PRE-TRIMMED no se compara" "pre-trimmed"                  "$S"
notiene "y no dice PRE-TRIMMED%"    "PRE-TRIMMED%"                  "$S"

echo "== 11b. verificar compara contra la retencion_est de la CORRIDA"
# Sin esto, una corrida con su propia fila se juzga con la expectativa del
# proyecto — justo al reves de para que existe la fila.
printf 'aa\tPRJ_A\tSRR_A1\tIllumina_universal\t%s\t91\t22\t40\tPARECE sRNA-seq\t2026-09-24\n' \
  "$UNIV" >> "$TABLA"
S=$(corre verificar); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0: 100 medido contra 40 esperado" || mal "exit != 0 (rc=$RC)"
tiene "usa el 40 de la corrida"  "40%"      "$S"
tiene "y nombra a SRR_A1"        "SRR_A1"   "$S"
# Las otras del mismo proyecto siguen con el 95 del proyecto.
tiene "SRR_A2 sigue en 95"       "95%"      "$S"
tabla_completa

echo "== 12. verificar atrapa el adaptador equivocado (el fallo que no hace ruido)"
# cutadapt corre con --trimmed-only: una secuencia que no corresponde no da
# error, deja un .t.fq.gz casi vacio y el pipeline sigue. Es el caso de maggi.
rm -rf "$TMP/trim"
S=$(SEC_BUENA="$RA3" corre correr); RC=$?
[[ $RC -eq 0 ]] && ok "el recorte no falla (por eso hace falta verificar)" || mal "el recorte no falla (rc=$RC)"
S=$(corre verificar); RC=$?
[[ $RC -ne 0 ]] && ok "pero verificar sale != 0" || mal "pero verificar sale != 0 (rc=$RC)"
tiene "nombra la corrida vacia"     "SRR_A3"                        "$S"
tiene "y dice que es"               "VACIA"                         "$S"
tiene "y a donde ir"                "perfil --proyectos"            "$S"
tiene "las buenas siguen en ok"     "PRJ_A "                        "$S"

echo "== 13. con SOLAPAR=1 se vuelca la tanda siguiente durante el recorte"
rm -rf "$TMP/trim"; : > "$LOG_FQ"
S=$(SOL=1 corre correr); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
[[ $(grep -c '^SALIDA=' "$LOG_FQ") -eq 5 ]] \
  && ok "vuelca cada corrida una sola vez" \
  || mal "vuelca cada corrida una sola vez (fueron $(grep -c '^SALIDA=' "$LOG_FQ"))"
S=$(corre estado)
tiene "y llega a todas"             "total: 5 recortadas, 0 faltan" "$S"

echo "== 14. sin el .sra lo dice y no inventa"
rm -rf "$TMP/trim" "$DEST/aa/SRR_A2.sra"
S=$(corre plan)
tiene "marca la que falta"          "FALTA el .sra"                "$S"
tiene "y dice como traerla"         "drive_pull.sh sra"            "$S"

echo "== 14b. rehacer saca del registro para que correr las vuelva a hacer"
# El recorte es idempotente por diseno, y eso estorba justo cuando lo que hay
# que rehacer es un recorte MALO: en gadmo/duplicado 6 de 12 salieron con la
# secuencia equivocada y figuraban como hechas.
rm -rf "$TMP/trim"; tabla_completa
corre correr >/dev/null 2>&1
P="$TMP/trim/aa_primario"
[[ -f "$P/trim/SRR_A1.t.fq.gz" ]] && ok "arranca con SRR_A1 recortada" \
  || mal "arranca con SRR_A1 recortada"

S=$(corre rehacer aa/primario SRR_A1); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
tiene "dice cuantas saco"  "sacadas del registro: 1"  "$S"
[[ -f "$P/trim/SRR_A1.t.fq.gz" ]] && mal "borra su .t.fq.gz" || ok "borra su .t.fq.gz"
# SRR_A3 es la otra del mismo proyecto que llego a recortarse: rehacer una no
# puede tocarla.
[[ -f "$P/trim/SRR_A3.t.fq.gz" ]] && ok "y no toca a su hermana SRR_A3" \
  || mal "y no toca a su hermana SRR_A3"
tiene "que sigue en el ledger" "SRR_A3" "$(cat "$P/recortadas.tsv")"
notiene "sale del ledger"     "SRR_A1"  "$(cat "$P/recortadas.tsv")"
# inputs.json tambien: es de donde leen los comandos de aguas abajo.
notiene "y de inputs.json"    "SRR_A1"  "$(cat "$P/inputs.json")"
# Y volver a correr la rehace.
corre correr >/dev/null 2>&1
[[ -f "$P/trim/SRR_A1.t.fq.gz" ]] && ok "y correr la vuelve a recortar" \
  || mal "y correr la vuelve a recortar"

echo "== 14c. rehacer no borra la salida de una PRE-TRIMMED"
# Su `fichero` apunta a untrimmed/, o sea a su ENTRADA: borrarla seria tirar el
# fastq original y dejar la corrida sin forma de rehacerse.
B="$TMP/trim/bb_primario"
_pre=$(awk -F'\t' 'NR>1 {print $3; exit}' "$B/recortadas.tsv" 2>/dev/null || true)
corre rehacer bb/primario SRR_B1 >/dev/null 2>&1
if [[ -n "${_pre:-}" ]]; then
  _q="$_pre"; [[ "$_q" == /* ]] || _q="$B/$_q"
  [[ -f "$_q" ]] && ok "el fastq original sigue en disco" \
    || mal "el fastq original sigue en disco ($_q)"
else
  mal "el ledger de bb_primario no tiene filas"
fi

echo "== 14d. rehacer exige las corridas, no rehace un proyecto entero"
S=$(corre rehacer aa/primario 2>&1 || true)
tiene "pide los RUN"          "decime QUE corridas"  "$S"
S=$(corre rehacer aa 2>&1 || true)
tiene "y exige <org>/<rol>"   "uso:"                 "$S"
S=$(corre rehacer aa/primario SRR_NOEXISTE 2>&1 || true)
tiene "y avisa si no estaba"  "no estaban en el registro"  "$S"

echo "== 15. errores"
tiene "modo desconocido"            "modo desconocido"             "$(corre nosequé 2>&1 || true)"
S=$(MANIFEST="$TMP/noexiste.tsv" ADAPTADORES_TSV="$TABLA" SRA_DEST="$DEST" \
    PROY_DIR="$TMP/trim" bash "$T" plan 2>&1 || true)
tiene "sin manifiesto no sigue"     "no existe"                    "$S"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
