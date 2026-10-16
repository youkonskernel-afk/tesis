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
tiene(){ [[ -z "$2" || "$2" == *$'\n'* ]] && { mal "$1 — patrón vacío o con salto de línea: grep lo parte en varios y uno vacío matchea todo"; return; }; grep -qF -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
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
# FQ_SEQ cambia la lectura: para los kits 4N hace falta una de mas de 23 nt, o
# despues de cortar 8 queda debajo de la ventana.
# FQ_N lecturas iguales (1 por defecto).
q="${FQ_SEQ:-ACGTACGTACGTACGTACGTAC}"
mkdir -p "$O"; : > "$O/$N"
for ((k = 1; k <= ${FQ_N:-1}; k++)); do
  printf '@r%d\n%s\n+\n%s\n' "$k" "$q" "$(printf 'I%.0s' $(seq ${#q}))" >> "$O/$N"
done
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
# Sin nada que recortar solo dice la version. Con -a (lo que llaman `muestra` y
# la primera pasada de un kit 4N) deja la mitad de las lecturas, como si la otra
# mitad no tuviera adaptador (--trimmed-only). -u N corta N del 5' y -u -N del
# 3', y --minimum-length filtra, como el de verdad: es lo que se prueba en la
# segunda pasada. Sin -o escribe a stdout y el reporte va a stderr, tambien
# como el de verdad (medido con cutadapt 5.2); con -o, el reporte a stdout.
[[ " $* " == *" -a "* || " $* " == *" -u "* || " $* " == *" -o "* ]] || { echo "cutadapt 4.0"; exit 0; }
echo "CA $*" >> "$LOG_CA"
# El programa de python le llega por stdin (el heredoc), asi que la entrada '-'
# de la segunda pasada se guarda antes en un fichero.
args=("$@")
if [[ "${args[-1]}" == "-" ]]; then t=$(mktemp); cat > "$t"; args[-1]="$t"; fi
python3 - "${args[@]}" <<'PY'
import gzip, sys
a = sys.argv[1:]
out = a[a.index('-o') + 1] if '-o' in a else None
inp = a[-1]
if inp == '-':
    txt = sys.stdin.read()
elif inp.endswith('.gz'):
    txt = gzip.open(inp, 'rt').read()
else:
    txt = open(inp).read()
lineas = txt.splitlines()
recs = [lineas[i:i + 4] for i in range(0, len(lineas), 4)]
n_in = len(recs)
if '-a' in a:
    recs = recs[::2]
minimo = int(a[a.index('--minimum-length') + 1]) if '--minimum-length' in a else 0
for i, x in enumerate(a):
    if x == '-u':
        c = int(a[i + 1])
        recs = [[r[0], r[1][c:], r[2], r[3][c:]] if c > 0 else [r[0], r[1][:c], r[2], r[3][:c]]
                for r in recs]
recs = [r for r in recs if len(r[1]) >= minimo]
texto = ''.join('\n'.join(r) + '\n' for r in recs)
if out is None:
    sys.stdout.write(texto)
else:
    with (gzip.open(out, 'wt') if out.endswith('.gz') else open(out, 'w')) as h:
        h.write(texto)
rep = (f"Total reads processed:                   {n_in:,}\n"
       f"Reads written (passing filters):         {len(recs):,}\n")
(sys.stderr if out is None else sys.stdout).write(rep)
PY
CA
cat > "$TMP/bin/fastq-dump" <<'FQD'
#!/usr/bin/env bash
# -X N lecturas a stdout (-Z). Anota como fue llamado.
echo "FQD $*" >> "$LOG_FQD"
n=0; while [[ $# -gt 0 ]]; do case "$1" in -X) n="$2"; shift 2 ;; *) shift ;; esac; done
for ((i = 0; i < n; i++)); do printf '@r%d\nACGTACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIIIIIII\n' "$i"; done
FQD
chmod +x "$TMP/bin"/*
export PATH="$TMP/bin:$PATH"
export LOG_FQ="$TMP/fq.log" LOG_YA="$TMP/ya.log" LOG_CA="$TMP/ca.log" LOG_FQD="$TMP/fqd.log"
: > "$LOG_FQ"; : > "$LOG_YA"; : > "$LOG_CA"; : > "$LOG_FQD"

# PRESUPUESTO_GB=2 con corridas de ~2 GB estimados -> una tanda por corrida, que
# es lo que ejercita el ledger. SOLAPAR=0 para que el log de fasterq sea
# determinista; el solapado tiene su propio chequeo.
EXT="$TMP/extremos.tsv"     # no existe salvo en el 17: no cortar es lo normal
corre() {
  MANIFEST="$MAN" ADAPTADORES_TSV="$TABLA" EXTREMOS_TSV="$EXT" SRA_DEST="$DEST" PROY_DIR="$TMP/trim" \
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
# Con 3 tandas de una corrida cada una, "el primero de la tanda" ES la corrida:
# el caso de arriba no distingue y la mutacion pasaba. Hace falta UNA tanda que
# cruce los dos BioProjects, que es justo maggi_primario.
mv "$TMP/trim/aa_primario" "$TMP/aa_primario.bak"
S7=$(PRES=100 corre correr aa/primario)
tiene "una sola tanda con los dos BioProjects" "en 1 tanda(s)" "$S7"
L7=$(cat "$TMP/trim/aa_primario/recortadas.tsv" 2>/dev/null)
tiene "y aun así SRR_A3 -> PRJ_A2" "$(printf 'SRR_A3\tPRJ_A2\t')" "$L7"
rm -rf "$TMP/trim/aa_primario"; mv "$TMP/aa_primario.bak" "$TMP/trim/aa_primario"

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

echo "== 14c2. rehacer solo borra lo que está en trim/, aunque la fila no diga PRE-TRIMMED"
# La segunda defensa: una fila con conteos cuyo fichero apunta afuera de trim/
# (una ruta absoluta, un ledger viejo) tampoco es salida de YASMA. Sin este
# caso, la mutacion que la sacaba pasaba: 14c solo cubre la primera defensa.
A="$TMP/trim/aa_primario"
mkdir -p "$A/untrimmed"; printf 'ORIGINAL' > "$A/untrimmed/SRR_FUERA.fastq.gz"
printf 'SRR_FUERA\tPRJ\tuntrimmed/SRR_FUERA.fastq.gz\t1000\t500\t50.0\t2026-10-02T00:00:00Z\n' \
  >> "$A/recortadas.tsv"
corre rehacer aa/primario SRR_FUERA >/dev/null 2>&1
[[ -f "$A/untrimmed/SRR_FUERA.fastq.gz" ]] && ok "lo de afuera de trim/ sigue en disco" \
  || mal "lo de afuera de trim/ sigue en disco"
notiene "y sale del registro igual" "SRR_FUERA" "$(cat "$A/recortadas.tsv")"
rm -f "$A/untrimmed/SRR_FUERA.fastq.gz"

echo "== 14d. rehacer exige las corridas, no rehace un proyecto entero"
S=$(corre rehacer aa/primario 2>&1 || true)
tiene "pide los RUN"          "decime QUE corridas"  "$S"
S=$(corre rehacer aa 2>&1 || true)
tiene "y exige <org>/<rol>"   "uso:"                 "$S"
S=$(corre rehacer aa/primario SRR_NOEXISTE 2>&1 || true)
tiene "y avisa si no estaba"  "no estaban en el registro"  "$S"

echo "== 14e. un proyecto con kit RA3 y ALGUNAS corridas ya recortadas"
# Es gadmo PRJNA328800 medido: 6 de 12 corridas traen reads de 22 nt sin
# adaptador —vienen ya recortadas— y con la fila RA3 del proyecto cutadapt
# (--trimmed-only) las tiraba casi enteras. La fila PRE-TRIMMED de la corrida
# tiene que ganarle a la del proyecto, su fastq tiene que quedar (ES su
# salida), y verificar no la puede marcar VACIA.
tabla_completa
printf 'aa\tPRJ_A\tSRR_A2\t-\tPRE-TRIMMED\t0\t-\t100\tYA RECORTADA\t2026-09-30\n' >> "$TABLA"
# El escenario 14 borra SRR_A2.sra a proposito: sin reponerlo, esta corrida
# ni entra y el escenario pasaria por no probar nada.
head -c 2048 /dev/zero > "$DEST/aa/SRR_A2.sra"
rm -rf "$TMP/trim/aa_primario"; : > "$LOG_YA"
S=$(corre correr aa/primario); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
LY=$(cat "$LOG_YA")
tiene "la ya recortada va PRE-TRIMMED"   '"SRR_A2.fastq.gz": "PRE-TRIMMED"'   "$LY"
tiene "y la otra sigue con RA3"          "\"SRR_A1.fastq\": \"$RA3\""         "$LY"
[[ -s "$TMP/trim/aa_primario/untrimmed/SRR_A2.fastq.gz" ]] \
  && ok "su fastq queda, comprimido" || mal "su fastq queda, comprimido"
[[ ! -e "$TMP/trim/aa_primario/untrimmed/SRR_A1.fastq" ]] \
  && ok "el de la RA3 se borra" || mal "el de la RA3 se borra"
S=$(corre verificar aa/primario)
L2=$(grep 'SRR_A2' <<<"$S")
tiene "verificar la reconoce"            "pre-trimmed"   "$L2"
notiene "y no la marca VACIA"            "VACIA"         "$L2"

echo "== 16. muestra: N lecturas por corrida, el mismo recorte que yasma"
# Para la prueba 4N sin recortar y alinear el proyecto entero (~2.5 h en rhirr).
tabla_completa
: > "$LOG_CA"; : > "$LOG_FQD"
S=$(corre muestra aa/primario 6); RC=$?
M="$TMP/trim/_muestras/aa_primario"
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC): $S"
LED=$(cat "$M/recortadas.tsv" 2>/dev/null)
tiene "ledger con las columnas del de verdad" \
  "$(printf 'run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc')" "$LED"
tiene "las tres corridas, también la del otro BioProject" \
  "$(printf 'SRR_A3\tPRJ_A2\ttrim/SRR_A3.t.fq.gz\t6\t3\t50.0')" "$LED"
tiene "fastq-dump con -X N a stdout"    "-X 6 -Z $DEST/aa/SRR_A1.sra"  "$(cat "$LOG_FQD")"
CA=$(cat "$LOG_CA")
tiene "cutadapt con los argumentos de yasma trim" \
  "-a $RA3 --minimum-length 15 --maximum-length 50 -j 2 -O 4 --max-n 0 --trimmed-only" "$CA"
tiene "y cada BioProject con SU adaptador" "-a $UNIV" "$CA"
[[ -s "$M/trim/SRR_A1.t.fq.gz" && ! -e "$M/untrimmed/SRR_A1.fq" ]] \
  && ok "deja el .t.fq.gz y borra el crudo" || mal "deja el .t.fq.gz y borra el crudo"
printf 'aa\tPRJ_A\tSRR_A2\t-\tPRE-TRIMMED\t0\t-\t100\tYA RECORTADA\t2026-09-22\n' >> "$TABLA"
: > "$LOG_CA"
S=$(corre muestra aa/primario 6)
LED=$(cat "$M/recortadas.tsv")
tiene "la fila de la corrida le gana: SRR_A2 va sin recortar" \
  "$(printf 'SRR_A2\tPRJ_A\tuntrimmed/SRR_A2.fq.gz\t6\t6')" "$LED"
notiene "y cutadapt no la toca"          "SRR_A2"                       "$(cat "$LOG_CA")"
tabla_completa
S=$(corre muestra aa/primario 4 SRR_A3)
LED=$(cat "$M/recortadas.tsv")
tiene "con RUNs, solo esas"              "SRR_A3"                       "$LED"
notiene "y no las otras"                 "SRR_A1"                       "$LED"
S=$(corre muestra bb/primario 5)
tiene "una PRE-TRIMMED no se recorta" \
  "$(printf 'SRR_B1\tPRJ_B\tuntrimmed/SRR_B1.fq.gz\t5\t5')" "$(cat "$TMP/trim/_muestras/bb_primario/recortadas.tsv")"
mv "$DEST/aa/SRR_A1.sra" "$DEST/aa/SRR_A1.sra.no"
S=$(corre muestra aa/primario 6); RC=$?
mv "$DEST/aa/SRR_A1.sra.no" "$DEST/aa/SRR_A1.sra"
tiene "sin el .sra lo dice"              "FALTA el .sra: SRR_A1"        "$S"
[[ $RC -eq 0 ]] && ok "y sigue con las otras" || mal "y sigue con las otras (rc=$RC)"
S=$(corre muestra aa/primario cero 2>&1); RC=$?
[[ $RC -ne 0 ]] && ok "un N que no es número corta" || mal "un N que no es número corta"
S=$(corre muestra aa 2>&1); RC=$?
[[ $RC -ne 0 ]] && ok "sin rol corta" || mal "sin rol corta"

echo "== 17. kits 4N: data/extremos.tsv corta las puntas, y lo recortado antes se rehace"
# rhirr_duplicado: 0.17% alineado; cortando 4 nt de cada punta, 82.07%. La fila
# llega DESPUES de recortar (hay que alinear para verlo), asi que el caso normal
# es una corrida ya recortada sin cortar que hay que rehacer.
tabla_completa; rm -rf "$TMP/trim" "$EXT"; : > "$LOG_CA"; : > "$LOG_YA"
Q30=GGGGTTTTACGTACGTACGTACGTACCCCC     # 4 + 22 + 4
export FQ_SEQ=$Q30
fila_ext() {   # org proy run c5 c3
  { [[ -s "$EXT" ]] || printf 'org\tbioproject\trun\tcortar_5p\tcortar_3p\talin_sin_pct\talin_con_pct\tlecturas\tfuente\tfecha_utc\n'
    printf '%s\t%s\t%s\t%s\t%s\t0.17\t82.07\t900\tbanco\t2026-10-16\n' "$@"; } >> "$EXT"
}
corre correr aa/duplicado >/dev/null
D="$TMP/trim/aa_duplicado"
tiene "sin fila: el ledger dice que no se cortó" "$(printf 'SRR_D1\tPRJ_D\t')" "$(cat "$D/recortadas.tsv")"
[[ "$(awk -F'\t' 'NR==2 {print $8}' "$D/recortadas.tsv")" == "-" ]] \
  && ok "columna extremos = -" || mal "columna extremos = - ($(sed -n 2p "$D/recortadas.tsv"))"
fila_ext aa PRJ_D - 4 4

echo "  -- 17a. la fila nueva vuelve desactualizado lo recortado"
S=$(corre plan aa/duplicado)
tiene "plan muestra el corte"            "4/4"                                  "$S"
tiene "y que se rehace"                  "recortada con otro corte: se rehace"  "$S"
S=$(corre estado aa/duplicado)
grep -E '^aa_duplicado ' <<<"$S" | tr -s ' ' | grep -q '^aa_duplicado 0 1' \
  && ok "estado la cuenta como faltante" || mal "estado la cuenta como faltante ($S)"
S=$(corre verificar aa/duplicado); RC=$?
[[ $RC -ne 0 ]] && ok "verificar sale != 0" || mal "verificar sale != 0"
tiene "y dice por qué"  "EXTREMOS DESACTUALIZADOS: cortado -, la tabla pide 4/4"  "$S"

echo "  -- 17b. correr la rehace con cutadapt en dos pasadas, sin yasma"
: > "$LOG_CA"; : > "$LOG_YA"
S=$(FQ_N=4 corre correr aa/duplicado); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC): $S"
tiene "avisa qué rehace"             "SRR_D1: tiene -, la tabla pide 4/4"  "$S"
CA=$(cat "$LOG_CA")
tiene "1ra pasada: el 5' y lo de yasma trim, sin ventana" \
  "-u 4 -a $RA3 -j 2 -O 4 --max-n 0 --trimmed-only $D/untrimmed/SRR_D1.fastq" "$CA"
tiene "2da: el 3' y la ventana, desde stdin" \
  "-u -4 --minimum-length 15 --maximum-length 50 -j 2 -o $D/trim/.SRR_D1.extremos.fq.gz -" "$CA"
notiene "yasma no la ve"              "SRR_D1"   "$(cat "$LOG_YA")"
L=$(gzip -dc "$D/trim/SRR_D1.t.fq.gz" 2>/dev/null | sed -n 2p)
[[ "$L" == TTTTACGTACGTACGTACGTAC ]] && ok "queda el inserto sin las 8 bases" \
  || mal "queda el inserto sin las 8 bases (quedó '$L')"
[[ ! -e "$D/trim/.SRR_D1.extremos.fq.gz" ]] && ok "sin el temporal" || mal "sin el temporal"
[[ ! -e "$D/untrimmed/SRR_D1.fastq" ]] && ok "y borra el fastq crudo" || mal "y borra el fastq crudo"
R=$(grep '^SRR_D1' "$D/recortadas.tsv")
# 4 lecturas, la 1ra pasada tira la mitad (sin adaptador) y la 2da ninguna: lo
# que entra es lo de la 1ra y lo que sale lo de la 2da.
[[ "$(cut -f4-6,8 <<<"$R")" == "$(printf '4\t2\t50.0\t4/4')" ]] \
  && ok "ledger: entra 4, sale 2, cortada 4/4" || mal "ledger: entra 4, sale 2, cortada 4/4 ($R)"
tiene "inputs.json la tiene" "trim/SRR_D1.t.fq.gz" "$(cat "$D/inputs.json")"
S=$(corre verificar aa/duplicado); RC=$?
[[ $RC -eq 0 ]] && ok "verificar pasa" || mal "verificar pasa (rc=$RC): $S"
tiene "y dice que se cortó"  "ok, extremos 4/4 cortados"  "$S"
S=$(corre correr aa/duplicado)
tiene "idempotente"  "nada que recortar"  "$S"

echo "  -- 17c. la cuenta usa la 2da pasada: lo que no entra en la ventana no sobrevive"
# 4 + 14 + 4: con el adaptador sale de 22 y pasaria la ventana; el inserto de
# verdad son 14 nt, y la ventana va sobre el inserto.
corre rehacer aa/duplicado SRR_D1 >/dev/null
FQ_SEQ=GGGGACGTACGTACGTACCCCC corre correr aa/duplicado >/dev/null
R=$(grep '^SRR_D1' "$D/recortadas.tsv")
[[ "$(cut -f4,5 <<<"$R")" == "$(printf '1\t0')" ]] \
  && ok "entra 1, sale 0" || mal "entra 1, sale 0 ($R)"
S=$(corre verificar aa/duplicado)
tiene "y una corrida cortada que queda vacía es VACIA" "VACIA" "$S"

echo "  -- 17d. sacar la fila tambien desactualiza: vuelve a yasma"
rm -f "$EXT"; : > "$LOG_YA"
S=$(corre correr aa/duplicado)
tiene "avisa"            "SRR_D1: tiene 4/4, la tabla pide -"  "$S"
tiene "yasma la recorta" "TRIM SRR_D1"                          "$(cat "$LOG_YA")"
[[ "$(awk -F'\t' '$1=="SRR_D1" {print $8}' "$D/recortadas.tsv")" == "-" ]] \
  && ok "y el ledger vuelve a -" || mal "y el ledger vuelve a -"

echo "  -- 17e. una tanda mezclada: la 4N por cutadapt, las otras por yasma"
fila_ext aa PRJ_A SRR_A1 4 4
rm -rf "$TMP/trim/aa_primario"; : > "$LOG_YA"; : > "$LOG_CA"
S=$(PRES=100 corre correr aa/primario); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC): $S"
tiene "una tanda"                  "en 1 tanda(s)"           "$S"
LY=$(cat "$LOG_YA")
tiene "yasma recibe las otras dos" "LIBS=untrimmed/SRR_A2.fastq,untrimmed/SRR_A3.fastq" "$LY"
notiene "y no la 4N"               "SRR_A1"                  "$LY"
P="$TMP/trim/aa_primario"
J=$(python3 -c "import json;print(' '.join(json.load(open('$P/inputs.json'))['trimmed_libraries']))")
for r in SRR_A1 SRR_A2 SRR_A3; do tiene "inputs.json tiene $r" "$r.t.fq.gz" "$J"; done
[[ "$(awk -F'\t' 'NR>1 {print $1 ":" $8}' "$P/recortadas.tsv" | sort | tr '\n' ' ')" == "SRR_A1:4/4 SRR_A2:- SRR_A3:- " ]] \
  && ok "ledger: A1 cortada, las otras no" \
  || mal "ledger: A1 cortada, las otras no ($(awk -F'\t' 'NR>1 {print $1 ":" $8}' "$P/recortadas.tsv" | tr '\n' ' '))"

echo "  -- 17f. una tanda SOLO 4N despues de tandas normales no pisa a las anteriores"
# yasma no corre, asi que nadie deja trimmed_libraries con lo de esta tanda; si
# se queda con la lista entera, cosechar reescribe A1 y A2 como PRE-TRIMMED.
rm -f "$EXT"; fila_ext aa PRJ_A2 SRR_A3 4 4
rm -rf "$P"
S=$(corre correr aa/primario); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC): $S"
tiene "3 tandas" "en 3 tanda(s)" "$S"
notiene "A1 y A2 conservan sus conteos" "PRE-TRIMMED" "$(cat "$P/recortadas.tsv")"
[[ "$(awk -F'\t' '$1=="SRR_A3" {print $8}' "$P/recortadas.tsv")" == "4/4" ]] \
  && ok "y A3 queda cortada" || mal "y A3 queda cortada"

echo "  -- 17g. un ledger de antes de la columna cuenta como sin cortar"
cut -f1-7 "$P/recortadas.tsv" > "$P/r7" && mv "$P/r7" "$P/recortadas.tsv"
S=$(corre plan aa/primario)
tiene "A3 vuelve a estar desactualizada" "recortada con otro corte: se rehace" "$S"
S=$(corre verificar aa/primario)
tiene "verificar lo lee igual"  "EXTREMOS DESACTUALIZADOS: cortado -, la tabla pide 4/4"  "$S"

echo "  -- 17h. una PRE-TRIMMED con corte no se inventa"
fila_ext bb PRJ_B - 4 4
S=$(corre plan); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0"
tiene "la nombra" "bb PRJ_B SRR_B1" "$S"
tiene "y por qué" "no pasa por cutadapt" "$S"

echo "  -- 17i. una tabla rota hace fallar, no recorta sin cortar"
rm -f "$EXT"; fila_ext aa PRJ_D - 4 cuarenta
S=$(corre correr aa/duplicado); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0"
tiene "nombra la tabla" "no pude leer $EXT" "$S"
rm -f "$EXT"; unset FQ_SEQ

echo "== 15. errores"
tiene "modo desconocido"            "modo desconocido"             "$(corre nosequé 2>&1 || true)"
S=$(MANIFEST="$TMP/noexiste.tsv" ADAPTADORES_TSV="$TABLA" SRA_DEST="$DEST" \
    PROY_DIR="$TMP/trim" bash "$T" plan 2>&1 || true)
tiene "sin manifiesto no sigue"     "no existe"                    "$S"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
