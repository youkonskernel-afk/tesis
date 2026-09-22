#!/usr/bin/env bash
# trim.sh con yasma/cutadapt/fasterq-dump falsos en el PATH.
#
# Lo que este banco protege, por orden de gravedad:
#  1. Que un proyecto sin adaptador en la tabla HAGA FALLAR el recorte, en vez
#     de recortar con una secuencia adivinada.
#  2. Que un adaptador 5p o sin identificar tambien lo haga fallar: no sirven
#     como -a de cutadapt.
#  3. Que inputs.json quede como YASMA lo espera —rutas relativas al
#     output_directory, adapters como dict— porque si queda mal `yasma trim`
#     descarta librerias sin error.
#  4. Que el fastq no se llame <RUN>_1.fastq: check_paired_end de YASMA lo
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

MAN="$TMP/man.tsv"
{ printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n'
  printf 'aa\tSRR_A1\tPRJ_A\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRR_A2\tPRJ_A\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'bb\tSRR_B1\tPRJ_B\tduplicado\tapl\t100\t100\t35\tncRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
} > "$MAN"

TABLA="$TMP/adaptadores.tsv"
tabla_completa() {
  { printf '# comentario que se tiene que saltear\n'
    printf 'org\tbioproject\tfamilia\tsecuencia\tadapt_pct\tinserto_modal\tveredicto\tfecha_utc\n'
    printf 'aa\tPRJ_A\tRA3\tTGGAATTCTCGGGTGCCAAGG\t98\t22\tPARECE sRNA-seq\t2026-09-22\n'
    printf 'bb\tPRJ_B\t-\tPRE-TRIMMED\t0\t-\tYA RECORTADA\t2026-09-22\n'
  } > "$TABLA"
}

DEST="$TMP/sra"; mkdir -p "$DEST/aa" "$DEST/bb"
for r in SRR_A1 SRR_A2; do head -c 2048 /dev/zero > "$DEST/aa/$r.sra"; done
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
    *) shift ;;
  esac
done
mkdir -p "$O"; printf '@r1\nACGTACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIIIIIII\n' > "$O/$N"
echo "SALIDA=$O/$N" >> "$LOG_FQ"
FD
cat > "$TMP/bin/yasma" <<'YA'
#!/usr/bin/env bash
# Registra el call y simula la salida de yasma trim leyendo inputs.json.
echo "CWD=$PWD" >> "$LOG_YA"
echo "ARGS=$*" >> "$LOG_YA"
python3 - "$PWD" >> "$LOG_YA" <<'PY'
import json, pathlib, sys
d = pathlib.Path(sys.argv[1])
j = json.loads((d / 'inputs.json').read_text())
print("LIBS=" + ",".join(j['untrimmed_libraries']))
print("ADAPTERS=" + json.dumps(j['adapters'], sort_keys=True))
print("VENTANA=%s-%s" % (j.get('min_length'), j.get('max_length')))
(d / 'trim').mkdir(exist_ok=True)
# Imita al yasma real, verificado contra v1.1.1:
#  - el nombre de salida es <RUN>.t.fq.gz ('.t' + el formato, que trae punto)
#  - una libreria PRE-TRIMMED NO produce fichero nuevo: se anota la ruta original
#  - y escribe trimmed_libraries en inputs.json, que es de donde trim.sh lee
trimmed = []
for lib in j['untrimmed_libraries']:
    nombre = pathlib.Path(lib).name
    sec = j['adapters'][nombre]
    base = nombre.split('.')[0]
    print(f"TRIM {base} adapter={sec}")
    if sec == 'PRE-TRIMMED':
        trimmed.append(str((d / lib).resolve()))
        continue
    out = d / 'trim' / (base + '.t.fq.gz')
    out.write_bytes(b'x' * 64)
    trimmed.append(f"trim/{out.name}")
j['trimmed_libraries'] = trimmed
(d / 'inputs.json').write_text(json.dumps(j, indent=1))
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

corre() {
  MANIFEST="$MAN" ADAPTADORES_TSV="$TABLA" SRA_DEST="$DEST" \
  TRIM_DIR="$TMP/trim" bash "$T" "$@" 2>&1
}

echo "== 1. un proyecto sin fila en la tabla hace fallar, no adivina"
{ printf 'org\tbioproject\tfamilia\tsecuencia\tadapt_pct\tinserto_modal\tveredicto\tfecha_utc\n'
  printf 'aa\tPRJ_A\tRA3\tTGGAATTCTCGGGTGCCAAGG\t98\t22\tPARECE sRNA-seq\t2026-09-22\n'
} > "$TABLA"
S=$(corre plan); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "nombra el proyecto que falta" "bb PRJ_B"                       "$S"
tiene "y dice como medirlo"          "perfil --proyectos"             "$S"
notiene "no lista corridas"          "SRR_A1"                         "$S"

echo "== 2. un adaptador 5p tampoco: no sirve como -a"
tabla_completa
sed -i 's/^bb\tPRJ_B\t-\tPRE-TRIMMED/bb\tPRJ_B\t5p:RA5\tGATCGTCGGACTGTAGAACTCTGAAC/' "$TABLA"
S=$(corre plan); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "explica por que"              "dimero o quimera"               "$S"

echo "== 2b. y un sin_identificar igual"
sed -i 's/^bb\tPRJ_B\t5p:RA5\t[A-Z]*/bb\tPRJ_B\t??:sin_identificar\tCGCCTTGGCCGT/' "$TABLA"
S=$(corre plan); RC=$?
[[ $RC -ne 0 ]] && ok "exit != 0" || mal "exit != 0 (rc=$RC)"
tiene "lo nombra"                    "sin_identificar"                "$S"

echo "== 3. plan con la tabla completa"
tabla_completa
S=$(corre plan); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
tiene "lista las 3 corridas"         "SRR_B1"                         "$S"
tiene "muestra la secuencia"         "TGGAATTCTCGGGTGCCAAGG"          "$S"
tiene "el PRE-TRIMMED aparece"       "PRE-TRIMMED"                    "$S"
tiene "cuenta bien"                  "3 corridas: 0 ya recortadas, 3 por recortar, 0 sin .sra" "$S"

echo "== 4. correr: como se llama a yasma y que ve en inputs.json"
S=$(corre correr); RC=$?
[[ $RC -eq 0 ]] && ok "exit 0" || mal "exit 0 (rc=$RC)"
LY=$(cat "$LOG_YA")
tiene "corre parado en el dir del org" "CWD=$TMP/trim/aa"             "$LY"
tiene "y le pasa -o absoluto"          "-o $TMP/trim/aa"              "$LY"
tiene "la ventana va explicita"        "--min_length 15 --max_length 50" "$LY"
tiene "inputs.json trae la ventana"    "VENTANA=15-50"                "$LY"
tiene "rutas relativas al outdir"      "LIBS=untrimmed/SRR_A1.fastq"  "$LY"
notiene "y NO absolutas"               "LIBS=$TMP"                    "$LY"
tiene "adapters es un dict por fichero" '"SRR_A1.fastq": "TGGAATTCTCGGGTGCCAAGG"' "$LY"
tiene "el PRE-TRIMMED llega tal cual"  '"SRR_B1.fastq": "PRE-TRIMMED"' "$LY"

echo "== 5. el fastq no lleva sufijo _1 (check_paired_end lo tomaria por par)"
LF=$(cat "$LOG_FQ")
tiene "se pide como <RUN>.fastq"       "SALIDA=$TMP/trim/aa/untrimmed/SRR_A1.fastq" "$LF"
notiene "nada de _1.fastq"             "_1.fastq"                     "$LF"

echo "== 6. fasterq-dump con -t, o los temporales van al CWD (109 GB una vez)"
tiene "le pasa -t"                     "TMP=$TMP/trim/.tmp"           "$LF"

echo "== 6b. lo recortado se lee de trimmed_libraries, no de un nombre adivinado"
# Adivinar el nombre fallo de dos formas contra el yasma real: .t.fq.gz en vez
# de .tfq.gz, y una PRE-TRIMMED que no deja fichero nuevo. Con el nombre mal,
# nada contaba como recortado: ni idempotencia ni estado.
tiene "el nombre es .t.fq.gz"          "SRR_A1.t.fq.gz"  "$(ls "$TMP/trim/aa/trim")"
notiene "y no .tfq.gz"                 "SRR_A1.tfq.gz"   "$(ls "$TMP/trim/aa/trim")"
[[ ! -e "$TMP/trim/bb/trim/SRR_B1.t.fq.gz" ]] \
  && ok "la PRE-TRIMMED no produce fichero nuevo" \
  || mal "la PRE-TRIMMED no produce fichero nuevo"
S6=$(corre estado)
grep -E '^bb ' <<<"$S6" | tr -s ' ' | grep -q '^bb 1 0' \
  && ok "y aun asi cuenta como recortada" \
  || mal "y aun asi cuenta como recortada (dijo '$(grep -E '^bb ' <<<"$S6" | tr -s ' ')')"

echo "== 7. es idempotente: relanzar no rehace lo hecho"
: > "$LOG_YA"; : > "$LOG_FQ"
S=$(corre correr)
tiene "lo dice"                        "nada que recortar"            "$S"
notiene "no vuelve a llamar fasterq"   "SALIDA="                      "$(cat "$LOG_FQ")"

echo "== 8. estado cuenta lo recortado"
S=$(corre estado)
grep -E '^aa ' <<<"$S" | tr -s ' ' | grep -q '^aa 2 0' \
  && ok "aa: 2 recortadas" || mal "aa: 2 recortadas (dijo '$(grep -E '^aa ' <<<"$S" | tr -s ' ')')"
tiene "total"                          "total: 3 recortadas, 0 faltan" "$S"

echo "== 9. sin el .sra lo dice y no inventa"
rm -rf "$TMP/trim" "$DEST/aa/SRR_A2.sra"
S=$(corre plan)
tiene "marca la que falta"             "FALTA el .sra"                "$S"
tiene "y dice como traerla"            "drive_pull.sh sra"            "$S"

echo "== 10. errores"
tiene "modo desconocido"               "modo desconocido"             "$(corre nosequé 2>&1 || true)"
S=$(MANIFEST="$TMP/noexiste.tsv" ADAPTADORES_TSV="$TABLA" SRA_DEST="$DEST" \
    TRIM_DIR="$TMP/trim" bash "$T" plan 2>&1 || true)
tiene "sin manifiesto no sigue"        "no existe"                    "$S"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
