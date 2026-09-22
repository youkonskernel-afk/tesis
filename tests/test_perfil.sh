#!/usr/bin/env bash
# 'perfil' con un fastq-dump falso que emite librerias sinteticas de tres tipos.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R="$RAIZ/scripts/fetch_runs.sh"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FALLAS=0
ok(){ printf '  ok   %s\n' "$1"; }; mal(){ printf '  MAL  %s\n' "$1"; FALLAS=$((FALLAS+1)); }
tiene(){ grep -qF -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
notiene(){ grep -qF -- "$2" <<<"$3" && mal "$1 — no debia: $2" || ok "$1"; }

mkdir -p "$TMP/bin"
cat > "$TMP/bin/fastq-dump" <<'FD'
#!/usr/bin/env python3
import random, sys
random.seed(7)
run = sys.argv[-1]
AD = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
B = "ACGT"
def rnd(n): return ''.join(random.choice(B) for _ in range(n))
VIEJO = "TCGTATGCCGTCTTCTGCTTG"   # adaptador Illumina de 2011
n = 2000
for i in range(n):
    if 'SRNA' in run:          # sRNA: inserto 21-24 + adaptador, read 150
        ins = random.choice([21, 22, 22, 22, 23, 24])
        s = (rnd(ins) + AD + rnd(150))[:150]
    elif 'MRNA' in run:        # mRNA: inserto mas largo que el read, sin adaptador
        s = rnd(150)
    elif 'TRIM' in run:        # YA RECORTADA: el read ES el inserto, sin adaptador
        s = rnd(random.choice([21, 22, 22, 23, 24, 30]))
    elif 'TRF' in run:         # tRF: inserto 38 nt, dentro de la ventana 15-50
        s = (rnd(random.choice([37, 38, 38, 39])) + AD + rnd(150))[:150]
    elif 'RA5' in run:         # adaptador 5p: dimero, no read-through 3'
        s = (rnd(22) + "GATCGTCGGACTGTAGAACTCTGAAC" + rnd(150))[:150]
    elif 'VIEJO' in run:       # adaptador de 2011
        s = (rnd(22) + VIEJO + rnd(150))[:150]
    else:                      # dudosa: adaptador pero inserto fuera de 15-50
        ins = random.choice([80, 95, 110])
        s = (rnd(ins) + AD + rnd(150))[:150]
    print(f"@r{i}\n{s}\n+\n{'I'*len(s)}")
FD
cat > "$TMP/bin/curl" <<'CURL'
#!/usr/bin/env bash
url="${@: -1}"
printf 'run_accession\tread_count\tbase_count\tlibrary_strategy\tlibrary_layout\tlibrary_source\n'
case "$url" in
  *PRJ_BUENO*)
    # la primera NO pasa el filtro (pocas reads); la segunda si
    printf 'SRNA_MALA\t1000\t20000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
    printf 'SRNA_OK\t9000000\t200000000\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n' ;;
  *PRJ_VACIO*)
    printf 'MRNA_X\t60000000\t9000000000\tRNA-Seq\tPAIRED\tTRANSCRIPTOMIC\n' ;;
esac
exit 0
CURL
chmod +x "$TMP/bin/curl"
chmod +x "$TMP/bin/fastq-dump"
export PATH="$TMP/bin:$PATH" MANIFEST="$TMP/noexiste.tsv"

echo "== libreria de sRNA real (inserto 21-24 nt)"
S=$(bash "$R" perfil SRNA1 -n 2000 2>&1); echo "$S" | sed -n '3,12p'
tiene "detecta adaptador en casi todos"  "con adaptador: 2000 (100%)" "$S"
tiene "el inserto modal es de sRNA"      "22 nt"                      "$S"
tiene "veredicto correcto"               ">>> PARECE sRNA-seq"        "$S"
# CUAL adaptador, no solo cuanto: es el dato que el paso de recorte le pasa a
# fastp, y "98% de adaptador" no sirve para eso.
tiene "nombra el adaptador"              "adaptador: RA3"            "$S"

echo "== libreria de mRNA (sin adaptador visible)"
S=$(bash "$R" perfil MRNA1 -n 2000 2>&1)
tiene "casi ningun adaptador"            "con adaptador: 0 (0%)"      "$S"
tiene "veredicto correcto"               ">>> NO PARECE sRNA-seq"     "$S"
tiene "y explica por que"                "inserto es mas largo que el read" "$S"

echo "== dudosa: hay adaptador pero el inserto es largo"
S=$(bash "$R" perfil OTRA1 -n 2000 2>&1)
tiene "veredicto correcto"               ">>> DUDOSA"                 "$S"

echo "== REGRESION: ya recortada no es mRNA"
S=$(bash "$R" perfil TRIM1 -n 2000 2>&1); echo "$S" | sed -n '3,6p'
tiene "verdicto YA RECORTADA"      ">>> YA RECORTADA"          "$S"
rex2(){ grep -qE -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
rex2  "dice el largo del read"     "largo mediano: 2[0-9] nt"  "$S"
tiene "aclara que no hay que recortar" "no necesita recorte"   "$S"

echo "== mRNA sigue dando NO PARECE (no rompi el caso de SRR23277331)"
S=$(bash "$R" perfil MRNA1 -n 2000 2>&1)
tiene "sigue NO PARECE"            ">>> NO PARECE sRNA-seq"    "$S"
tiene "y nombra el largo"          "reads de 150 nt"           "$S"

echo "== tRF de 38 nt entra: esta dentro de la ventana 15-50 del proyecto"
S=$(bash "$R" perfil TRF1 -n 2000 2>&1)
tiene "PARECE, no DUDOSA"          ">>> PARECE sRNA-seq"       "$S"
tiene "nombra la ventana"          "ventana de fastp (15-50 nt)" "$S"

echo "== inserto realmente largo sigue DUDOSA"
S=$(bash "$R" perfil OTRA1 -n 2000 2>&1)
tiene "DUDOSA"                     ">>> DUDOSA"                "$S"

echo "== adaptador viejo de 2011 se detecta"
S=$(bash "$R" perfil VIEJO1 -n 2000 2>&1)
tiene "lo encuentra"               "con adaptador: 2000 (100%)" "$S"
tiene "y da PARECE"                ">>> PARECE sRNA-seq"        "$S"
tiene "y lo nombra bien"           "adaptador: smallRNA_2011"   "$S"
grep -q "adaptador: RA3" <<<"$S" && mal "no lo confunde con el moderno" \
  || ok "no lo confunde con el moderno"

echo "== un adaptador 5p se detecta pero NO se puede recortar con el"
S=$(bash "$R" perfil RA5_1 -n 2000 2>&1)
tiene "lo nombra como 5p"          "adaptador: 5p:RA5"          "$S"
tiene "y avisa que no sirve"       "No sirve para recortar"     "$S"

echo "== sin adaptador no inventa uno"
for _r in TRIM1 MRNA1; do
  _S=$(bash "$R" perfil $_r -n 2000 2>&1)
  grep -qE "^   adaptador:" <<<"$_S" && mal "$_r: no reporta adaptador" \
    || ok "$_r: no reporta adaptador"
done

echo "== perfilar un BioProject candidato (no esta en el manifiesto)"
S=$(bash "$R" perfil PRJ_BUENO -n 2000 2>&1); echo "$S" | sed -n '1,4p'
tiene "resuelve el proyecto a una corrida" "PRJ_BUENO -> corrida representativa: SRNA_OK" "$S"
tiene "y salta la que no pasa el filtro"   "SRNA_OK"        "$S"
grep -q "SRNA_MALA" <<<"$S" && mal "no usa la que no pasa el filtro" || ok "no usa la que no pasa el filtro"
tiene "perfila de verdad"                  ">>> PARECE sRNA-seq" "$S"

echo "== un proyecto sin corridas que pasen el filtro"
S=$(bash "$R" perfil PRJ_VACIO -n 2000 2>&1 || true)
tiene "lo dice y no inventa" "no encontré en PRJ_VACIO ninguna corrida" "$S"

echo "== errores"
tiene "exige un RUN"  "uso:"  "$(bash "$R" perfil 2>&1 || true)"


# ---- modo --proyectos ----
echo
echo "=== perfil --proyectos"
MAN="$TMP/man.tsv"
{ printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n'
  # 2 corridas del mismo proyecto: solo una tiene que perfilarse
  printf 'aa\tSRNA1\tPRJ_A\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'aa\tSRNA2\tPRJ_A\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'bb\tMRNA1\tPRJ_B\tprimario\tentr\t100\t100\t150\tRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'cc\tOTRA1\tPRJ_C\tduplicado\tapl\t100\t100\t150\tncRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  # PRJ_D es el caso cloro: ya recortada. Va en la tabla porque el veredicto
  # que se escribe ahi es una variable DISTINTA del mensaje '>>>' de arriba, y
  # sin esta fila se puede romper uno sin que el otro se entere.
  printf 'dd\tTRIM1\tPRJ_D\tprimario\tapl\t100\t100\t30\tncRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
} > "$MAN"

S=$(MANIFEST="$MAN" bash "$R" perfil --proyectos -n 2000 2>&1); RC=$?
echo "$S" | sed -n '/RESUMEN$/,$p'

[[ $(grep -c '^== ' <<<"$S") -eq 4 ]] && ok "perfila 4 proyectos, no 5 corridas" \
  || mal "perfila 4 proyectos (vio $(grep -c '^== ' <<<"$S"))"
grep -q "SRNA2" <<<"$S" && mal "no repite proyecto" || ok "no repite proyecto"
rex(){ grep -qE -- "$2" <<<"$3" && ok "$1" || mal "$1 — falta: $2"; }
rex   "PRJ_A pasa"            "PRJ_A .*PARECE sRNA-seq"  "$S"
rex   "PRJ_B marcado"         "PRJ_B .*NO PARECE"        "$S"
rex   "PRJ_C marcado"         "PRJ_C .*DUDOSA"           "$S"
# La columna VEREDICTO de la tabla, no el mensaje '>>>'. Es lo que decide el
# exit code, y era lo unico del arreglo de cloro que el banco no miraba.
rex   "PRJ_D: YA RECORTADA en la tabla" "PRJ_D .*YA RECORTADA" "$S"
rex   "cuenta las corridas"   "PRJ_B .*1 corridas"       "$S"
tiene "muestra el % adaptador" "100%"                     "$S"
tiene "la tabla trae el largo de read" "READ"                "$S"
tiene "y la columna de adaptador"      "ADAPTADOR"           "$S"
rex   "PRJ_A dice cual"       "PRJ_A .*RA3"              "$S"
rex   "PRJ_D (recortada) dice -" "PRJ_D .*YA RECORTADA +-" "$S"
tiene "avisa cuantos fallan"  "2 proyecto(s) sin veredicto favorable"  "$S"
[[ $RC -ne 0 ]] && ok "exit != 0 si alguno falla" || mal "exit != 0 si alguno falla (rc=$RC)"

# todos buenos -> exit 0. Con un YA RECORTADA adentro: es el caso cloro, y si
# ese veredicto dejara de contar como favorable esto tiene que ponerse rojo.
{ printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n'
  printf 'aa\tSRNA1\tPRJ_A\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
  printf 'dd\tTRIM1\tPRJ_D\tprimario\tapl\t100\t100\t30\tncRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
} > "$MAN"
S2=$(MANIFEST="$MAN" bash "$R" perfil --proyectos -n 2000 2>&1); RC2=$?
[[ $RC2 -eq 0 ]] && ok "exit 0 si todos pasan" || mal "exit 0 si todos pasan (rc=$RC2)"
tiene "lo dice" "son sRNA-seq (PARECE o YA RECORTADA)" "$S2"

echo "== --tsv: filas listas para data/adaptadores.tsv"
# Existe para no pasar a mano lo que la herramienta ya midio: copiar 19 filas de
# una tabla formateada es donde se cuela un error que despues no falla.
S3=$(MANIFEST="$MAN" bash "$R" perfil --proyectos --tsv -n 2000 2>&1)
tiene "tiene el bloque"        "PARA data/adaptadores.tsv"   "$S3"
# 8 columnas, las de adaptadores.tsv
_fila=$(sed -n "/PARA data/,\$p" <<<"$S3" | grep -P "^aa\tPRJ_A\t")
[[ $(awk -F'\t' '{print NF}' <<<"$_fila") -eq 8 ]] \
  && ok "8 columnas" || mal "8 columnas (tiene $(awk -F'\t' '{print NF}' <<<"$_fila"))"
[[ $(cut -f3 <<<"$_fila") == "RA3" ]] && ok "columna familia" \
  || mal "columna familia (dio '$(cut -f3 <<<"$_fila")')"
# La secuencia COMPLETA, no el prefijo que se usa para detectar
[[ $(cut -f4 <<<"$_fila") == "TGGAATTCTCGGGTGCCAAGG" ]] \
  && ok "la secuencia completa, no el prefijo" \
  || mal "la secuencia completa (dio '$(cut -f4 <<<"$_fila")')"
[[ $(cut -f6 <<<"$_fila") == "22" ]] && ok "inserto sin la unidad" \
  || mal "inserto sin la unidad (dio '$(cut -f6 <<<"$_fila")')"
# Una ya recortada va PRE-TRIMMED: yasma la pasa de largo sin llamar a cutadapt
_trim=$(sed -n "/PARA data/,\$p" <<<"$S3" | grep -P "^dd\tPRJ_D\t")
[[ $(cut -f4 <<<"$_trim") == "PRE-TRIMMED" ]] \
  && ok "la ya recortada va PRE-TRIMMED" \
  || mal "la ya recortada va PRE-TRIMMED (dio '$(cut -f4 <<<"$_trim")')"
# Su familia y su inserto salen de los poquisimos reads que igual matchearon
# (5 de 20 000 en cloro): no es una medicion, asi que van a '-'.
[[ $(cut -f3 <<<"$_trim") == "-" ]] && ok "y sin familia inventada" \
  || mal "y sin familia inventada (dio '$(cut -f3 <<<"$_trim")')"
[[ $(cut -f6 <<<"$_trim") == "-" ]] && ok "ni inserto inventado" \
  || mal "ni inserto inventado (dio '$(cut -f6 <<<"$_trim")')"
tiene "avisa de las que no se pueden recortar" "NO se pueden recortar" "$S3"

echo "== sin --tsv no imprime el bloque"
notiene "no esta"  "PARA data/adaptadores.tsv"  "$S2"
tiene "y cuenta los 2" "Los 2 proyectos son sRNA-seq" "$S2"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
