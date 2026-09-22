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
        # El caso real (cloro PRJEB43636) tuvo 5 matches de 20 000: coincidencias
        # al azar, no adaptador. Con 0 exactos el banco no reproduce el problema
        # de que la tabla reporte una familia calculada sobre ese punado.
        if i % 400 == 0:
            s = rnd(19) + AD[:13] + rnd(4)   # 13 nt: el prefijo que busca perfil
    elif 'DIMERO' in run:      # el caso gadmo PRJNA328800: inserto modal 10 nt
        ins = 10 if i % 100 < 45 else random.choice([22, 23, 32, 33])
        s = (rnd(ins) + AD + rnd(150))[:150]
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

echo "== REGRESION: una familia no se reporta sobre un punado de reads"
# cloro PRJEB43636 dio 5 matches de 20 000 —al azar— y con eso la tabla decia
# "RA3 (100% de los que tienen)" mientras el TSV ponia '-'. El mismo comando,
# dos respuestas. El umbral ahora es el mismo que usa el veredicto para decir
# que no hay adaptador, asi que no pueden contradecirse.
S=$(bash "$R" perfil TRIM1 -n 2000 2>&1)
tiene "detecta los pocos matches"  "con adaptador: 5 (0%)"  "$S"
notiene "pero NO reporta familia"  "   adaptador: "         "$S"

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

echo "== inserto modal FUERA de la ventana: pasa, pero lo dice"
# El umbral de PARECE es 30%, asi que un proyecto puede pasar con la mayoria de
# los insertos fuera. gadmo PRJNA328800 pasa con 51% y su inserto modal es de
# 10 nt: "el inserto cae dentro de la ventana" ahi era falso, y al lado de un
# INSERTO de 10 nt en la tabla se leia como contradiccion.
S=$(bash "$R" perfil DIMERO1 -n 2000 2>&1)
tiene "sigue siendo PARECE"        ">>> PARECE sRNA-seq"       "$S"
rex2  "da la fraccion medida"      "el [0-9]+% de"             "$S"
tiene "y avisa del inserto modal"  "el inserto MODAL es de 10 nt, fuera" "$S"
tiene "y que eso se descarta"      "a proposito"               "$S"

echo "== con el inserto modal DENTRO, no avisa de mas"
S=$(bash "$R" perfil SRNA1 -n 2000 2>&1)
notiene "sin el aviso"             "inserto MODAL"             "$S"

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
  # PRJ_E es el caso galga PRJEB12164: casi todo tiene adaptador, pero el
  # inserto modal es un dimero y muere en el piso de 15 nt. La retencion es
  # mucho menor que el adapt_pct.
  printf 'ee\tDIMERO1\tPRJ_E\tprimario\tapl\t100\t100\t150\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
} > "$MAN"

S=$(MANIFEST="$MAN" bash "$R" perfil --proyectos -n 2000 2>&1); RC=$?
echo "$S" | sed -n '/RESUMEN$/,$p'

[[ $(grep -c '^== ' <<<"$S") -eq 5 ]] && ok "perfila 5 proyectos, no 6 corridas" \
  || mal "perfila 5 proyectos (vio $(grep -c '^== ' <<<"$S"))"
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
  # El dimero pasa el veredicto igual, y es el unico con la retencion muy por
  # debajo del adapt_pct: lo necesita el bloque --tsv de mas abajo.
  printf 'ee\tDIMERO1\tPRJ_E\tprimario\tapl\t100\t100\t150\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC\n'
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
[[ $(awk -F'\t' '{print NF}' <<<"$_fila") -eq 9 ]] \
  && ok "9 columnas" || mal "9 columnas (tiene $(awk -F'\t' '{print NF}' <<<"$_fila"))"
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

# La retencion es la INTERSECCION de los dos filtros de cutadapt —tener
# adaptador Y caer en 15-50—, no el adapt_pct. galga PRJEB12164 tiene 95% de
# adaptador y retiene 52%, porque el 20% de sus insertos mide 6-7 nt.
_dim=$(sed -n "/PARA data/,\$p" <<<"$S3" | grep -P "^ee\tPRJ_E\t")
_ap=$(cut -f5 <<<"$_dim"); _ret=$(cut -f7 <<<"$_dim")
[[ "$_ap" == "100" && "$_ret" -lt 70 ]] \
  && ok "retencion < adapt_pct cuando el inserto cae fuera ($_ap% vs $_ret%)" \
  || mal "retencion < adapt_pct (adapt=$_ap ret=$_ret)"
tiene "y la tabla trae la columna" "RETIENE" "$S3"
# La PRE-TRIMMED retiene 100 porque NO se le aplica ningun filtro, ni el de
# longitud. El 0% que salia era la fraccion en ventana, que sin adaptador es 0
# por definicion y no dice nada.
[[ $(cut -f7 <<<"$_trim") == "100" ]] \
  && ok "la PRE-TRIMMED retiene 100, no 0" \
  || mal "la PRE-TRIMMED retiene 100 (dio $(cut -f7 <<<"$_trim"))"
# La tabla y el TSV salen del MISMO comando: no pueden decir cosas distintas.
_tab=$(sed -n "/RESUMEN/,/PARA data/p" <<<"$S3" | grep -E "^dd +PRJ_D")
_tsv=$(sed -n "/PARA data/,\$p" <<<"$S3" | grep -P "^dd\tPRJ_D\t")
grep -qE " -$" <<<"$_tab" && ok "la tabla tampoco inventa familia" \
  || mal "la tabla tampoco inventa familia (dijo '$_tab')"
[[ $(cut -f3 <<<"$_tsv") == "-" ]] && ok "y el TSV dice lo mismo" \
  || mal "y el TSV dice lo mismo (dio '$(cut -f3 <<<"$_tsv")')"

echo "== sin --tsv no imprime el bloque"
notiene "no esta"  "PARA data/adaptadores.tsv"  "$S2"
tiene "y cuenta los 3" "Los 3 proyectos son sRNA-seq" "$S2"

echo; [[ $FALLAS -eq 0 ]] && echo "TODO OK" || { echo "$FALLAS fallas"; exit 1; }
