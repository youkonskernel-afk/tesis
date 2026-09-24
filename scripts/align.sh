#!/usr/bin/env bash
# Alineamiento con `yasma align` (bowtie1 nativo, pesado estilo ShortStack3).
#
#   ./scripts/align.sh genoma    [<org>]            verifica y descomprime el FASTA
#   ./scripts/align.sh plan      [<org>[/<rol>]]    que se haria, sin hacerlo
#   ./scripts/align.sh correr    [<org>[/<rol>]]    alinea
#   ./scripts/align.sh estado    [<org>[/<rol>]]    que esta alineado y que falta
#   ./scripts/align.sh verificar [<org>[/<rol>]]    el BAM contra el manifiesto
#   ./scripts/align.sh ledger                       junta los alineado.tsv en data/
#
# El filtro es un organismo (`galga`) o un organismo y un rol (`galga/primario`),
# igual que en trim.sh. La unidad es el proyecto YASMA <org>_<rol>: un BAM por
# organismo y rol, con un @RG por corrida.
#
# Variables: PROY_DIR, GENOMES_DIR, BAM_DIR, MANIFEST, CORES, MAX_MULTI,
# MAX_RANDOM, UNIQUE_LOCALITY, OFFRATE, FORZAR.
#
# POR QUE `yasma align` Y NO UN BOWTIE PROPIO
# Hasta que se leyo el codigo de v1.1.1, este repo daba por hecho que `yasma
# align` envolvia a ShortStack y que usarlo costaria el `-m 50` medido en danre.
# Las dos cosas son falsas: el wrapper de ShortStack (src/yasma/align.py) esta
# COMENTADO en __init__.py, y el comando que se registra es nativealign.py, un
# bowtie1 nativo cuyo `--max_multi` vale 50 por defecto. O sea que el -m 50 no
# se pierde: es el default, y aca va explicito igual.
#
# Lo que hace, por libreria y contra el mismo indice:
#   1. unique  bowtie -v 1 -S -m 1 --best --strata --offrate N --max <f>
#   2. multi   sobre ese fichero, -m <max_multi> -a --best --strata; pesa cada
#              posicion por la cobertura unica en una ventana de
#              unique_locality/2 y elige por sorteo ponderado. Empate con mas de
#              max_random sitios -> sin mapear (XY:Z:Q)
#   3. over    los que pasan max_multi NO se tiran: van al BAM como no mapeados
#              con XY:Z:H
#
# Tres cosas que esto nos da y que habia que construir a mano:
#   - `-v 1` cuenta mismatches e IGNORA las calidades, asi que la calidad
#     sintetica unica de las 2 corridas SRA Lite (SRR317135, SRR1066790) no
#     cambia el alineamiento. Cierra esa nota de metodos.
#   - Escribe @RG por libreria en la cabecera y por read. `yasma tradeoff` hace
#     header['RG'] sin .get(): un BAM sin read groups no degrada, tira KeyError.
#     El nombre sale de get_rg(), que pela .gz/.t/.fq -> SRR123.t.fq.gz da
#     SRR123. Un @RG por corrida, ya resuelto.
#   - Aplica la ventana 15-50 y descarta los reads con N AL ALINEAR (XY:Z:F).
#     Eso empareja a cloro PRJEB43636 —la unica PRE-TRIMMED, que sale del recorte
#     sin pasar por ninguna ventana— con las otras 18.
#
# DOS COSAS DEL GENOMA, LAS DOS MEDIDAS CONTRA EL BINARIO
#
# 1. TIENE QUE ESTAR ADENTRO DEL -o. `ic.check()` hace
#    value.relative_to(output_directory) sin protegerlo: un genoma compartido
#    fuera del proyecto no da un mensaje, tira ValueError. Se resuelve con un
#    symlink <proyecto>/genome -> <genomes>/<org>, que deja la ruta
#    lexicalmente adentro (validate_path usa .absolute(), no .resolve()) y el
#    indice .ebwt compartido entre el primario y el duplicado.
#
# 2. NO PUEDE IR COMPRIMIDO CON gzip
# make_bam_header() hace pysam.FastaFile(genome_file), y bowtie-build recibe el
# mismo fichero. Nuestros ensamblados estan en <genomes>/<org>/<acc>.fna.gz
# hechos con `gzip -c`, y pysam sobre eso tira OSError (medido). Por eso existe
# el modo `genoma`: verifica el sha256 del .gz contra data/genomas.sha256 —que
# es el ancla de reproducibilidad, el unico registro de contra que se alineo— y
# recien entonces descomprime. El .fna y el indice .ebwt son derivados: no se
# respaldan y estan en .gitignore.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${MANIFEST:-$ROOT/data/srr_manifest.tsv}"
GENOMAS_LEDGER="${GENOMAS_LEDGER:-$ROOT/data/genomas.sha256}"
# shellcheck source=_drive_lib.sh
. "$ROOT/scripts/_drive_lib.sh"
PROY_DIR="${PROY_DIR:-$(ruta_proyectos)}"
GENOMES_DIR="${GENOMES_DIR:-$(ruta_local genomas)}"
# El BAM queda donde YASMA lo deja y ademas se enlaza aca, que es de donde
# drive_push.sh lo sube. Ver enlazar_bam().
BAM_DIR="${BAM_DIR:-$(ruta_local bam)}"
CORES="${CORES:-$(nproc 2>/dev/null || echo 1)}"
# Explicitos aunque coincidan con los defaults de yasma align: un default que
# cambie en una version nueva no avisa, y el -m 50 esta medido en danre — subirlo
# inunda la anotacion de fragmentos de tRNA.
MAX_MULTI="${MAX_MULTI:-50}"
MAX_RANDOM="${MAX_RANDOM:-3}"
UNIQUE_LOCALITY="${UNIQUE_LOCALITY:-50}"
OFFRATE="${OFFRATE:-3}"
FORZAR="${FORZAR:-0}"

LEDGER_ALIN=alineado.tsv

die() { echo "ERROR: $*" >&2; exit 1; }

# --- manifiesto (mismo criterio que trim.sh) ----------------------------------

proyectos() {
  local filtro="${1:-}" org rol=""
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: fetch_runs.sh manifest"
  if [[ "$filtro" == */* ]]; then org="${filtro%%/*}"; rol="${filtro##*/}"; else org="$filtro"; fi
  awk -F'\t' -v o="$org" -v r="$rol" '
    NR>1 && (o=="" || $1==o) && (r=="" || $4==r) {
      k = $1 "\t" $4
      if (!(k in visto)) { visto[k]=1; print k }
    }' "$MANIFEST"
}

corridas() {
  awk -F'\t' -v o="$1" -v r="$2" 'NR>1 && $1==o && $4==r {print $2}' "$MANIFEST"
}

organismos() {
  local filtro="${1:-}"
  proyectos "$filtro" | cut -f1 | awk '!v[$0]++'
}

# --- el genoma ----------------------------------------------------------------

# org -> "<accession>\t<sha256 esperado>"
genoma_de() {
  awk -F'\t' -v o="$1" 'NR>1 && $1==o {print $2"\t"$4; exit}' "$GENOMAS_LEDGER"
}

# Deja en la salida la ruta del FASTA descomprimido, listo para bowtie.
# Verifica el sha256 del .gz ANTES de descomprimir: es lo unico que permite
# decir, dentro de dos años, contra que se alineo.
# Lo unico que va a stdout es la ruta: quien la llama hace fna=$(preparar_genoma
# ...) y cualquier mensaje que se escape ahi se le pega a la ruta. Por eso todo
# lo informativo va a stderr — se ve igual en la terminal y no contamina.
preparar_genoma() {
  local org="$1" acc esp dir gz fna real
  IFS=$'\t' read -r acc esp < <(genoma_de "$org")
  [[ -n "${acc:-}" ]] || die "$org no está en $(basename "$GENOMAS_LEDGER")"
  dir="$GENOMES_DIR/$org"; gz="$dir/$acc.fna.gz"; fna="$dir/$acc.fna"

  if [[ ! -s "$gz" && ! -s "$fna" ]]; then
    die "falta el genoma de $org: $gz
  Traelo con: ./scripts/drive_pull.sh genomas $org --go"
  fi

  if [[ -s "$gz" ]]; then
    real=$(sha256sum "$gz" | cut -d' ' -f1)
    [[ "$real" == "$esp" ]] || die "el sha256 de $gz no coincide con el ledger
  esperado: $esp
  medido  : $real
  No alineo contra un genoma que no es el que quedó registrado."
    echo "   $org  $acc  sha256 ok" >&2
  else
    echo "   $org  $acc  AVISO: solo está el .fna; sin el .gz no puedo comprobar el sha256" >&2
  fi

  if [[ ! -s "$fna" ]]; then
    echo "   descomprimiendo $acc.fna.gz (pysam no lee gzip plano)" >&2
    gzip -cd "$gz" > "$fna.parcial" || { rm -f "$fna.parcial"; die "no pude descomprimir $gz"; }
    mv "$fna.parcial" "$fna"
  fi
  echo "$fna"
}

cmd_genoma() {
  local filtro="${1:-}" org fna
  echo "genomas    : $GENOMES_DIR"
  echo "ledger     : $GENOMAS_LEDGER"
  echo
  for org in $(organismos "$filtro"); do
    fna=$(preparar_genoma "$org")
    printf '   %-7s %s  (%s)\n' "$org" "$(basename "$fna")" \
      "$(du -h "$fna" 2>/dev/null | cut -f1)"
    # El indice se construye solo dentro de yasma align, pero decirlo aca evita
    # la sorpresa de una primera corrida que tarda media hora "sin hacer nada".
    if [[ -f "${fna%.fna}.rev.1.ebwt" ]]; then
      echo "           índice bowtie: ya está"
    else
      echo "           índice bowtie: falta — yasma lo construye con bowtie-build --offrate $OFFRATE"
    fi
  done
}

# --- estado del recorte y del alineamiento ------------------------------------

recortadas_de() {
  local dir="$1"
  [[ -d "$dir" ]] || return 0
  python3 - "$dir" <<'PY'
import json, pathlib, sys
d = pathlib.Path(sys.argv[1])
rutas = []
led = d / 'recortadas.tsv'
if led.is_file():
    for ln in led.read_text().splitlines()[1:]:
        c = ln.split('\t')
        if len(c) > 2:
            rutas.append(c[2])
f = d / 'inputs.json'
if f.is_file():
    try:
        rutas += json.loads(f.read_text()).get('trimmed_libraries') or []
    except (json.JSONDecodeError, OSError):
        pass
visto = set()
for r in rutas:
    q = pathlib.Path(r)
    if not q.is_absolute():
        q = d / q
    if q.is_file() and q.stat().st_size > 0:
        run = q.name.split('.')[0]
        if run not in visto:
            visto.add(run); print(run)
PY
}

bam_de() { echo "$1/align/alignment.bam"; }

# Un BAM viejo respecto de lo recortado es un BAM que no tiene todo adentro. Es
# el caso de "recorte otra tanda y me olvide de re-alinear", que no falla: el
# alineamiento simplemente no incluye las corridas nuevas.
bam_al_dia() {
  local dir="$1" bam nuevo
  bam=$(bam_de "$dir")
  [[ -s "$bam" ]] || return 1
  nuevo=$(find "$dir/trim" "$dir/untrimmed" -newer "$bam" -type f 2>/dev/null | head -1)
  [[ -z "$nuevo" ]]
}

cmd_estado() {
  local filtro="${1:-}" org rol dir n_rec n_man bam est
  printf '%-18s %10s %10s  %s\n' PROYECTO CORRIDAS RECORTADAS BAM
  local t_bam=0 t_proy=0
  while IFS=$'\t' read -r org rol; do
    dir="$PROY_DIR/${org}_${rol}"
    n_man=$(corridas "$org" "$rol" | grep -c .)
    n_rec=$(recortadas_de "$dir" | grep -c . || true)
    bam=$(bam_de "$dir")
    if [[ ! -s "$bam" ]]; then est="falta"
    elif bam_al_dia "$dir"; then est="$(du -h "$bam" | cut -f1)"; t_bam=$((t_bam+1))
    else est="DESACTUALIZADO"; fi
    printf '%-18s %10d %10d  %s\n' "${org}_${rol}" "$n_man" "$n_rec" "$est"
    t_proy=$((t_proy+1))
  done < <(proyectos "$filtro")
  echo
  echo "total: $t_bam de $t_proy proyectos con BAM al día"
}

cmd_plan() {
  local filtro="${1:-}"
  echo "manifiesto : $MANIFEST"
  echo "proyectos  : $PROY_DIR"
  echo "genomas    : $GENOMES_DIR"
  echo "BAMs a     : $BAM_DIR  (para drive_push.sh bam <org>)"
  echo "bowtie     : -v 1, -m $MAX_MULTI, max_random $MAX_RANDOM, locality $UNIQUE_LOCALITY, offrate $OFFRATE, $CORES cores"
  echo
  local org rol dir n_man n_rec acc _esp listos=0 sin_trim=0 hechos=0 total=0
  printf '%-18s %-16s %10s %10s  %s\n' PROYECTO GENOMA CORRIDAS RECORTADAS ESTADO
  while IFS=$'\t' read -r org rol; do
    dir="$PROY_DIR/${org}_${rol}"
    n_man=$(corridas "$org" "$rol" | grep -c .)
    n_rec=$(recortadas_de "$dir" | grep -c . || true)
    IFS=$'\t' read -r acc _esp < <(genoma_de "$org")
    local est
    if [[ -s "$(bam_de "$dir")" ]] && bam_al_dia "$dir" && [[ "$FORZAR" != "1" ]]; then
      est="ya alineado"; hechos=$((hechos+1))
    elif [[ "$n_rec" -eq 0 ]]; then
      est="NADA RECORTADO"; sin_trim=$((sin_trim+1))
    elif [[ "$n_rec" -lt "$n_man" ]]; then
      est="RECORTE INCOMPLETO ($((n_man-n_rec)) faltan)"; sin_trim=$((sin_trim+1))
    else
      est="por alinear"; listos=$((listos+1))
    fi
    printf '%-18s %-16s %10d %10d  %s\n' "${org}_${rol}" "${acc:-SIN GENOMA}" "$n_man" "$n_rec" "$est"
    total=$((total+1))
  done < <(proyectos "$filtro")
  [[ $total -gt 0 ]] || die "el filtro '${filtro:-(todo)}' no encontró ningún proyecto"
  echo
  echo "$total proyectos: $hechos ya alineados, $listos por alinear, $sin_trim sin recorte completo"
  [[ $sin_trim -eq 0 ]] || echo "Recortá primero con: ./scripts/trim.sh correr <org>[/<rol>]" >&2
}

# El BAM se queda donde YASMA lo puso —inputs.json guarda esa ruta absoluta y
# `tradeoff` la lee— y ademas se enlaza a <bams>/<org>/<rol>.bam, que es de donde
# drive_push.sh sube. Un hard link no cuesta disco y las dos rutas quedan
# validas; si el enlace no se puede (otro filesystem) se copia y se avisa,
# porque 340 GB duplicados son una decision y no un detalle.
enlazar_bam() {
  local org="$1" rol="$2" dir="$3" bam destino f
  bam=$(bam_de "$dir")
  mkdir -p "$BAM_DIR/$org"
  for f in "$bam" "$bam.bai"; do
    [[ -s "$f" ]] || continue
    destino="$BAM_DIR/$org/$rol.bam${f##*.bam}"
    rm -f "$destino"
    if ! ln "$f" "$destino" 2>/dev/null; then
      echo "   aviso: no pude hacer hard link, copio ($(du -h "$f" | cut -f1))" >&2
      cp "$f" "$destino"
    fi
  done
  echo "   BAM en $BAM_DIR/$org/$rol.bam"
}

# Que se alineo y contra que. inputs.json guarda la RUTA del genoma, que no dice
# nada dentro de dos años; esto guarda el accession y el sha256 del .gz, que es
# lo que va en metodos, mas los parametros de bowtie que se usaron.
registrar() {
  local dir="$1" org="$2" rol="$3" acc="$4" sha="$5"
  local f="$dir/$LEDGER_ALIN"
  [[ -f "$f" ]] || printf 'proyecto\torg\trol\taccession\tsha256_gz\tmax_multi\tmax_random\tunique_locality\toffrate\tcorridas\tfecha_utc\n' > "$f"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${org}_${rol}" "$org" "$rol" "$acc" "$sha" \
    "$MAX_MULTI" "$MAX_RANDOM" "$UNIQUE_LOCALITY" "$OFFRATE" \
    "$(recortadas_de "$dir" | grep -c . || true)" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$f"
}

cmd_correr() {
  local filtro="${1:-}"
  command -v yasma >/dev/null || die "no está yasma en el PATH (ver docs/yasma.md)"
  command -v bowtie >/dev/null || die "no está bowtie: yasma align lo llama directo"
  command -v bowtie-build >/dev/null || die "no está bowtie-build: hace falta para el índice"

  local org rol dir n_man n_rec fna acc esp n_proy=0
  while IFS=$'\t' read -r org rol; do
    n_proy=$((n_proy+1))
    dir="$PROY_DIR/${org}_${rol}"
    echo "== ${org}_${rol}"

    if [[ -s "$(bam_de "$dir")" ]] && bam_al_dia "$dir" && [[ "$FORZAR" != "1" ]]; then
      echo "   ya alineado (FORZAR=1 para rehacerlo)"
      continue
    fi

    # Alinear un recorte incompleto no falla: sale un BAM al que le faltan
    # corridas, y nada lo dice. Es el mismo modo de fallo que el adaptador
    # equivocado, un paso mas abajo.
    n_man=$(corridas "$org" "$rol" | grep -c .)
    n_rec=$(recortadas_de "$dir" | grep -c . || true)
    if [[ "$n_rec" -lt "$n_man" ]]; then
      echo "   recorte incompleto: $n_rec de $n_man. Corré primero:" >&2
      echo "     ./scripts/trim.sh correr $org/$rol" >&2
      die "no alineo un proyecto a medio recortar"
    fi

    IFS=$'\t' read -r acc esp < <(genoma_de "$org")
    fna=$(preparar_genoma "$org")

    # EL GENOMA TIENE QUE ESTAR ADENTRO DEL -o. `ic.check()` hace
    # value.relative_to(output_directory) SIN protegerlo, asi que un genoma
    # compartido fuera del proyecto tira ValueError y no un mensaje. Es la misma
    # trampa que las librerias de `yasma adapter`, una fase mas abajo.
    #
    # La salida es un symlink al directorio de genomas del organismo, no una
    # copia: `validate_path` hace Path(v).absolute() y NO .resolve(), asi que la
    # ruta queda lexicalmente adentro del proyecto y pasa el chequeo, mientras
    # que el fichero real —y sobre todo el indice .ebwt que yasma construye al
    # lado— siguen siendo UNO por organismo. Con un symlink por fichero en vez
    # de por directorio, bowtie-build correria dos veces por organismo (horas,
    # en un genoma de 1 Gb) y habria 18 indices en vez de 9.
    ln -sfn "$GENOMES_DIR/$org" "$dir/genome"

    # Parado adentro del directorio del proyecto, con -o absoluto: YASMA guarda
    # las rutas relativas al output_directory y despues las abre desde el CWD.
    ( cd "$dir" && yasma align -o "$dir" -g "$dir/genome/$(basename "$fna")" \
        --cores "$CORES" --max_multi "$MAX_MULTI" --max_random "$MAX_RANDOM" \
        --unique_locality "$UNIQUE_LOCALITY" --offrate "$OFFRATE" \
        --min_length 15 --max_length 50 --override </dev/null ) \
      || die "yasma align falló en ${org}_${rol} (log: $dir/align/log.txt)"

    [[ -s "$(bam_de "$dir")" ]] || die "yasma align terminó sin dejar $(bam_de "$dir")"
    registrar "$dir" "$org" "$rol" "$acc" "$esp"
    enlazar_bam "$org" "$rol" "$dir"
  done < <(proyectos "$filtro")
  [[ $n_proy -gt 0 ]] || die "el filtro '${filtro:-(todo)}' no encontró ningún proyecto"
  echo
  cmd_estado "$filtro"
}

# El BAM contra el manifiesto. Alinear contra el genoma equivocado NO falla:
# bowtie alinea mal y sigue, igual que cutadapt con la secuencia equivocada un
# paso mas arriba. Lo que lo delata es la fraccion alineada.
#
# Las columnas salen de align/library_stats.txt, que YASMA escribe con los
# conteos por read group: umap(U) mmap_wg(P) mmap_nw(R) xmap_nw(Q) xmap_ma(H)
# xmap_nv(N) xmap_fr(F).
#
# ALIN y SIN_AL no son complementarias y por eso van las dos. ALIN=(U+P+R) es lo
# que quedo COLOCADO en un locus; SIN_AL=N es lo que no alineo en NINGUNA parte
# del genoma. Lo del medio —Q por encima de max_random, H por encima de -m— si
# alineo, solo que no se coloco. Son diagnosticos distintos: N alto significa
# que los reads no son de este genoma (ensamblado equivocado, contaminacion, o
# el huesped en un experimento de infeccion), mientras que H alto significa un
# genoma repetitivo, que es el fenomeno de los tRF medido en danre. Con una sola
# columna los dos casos se ven igual: "poco alineado".
cmd_verificar() {
  local filtro="${1:-}" org rol dir stats fallas=0 filas=0
  printf '%-18s %-12s %10s %7s %7s %7s %7s  %s\n' \
    PROYECTO CORRIDA READS ALIN SIN_AL ">m$MAX_MULTI" FILTR VEREDICTO
  while IFS=$'\t' read -r org rol; do
    dir="$PROY_DIR/${org}_${rol}"
    stats="$dir/align/library_stats.txt"
    if [[ ! -f "$stats" ]]; then
      echo "   sin alinear: ${org}_${rol}" >&2; fallas=$((fallas+1)); continue
    fi
    # Las corridas del manifiesto que NO aparecen en el BAM. Un @RG que falta es
    # una libreria que se perdio en el camino y el BAM no lo dice.
    local run
    while read -r run; do
      grep -q -- "	$run	" "$stats" \
        || { printf '%-18s %-12s %10s %7s %7s %7s %7s  %s\n' "${org}_${rol}" "$run" - - - - - \
               "FALTA: no tiene @RG en el BAM"; fallas=$((fallas+1)); }
    done < <(corridas "$org" "$rol")

    while IFS=$'\t' read -r _p run u p r q h n f; do
      [[ "$run" == "library" || -z "${f:-}" ]] && continue
      filas=$((filas+1))
      local linea
      linea=$(awk -v u="$u" -v p="$p" -v r="$r" -v q="$q" -v h="$h" -v n="$n" -v f="$f" \
                  -v m="$MAX_MULTI" 'BEGIN{
        tot = u+p+r+q+h+n+f
        if (tot == 0) { print "0\t-\t-\t-\t-\tVACIA — 0 reads en el BAM"; exit }
        al = 100*(u+p+r)/tot; sa = 100*n/tot; ov = 100*h/tot; fr = 100*f/tot
        if (al < 10)       v = "MUY BAJA — ¿el genoma correcto?"
        else if (sa > 50)  v = "ok, pero " int(sa) "% no alinea en ninguna parte (¿reads de otro organismo?)"
        else if (ov > 50)  v = "ok, pero " int(ov) "% se pasa de -m " m " (mirar largos antes de tocarlo)"
        else               v = "ok"
        printf "%d\t%.1f%%\t%.1f%%\t%.1f%%\t%.1f%%\t%s", tot, al, sa, ov, fr, v }')
      IFS=$'\t' read -r tot al sa ov fr ver <<<"$linea"
      [[ "$ver" == ok* ]] || fallas=$((fallas+1))
      printf '%-18s %-12s %10s %7s %7s %7s %7s  %s\n' \
        "${org}_${rol}" "$run" "$tot" "$al" "$sa" "$ov" "$fr" "$ver"
    done < "$stats"
  done < <(proyectos "$filtro")
  echo
  if [[ $fallas -eq 0 ]]; then
    echo "$filas librerías verificadas, ninguna fuera de lo esperado"
    return 0
  fi
  echo "$filas librerías verificadas, $fallas con problemas" >&2
  echo "Una fracción alineada muy baja es, casi siempre, el genoma equivocado:" >&2
  echo "bowtie no falla por eso, alinea mal y sigue. Comprobá el accession con" >&2
  echo "  ./scripts/align.sh genoma $filtro   y  ./scripts/fetch_genomes.sh verificar" >&2
  return 1
}

# Junta los alineado.tsv de cada proyecto en data/alineamientos.tsv.
#
# El per-proyecto vive en proyectos/<org>_<rol>/, que esta en .gitignore —es un
# directorio de trabajo— asi que no hay nada que commitear ahi. Este es el que
# va a git: es texto, es chico, y es el unico registro de contra que ensamblado
# y con que parametros se alineo cada proyecto. Mismo criterio que
# data/genomas.sha256 y data/sra_md5.tsv.
cmd_ledger() {
  local salida="${1:-$ROOT/data/alineamientos.tsv}" org rol dir f n=0
  local tmp; tmp=$(mktemp); trap 'rm -f "$tmp"' RETURN
  while IFS=$'\t' read -r org rol; do
    f="$PROY_DIR/${org}_${rol}/$LEDGER_ALIN"
    [[ -f "$f" ]] || continue
    # Solo la ULTIMA fila de cada proyecto: el per-proyecto acumula una por
    # corrida de align, y lo que vale es con que se alineo el BAM que quedo.
    tail -n +2 "$f" | tail -1 >> "$tmp"
    n=$((n+1))
  done < <(proyectos "")
  # El bloque de comentarios de arriba del fichero se conserva: documenta las
  # columnas y por que existe, y regenerar el ledger no tiene por que borrarlo.
  {
    [[ -f "$salida" ]] && sed -n '/^#/p' "$salida"
    printf 'proyecto\torg\trol\taccession\tsha256_gz\tmax_multi\tmax_random\tunique_locality\toffrate\tcorridas\tfecha_utc\n'
    sort "$tmp"
  } > "$salida.nuevo" && mv "$salida.nuevo" "$salida"
  echo "$n proyecto(s) -> $salida"
  [[ $n -gt 0 ]] || echo "   (ninguno alineado todavía)" >&2
}

[[ $# -ge 1 ]] || { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
case "$1" in
  genoma)    shift; cmd_genoma    "${1:-}" ;;
  plan)      shift; cmd_plan      "${1:-}" ;;
  estado)    shift; cmd_estado    "${1:-}" ;;
  correr)    shift; cmd_correr    "${1:-}" ;;
  verificar) shift; cmd_verificar "${1:-}" ;;
  ledger)    shift; cmd_ledger    "${1:-}" ;;
  -h|--help) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//' ;;
  *) die "modo desconocido: $1 (genoma|plan|correr|estado|verificar|ledger)" ;;
esac
