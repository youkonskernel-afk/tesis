#!/usr/bin/env bash
# Resuelve los BioProjects de data/organismos.tsv a corridas y las descarga.
# Reemplaza al gen_manifest.sh de R1, que asumia un proyecto por organismo.
#
#   ./scripts/fetch_runs.sh manifest             # consulta ENA -> data/srr_manifest.tsv
#   ./scripts/fetch_runs.sh estado               # que hay bajado y que falta
#   ./scripts/fetch_runs.sh prefetch [ORG] [-n N]  # descarga los .sra
#   ./scripts/fetch_runs.sh diag RUN [RUN...]    # por que falla una corrida
#   ./scripts/fetch_runs.sh ledger [ORG]         # rehace md5 de lo que ya esta bajado
#   ./scripts/fetch_runs.sh buscar 'Especie'     # proyectos de sRNA-seq de una especie
#   ./scripts/fetch_runs.sh perfil RUN|PRJ [-n N]  # es sRNA-seq de verdad?
#   ./scripts/fetch_runs.sh perfil --proyectos   # una corrida de cada proyecto
#
# El destino de los .sra sale de SRA_DEST (por defecto SRA_CACHE). En Colab se
# apunta al mount de Drive y SRA_STAGING al disco efimero de la VM: prefetch
# escribe en la VM, se valida, y recien ahi se mueve a Drive. Escribir GB
# directo al FUSE de Drive es lento e inestable.
#
#   -n N   baja como mucho N corridas y termina. Para sesiones de Colab, que
#          se mueren solas: la siguiente retoma donde quedo.
#
# 9 organismos x (primario + duplicado) = 18 proyectos, 19 accessions: el
# primario de maggi son dos BioProjects combinados (PRJNA154615 + PRJNA232734).
#
# Requiere: curl, jq, y prefetch (sra-tools) para la fase de descarga.
# La sesion cloud NO tiene red hacia la ENA; esto corre en la maquina local.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Sobreescribible para poder probar contra una spec sintetica, igual que
# MANIFEST y que GENOMES_SPEC en fetch_genomes.sh.
SPEC="${ORGANISMOS:-$ROOT/data/organismos.tsv}"
# Corridas excluidas a mano, con motivo medido. Ver el encabezado del archivo.
EXCLUIDAS="${EXCLUIDAS:-$ROOT/data/excluidas.tsv}"
MANIFEST="${MANIFEST:-$ROOT/data/srr_manifest.tsv}"
CACHE="${SRA_CACHE:-/home/dev/sra_cache}"
# Destino final de los .sra. En Colab: /content/drive/MyDrive/tesis/80_sra
DEST="${SRA_DEST:-$CACHE}"
# Donde prefetch escribe primero. En Colab: disco efimero de la VM.
STAGING="${SRA_STAGING:-$DEST}"
LEDGER="${SRA_LEDGER:-$ROOT/data/sra_md5.tsv}"
# Margen de disco libre exigido antes de cada descarga, en GB.
MIN_LIBRE_GB="${MIN_LIBRE_GB:-8}"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_BUSCAR="https://www.ebi.ac.uk/ena/portal/api/search"

# Criterio de seleccion. Cambiarlo cambia el dataset: dejarlo explicito y
# versionado es lo que hace el manifiesto reproducible.
MIN_READS="${MIN_READS:-3000000}"
# Solo datos de RNA. library_source=TRANSCRIPTOMIC es el filtro duro: excluye
# corridas GENOMIC que varios BioProjects mezclan con las de RNA.
FUENTE_OK="TRANSCRIPTOMIC"
# RNA-Seq solo SINGLE: el PAIRED de un proyecto de RNA-Seq no es sRNA-seq y
# contaminaria la anotacion. miRNA-Seq y ncRNA-Seq entran con cualquier layout.
ESTRATEGIAS="miRNA-Seq ncRNA-Seq RNA-Seq"

# El criterio de seleccion en UNA funcion. Si se duplica entre el manifiesto y
# la busqueda, se termina eligiendo un proyecto que despues no entra al
# manifiesto — que es exactamente como sclsc quedo sin duplicado.
# Sin grep: esto corre una vez por corrida y pueden ser miles.
pasa_filtro() {
  local src="$1" strat="$2" layout="$3" rc="$4"
  [[ "$src" == "$FUENTE_OK" ]] || return 1
  [[ " $ESTRATEGIAS " == *" $strat "* ]] || return 1
  [[ "$strat" == "RNA-Seq" && "$layout" != "SINGLE" ]] && return 1
  [[ -n "$rc" && "$rc" =~ ^[0-9]+$ && "$rc" -ge "$MIN_READS" ]] || return 1
  return 0
}

SEP=$'\x1f'
usage() { sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }
die() { echo "error: $*" >&2; exit 1; }

# Campos de organismos.tsv: org especie reino clado rol bioproject runs
# spots_M strategy assembly nota set_modelo
proyectos() {
  grep -v '^[[:space:]]*#' "$SPEC" | tail -n +2 | awk -F'\t' -v s="$SEP" \
    'NF>=12 { printf "%s%s%s%s%s%s%s\n", $1, s, $6, s, $5, s, $12 }'
}

org_valido_manifest() {
  grep -v '^[[:space:]]*#' "$SPEC" | tail -n +2 | cut -f1 | grep -qx "$1"
}

cmd_manifest() {
  command -v curl >/dev/null || die "falta curl"
  command -v jq   >/dev/null || die "falta jq"

  local tmp; tmp=$(mktemp)
  printf 'org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\tavg_len\tstrategy\tlayout\tsource\n' > "$tmp"

  # Las exclusiones se aplican ACA, no editando el manifiesto despues: una
  # edicion a mano la deshace la proxima regeneracion, sin avisar.
  local excl; excl=$(mktemp)
  if [[ -f "$EXCLUIDAS" ]]; then
    grep -v '^[[:space:]]*#' "$EXCLUIDAS" 2>/dev/null | tail -n +2 \
      | awk -F'\t' 'NF>=1 && $1!="" {print $1}' > "$excl" || true
  else
    : > "$excl"
  fi
  local n_excl; n_excl=$(wc -l < "$excl")

  local total=0 desc=0 fuera=0
  while IFS="$SEP" read -r org acc rol setm; do
    [[ "$acc" =~ ^PRJ ]] || { echo "salto $org: accession invalido '$acc'" >&2; continue; }
    echo "== $org $rol $acc" >&2

    local resp
    resp=$(curl -sS --fail --max-time 180 --retry 3 --retry-delay 2 \
      "$ENA?accession=$acc&result=read_run&format=tsv&fields=run_accession,read_count,base_count,library_strategy,library_layout,library_source" \
      ) || die "la ENA no respondio para $acc"

    local n_crudo=0 n_ok=0
    while IFS=$'\t' read -r run rc bc strat layout src; do
      [[ "$run" == "run_accession" || -z "$run" ]] && continue
      n_crudo=$((n_crudo+1))

      pasa_filtro "$src" "$strat" "$layout" "$rc" || continue
      if grep -qx -- "$run" "$excl"; then
        echo "   EXCLUIDA $run (ver $(basename "$EXCLUIDAS"))" >&2
        fuera=$((fuera+1)); continue
      fi

      local avg="NA"
      [[ -n "$bc" && "$bc" =~ ^[0-9]+$ && "$rc" -gt 0 ]] && avg=$(( bc / rc ))

      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$org" "$run" "$acc" "$rol" "$setm" "$rc" "${bc:-NA}" "$avg" \
        "$strat" "$layout" "$src" >> "$tmp"
      n_ok=$((n_ok+1))
    done <<< "$resp"

    echo "   $n_ok / $n_crudo corridas pasan el filtro" >&2
    total=$((total+n_crudo)); desc=$((desc+n_crudo-n_ok))
  done < <(proyectos)

  mv "$tmp" "$MANIFEST"
  rm -f "$excl"
  echo
  echo "manifiesto: $MANIFEST"
  echo "corridas vistas=$total  descartadas=$(( desc - fuera ))  excluidas=$fuera  retenidas=$(( $(wc -l < "$MANIFEST") - 1 ))"
  [[ $n_excl -gt 0 ]] && echo "($n_excl corrida(s) en $(basename "$EXCLUIDAS"), con motivo medido)"
  # Antes esto listaba corrida por corrida las que pasan de 50 nt y cortaba en
  # 20. Con 144 de 416 por encima de ese umbral, el listado se llenaba de
  # librerias normales sin recortar (51 nt = 50 ciclos, 65-75 = 75 ciclos, que
  # fastp resuelve recortando adaptador) y enterraba el unico caso que de
  # verdad no parece sRNA-seq. Agrupar y separar lo sospechoso lo hace visible.
  echo "Longitud de read por proyecto (avg_len = base_count/read_count):"
  awk -F'\t' 'NR>1 && $8!="NA" {
      k=$1"\t"$3"\t"$4"\t"$10; n[k]++
      if ($8+0>mx[k]) mx[k]=$8+0
      if (mn[k]==0 || $8+0<mn[k]) mn[k]=$8+0
    }
    END {
      for (k in n) {
        split(k, a, "\t")
        printf "  %-7s %-14s %-10s %-7s %4d  %s\n", a[1], a[2], a[3], a[4], n[k],
               (mn[k]==mx[k] ? mn[k]" nt" : mn[k]"-"mx[k]" nt")
      }
    }' "$MANIFEST" | sort
  echo

  # PAIRED o reads muy largos: no parece sRNA-seq aunque la ENA lo etiquete
  # miRNA-Seq. El filtro solo exige SINGLE para RNA-Seq, asi que esto pasa.
  local raras; raras=$(awk -F'\t' 'NR>1 && ($10=="PAIRED" || ($8!="NA" && $8+0>200)) {
      printf "  %s %s  %s nt  %s  %s  (%s)\n", $1, $2, $8, $10, $9, $3 }' "$MANIFEST")
  if [[ -n "$raras" ]]; then
    echo "SOSPECHOSAS — PAIRED o reads muy largos, revisar antes de alinear:"
    echo "$raras"
    echo
  fi
  echo "Por organismo y rol:"
  awk -F'\t' 'NR>1 {c[$1" "$4]++} END {for (k in c) printf "  %-22s %4d\n", k, c[k]}' "$MANIFEST" | sort
}

# Un .sra puede estar en el layout de prefetch (<dir>/<RUN>/<RUN>.sra), en el
# plano (<dir>/<RUN>.sra) o en el nuestro por organismo (<dir>/<org>/<RUN>.sra).
# El cache heredado de R1 usa el primero.
ruta_sra() {
  local dir="$1" org="$2" run="$3" c
  for c in "$dir/$org/$run.sra" "$dir/$run/$run.sra" "$dir/$run.sra"; do
    [[ -s "$c" ]] && { echo "$c"; return 0; }
  done
  return 1
}

# Despues de prefetch, sobre STAGING. Es a proposito mas permisiva que
# ruta_sra: ruta_sra la usan cmd_estado y cola contra DEST, que en Colab es el
# mount de Drive, y un find recursivo por cada una de las 416 corridas sobre
# FUSE seria lentisimo. Aca es una sola corrida recien bajada al disco local.
#
# Acepta .sralite, que es un formato valido de SRA y fasterq-dump lee. NO
# acepta los ficheros originales del envio (.fastq.gz, .bam, .sff): el pipeline
# espera un .sra y tratarlos como si lo fueran romperia el contrato en silencio,
# que es justo el modo de falla que estos cambios vienen a sacar.
hallar_descarga() {
  local dir="$1" org="$2" run="$3" c
  for c in "$dir/$org/$run.sra"     "$dir/$run/$run.sra"     "$dir/$run.sra" \
           "$dir/$org/$run.sralite" "$dir/$run/$run.sralite" "$dir/$run.sralite"; do
    [[ -s "$c" ]] && { echo "$c"; return 0; }
  done
  c=$(find "$dir/$run" -maxdepth 2 -type f \
        \( -name '*.sra' -o -name '*.sralite' \) -size +0 -print -quit 2>/dev/null) || true
  [[ -n "$c" ]] && { echo "$c"; return 0; }
  return 1
}

# Ultimas lineas de un log, indentadas. prefetch es ruidoso y lo que importa
# siempre esta al final.
eco_log() {
  local t; t=$(tail -n 15 "$1" 2>/dev/null || true)
  if [[ -n "$t" ]]; then sed 's/^/     | /' <<<"$t"
  else echo "     | (prefetch no imprimió nada)"; fi
}

# Que dejo prefetch en el staging, para una corrida. Cubre el layout de
# directorio y el plano.
listar_staging() {
  local dir="$1" run="$2" dejo
  dejo=$(find "$dir" -maxdepth 2 \( -type f -o -type l \) -name "$run*" \
           -printf '     %10s  %p\n' 2>/dev/null || true)
  echo "${dejo:-     (nada)}"
}

cmd_estado() {
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: $0 manifest"
  echo "manifiesto: $MANIFEST"
  echo "destino   : $DEST"
  echo

  local hay=0 falta=0
  declare -A h f
  while IFS=$'\t' read -r org run resto; do
    [[ "$org" == "org" ]] && continue
    if ruta_sra "$DEST" "$org" "$run" >/dev/null; then
      hay=$((hay+1)); h[$org]=$(( ${h[$org]:-0} + 1 ))
    else
      falta=$((falta+1)); f[$org]=$(( ${f[$org]:-0} + 1 ))
    fi
  done < "$MANIFEST"

  printf "%-8s %6s %6s\n" ORG TIENE FALTA
  local o
  for o in $(awk -F'\t' 'NR>1{print $1}' "$MANIFEST" | sort -u); do
    printf "%-8s %6d %6d\n" "$o" "${h[$o]:-0}" "${f[$o]:-0}"
  done
  echo
  echo "total: $hay bajadas, $falta faltan"
}

libre_gb() { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

# La columna formato distingue sra de sralite. Importa porque el .sralite se
# mueve a DEST como <RUN>.sra igual que los demas —fasterq-dump lee los dos— y
# ahi la extension, que era lo unico que lo decia, desaparece. Recuperarlo
# despues costaria un vdb-dump --info por archivo sobre el mount de Drive.
registrar_md5() {
  local org="$1" run="$2" sum="$3" fmt="${4:-sra}" tab
  tab=$(printf '\t')
  mkdir -p "$(dirname "$LEDGER")"
  touch "$LEDGER"
  {
    grep -v "^org${tab}" "$LEDGER" 2>/dev/null | grep -v "^${org}${tab}${run}${tab}" || true
    printf '%s\t%s\t%s\t%s\t%s\n' "$org" "$run" "$sum" "$fmt" "$(date -u +%Y-%m-%d)"
  } | sort > "$LEDGER.tmp"
  { printf 'org\trun\tmd5\tformato\tfecha_utc\n'; cat "$LEDGER.tmp"; } > "$LEDGER"
  rm -f "$LEDGER.tmp"
}

# Organismos de entrenamiento: son los que definen si el modelo sirve, asi que
# hay un modo de orden que los pone primero. Ver set_modelo en organismos.tsv.
ORGS_ENTRENAMIENTO="gadmo galga maggi"

# Construye la cola de pendientes (lo que NO esta ya en DEST) y la ordena.
#   alfabetico   : por organismo, y dentro del organismo los grandes primero
#   chico        : organismos con menos pendientes primero -> completa antes
#   entrenamiento: gadmo/galga/maggi primero, despues el resto
cola() {
  local filtro="$1" orden="$2" pend="$3"
  : > "$pend"
  local org run rc
  while IFS=$'\t' read -r org run rc; do
    ruta_sra "$DEST" "$org" "$run" >/dev/null && continue
    printf '%s\t%s\t%s\n' "$org" "$run" "$rc" >> "$pend"
  done < <(awk -F'\t' -v o="$filtro" 'NR>1 && (o=="" || $1==o) {print $1"\t"$2"\t"$6}' "$MANIFEST")

  case "$orden" in
    chico)
      awk -F'\t' 'NR==FNR {c[$1]++; next} {print c[$1]"\t"$0}' "$pend" "$pend" \
        | sort -t$'\t' -k1,1n -k2,2 -k4,4nr | cut -f2-
      ;;
    entrenamiento)
      awk -F'\t' -v e="$ORGS_ENTRENAMIENTO" '
        BEGIN { split(e, a, " "); for (i in a) tr[a[i]] = 1 }
        { print (($1 in tr) ? 0 : 1) "\t" $0 }' "$pend" \
        | sort -t$'\t' -k1,1n -k2,2 -k4,4nr | cut -f2-
      ;;
    *)
      sort -t$'\t' -k1,1 -k3,3nr "$pend"
      ;;
  esac
}

cmd_prefetch() {
  command -v prefetch >/dev/null || die "falta prefetch (sra-tools). Activá el entorno srna2."
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: $0 manifest"
  local filtro="$1" limite="$2" horas="${3:-}" orden="${4:-alfabetico}"
  mkdir -p "$DEST" "$STAGING"

  local pend; pend=$(mktemp)
  local ordenada; ordenada=$(mktemp)
  cola "$filtro" "$orden" "$pend" > "$ordenada"
  local total; total=$(wc -l < "$ordenada")
  rm -f "$pend"

  if [[ "$total" -eq 0 ]]; then
    echo "no hay nada pendiente${filtro:+ para $filtro}."
    rm -f "$ordenada"; cmd_estado; return 0
  fi

  local t0; t0=$(date +%s)
  local corte=""
  [[ -n "$horas" ]] && corte=$(awk -v h="$horas" 'BEGIN{printf "%d", h*3600}')

  echo "cola: $total corridas pendientes, orden=$orden${limite:+, tope $limite}${horas:+, corte a las ${horas}h}"
  echo

  local n=0 ok=0 fallos=0 bytes=0
  local org run rc
  while IFS=$'\t' read -r org run rc; do
    [[ -n "$limite" && $n -ge $limite ]] && { echo; echo "corte por -n $limite"; break; }
    local trans=$(( $(date +%s) - t0 ))
    if [[ -n "$corte" && $trans -ge $corte ]]; then
      echo; echo "corte por tiempo (${horas}h). Re-ejecutá para seguir donde quedó."
      break
    fi

    local libre; libre=$(libre_gb "$STAGING")
    if [[ -n "$libre" && "$libre" -lt "$MIN_LIBRE_GB" ]]; then
      echo "ALTO: quedan ${libre}G libres en $STAGING (mínimo $MIN_LIBRE_GB). Cortando." >&2
      break
    fi

    n=$((n+1))
    printf '== [%d/%d] %s %s (%s reads)\n' "$n" "$total" "$org" "$run" "$rc"
    rm -rf "${STAGING:?}/$run"
    # La salida de prefetch se captura, no se tira: un fallo mudo cuesta una
    # corrida entera para diagnosticarse, y las de maggi costaron dos.
    local log; log=$(mktemp)
    if ! prefetch --output-directory "$STAGING" --max-size u "$run" >"$log" 2>&1; then
      echo "   FALLO prefetch (exit != 0)" >&2; eco_log "$log" >&2
      rm -f "$log"; fallos=$((fallos+1)); rm -rf "${STAGING:?}/$run"; continue
    fi

    local src
    if ! src=$(hallar_descarga "$STAGING" "$org" "$run"); then
      # prefetch salio 0 y no dejo un .sra reconocible. Las dos causas
      # plausibles: la corrida solo existe en formato original (.fastq.gz,
      # .bam, .sff), o el resolver no devolvio nada y sra-tools igual salio 0.
      # Para distinguirlas: $0 diag "$run".
      echo "   FALLO: prefetch salió 0 pero no dejó .sra" >&2
      eco_log "$log" >&2
      echo "   quedó en el staging:" >&2
      listar_staging "$STAGING" "$run" >&2
      echo "   para ver por qué:  $0 diag $run" >&2
      rm -f "$log"; fallos=$((fallos+1)); rm -rf "${STAGING:?}/$run"; continue
    fi
    rm -f "$log"

    local fmt=sra
    case "$src" in
      *.sralite)
        fmt=sralite
        echo "   AVISO: $run vino en formato lite (.sralite). Las calidades son" >&2
        echo "          sintéticas, asi que el filtro de calidad de fastp ve una" >&2
        echo "          constante en esta corrida y no en las demás. Queda" >&2
        echo "          anotado en el ledger; declararlo en métodos." >&2 ;;
    esac

    # Un .sra truncado NO falla ruidosamente: alinea de menos. Validar antes de
    # darlo por bueno, y nunca mover al destino uno que no valida.
    if command -v vdb-validate >/dev/null; then
      if ! vdb-validate "$src" >/dev/null 2>&1; then
        echo "   FALLO vdb-validate — descartado, se reintenta después" >&2
        rm -rf "${STAGING:?}/$run" "$src"; fallos=$((fallos+1)); continue
      fi
    fi

    mkdir -p "$DEST/$org"
    local final="$DEST/$org/$run.sra"
    local sz; sz=$(stat -c%s "$src" 2>/dev/null || echo 0)
    if [[ "$src" != "$final" ]]; then
      mv "$src" "$final" || { echo "   FALLO al mover" >&2; fallos=$((fallos+1)); continue; }
      rm -rf "${STAGING:?}/$run"
    fi
    bytes=$((bytes + sz))

    registrar_md5 "$org" "$run" "$(md5sum "$final" | cut -d' ' -f1)" "$fmt"

    # Progreso con ETA: sirve para decidir si conviene subir --horas o cortar.
    trans=$(( $(date +%s) - t0 ))
    awk -v sz="$sz" -v n="$n" -v tot="$total" -v b="$bytes" -v s="$trans" 'BEGIN {
      eta = (n > 0 && s > 0) ? (s/n) * (tot-n) / 60 : 0
      printf "   ok  %.0f MB   acumulado %.1f GB   %.0f min   ETA ~%.0f min\n",
             sz/1e6, b/1e9, s/60, eta
    }'
    ok=$((ok+1))
  done < "$ordenada"
  rm -f "$ordenada"

  echo
  echo "bajadas=$ok fallos=$fallos"
  [[ $ok -gt 0 ]] && echo "Commitear el ledger: git add $LEDGER"
  cmd_estado
}

# Adaptadores 3\' de los kits de sRNA. Se busca el prefijo, no la secuencia
# entera: el read se corta antes de terminarla. La lista incluye los Illumina
# VIEJOS a proposito — hay librerias de 2011 en el dataset, y con solo los
# modernos un 0% puede significar "no esta en la lista" en vez de "no hay".
#   TruSeq Small RNA / generico / NEBNext / Qiagen  +  Illumina v1.5 y GEX
# Prefijos de adaptador 3' con nombre. El nombre no es cosmetico: el paso de
# recorte necesita saber CUAL se encontro para pasarselo a fastp, y "hay
# adaptador en el 98%" no alcanza para eso. Los tres ultimos son librerias
# viejas: agregarlos hizo que maggi PRJNA154615 pasara de 0% a 92%.
# El nombre de cada prefijo esta cotejado contra la tabla de 161 adaptadores de
# YASMA (src/yasma/adapter.py), no puesto de memoria. Dos cosas que salieron de
# ese cotejo y hay que respetar:
#
#  - RA5 es un adaptador 5'. Encontrarlo en un read significa dimero o quimera,
#    NO read-through 3', asi que NO sirve como `-a` de cutadapt. Por eso lleva
#    el prefijo `5p:`, que el recorte usa para no intentar recortar con el.
#  - CGCCTTGGCCGT no aparece en ninguno de los 161. Queda como sonda porque
#    detecta algo, pero marcado `??:` — no se recorta a ciegas con una secuencia
#    cuya procedencia no se puede decir.
#
# El prefijo sirve para DETECTAR; la secuencia que va a cutadapt es la del
# adaptador completo y esta en data/adaptadores.tsv. No son lo mismo: TGGAATTCTCGGG
# es el prefijo compartido por 50 entradas de la familia RPI.
# De familia detectada a la secuencia COMPLETA que se le pasa a cutadapt. No es
# lo mismo que el prefijo de ADAPTADORES: ese detecta, este recorta. Un 5p o un
# sin_identificar no tienen entrada a proposito — no se recorta con ellos, y
# trim.sh lo rechaza mirando la familia.
SECUENCIAS="RA3:TGGAATTCTCGGGTGCCAAGG TruSeq_universal:AGATCGGAAGAGCACACGTCTGAACTCCAGTCA smallRNA_2011:TCGTATGCCGTCTTCTGCTTG"

ADAPTADORES="TGGAATTCTCGGG:RA3 AGATCGGAAGAGC:TruSeq_universal GATCGTCGGACTG:5p:RA5 ATCTCGTATGCCG:smallRNA_2011 TCGTATGCCGTCTTCTGCTTG:smallRNA_2011 CGCCTTGGCCGT:??:sin_identificar"

# Ventana de fastp, que es la que decide que entra al pipeline. El veredicto la
# usa en vez de un 18-30 propio: el proyecto eligio 15-50 a proposito para no
# truncar tRFs (30-40 nt) ni dejar los piRNAs (24-32) sin margen.
VENT_MIN="${VENT_MIN:-15}"
VENT_MAX="${VENT_MAX:-50}"

# Una corrida es sRNA-seq si el adaptador 3\' aparece temprano: el inserto es
# corto y el secuenciador siguio leyendo. Lo que NO se puede saber mirando
# avg_len es justamente esto — un sRNA de 22 nt corrido en 2x150 da avg_len 273
# igual que un mRNA. Pasa con SRR23277331 de prupe, la unica PAIRED de las 416.
#
# Corre donde esten los .sra: en Colab con SRA_DEST al mount de Drive.
# Donde cmd_perfil deja el resumen de una linea, para que --proyectos lo junte.
RESUMEN_PERFIL=""

# De un BioProject a una corrida representativa, aplicando el mismo filtro que
# el manifiesto. Hace falta para evaluar un proyecto CANDIDATO, que por
# definicion todavia no esta en el manifiesto — el caso de sclsc, que necesita
# duplicado y no se puede adoptar sin mirarlo primero.
corrida_de_proyecto() {
  local acc="$1" resp run rc bc strat layout src
  resp=$(curl -sS --fail --max-time 180 --retry 3 --retry-delay 2 \
    "$ENA?accession=$acc&result=read_run&format=tsv&fields=run_accession,read_count,base_count,library_strategy,library_layout,library_source") \
    || return 1
  while IFS=$'\t' read -r run rc bc strat layout src; do
    [[ "$run" == "run_accession" || -z "$run" ]] && continue
    pasa_filtro "$src" "$strat" "$layout" "$rc" || continue
    echo "$run"
    return 0
  done <<< "$resp"
  return 1
}

cmd_perfil() {
  command -v fastq-dump >/dev/null || die "falta fastq-dump (sra-tools)"
  local run="${1:-}" n="${2:-20000}"
  [[ -n "$run" ]] || die "uso: $0 perfil RUN|PRJ [-n N]"

  if [[ "$run" =~ ^PRJ ]]; then
    command -v curl >/dev/null || die "falta curl"
    local proy="$run"
    run=$(corrida_de_proyecto "$proy") \
      || die "no encontré en $proy ninguna corrida que pase el filtro"
    echo "== $proy -> corrida representativa: $run"
  fi

  # Si ya esta bajado se usa el archivo; si no, fastq-dump lo resuelve por red.
  local org src=""
  org=$(awk -F'\t' -v r="$run" 'NR>1 && $2==r {print $1; exit}' "$MANIFEST" 2>/dev/null || true)
  [[ -n "$org" ]] && src=$(ruta_sra "$DEST" "$org" "$run" 2>/dev/null || true)

  echo "== $run${org:+  ($org)}   primeros $n spots"
  [[ -n "$src" ]] && echo "   fuente: $src" || echo "   fuente: la red (no esta bajada)"

  local fq; fq=$(mktemp)
  fastq-dump --split-spot -X "$n" -Z "${src:-$run}" 2>/dev/null > "$fq" || {
    rm -f "$fq"; die "fastq-dump no pudo leer $run"; }

  awk -v ads="$ADAPTADORES" -v vmin="$VENT_MIN" -v vmax="$VENT_MAX" '
    BEGIN {
      na = split(ads, campo, " ")
      for (i = 1; i <= na; i++) {
        # El nombre es TODO lo que sigue al primer ':', porque hay nombres con
        # ':' adentro (5p:RA5, ??:sin_identificar) y eso es el marcador de que
        # no se puede recortar con ellos.
        j = index(campo[i], ":")
        A[i] = substr(campo[i], 1, j - 1); NOM[i] = substr(campo[i], j + 1)
      }
    }
    NR % 4 == 2 {
      total++
      lr[length($0)]++           # longitud del read: sin esto, un 0% no se
      suma_lr += length($0)      # puede interpretar (ver abajo)
      mejor = 0; cual = ""
      for (i = 1; i <= na; i++) {
        p = index($0, A[i])
        if (p > 0 && (mejor == 0 || p < mejor)) { mejor = p; cual = NOM[i] }
      }
      if (mejor > 0) {
        con++; ins = mejor - 1; h[ins]++; ad[cual]++
        if (ins >= vmin && ins <= vmax) dentro++
      }
    }
    END {
      if (total == 0) { print "   sin reads"; exit }

      # Mediana de la longitud del read.
      n = 0; for (k in lr) { largos[n++] = k + 0 }
      for (a = 0; a < n; a++) for (b = a+1; b < n; b++)
        if (largos[b] < largos[a]) { t = largos[a]; largos[a] = largos[b]; largos[b] = t }
      acum = 0; med_lr = largos[0]
      for (a = 0; a < n; a++) { acum += lr[largos[a]]; if (acum >= total/2) { med_lr = largos[a]; break } }

      printf "   reads: %d   largo mediano: %d nt   con adaptador: %d (%.0f%%)\n",
             total, med_lr, con, 100*con/total

      # Cual adaptador, no solo cuanto. Es el dato que el recorte necesita.
      ad_top = "-"; ad_mx = 0
      for (k in ad) if (ad[k] + 0 > ad_mx) { ad_mx = ad[k] + 0; ad_top = k }
      if (con > 0) {
        printf "   adaptador: %s (%.0f%% de los que tienen)\n", ad_top, 100*ad_mx/con
        for (k in ad) if (k != ad_top) printf "     tambien: %s  %.0f%%\n", k, 100*ad[k]/con
        # Un 5p o un sin_identificar no se puede usar como `-a`: decirlo aca en
        # vez de dejar que el recorte lo descubra con una secuencia equivocada.
        if (ad_top ~ /^5p:/) {
          print "   OJO: el mayoritario es un adaptador 5p: es dimero o quimera,"
          print "        no read-through. No sirve para recortar."
        }
        if (ad_top ~ /^[?][?]:/) {
          print "   OJO: el mayoritario no esta identificado. No recortar con el."
        }
      }
      if (con > 0) {
        printf "   inserto (largo antes del adaptador), los mas frecuentes:\n"
        cn = 0
        for (k in h) { ord[cn++] = k }
        for (a = 0; a < cn; a++) for (b = a+1; b < cn; b++)
          if (h[ord[b]] + 0 > h[ord[a]] + 0) { t = ord[a]; ord[a] = ord[b]; ord[b] = t }
        for (a = 0; a < cn && a < 8; a++)
          printf "     %3d nt  %6d  %5.1f%%\n", ord[a], h[ord[a]], 100*h[ord[a]]/total
        printf "   inserto dentro de la ventana de fastp (%d-%d nt): %.0f%%\n",
               vmin, vmax, 100*dentro/total
      }
      print ""

      modal = 0; mx = 0
      for (k in h) if (h[k] + 0 > mx) { mx = h[k] + 0; modal = k }

      if (con >= 0.5*total && dentro >= 0.3*total) {
        ver = "PARECE sRNA-seq"
        print "   >>> PARECE sRNA-seq: el adaptador aparece temprano y el inserto cae"
        printf "       dentro de la ventana %d-%d nt.\n", vmin, vmax
      } else if (con < 0.2*total && med_lr <= vmax) {
        # Sin adaptador PERO reads cortos: el read ES el inserto. Confundir esto
        # con mRNA casi hace tirar las 34 corridas de cloro PRJEB43636.
        ver = "YA RECORTADA"
        printf "   >>> YA RECORTADA: no hay adaptador porque ya se lo sacaron — los reads\n"
        printf "       miden %d nt, que es el inserto. Es sRNA-seq y no necesita recorte.\n", med_lr
      } else if (con < 0.2*total) {
        ver = "NO PARECE"
        printf "   >>> NO PARECE sRNA-seq: sin adaptador y con reads de %d nt, o sea que el\n", med_lr
        print  "       inserto es mas largo que el read. Es lo que se espera de mRNA."
      } else {
        ver = "DUDOSA"
        printf "   >>> DUDOSA: hay adaptador pero el inserto cae fuera de %d-%d nt.\n", vmin, vmax
      }
      printf "RESUMEN\t%.0f\t%s\t%d nt\t%s\t%s\n", 100*con/total,
             (con>0 ? modal" nt" : "-"), med_lr, ver, ad_top > "/dev/stderr"
    }' "$fq" 2> "$fq.res"
  RESUMEN_PERFIL=$(grep '^RESUMEN' "$fq.res" 2>/dev/null | cut -f2- || true)
  cat "$fq.res" | grep -v '^RESUMEN' >&2 || true
  rm -f "$fq" "$fq.res"
}

# Una corrida representativa de CADA proyecto. La etiqueta de estrategia de la
# ENA no establece que tipo de libreria es —SRR23277331 decia miRNA-Seq y no
# tenia un solo read con adaptador—, y en el manifiesto hay 3 proyectos
# primarios enteros etiquetados RNA-Seq, dos de ellos de organismos de
# entrenamiento. Son minutos; alinear a ciegas son 30-40 h.
cmd_perfil_proyectos() {
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: $0 manifest"
  local n="${1:-20000}" tsv="${2:-0}"
  local tabla; tabla=$(mktemp)
  local org proy rol strat run
  while IFS=$'\t' read -r org proy rol strat run; do
    cmd_perfil "$run" "$n" || true
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$org" "$proy" "$rol" "$strat" "$run" \
      "${RESUMEN_PERFIL:-?\t?\t?\tSIN DATO\t?}" >> "$tabla"
    echo
  done < <(awk -F'\t' 'NR>1 && !(($1 FS $3) in v) {v[$1 FS $3]=1;
             print $1"\t"$3"\t"$4"\t"$9"\t"$2}' "$MANIFEST" | sort)

  echo "==================== RESUMEN"
  printf '%-7s %-14s %-10s %-11s %6s %8s %7s  %-17s %s\n' \
    ORG PROYECTO ROL ETIQUETA ADAPT INSERTO READ VEREDICTO ADAPTADOR
  awk -F'\t' '{printf "%-7s %-14s %-10s %-11s %5s%% %8s %7s  %-17s %s\n", \
    $1,$2,$3,$4,$6,$7,$8,$9,$10}' "$tabla"
  echo
  # Filas listas para data/adaptadores.tsv. Existe para no pasar a mano lo que
  # la herramienta ya midio: copiar 19 filas de una tabla formateada es
  # exactamente donde se cuela un error que despues no falla ruidosamente.
  if [[ "$tsv" == "1" ]]; then
    echo
    echo "==================== PARA data/adaptadores.tsv"
    awk -F'\t' -v secs="$SECUENCIAS" -v hoy="$(date -u +%Y-%m-%d)" '
      BEGIN {
        ns = split(secs, S, " ")
        for (i = 1; i <= ns; i++) { j = index(S[i], ":");
          SEC[substr(S[i], 1, j-1)] = substr(S[i], j+1) }
      }
      {
        fam = $10; ver = $9; sec = "-"
        # La tabla para leer dice "22 nt"; un TSV lleva el numero solo.
        ins = $7; sub(/ nt$/, "", ins)
        # Una libreria ya recortada se marca PRE-TRIMMED: yasma la pasa de largo
        # sin llamar a cutadapt.
        if (ver == "YA RECORTADA") sec = "PRE-TRIMMED"
        else if (fam in SEC)       sec = SEC[fam]
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               $1, $2, fam, sec, $6, ins, ver, hoy
      }' "$tabla"
    echo
    echo "Las filas con secuencia '-' NO se pueden recortar: o la familia es 5p"
    echo "(dimero, no read-through) o no esta identificada. trim.sh las rechaza."
  fi

  local malos; malos=$(awk -F'\t' '$9!="PARECE sRNA-seq" && $9!="YA RECORTADA"' "$tabla" | wc -l)
  if [[ "$malos" -gt 0 ]]; then
    echo "$malos proyecto(s) sin veredicto favorable. Revisar antes de alinear:"
    awk -F'\t' -v man="$MANIFEST" '
      BEGIN { while ((getline l < man) > 0) { split(l, f, "\t"); n[f[1] "\t" f[3]]++ } }
      $9 != "PARECE sRNA-seq" && $9 != "YA RECORTADA" {
        printf "  %-7s %-14s %-10s %3d corridas etiquetadas %s -> %s\n",
               $1, $2, $3, n[$1 "\t" $2], $4, $9
      }' "$tabla"
  else
    echo "Los $(wc -l < "$tabla") proyectos son sRNA-seq (PARECE o YA RECORTADA)."
  fi
  rm -f "$tabla"
  [[ "$malos" -eq 0 ]]
}

# Busca proyectos de sRNA-seq de una especie en la ENA, con el MISMO filtro que
# arma el manifiesto. Para cuando un duplicado resulta no ser sRNA-seq y hay que
# reemplazarlo: es el caso de sclsc, cuyo PRJNA985401 es RNA-Seq y quedo afuera
# entero, dejandolo sin set de validacion.
cmd_buscar() {
  command -v curl >/dev/null || die "falta curl"
  local esp="${1:-}"
  [[ -n "$esp" ]] || die "uso: $0 buscar 'Nombre cientifico'"

  local tmp; tmp=$(mktemp)
  curl -sS --fail --max-time 300 --retry 3 --retry-delay 2 -G "$ENA_BUSCAR" \
    --data-urlencode "result=read_run" \
    --data-urlencode "query=scientific_name=\"$esp\" AND library_source=\"$FUENTE_OK\"" \
    --data-urlencode "fields=run_accession,study_accession,read_count,base_count,library_strategy,library_layout,library_source" \
    --data-urlencode "format=tsv" \
    --data-urlencode "limit=0" > "$tmp" \
    || { rm -f "$tmp"; die "la ENA no respondio para '$esp'"; }

  # Proyectos que ya estan en la spec, para no proponer uno que ya se usa.
  local ya; ya=$(mktemp)
  awk -F'\t' '!/^[[:space:]]*#/ && NR>1 && NF>=6 {print $6"\t"$1" "$5}' "$SPEC" > "$ya"

  local acum; acum=$(mktemp)
  local run est rc bc strat layout src n_crudo=0 n_ok=0
  while IFS=$'\t' read -r run est rc bc strat layout src; do
    [[ "$run" == "run_accession" || -z "$run" ]] && continue
    n_crudo=$((n_crudo+1))
    pasa_filtro "$src" "$strat" "$layout" "$rc" || continue
    printf '%s\t%s\t%s\n' "$est" "$rc" "$strat" >> "$acum"
    n_ok=$((n_ok+1))
  done < "$tmp"

  echo "== $esp"
  echo "   $n_crudo corridas TRANSCRIPTOMIC en la ENA; $n_ok pasan el filtro del proyecto"
  echo
  if [[ "$n_ok" -eq 0 ]]; then
    echo "   Ningun proyecto de esta especie tiene datos que entren al manifiesto."
    rm -f "$tmp" "$ya" "$acum"; return 0
  fi

  printf '   %-14s %8s %9s  %-22s %s\n' PROYECTO CORRIDAS 'SPOTS(M)' ESTRATEGIAS 'EN LA SPEC'
  awk -F'\t' -v yaf="$ya" '
    BEGIN { while ((getline l < yaf) > 0) { split(l, a, "\t"); spec[a[1]] = a[2] } }
    { n[$1]++; spots[$1] += $2; if (index(e[$1], $3) == 0) e[$1] = e[$1] (e[$1] ? "," : "") $3 }
    END {
      for (p in n)
        printf "   %-14s %8d %9.0f  %-22s %s\n", p, n[p], spots[p]/1e6, e[p], (p in spec ? spec[p] : "-")
    }' "$acum" | sort -k2,2nr
  echo
  echo "   Un candidato a duplicado tiene que decir '-' en la ultima columna:"
  echo "   un proyecto que ya esta en la spec no valida nada de forma independiente."
  rm -f "$tmp" "$ya" "$acum"
}

# Por que falla una corrida. Necesita red hacia NCBI: corre en Colab o en la
# maquina local, nunca en la sesion cloud (el gateway responde 403 al CONNECT).
# Ver la celda 5 de notebooks/10_descarga_runs.ipynb.
cmd_diag() {
  [[ $# -ge 1 ]] || die "uso: $0 diag RUN [RUN...]"
  local run tmp rc dejo
  for run in "$@"; do
    echo "=============== $run"

    echo "-- srapath: que URL resuelve el resolver de SRA"
    if command -v srapath >/dev/null; then
      srapath "$run" 2>&1 | sed 's/^/   /' || echo "   (srapath salió con error)"
    else
      echo "   (no hay srapath en el PATH)"
    fi

    echo "-- vdb-dump --info: que cree SRA que existe"
    if command -v vdb-dump >/dev/null; then
      vdb-dump --info "$run" 2>&1 | head -20 | sed 's/^/   /' || true
    else
      echo "   (no hay vdb-dump en el PATH)"
    fi

    echo "-- prefetch, con la salida a la vista"
    if ! command -v prefetch >/dev/null; then
      echo "   (no hay prefetch en el PATH)"; echo; continue
    fi
    tmp=$(mktemp -d)
    set +e
    prefetch --output-directory "$tmp" --max-size u "$run" 2>&1 | sed 's/^/   /'
    rc=${PIPESTATUS[0]}
    set -e
    echo "   exit=$rc"

    echo "-- qué quedó en el disco"
    dejo=$(find "$tmp" -type f -printf '   %10s  %P\n' 2>/dev/null || true)
    echo "${dejo:-   (nada)}"
    rm -rf "$tmp"
    echo
  done
}

# Repara el ledger: calcula el md5 de las corridas que estan en DEST pero no
# figuran en data/sra_md5.tsv. Pasa cuando la sesion de Colab que las bajo se
# muere antes de que el ledger llegue a git: el dato esta en Drive y el
# checksum se perdio. Solo toca lo que falta, asi que correrlo de mas no hace
# nada.
#
# El formato (sra|sralite) no se puede deducir del archivo: al moverlo a DEST se
# guarda como <RUN>.sra en los dos casos. Por eso --formato, que aplica a las
# que falten. Default sra.
cmd_ledger() {
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: $0 manifest"
  local filtro="${1:-}" fmt="${2:-sra}"
  local org run src tab n=0 ya=0 sin=0
  tab=$(printf '\t')
  mkdir -p "$(dirname "$LEDGER")"; touch "$LEDGER"
  while IFS=$'\t' read -r org run; do
    if grep -q "^${org}${tab}${run}${tab}" "$LEDGER" 2>/dev/null; then
      ya=$((ya+1)); continue
    fi
    if ! src=$(ruta_sra "$DEST" "$org" "$run"); then
      echo "   FALTA bajar: $org $run" >&2; sin=$((sin+1)); continue
    fi
    printf '== %s %s  (formato=%s, %s MB)\n' "$org" "$run" "$fmt" \
      "$(awk -v b="$(stat -c%s "$src")" 'BEGIN{printf "%.0f", b/1e6}')"
    registrar_md5 "$org" "$run" "$(md5sum "$src" | cut -d' ' -f1)" "$fmt"
    n=$((n+1))
  done < <(awk -F'\t' -v o="$filtro" 'NR>1 && (o=="" || $1==o) {print $1"\t"$2}' "$MANIFEST")

  # El ledger espeja el manifiesto. Una fila de una corrida que ya no esta en el
  # manifiesto —porque se excluyo— haria que la celda de "guardar el ledger"
  # proponga volver a commitearla, deshaciendo la exclusion. Ya paso una vez.
  local sobran; sobran=$(awk -F'\t' -v m="$MANIFEST" '
    BEGIN { while ((getline l < m) > 0) { split(l, f, "\t"); en[f[2]] = 1 } }
    NR > 1 && !($2 in en) { print $2 }' "$LEDGER")
  if [[ -n "$sobran" ]]; then
    local ns; ns=$(wc -l <<< "$sobran")
    echo
    echo "$ns fila(s) del ledger ya no estan en el manifiesto; las saco:"
    sed 's/^/  /' <<< "$sobran"
    local tmp_l; tmp_l=$(mktemp)
    awk -F'\t' -v m="$MANIFEST" '
      BEGIN { while ((getline l < m) > 0) { split(l, f, "\t"); en[f[2]] = 1 } }
      NR == 1 || ($2 in en)' "$LEDGER" > "$tmp_l"
    mv "$tmp_l" "$LEDGER"
  fi

  echo
  echo "agregadas=$n  ya estaban=$ya  sin bajar=$sin"
  [[ $n -gt 0 ]] && echo "Commitear el ledger: git add $LEDGER"
  return 0
}

[[ $# -ge 1 ]] || usage
case "$1" in
  -h|--help) usage 0 ;;
  manifest)  shift; cmd_manifest ;;
  estado)    shift; cmd_estado ;;
  diag)      shift; cmd_diag "$@" ;;
  buscar)    shift; cmd_buscar "${1:-}" ;;
  perfil)
    shift
    RUN_P=""; NP=""; TODOS=0; TSV=0
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --proyectos) TODOS=1; shift ;;
        --tsv) TSV=1; shift ;;
        -n) NP="${2:-}"; shift 2 ;;
        -n*) NP="${1#-n}"; shift ;;
        *) RUN_P="$1"; shift ;;
      esac
    done
    if [[ $TODOS -eq 1 ]]; then cmd_perfil_proyectos "${NP:-20000}" "$TSV"
    else cmd_perfil "$RUN_P" "${NP:-20000}"; fi ;;
  ledger)
    shift
    ORG_F=""; FMT=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --formato) FMT="${2:-}"; shift 2 ;;
        *) ORG_F="$1"; shift ;;
      esac
    done
    [[ -z "$ORG_F" ]] || org_valido_manifest "$ORG_F" || die "organismo desconocido: $ORG_F"
    case "${FMT:-sra}" in
      sra|sralite) ;;
      *) die "formato desconocido: $FMT (sra|sralite)" ;;
    esac
    cmd_ledger "$ORG_F" "${FMT:-sra}" ;;
  prefetch)
    shift
    ORG_F=""; LIMITE=""; HORAS=""; ORDEN=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -n) LIMITE="${2:-}"; shift 2 ;;
        -n*) LIMITE="${1#-n}"; shift ;;
        --horas) HORAS="${2:-}"; shift 2 ;;
        --orden) ORDEN="${2:-}"; shift 2 ;;
        *) ORG_F="$1"; shift ;;
      esac
    done
    [[ -z "$ORG_F" ]] || org_valido_manifest "$ORG_F" || die "organismo desconocido: $ORG_F"
    case "${ORDEN:-alfabetico}" in
      alfabetico|chico|entrenamiento) ;;
      *) die "orden desconocido: $ORDEN (alfabetico|chico|entrenamiento)" ;;
    esac
    cmd_prefetch "$ORG_F" "$LIMITE" "${HORAS:-}" "${ORDEN:-alfabetico}" ;;
  *) echo "comando desconocido: $1" >&2; usage ;;
esac
