#!/usr/bin/env bash
# Resuelve los BioProjects de data/organismos.tsv a corridas y las descarga.
# Reemplaza al gen_manifest.sh de R1, que asumia un proyecto por organismo.
#
#   ./scripts/fetch_runs.sh manifest             # consulta ENA -> data/srr_manifest.tsv
#   ./scripts/fetch_runs.sh estado               # que hay bajado y que falta
#   ./scripts/fetch_runs.sh prefetch [ORG] [-n N]  # descarga los .sra
#   ./scripts/fetch_runs.sh diag RUN [RUN...]    # por que falla una corrida
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
SPEC="$ROOT/data/organismos.tsv"
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

# Criterio de seleccion. Cambiarlo cambia el dataset: dejarlo explicito y
# versionado es lo que hace el manifiesto reproducible.
MIN_READS="${MIN_READS:-3000000}"
# Solo datos de RNA. library_source=TRANSCRIPTOMIC es el filtro duro: excluye
# corridas GENOMIC que varios BioProjects mezclan con las de RNA.
FUENTE_OK="TRANSCRIPTOMIC"
# RNA-Seq solo SINGLE: el PAIRED de un proyecto de RNA-Seq no es sRNA-seq y
# contaminaria la anotacion. miRNA-Seq y ncRNA-Seq entran con cualquier layout.
ESTRATEGIAS="miRNA-Seq ncRNA-Seq RNA-Seq"

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

  local total=0 desc=0
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

      [[ "$src" == "$FUENTE_OK" ]] || continue
      grep -qw -- "$strat" <<<"$ESTRATEGIAS" || continue
      [[ "$strat" == "RNA-Seq" && "$layout" != "SINGLE" ]] && continue
      [[ -n "$rc" && "$rc" =~ ^[0-9]+$ && "$rc" -ge "$MIN_READS" ]] || continue

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
  echo
  echo "manifiesto: $MANIFEST"
  echo "corridas vistas=$total  descartadas=$desc  retenidas=$(( $(wc -l < "$MANIFEST") - 1 ))"
  echo
  echo "Reads largos (avg_len > 50): necesitan pre-trim antes de fastp."
  awk -F'\t' 'NR>1 && $8!="NA" && $8+0>50 {print "  " $1" "$2"  "$8" nt  ("$3")"}' "$MANIFEST" | head -20
  echo
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
  awk -F'\t' -v dest="$DEST" '
    NR>1 { total[$1]++ }
    END { }' "$MANIFEST" >/dev/null

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

[[ $# -ge 1 ]] || usage
case "$1" in
  -h|--help) usage 0 ;;
  manifest)  shift; cmd_manifest ;;
  estado)    shift; cmd_estado ;;
  diag)      shift; cmd_diag "$@" ;;
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
