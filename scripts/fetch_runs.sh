#!/usr/bin/env bash
# Resuelve los BioProjects de data/organismos.tsv a corridas y las descarga.
# Reemplaza al gen_manifest.sh de R1, que asumia un proyecto por organismo.
#
#   ./scripts/fetch_runs.sh manifest             # consulta ENA -> data/srr_manifest.tsv
#   ./scripts/fetch_runs.sh estado               # que hay bajado y que falta
#   ./scripts/fetch_runs.sh prefetch [ORG] [-n N]  # descarga los .sra
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
usage() { sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }
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

registrar_md5() {
  local org="$1" run="$2" sum="$3" tab
  tab=$(printf '\t')
  mkdir -p "$(dirname "$LEDGER")"
  touch "$LEDGER"
  {
    grep -v "^org${tab}" "$LEDGER" 2>/dev/null | grep -v "^${org}${tab}${run}${tab}" || true
    printf '%s\t%s\t%s\t%s\n' "$org" "$run" "$sum" "$(date -u +%Y-%m-%d)"
  } | sort > "$LEDGER.tmp"
  { printf 'org\trun\tmd5\tfecha_utc\n'; cat "$LEDGER.tmp"; } > "$LEDGER"
  rm -f "$LEDGER.tmp"
}

cmd_prefetch() {
  command -v prefetch >/dev/null || die "falta prefetch (sra-tools). Activá el entorno srna2."
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: $0 manifest"
  local filtro="$1" limite="$2"
  mkdir -p "$DEST" "$STAGING"

  local n=0 ok=0 fallos=0
  # Orden: primero los organismos con menos faltantes, para completar organismos
  # enteros antes de empezar otros. Un organismo completo se puede alinear; uno
  # a medias no sirve para nada. Dentro del organismo, los grandes primero.
  while IFS=$'\t' read -r org run rc; do
    [[ -n "$limite" && $n -ge $limite ]] && { echo; echo "corte por -n $limite"; break; }

    if ruta_sra "$DEST" "$org" "$run" >/dev/null; then continue; fi

    local libre; libre=$(libre_gb "$STAGING")
    if [[ -n "$libre" && "$libre" -lt "$MIN_LIBRE_GB" ]]; then
      echo "ALTO: quedan ${libre}G libres en $STAGING (mínimo $MIN_LIBRE_GB). Cortando." >&2
      break
    fi

    n=$((n+1))
    echo "== [$n] $org $run ($rc reads)  libre=${libre:-?}G"
    rm -rf "${STAGING:?}/$run"
    if ! prefetch --output-directory "$STAGING" --max-size u "$run"; then
      echo "   FALLO prefetch $run" >&2; fallos=$((fallos+1)); rm -rf "${STAGING:?}/$run"; continue
    fi

    local src; src=$(ruta_sra "$STAGING" "$org" "$run") || {
      echo "   FALLO: prefetch no dejó .sra para $run" >&2; fallos=$((fallos+1)); continue; }

    # Un .sra truncado NO falla ruidosamente: alinea de menos. Validar antes de
    # darlo por bueno, y nunca mover al destino uno que no valida.
    if command -v vdb-validate >/dev/null; then
      if ! vdb-validate "$src" >/dev/null 2>&1; then
        echo "   FALLO vdb-validate $run — descartado, se reintenta después" >&2
        rm -rf "${STAGING:?}/$run" "$src"; fallos=$((fallos+1)); continue
      fi
    else
      echo "   aviso: sin vdb-validate, no puedo verificar integridad" >&2
    fi

    mkdir -p "$DEST/$org"
    local final="$DEST/$org/$run.sra"
    if [[ "$src" != "$final" ]]; then
      mv "$src" "$final" || { echo "   FALLO al mover $run" >&2; fallos=$((fallos+1)); continue; }
      rm -rf "${STAGING:?}/$run"
    fi

    registrar_md5 "$org" "$run" "$(md5sum "$final" | cut -d' ' -f1)"
    echo "   ok  $(du -h "$final" | cut -f1)"
    ok=$((ok+1))
  done < <(awk -F'\t' -v o="$filtro" 'NR>1 && (o=="" || $1==o) {print $1"\t"$2"\t"$6}' "$MANIFEST" \
           | sort -t$'\t' -k1,1 -k3,3nr)

  echo
  echo "bajadas=$ok fallos=$fallos"
  [[ $ok -gt 0 ]] && echo "Commitear el ledger: git add $LEDGER"
  cmd_estado
}

[[ $# -ge 1 ]] || usage
case "$1" in
  -h|--help) usage 0 ;;
  manifest)  shift; cmd_manifest ;;
  estado)    shift; cmd_estado ;;
  prefetch)
    shift
    ORG_F=""; LIMITE=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -n) LIMITE="${2:-}"; shift 2 ;;
        -n*) LIMITE="${1#-n}"; shift ;;
        *) ORG_F="$1"; shift ;;
      esac
    done
    [[ -z "$ORG_F" ]] || org_valido_manifest "$ORG_F" || die "organismo desconocido: $ORG_F"
    cmd_prefetch "$ORG_F" "$LIMITE" ;;
  *) echo "comando desconocido: $1" >&2; usage ;;
esac
