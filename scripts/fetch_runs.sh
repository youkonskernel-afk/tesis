#!/usr/bin/env bash
# Resuelve los BioProjects de data/organismos.tsv a corridas y las descarga.
# Reemplaza al gen_manifest.sh de R1, que asumia un proyecto por organismo.
#
#   ./scripts/fetch_runs.sh manifest          # consulta ENA -> data/srr_manifest.tsv
#   ./scripts/fetch_runs.sh estado            # que hay bajado y que falta
#   ./scripts/fetch_runs.sh prefetch [ORG]    # descarga los .sra al cache
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
usage() { sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }
die() { echo "error: $*" >&2; exit 1; }

# Campos de organismos.tsv: org especie reino clado rol bioproject runs
# spots_M strategy assembly nota set_modelo
proyectos() {
  grep -v '^[[:space:]]*#' "$SPEC" | tail -n +2 | awk -F'\t' -v s="$SEP" \
    'NF>=12 { printf "%s%s%s%s%s%s%s\n", $1, s, $6, s, $5, s, $12 }'
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

cmd_estado() {
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: $0 manifest"
  local n; n=$(( $(wc -l < "$MANIFEST") - 1 ))
  echo "manifiesto: $n corridas"
  local hay=0 falta=0
  while IFS=$'\t' read -r org run resto; do
    [[ "$org" == "org" ]] && continue
    if [[ -s "$CACHE/$run/$run.sra" || -s "$CACHE/$run.sra" ]]; then
      hay=$((hay+1)); else falta=$((falta+1)); fi
  done < "$MANIFEST"
  echo "en cache ($CACHE): $hay    faltan: $falta"
}

cmd_prefetch() {
  command -v prefetch >/dev/null || die "falta prefetch (sra-tools). Activá el entorno srna2."
  [[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: $0 manifest"
  local filtro="${1:-}"
  mkdir -p "$CACHE"
  # Orden por profundidad descendente (LPT): las corridas grandes primero para
  # que la cola no termine con un job largo solo al final.
  awk -F'\t' -v o="$filtro" 'NR>1 && (o=="" || $1==o) {print $6"\t"$1"\t"$2}' "$MANIFEST" \
  | sort -k1,1nr \
  | while IFS=$'\t' read -r rc org run; do
      if [[ -s "$CACHE/$run/$run.sra" || -s "$CACHE/$run.sra" ]]; then
        echo "== $org $run: ya está"; continue
      fi
      echo "== $org $run ($rc reads)"
      prefetch --output-directory "$CACHE" --max-size u "$run" \
        || echo "FALLO $run — seguir y reintentar después" >&2
    done
  echo; cmd_estado
}

[[ $# -ge 1 ]] || usage
case "$1" in
  -h|--help) usage 0 ;;
  manifest)  shift; cmd_manifest ;;
  estado)    shift; cmd_estado ;;
  prefetch)  shift; cmd_prefetch "${1:-}" ;;
  *) echo "comando desconocido: $1" >&2; usage ;;
esac
