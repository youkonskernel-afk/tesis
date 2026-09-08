#!/usr/bin/env bash
# Resuelve y descarga los ensamblados de referencia. Ver data/genomas.tsv.
#
#   ./scripts/fetch_genomes.sh resolve [ORG]   # consulta NCBI, no descarga
#   ./scripts/fetch_genomes.sh fetch   [ORG]   # descarga los 'verificado'
#   ./scripts/fetch_genomes.sh estado          # qué falta
#
# Flujo: 'resolve' te dice qué ensamblado es el vigente y si el candidato de
# data/genomas.tsv coincide. Confirmás a mano, cambiás estado a 'verificado',
# y recién ahí 'fetch' lo baja. Un ensamblado equivocado no falla ruidosamente:
# alinea peor y contamina la anotación, así que el paso manual es a propósito.
#
# Para subir a Drive después:  ./scripts/drive_push.sh genomas --go
#
# Requiere: curl, jq. La sesión cloud NO tiene red hacia NCBI; esto corre local.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SPEC="$ROOT/data/genomas.tsv"
DEST="${GENOMES_DIR:-$ROOT/genomes}"
LEDGER="$ROOT/data/genomas.sha256"
API="https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

usage() { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }
die()   { echo "error: $*" >&2; exit 1; }

for bin in curl jq; do
  command -v "$bin" >/dev/null || die "falta $bin"
done
[[ -f "$SPEC" ]] || die "no existe $SPEC"

# Filas de la spec, sin comentarios ni cabecera. Campos: org especie fuente
# assembly accession estado nota
# OJO: no usar IFS=$'\t' para leer estas filas. El tab es un carácter de
# espacio de IFS, así que bash colapsa tabs consecutivos y las filas con
# accession vacío se leen corridas un campo. Se re-separa con \x1f, que no es
# whitespace y por lo tanto preserva los campos vacíos.
SEP=$'\x1f'
filas() {
  grep -v '^[[:space:]]*#' "$SPEC" | tail -n +2 | awk -F'\t' -v o="${1:-}" -v s="$SEP" \
    'NF>=6 && (o=="" || $1==o) { for(i=1;i<=7;i++) printf "%s%s", (i>1?s:""), $i; print "" }'
}

# El checksum también va a git. Guardarlo solo junto al FASTA en Drive no prueba
# nada: quien reemplace el genoma reemplaza el checksum con él. El registro
# versionado es lo que permite detectar después que el ensamblado cambió.
registrar() {
  local org="$1" acc="$2" asm="$3" sum="$4" tab
  tab=$(printf '\t')
  touch "$LEDGER"
  {
    grep -v "^org${tab}" "$LEDGER" 2>/dev/null | grep -v "^${org}${tab}" || true
    printf '%s\t%s\t%s\t%s\t%s\n' "$org" "$acc" "$asm" "$sum" "$(date -u +%Y-%m-%d)"
  } | sort > "$LEDGER.tmp"
  {
    printf 'org\taccession\tassembly\tsha256\tfecha_utc\n'
    cat "$LEDGER.tmp"
  } > "$LEDGER"
  rm -f "$LEDGER.tmp"
}

cmd_estado() {
  printf "%-8s %-12s %-32s %s\n" ORG ESTADO ASSEMBLY ACCESSION
  filas "${1:-}" | while IFS="$SEP" read -r org esp fuente asm acc estado nota; do
    printf "%-8s %-12s %-32s %s\n" "$org" "$estado" "$asm" "${acc:--}"
  done
  echo
  local n; n=$(filas | awk -F'\t' '$6=="candidato"' | wc -l)
  [[ $n -gt 0 ]] && echo "$n sin verificar. Corré: $0 resolve"
  return 0
}

cmd_resolve() {
  filas "${1:-}" | while IFS="$SEP" read -r org esp fuente asm acc estado nota; do
    [[ "$estado" == "heredado" ]] && { echo "== $org: heredado, se salta"; continue; }
    echo "== $org — $esp"

    if [[ -n "$acc" ]]; then
      echo "   candidato en la spec: $acc ($asm)"
      local_json=$(curl -sS --max-time 60 \
        "$API/genome/accession/$acc/dataset_report" 2>/dev/null || echo '{}')
      echo "$local_json" | jq -r '
        if (.reports|length) > 0 then
          .reports[0] |
          "   NCBI dice: \(.accession)  \(.assembly_info.assembly_name)  " +
          "\(.assembly_info.assembly_status // "?")  " +
          "org=\(.organism.organism_name)"
        else "   NCBI: accession NO encontrado — el candidato es incorrecto"
        end'
    fi

    # Cuál es el ensamblado de referencia vigente para la especie, sea cual sea
    # el candidato. Es la respuesta que realmente importa.
    esp_url=${esp// /%20}
    curl -sS --max-time 60 \
      "$API/genome/taxon/$esp_url/dataset_report?filters.reference_only=true&page_size=3" \
      2>/dev/null | jq -r '
        if (.reports|length) > 0 then
          .reports[] |
          "   referencia vigente: \(.accession)  \(.assembly_info.assembly_name)  " +
          "nivel=\(.assembly_info.assembly_level)  org=\(.organism.organism_name)"
        else "   sin ensamblado de referencia para esta especie — buscar por cepa"
        end'
    echo
  done
  echo "Si coincide, cambiá estado a 'verificado' en data/genomas.tsv y corré: $0 fetch"
}

cmd_fetch() {
  local bajados=0 saltados=0
  while IFS="$SEP" read -r org esp fuente asm acc estado nota; do
    case "$estado" in
      verificado) ;;
      candidato)
        echo "SALTO $org: estado 'candidato'. Verificá con: $0 resolve $org" >&2
        saltados=$((saltados+1)); continue ;;
      heredado)
        echo "SALTO $org: heredado, se baja con el pipeline viejo (config.sh)" >&2
        saltados=$((saltados+1)); continue ;;
      *) echo "SALTO $org: estado desconocido '$estado'" >&2
         saltados=$((saltados+1)); continue ;;
    esac
    [[ -n "$acc" ]] || { echo "SALTO $org: verificado pero sin accession" >&2; continue; }

    out="$DEST/$org"
    if [[ -f "$out/$acc.fna.gz" ]]; then
      echo "== $org: ya está en $out/$acc.fna.gz"; continue
    fi
    mkdir -p "$out"
    echo "== $org: bajando $acc ($asm)"
    tmp="$out/.$acc.zip.partial"
    curl -sS --fail --max-time 3600 -o "$tmp" \
      "$API/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" \
      || { rm -f "$tmp"; die "falló la descarga de $org ($acc)"; }
    # El zip de datasets trae ncbi_dataset/data/<acc>/<acc>_<asm>_genomic.fna
    unzip -p "$tmp" "ncbi_dataset/data/$acc/"'*_genomic.fna' | gzip -c > "$out/$acc.fna.gz" \
      || { rm -f "$tmp" "$out/$acc.fna.gz"; die "no pude extraer el FASTA de $org"; }
    rm -f "$tmp"
    ( cd "$out" && sha256sum "$acc.fna.gz" > "$acc.fna.gz.sha256" )
    registrar "$org" "$acc" "$asm" "$(cut -d" " -f1 < "$out/$acc.fna.gz.sha256")"
    echo "   $(du -h "$out/$acc.fna.gz" | cut -f1)  sha256 -> data/genomas.sha256"
    bajados=$((bajados+1))
  done < <(filas "${1:-}")

  echo
  echo "bajados=$bajados saltados=$saltados"
  if [[ $bajados -gt 0 ]]; then
    echo "Subir a Drive:  ./scripts/drive_push.sh genomas --go"
    echo "Y commitear:    git add data/genomas.sha256"
  fi
  return 0
}

[[ $# -ge 1 ]] || usage
case "$1" in
  -h|--help) usage 0 ;;
  estado)  shift; cmd_estado  "${1:-}" ;;
  resolve) shift; cmd_resolve "${1:-}" ;;
  fetch)   shift; cmd_fetch   "${1:-}" ;;
  *) echo "comando desconocido: $1" >&2; usage ;;
esac
