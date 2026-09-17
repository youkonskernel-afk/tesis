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
# Requiere: curl, jq, unzip. La sesión cloud NO tiene red hacia NCBI: esto corre
# en Colab (notebooks/descarga_genomas.ipynb) o en la máquina local.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SPEC="$ROOT/data/genomas.tsv"
DEST="${GENOMES_DIR:-$ROOT/genomes}"
LEDGER="$ROOT/data/genomas.sha256"
API="https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

usage() { sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }
die()   { echo "error: $*" >&2; exit 1; }

for bin in curl jq unzip; do
  command -v "$bin" >/dev/null || die "falta $bin"
done
[[ -f "$SPEC" ]] || die "no existe $SPEC"

# Filas de la spec, sin comentarios ni cabecera. Campos: org especie fuente
# assembly accession estado confianza nota
# OJO: no usar IFS=$'\t' para leer estas filas. El tab es un carácter de
# espacio de IFS, así que bash colapsa tabs consecutivos y las filas con
# accession vacío se leen corridas un campo. Se re-separa con \x1f, que no es
# whitespace y por lo tanto preserva los campos vacíos.
SEP=$'\x1f'
filas() {
  grep -v '^[[:space:]]*#' "$SPEC" | tail -n +2 | awk -F'\t' -v o="${1:-}" -v s="$SEP" \
    'NF>=7 && (o=="" || $1==o) { for(i=1;i<=8;i++) printf "%s%s", (i>1?s:""), $i; print "" }'
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

# Consulta a la API. Separa "NCBI dice que no existe" de "no pude preguntar":
# son cosas distintas y confundirlas es el mismo modo de falla que el prefetch
# que salia 0 sin bajar nada. rc: 0 ok, 1 red, 2 respuesta que no es JSON.
api() {
  local url="$1" dst="$2"
  : > "$dst.err"
  curl -sS --max-time 60 -o "$dst" "$url" 2>"$dst.err" || return 1
  jq -e . "$dst" >/dev/null 2>&1 || return 2
  return 0
}

# Formato compacto de un report de NCBI. Incluye nivel, N50 y tamano porque el
# criterio de eleccion es contiguidad y completitud (ver CLAUDE.md), y eso no se
# puede juzgar solo con el nombre del ensamblado.
JQ_FILA='
  def mb: if . == null or . == "" then "?"
          else ((tonumber? // 0) / 1000000 * 10 | round / 10 | tostring) + " Mb" end;
  def fila:
    "\(.accession)  \(.assembly_info.assembly_name // "?")"
    + "  nivel=\(.assembly_info.assembly_level // "?")"
    + "  N50=\((.assembly_stats.scaffold_n50 // .assembly_stats.contig_n50) | mb)"
    + "  total=\(.assembly_stats.total_sequence_length | mb)"
    + (if (.organism.infraspecific_names.strain // "") != ""
       then "  cepa=\(.organism.infraspecific_names.strain)" else "" end);
'

cmd_estado() {
  printf "%-8s %-11s %-10s %-32s %s\n" ORG ESTADO CONFIANZA ASSEMBLY ACCESSION
  filas "${1:-}" | while IFS="$SEP" read -r org esp fuente asm acc estado conf nota; do
    printf "%-8s %-11s %-10s %-32s %s\n" "$org" "$estado" "$conf" "$asm" "${acc:--}"
  done
  echo
  local n; n=$(filas | awk -F"$SEP" '$6=="candidato"' | wc -l)
  [[ $n -gt 0 ]] && echo "$n sin verificar. Corré: $0 resolve"
  return 0
}

# Cuando la especie no tiene ensamblado de referencia —el caso de cloro— hay que
# elegir por cepa. Sin esta lista, resolve dejaba un callejon sin salida: decia
# "buscar por cepa" y no daba con que buscarla.
listar_cepas() {
  local esp="$1" tmp="$2" rc=0 n
  api "$API/genome/taxon/${esp// /%20}/dataset_report?page_size=20" "$tmp/s.json" || rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "     (no pude listar los ensamblados: $([[ $rc -eq 1 ]] && echo red || echo "respuesta invalida"))"
    return 0
  fi
  n=$(jq '.reports | length' "$tmp/s.json")
  if [[ "$n" -eq 0 ]]; then
    echo "     NCBI no tiene ningun ensamblado para esta especie"
    return 0
  fi
  echo "     ensamblados disponibles ($n, los mas contiguos primero):"
  jq -r "$JQ_FILA"'
    [.reports[]]
    | sort_by(-((.assembly_stats.scaffold_n50 // .assembly_stats.contig_n50 // 0)
                | tonumber? // 0))
    | .[] | "     " + fila' "$tmp/s.json"
}

cmd_resolve() {
  local tmp; tmp=$(mktemp -d)
  local org esp fuente asm acc estado conf nota rc cand_ok ref_acc
  while IFS="$SEP" read -r org esp fuente asm acc estado conf nota; do
    [[ "$estado" == "heredado" ]] && { echo "== $org: heredado, se salta"; continue; }
    echo "== $org — $esp"
    echo "   spec       : ${acc:-(sin accession)}  $asm  confianza=$conf"

    # 1. El accession propuesto, existe?
    cand_ok=0
    if [[ -n "$acc" && "$acc" != "?" ]]; then
      rc=0; api "$API/genome/accession/$acc/dataset_report" "$tmp/c.json" || rc=$?
      if [[ $rc -eq 1 ]]; then
        echo "   candidato  : ERROR DE RED — $(head -c 120 "$tmp/c.json.err" | tr '\n' ' ')"
      elif [[ $rc -eq 2 ]]; then
        echo "   candidato  : ERROR — NCBI no devolvio JSON"
      elif [[ "$(jq '.reports | length' "$tmp/c.json")" -eq 0 ]]; then
        echo "   candidato  : NO EXISTE en NCBI"
      else
        cand_ok=1
        jq -r "$JQ_FILA"' .reports[0] | "   candidato  : " + fila
               + "  estado=\(.assembly_info.assembly_status // "?")"' "$tmp/c.json"
      fi
    fi

    # 2. La referencia vigente de la especie. Es la respuesta que mas importa:
    #    un accession puede existir y no ser el que corresponde.
    ref_acc=""; rc=0
    api "$API/genome/taxon/${esp// /%20}/dataset_report?filters.reference_only=true&page_size=3" \
        "$tmp/r.json" || rc=$?
    if [[ $rc -eq 1 ]]; then
      echo "   referencia : ERROR DE RED — $(head -c 120 "$tmp/r.json.err" | tr '\n' ' ')"
    elif [[ $rc -eq 2 ]]; then
      echo "   referencia : ERROR — NCBI no devolvio JSON"
    elif [[ "$(jq '.reports | length' "$tmp/r.json")" -eq 0 ]]; then
      echo "   referencia : la especie NO tiene ensamblado de referencia"
      listar_cepas "$esp" "$tmp"
    else
      ref_acc=$(jq -r '.reports[0].accession' "$tmp/r.json")
      jq -r "$JQ_FILA"' .reports[] | "   referencia : " + fila' "$tmp/r.json"
    fi

    # 3. El veredicto. Comparar lo puede hacer la maquina; decidir no, y por eso
    #    el estado lo sigue cambiando una persona.
    if [[ -z "$acc" || "$acc" == "?" ]]; then
      echo "   >>> SIN CANDIDATO — hay que elegir uno de la lista de arriba"
    elif [[ -n "$ref_acc" && "$acc" == "$ref_acc" ]]; then
      echo "   >>> COINCIDE — el candidato ES la referencia vigente"
    elif [[ -n "$ref_acc" ]]; then
      echo "   >>> DIFIERE — la referencia vigente es $ref_acc, no $acc. Gana NCBI."
    elif [[ $cand_ok -eq 1 ]]; then
      echo "   >>> El candidato existe pero la especie no tiene referencia vigente."
      echo "       Decidi a mano con los numeros de arriba."
    else
      echo "   >>> SIN RESPUESTA UTIL — no confirmes nada con esto"
    fi
    echo
  done < <(filas "${1:-}")
  rm -rf "$tmp"
  echo "COINCIDE      -> pone 'verificado' en data/genomas.tsv y corre: $0 fetch"
  echo "DIFIERE       -> corregi el accession Y el nombre del assembly, despues verifica"
  echo "ERROR DE RED  -> no es un veredicto: volve a correr resolve"
}

cmd_fetch() {
  local bajados=0 saltados=0
  while IFS="$SEP" read -r org esp fuente asm acc estado conf nota; do
    case "$estado" in
      verificado) ;;
      candidato)
        echo "SALTO $org: 'candidato' (confianza $conf). Verificá: $0 resolve $org" >&2
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
