#!/usr/bin/env bash
# Resuelve y descarga los ensamblados de referencia. Ver data/genomas.tsv.
#
#   ./scripts/fetch_genomes.sh resolve [ORG]   # consulta NCBI, no descarga
#   ./scripts/fetch_genomes.sh fetch   [ORG]   # descarga los 'verificado'
#   ./scripts/fetch_genomes.sh estado          # qué falta, según la spec
#   ./scripts/fetch_genomes.sh verificar [ORG] [--rapido]
#                                              # qué hay bajado de verdad, e íntegro
#   ./scripts/fetch_genomes.sh cepas ORG [--taxon NOMBRE] [--grep TEXTO]
#                                              # ensamblados de la especie, por cepa
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
# Sobreescribible para poder probar el script contra una spec sintetica, igual
# que MANIFEST en fetch_runs.sh.
SPEC="${GENOMES_SPEC:-$ROOT/data/genomas.tsv}"
DEST="${GENOMES_DIR:-$ROOT/genomes}"
LEDGER="${GENOMES_LEDGER:-$ROOT/data/genomas.sha256}"
API="https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-2}"; }
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
  local esp="$1" tmp="$2" filtro="${3:-}" rc=0 n total tok pag=0 listados
  local todos="$tmp/cepas_todos.json"
  : > "$todos"
  total=""
  tok=""

  # Pagina hasta agotar. Antes pedia page_size=20 de una: para cloro devolvio
  # exactamente 20, o sea que se topo con el limite, y el listado parcial salia
  # impreso como si fuera todo. Un truncamiento que no se anuncia es un dato
  # falso disfrazado de dato. El tope de 20 paginas es defensivo: si la API
  # devolviera siempre el mismo token, esto cortaria en vez de colgarse.
  while :; do
    pag=$((pag+1))
    [[ $pag -gt 20 ]] && { echo "     (corte tras 20 paginas: la API no termina de paginar)"; break; }
    api "$API/genome/taxon/${esp// /%20}/dataset_report?page_size=100${tok:+&page_token=$tok}" \
        "$tmp/s.json" || rc=$?
    if [[ $rc -ne 0 ]]; then
      echo "     (no pude listar los ensamblados: $([[ $rc -eq 1 ]] && echo red || echo "respuesta invalida"))"
      return 0
    fi
    [[ -z "$total" ]] && total=$(jq -r '.total_count // empty' "$tmp/s.json")
    jq -c '.reports[]?' "$tmp/s.json" >> "$todos"
    tok=$(jq -r '.next_page_token // empty' "$tmp/s.json")
    [[ -z "$tok" ]] && break
  done

  n=$(wc -l < "$todos")
  if [[ "$n" -eq 0 ]]; then
    echo "     NCBI no tiene ningun ensamblado para esta especie"
    return 0
  fi

  # El orden es sobre el conjunto completo, no por pagina.
  local orden="$tmp/cepas_orden.txt"
  jq -rs "$JQ_FILA"'
    sort_by(-((.assembly_stats.scaffold_n50 // .assembly_stats.contig_n50 // 0)
              | tonumber? // 0))
    | .[] | fila' "$todos" > "$orden"

  local aviso=""
  [[ -n "$total" && "$total" != "$n" ]] && aviso="  (la API dice $total en total — listado incompleto)"

  if [[ -n "$filtro" ]]; then
    listados=$(grep -ic -- "$filtro" "$orden" || true)
    echo "     $n ensamblados; $listados coinciden con '$filtro'$aviso"
    if [[ "$listados" -eq 0 ]]; then
      echo "     NINGUNO coincide con '$filtro'"
    else
      grep -i -- "$filtro" "$orden" | sed 's/^/     /'
    fi
  else
    echo "     ensamblados disponibles ($n, los mas contiguos primero)$aviso:"
    sed 's/^/     /' "$orden"
  fi
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
        # 'suppressed' es el dato mas decisivo de toda la salida y al final de una
        # linea larga se saltea. NCBI retira un ensamblado por algo: contaminacion,
        # o porque el que lo deposito lo reemplazo. No se usa, y punto.
        if [[ "$(jq -r '.reports[0].assembly_info.assembly_status // ""' "$tmp/c.json")" == "suppressed" ]]; then
          echo "   !!! RETIRADO por NCBI (suppressed) — no usar este accession"
        fi
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
      # Que la especie tenga referencia no quiere decir que la referencia sirva:
      # puede ser de otra cepa que la de los datos. Sin candidato en la spec no
      # hay nada contra que comparar, asi que hace falta ver las cepas igual.
      # Es el caso de cloro: NF-06 es la referencia y los datos son de IK726.
      if [[ -z "$acc" || "$acc" == "?" ]]; then
        echo "   (sin candidato en la spec: la referencia puede ser de otra cepa)"
        listar_cepas "$esp" "$tmp"
      fi
    fi

    # 3. El veredicto. Comparar lo puede hacer la maquina; decidir no, y por eso
    #    el estado lo sigue cambiando una persona.
    if [[ -z "$acc" || "$acc" == "?" ]]; then
      echo "   >>> SIN CANDIDATO — elegi de la lista de arriba, mirando la cepa"
    elif [[ -n "$ref_acc" && "$acc" == "$ref_acc" ]]; then
      echo "   >>> COINCIDE — el candidato ES la referencia vigente"
    elif [[ -n "$ref_acc" ]]; then
      echo "   >>> DIFIERE — la referencia vigente es $ref_acc, no $acc."
      echo "       Gana NCBI, SALVO que el candidato este elegido a proposito —por"
      echo "       ejemplo para que la cepa coincida con la de los datos—. En ese"
      echo "       caso el motivo va escrito en la nota de data/genomas.tsv."
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
  echo "DIFIERE       -> corregi el accession Y el nombre del assembly, despues verifica;"
  echo "                 o dejalo, si elegiste otro a proposito y esta en la nota"
  echo "ERROR DE RED  -> no es un veredicto: volve a correr resolve"
}

# Chequea que los ensamblados esten realmente bajados y sanos. No confia en la
# spec ni en el ledger: mira los archivos. No necesita red, asi que corre igual
# en Colab (GENOMES_DIR al mount de Drive) que en la maquina local.
#
# El criterio fuerte es el sha256 contra data/genomas.sha256: si coincide, el
# archivo es identico byte a byte al que se bajo, y como fetch lo escribio con
# gzip, la integridad del stream va implicita. Solo cuando NO hay entrada en el
# ledger hace falta chequear el gzip aparte, que es lo que detecta un archivo
# cortado (falla el CRC). Por eso no se lee el archivo dos veces.
#
# --rapido: solo existencia y la primera linea. No lee el archivo completo, que
# sobre el FUSE de Drive son ~1 GB.
cmd_verificar() {
  local filtro="${1:-}" rapido="${2:-0}"
  local org esp fuente asm acc estado conf nota
  local f sha_led sha_real pri probs mb
  local n_ok=0 n_mal=0 n_falta=0 n_sin=0

  echo "arbol de genomas: $DEST"
  [[ "$rapido" == "1" ]] && echo "(modo rapido: no se verifican checksums)"
  echo
  printf '%-7s %-11s %-18s %-9s %s\n' ORG ESTADO ACCESSION TAMANO CHEQUEO

  while IFS="$SEP" read -r org esp fuente asm acc estado conf nota; do
    # Sin accession no hay nada que chequear, y eso ES el problema: el proyecto
    # respalda los genomas justamente para poder decir despues contra que se
    # alineo. Los heredados de R1 no tienen ni accession ni FASTA en DEST: su
    # URL vive en config.sh y nada en el repo registra cual fue.
    if [[ -z "$acc" || "$acc" == "?" ]]; then
      printf '%-7s %-11s %-18s %-9s %s\n' "$org" "$estado" "-" "-" \
        "SIN RESPALDO — no hay accession en la spec"
      n_sin=$((n_sin+1)); continue
    fi

    f="$DEST/$org/$acc.fna.gz"
    if [[ ! -s "$f" ]]; then
      printf '%-7s %-11s %-18s %-9s %s\n' "$org" "$estado" "$acc" "-" \
        "FALTA — corre: $0 fetch $org"
      n_falta=$((n_falta+1)); continue
    fi

    mb=$(awk -v b="$(stat -c%s "$f")" 'BEGIN{printf "%.1f MB", b/1e6}')
    probs=""

    # Arranca con '>'? Barato, y descarta que el .gz tenga otra cosa adentro.
    pri=$(gzip -cd "$f" 2>/dev/null | head -1 || true)
    [[ "$pri" == ">"* ]] || probs="${probs}no arranca con '>' (no es FASTA); "

    sha_led=""
    [[ -f "$LEDGER" ]] && sha_led=$(awk -F'\t' -v o="$org" '$1==o {print $4}' "$LEDGER")

    if [[ "$rapido" != "1" ]]; then
      if [[ -n "$sha_led" ]]; then
        sha_real=$(sha256sum "$f" | cut -d" " -f1)
        [[ "$sha_real" == "$sha_led" ]] || probs="${probs}sha256 NO coincide con el ledger; "
      else
        probs="${probs}sin entrada en $(basename "$LEDGER"); "
        gzip -t "$f" 2>/dev/null || probs="${probs}gzip corrupto o truncado; "
      fi
    fi

    [[ -s "$f.sha256" ]] || probs="${probs}falta el .sha256 al lado; "

    if [[ -z "$probs" ]]; then
      printf '%-7s %-11s %-18s %-9s %s\n' "$org" "$estado" "$acc" "$mb" "OK"
      n_ok=$((n_ok+1))
    else
      printf '%-7s %-11s %-18s %-9s %s\n' "$org" "$estado" "$acc" "$mb" "${probs%; }"
      n_mal=$((n_mal+1))
    fi
  done < <(filas "$filtro")

  echo
  echo "ok=$n_ok  con problemas=$n_mal  sin bajar=$n_falta  sin respaldo=$n_sin"
  if [[ $n_sin -gt 0 ]]; then
    echo
    echo "Los 'sin respaldo' son los heredados de R1. No es un fallo de descarga:"
    echo "nunca se registro su accession. Mientras siga asi no se puede decir"
    echo "contra que ensamblado se alineo, que es lo que data/DRIVE.md dice que"
    echo "los genomas se respaldan para poder decir."
  fi
  [[ $n_mal -eq 0 && $n_falta -eq 0 ]] || return 1
  return 0
}

cmd_cepas() {
  local org="${1:-}" taxon="${2:-}" filtro="${3:-}" tmp esp n
  [[ -n "$org" ]] || die "uso: $0 cepas ORG [--taxon NOMBRE] [--grep TEXTO]"
  n=$(filas "$org" | wc -l)
  [[ "$n" -gt 0 ]] || die "organismo desconocido: $org"
  esp=$(filas "$org" | awk -F"$SEP" '{print $2}')
  # --taxon permite probar los sinonimos sin tocar la spec. Hace falta porque
  # NCBI agrupa por el nombre que uso quien deposito: la misma trampa que
  # Magallana / Crassostrea, que ya esta documentada en CLAUDE.md.
  [[ -n "$taxon" ]] && esp="$taxon"
  tmp=$(mktemp -d)
  echo "== $org — $esp${taxon:+  (taxon forzado)}"
  listar_cepas "$esp" "$tmp" "$filtro"
  rm -rf "$tmp"
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
  verificar)
    shift
    ORG_V=""; RAPIDO=0
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --rapido) RAPIDO=1; shift ;;
        *) ORG_V="$1"; shift ;;
      esac
    done
    cmd_verificar "$ORG_V" "$RAPIDO" ;;
  cepas)
    shift
    ORG_C=""; TAXON=""; GREP_C=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --taxon) TAXON="${2:-}"; shift 2 ;;
        --grep)  GREP_C="${2:-}"; shift 2 ;;
        *) ORG_C="$1"; shift ;;
      esac
    done
    cmd_cepas "$ORG_C" "$TAXON" "$GREP_C" ;;
  resolve) shift; cmd_resolve "${1:-}" ;;
  fetch)   shift; cmd_fetch   "${1:-}" ;;
  *) echo "comando desconocido: $1" >&2; usage ;;
esac
