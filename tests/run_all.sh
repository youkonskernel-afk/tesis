#!/usr/bin/env bash
# Corre todos los bancos y suma. Sin red: cada banco pone binarios falsos en el
# PATH (prefetch, curl, git, fastq-dump) para ejercitar el codigo de verdad.
#
#     tests/run_all.sh          # solo el resumen
#     tests/run_all.sh -v       # con la salida de cada banco
#
# Por que existen: estos bancos encontraron nueve bugs que ninguna lectura del
# codigo habia visto —la paginacion que cortaba en 20 de 78 ensamblados, la cepa
# que no se leia del biosample, el `grep -v` que mataba `manifest` bajo
# pipefail, el doble conteo de descartadas, los tres defectos del veredicto de
# `perfil`, el ledger que no espejaba el manifiesto—. Correrlos antes de cada
# push es mas barato que cualquiera de esos diagnosticos.
set -uo pipefail
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$RAIZ"

VERBOSE=0
[[ "${1:-}" == "-v" || "${1:-}" == "--verbose" ]] && VERBOSE=1

BANCOS=(
  tests/test_manifest.sh
  tests/test_buscar.sh
  tests/test_perfil.sh
  tests/test_prefetch.sh
  tests/test_estado.sh
  tests/test_trim.sh
  tests/test_align.sh
  tests/test_rutas.sh
  tests/test_drive.sh
  tests/test_drive_check.sh
  tests/test_genomes.sh
  tests/test_cepas.sh
  tests/test_verificar.sh
  tests/test_clon.py
  tests/test_colab_git.py
  tests/test_celda1.py
  tests/test_celda4.py
  tests/test_celda_alinear.py
  tests/test_celda_trim.py
)

malos=0 total_ok=0 total_mal=0
for b in "${BANCOS[@]}"; do
  if [[ ! -x "$b" ]]; then
    printf '[MAL] %-24s no existe o no es ejecutable\n' "$(basename "$b")"
    malos=$((malos + 1)); continue
  fi

  salida=$("$b" 2>&1); rc=$?
  # Cada banco imprime '  ok   <que>' y '  MAL  <que>'. Se cuentan aca para que
  # el resumen diga cuantos chequeos corrieron, no solo cuantos bancos.
  n_ok=$(grep -c '^  ok  ' <<<"$salida")
  n_mal=$(grep -c '^  MAL ' <<<"$salida")
  total_ok=$((total_ok + n_ok)); total_mal=$((total_mal + n_mal))

  if [[ $rc -eq 0 && $n_mal -eq 0 ]]; then
    printf '[OK ] %-24s %3d chequeos\n' "$(basename "$b")" "$n_ok"
    [[ $VERBOSE -eq 1 ]] && sed 's/^/       /' <<<"$salida"
  else
    malos=$((malos + 1))
    printf '[MAL] %-24s %3d ok, %d fallas (rc=%d)\n' \
      "$(basename "$b")" "$n_ok" "$n_mal" "$rc"
    # Las fallas se muestran siempre: un banco rojo sin decir que fallo no
    # sirve de nada.
    grep '^  MAL ' <<<"$salida" | sed 's/^/       /'
    [[ $VERBOSE -eq 1 ]] && sed 's/^/       /' <<<"$salida"
  fi
done

echo
echo "${#BANCOS[@]} bancos, $((total_ok + total_mal)) chequeos: $total_ok ok, $total_mal fallas"
[[ $malos -eq 0 ]] || echo "$malos banco(s) en rojo"
exit $(( malos == 0 ? 0 : 1 ))
