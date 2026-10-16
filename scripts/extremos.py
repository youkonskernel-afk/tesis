#!/usr/bin/env python3
"""Cuántos nt sacar de cada punta del inserto: data/extremos.tsv.

    ./scripts/extremos.py requeridos <org> <rol>           run, corte 5', corte 3'
    ./scripts/extremos.py desactualizadas <proy_dir> <org>  lo recortado con otro corte
    (--tabla y --manifest para no usar los de data/)

POR QUÉ EXISTE

rhirr_duplicado (PRJNA722321) alineó 0.17% contra un genoma verificado y con
un recorte que retuvo lo esperado. `de_quien.py` lo resolvió con el índice:
sacando 4 nt del 5' y 4 del 3' la misma muestra alinea 82.07%. Es la firma de un
kit con adaptadores 4N (NEXTflex y los "HD"): bases aleatorias en las dos
uniones, que quedan pegadas al inserto después de sacar el adaptador, y con
`bowtie -v 1` un read con 4 bases inventadas no alinea.

`yasma trim` no puede sacarlas: le pasa a cutadapt `-a`, la ventana y nada más.
Así que el corte es un dato aparte, por proyecto o por corrida igual que el
adaptador, y `trim.sh` lo aplica. Está en su propio fichero y no como columna de
`adaptadores.tsv` porque lo mide otra herramienta: el adaptador lo decide
`perfil`, mirando dónde empieza; el corte lo decide `de_quien.py`, alineando una
muestra con y sin cortar. Cada tabla tiene un solo dueño.

Un proyecto que no está acá no se corta. Al revés que con el adaptador —donde
faltar la fila hace fallar el recorte— acá faltar es lo normal, y si un
proyecto 4N se escapa no pasa en silencio: `align.sh verificar` lo marca
`MUY BAJA` y corre `de_quien.py` solo.

LA BÚSQUEDA

Igual que en `adaptadores.tsv`: una fila con el RUN exacto le gana a la fila
`-` del proyecto. Un BioProject puede mezclar kits —gadmo PRJNA328800 mezcló
corridas recortadas y sin recortar— y esto es lo mismo un nivel más abajo.

QUÉ SE HACE CON ESTO (en trim.sh)

    cutadapt -u <5'> -a <adaptador> -O 4 --max-n 0 --trimmed-only  (stdout)
      | cutadapt -u -<3'> --minimum-length 15 --maximum-length 50

Dos pasadas porque cutadapt aplica `-u` ANTES de buscar el adaptador: el 5'
se puede cortar ahí, pero las 4 bases del 3' están entre el inserto y el
adaptador, y solo quedan en la punta después de sacarlo. La ventana va en la
segunda pasada, sobre el inserto ya cortado: aplicarla antes dejaría afuera
insertos de 43-50 nt y adentro los de 7-14. Medido con cutadapt 5.2 sobre reads
sintéticos 4N: los 2233 insertos esperados, exactos, ni uno de más.
"""
import argparse
import csv
import os
import pathlib
import sys
import time

RAIZ = pathlib.Path(__file__).resolve().parent.parent
TABLA = pathlib.Path(os.environ.get("EXTREMOS_TSV") or RAIZ / "data" / "extremos.tsv")
MANIFEST = pathlib.Path(os.environ.get("MANIFEST") or RAIZ / "data" / "srr_manifest.tsv")
COLS = ["org", "bioproject", "run", "cortar_5p", "cortar_3p", "alin_sin_pct",
        "alin_con_pct", "lecturas", "fuente", "fecha_utc"]
# Más que esto no es un kit 4N, es otra cosa: mejor que falle a que corte a ciegas.
MAX_CORTE = 10

CABECERA = """\
# Cuantos nt sacar de cada punta del inserto despues de recortar el adaptador.
# Lo escribe `de_quien.py --registrar` (via `align.sh prueba`), no se edita a mano.
# Ver el docstring de scripts/extremos.py.
#
# COLUMNAS
#   org, bioproject, run   a que aplica; run `-` es todo el proyecto, y una fila
#                          con el RUN exacto le gana (un BioProject puede
#                          mezclar kits)
#   cortar_5p, cortar_3p   nt a sacar de cada punta (las bases aleatorias de un
#                          kit 4N)
#   alin_sin_pct           % de la muestra que alinea con bowtie -v 1 sin cortar
#   alin_con_pct           y cortando lo de las dos columnas anteriores: es la
#                          evidencia, y lo que justifica la fila
#   lecturas               cuantas lecturas de la muestra se probaron
#   fuente                 que lo midio
#   fecha_utc              cuando
"""


def _corte(v, fila):
    try:
        n = int(v)
    except (TypeError, ValueError):
        raise ValueError(f"corte no entero ({v!r}) en {fila}") from None
    if not 0 <= n <= MAX_CORTE:
        raise ValueError(f"corte fuera de 0-{MAX_CORTE} ({n}) en {fila}")
    return n


def leer(tabla=None):
    """-> lista de filas (dict). Sin fichero, lista vacía: no cortar es lo normal."""
    tabla = pathlib.Path(tabla or TABLA)
    if not tabla.is_file():
        return []
    lineas = [l for l in tabla.read_text().splitlines() if l.strip() and not l.startswith("#")]
    filas = list(csv.DictReader(lineas, delimiter="\t"))
    for f in filas:
        clave = f"{f.get('org')} {f.get('bioproject')} {f.get('run')}"
        f["cortar_5p"] = _corte(f.get("cortar_5p"), clave)
        f["cortar_3p"] = _corte(f.get("cortar_3p"), clave)
    return filas


def requerido(filas, org, proy, run):
    """(5', 3') para esa corrida: su fila si la tiene, si no la del proyecto."""
    gen = None
    for f in filas:
        if f["org"] != org or f["bioproject"] != proy:
            continue
        if f["run"] == run:
            return f["cortar_5p"], f["cortar_3p"]
        if f["run"] == "-":
            gen = (f["cortar_5p"], f["cortar_3p"])
    return gen or (0, 0)


def etiqueta(c5, c3):
    """Lo que va en la columna `extremos` del ledger del recorte."""
    return "-" if (c5, c3) == (0, 0) else f"{c5}/{c3}"


def requeridos(filas, org, rol, manifest=None):
    """[(run, 5', 3')] de cada corrida del proyecto, en el orden del manifiesto."""
    with open(manifest or MANIFEST) as h:
        man = list(csv.DictReader(h, delimiter="\t"))
    return [(r["run"], *requerido(filas, org, r["bioproject"], r["run"]))
            for r in man if r["org"] == org and r["rol"] == rol]


def aplicados(pdir):
    """{run: (bioproject, etiqueta)} del ledger del recorte. Un ledger de antes de
    esta columna no cortó nada, así que su fila vale '-'."""
    led = pathlib.Path(pdir) / "recortadas.tsv"
    if not led.is_file():
        return {}
    out = {}
    for ln in led.read_text().splitlines()[1:]:
        c = ln.split("\t")
        if len(c) >= 7:
            out[c[0]] = (c[1], c[7] if len(c) >= 8 and c[7] else "-")
    return out


def desactualizadas(pdir, org, filas):
    """[(run, aplicado, requerido)] de lo recortado con un corte distinto del que
    pide la tabla hoy. Es lo que hay que rehacer: la fila llegó después del
    recorte, o cambió."""
    out = []
    for run, (proy, ap) in aplicados(pdir).items():
        req = etiqueta(*requerido(filas, org, proy, run))
        if ap != req:
            out.append((run, ap, req))
    return out


def registrar(nuevas, borrar=(), tabla=None):
    """Escribe filas en la tabla.

    nuevas: dicts con COLS. Reemplazan la fila de su misma (org, bioproject, run).
    borrar: (org, bioproject) cuya medición cubrió TODAS las corridas: sus filas
        viejas se sacan antes, también las por corrida, que si no le seguirían
        ganando a la nueva fila del proyecto."""
    tabla = pathlib.Path(tabla or TABLA)
    previas = leer(tabla)
    if tabla.is_file():
        cab = "".join(l + "\n" for l in tabla.read_text().splitlines() if l.startswith("#"))
    else:
        cab = CABECERA
    borrar = set(borrar)
    claves = {(f["org"], f["bioproject"], f["run"]) for f in nuevas}
    quedan = [f for f in previas
              if (f["org"], f["bioproject"]) not in borrar
              and (f["org"], f["bioproject"], f["run"]) not in claves]
    todas = quedan + list(nuevas)
    todas.sort(key=lambda f: (f["org"], f["bioproject"], f["run"] != "-", f["run"]))
    with open(tabla, "w") as h:
        h.write(cab)
        h.write("\t".join(COLS) + "\n")
        for f in todas:
            h.write("\t".join(str(f.get(c, "-")) for c in COLS) + "\n")
    return todas


def fila(org, proy, run, c5, c3, sin, con, n, fuente):
    return {"org": org, "bioproject": proy, "run": run, "cortar_5p": c5, "cortar_3p": c3,
            "alin_sin_pct": f"{sin:.2f}", "alin_con_pct": f"{con:.2f}", "lecturas": n,
            "fuente": fuente, "fecha_utc": time.strftime("%Y-%m-%d", time.gmtime())}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tabla", type=pathlib.Path)
    ap.add_argument("--manifest", type=pathlib.Path)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("requeridos")
    r.add_argument("org")
    r.add_argument("rol")
    d = sub.add_parser("desactualizadas")
    d.add_argument("proyecto_dir", type=pathlib.Path)
    d.add_argument("org")
    a = ap.parse_args(argv)
    try:
        filas = leer(a.tabla)
    except ValueError as e:
        print(f"ERROR en {a.tabla or TABLA}: {e}", file=sys.stderr)
        return 2
    if a.cmd == "requeridos":
        for run, c5, c3 in requeridos(filas, a.org, a.rol, a.manifest):
            print(f"{run}\t{c5}\t{c3}")
    else:
        for run, apl, req in desactualizadas(a.proyecto_dir, a.org, filas):
            print(f"{run}\t{apl}\t{req}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
