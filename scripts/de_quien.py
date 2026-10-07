#!/usr/bin/env python3
"""¿De quién son los reads de un proyecto que casi no alinea?

    ./scripts/de_quien.py <proyecto_dir> <org> [--n 200000]

`align.sh verificar` dice CUÁNTO no alinea (SIN_AL), no POR QUÉ. Un SIN_AL alto
tiene dos lecturas opuestas: el pipeline está mal (genoma, recorte) o los reads
no son de este organismo. Paso con rhirr_duplicado: 99.9% sin alinear en las 13
corridas, con recorte y genoma verificados. Rhizophagus es simbionte obligado y
sclsc es necrotrofo: una librería de raíz colonizada o de tejido infectado es
casi toda del huésped vegetal, y eso no lo arregla ningún parámetro.

Lo que lo distingue en un minuto, sin red: los miRNAs de planta conservados
(miR166, miR156, miR159...) son de los sRNAs más abundantes de cualquier tejido
vegetal y no existen en hongos ni animales. Si un proyecto de hongo o de animal
los trae en cantidad, los reads son de una planta.

Muestrea las primeras --n lecturas de CADA librería del ledger (recortadas.tsv),
no solo de la primera: un BioProject puede mezclar cosas, y ya mezcló kits. Lee
las PRE-TRIMMED también, que en el ledger apuntan al fastq sin recortar.
"""
import argparse
import collections
import csv
import gzip
import pathlib
import sys

RAIZ = pathlib.Path(__file__).resolve().parent.parent

# Prefijo de 18 nt de miRNAs de planta muy conservados (miRBase, Arabidopsis y
# arroz comparten estos maduros). Prefijo y no secuencia entera: el 3' varía
# por isomiRs, el 5' no — el mismo criterio del etiquetado de positivos.
PLANTA = {
    "miR166": "TCGGACCAGGCTTCATTC",
    "miR156": "TGACAGAAGAGAGTGAGC",
    "miR159": "TTTGGATTGAAGGGAGCT",
    "miR168": "TCGCTTGGTGCAGATCGG",
    "miR167": "TGAAGCTGCCAGCATGAT",
    "miR396": "TTCCACAGCTTTCTTGAA",
    "miR319": "TTGGACTGAAGGGAGCTC",
    "miR160": "TGCCTGGCTCCCTGTATG",
}
# Por encima de esto, en un hongo o un animal, los reads son de una planta. En
# tejido vegetal estos ocho suman decenas de %; en un hongo puro deberían dar 0.
UMBRAL_PCT = 1.0


def reino_de(org, organismos=RAIZ / "data" / "organismos.tsv"):
    lineas = (l for l in organismos.read_text().splitlines() if not l.startswith("#"))
    for fila in csv.DictReader(lineas, delimiter="\t"):
        if fila["org"] == org:
            return fila["reino"]
    raise SystemExit(f"{org} no está en {organismos}")


def librerias(pdir):
    """Los ficheros que entraron al alineamiento, según el ledger del recorte."""
    led = pdir / "recortadas.tsv"
    if not led.exists():
        raise SystemExit(f"no hay {led}: ¿se recortó este proyecto?")
    filas = list(csv.DictReader(led.open(), delimiter="\t"))
    rutas = [(f["run"], pdir / f["fichero"]) for f in filas]
    faltan = [r for r, p in rutas if not p.exists()]
    if faltan:
        raise SystemExit(f"faltan en disco: {', '.join(faltan)}")
    return rutas


def muestrear(rutas, n):
    """Cuenta secuencias de las primeras n lecturas de cada librería."""
    total = collections.Counter()
    por_lib = {}
    for run, p in rutas:
        c = collections.Counter()
        abrir = gzip.open if p.suffix == ".gz" else open
        with abrir(p, "rt") as h:
            for i, linea in enumerate(h):
                if i % 4 == 1:
                    c[linea.strip()] += 1
                    if sum(c.values()) >= n:
                        break
        por_lib[run] = c
        total.update(c)
    return total, por_lib


def pct_planta(c):
    n = sum(c.values()) or 1
    return {k: 100 * sum(v for q, v in c.items() if q.startswith(s)) / n
            for k, s in PLANTA.items()}


def veredicto(reino, suma_pct):
    if reino == "Plantae":
        return "ESPERADO", "es una planta: estos miRNAs son suyos"
    if suma_pct >= UMBRAL_PCT:
        return ("HUÉSPED VEGETAL",
                "hay miRNAs de planta en cantidad: buena parte de los reads es de una "
                "planta, no de este organismo. No es el pipeline: es el experimento "
                "(tejido colonizado o infectado)")
    return ("SIN PLANTA",
            "no hay miRNAs de planta: el SIN_AL no es un huésped vegetal. "
            "Mirá el top de secuencias (rRNA, adaptador, otra especie)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("proyecto_dir", type=pathlib.Path)
    ap.add_argument("org")
    ap.add_argument("--n", type=int, default=200_000, help="lecturas por librería")
    ap.add_argument("--top", type=int, default=15)
    a = ap.parse_args(argv)

    reino = reino_de(a.org)
    rutas = librerias(a.proyecto_dir)
    total, por_lib = muestrear(rutas, a.n)
    n = sum(total.values())
    print(f"{a.proyecto_dir.name}: {len(rutas)} librerías, {n:,} lecturas "
          f"muestreadas ({a.n:,} por librería), reino {reino}\n")

    p = pct_planta(total)
    print("miRNAs de planta conservados (prefijo 5' de 18 nt):")
    for k, v in p.items():
        print(f"  {k:7s} {v:6.2f}%")
    suma = sum(p.values())
    print(f"  {'suma':7s} {suma:6.2f}%\n")

    print("por librería (suma de los ocho):")
    for run, c in por_lib.items():
        print(f"  {run:14s} {sum(pct_planta(c).values()):6.2f}%")

    print(f"\ntop {a.top} secuencias:")
    for q, v in total.most_common(a.top):
        marca = next((k for k, s in PLANTA.items() if q.startswith(s)), "")
        print(f"  {100 * v / n:6.2f}%  {len(q):2d}  {q}  {marca}")

    ver, porque = veredicto(reino, suma)
    print(f"\n>>> {ver}: {porque}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
