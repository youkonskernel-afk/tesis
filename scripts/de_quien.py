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

Y no era el huésped: rhirr_duplicado dio 0% de miRNAs de planta, y el top de
secuencias mostró otra cosa — las 15 más abundantes eran el MISMO inserto de 32
nt seguido de 4 nt distintos en cada una. Un genoma no produce eso; un kit con
adaptadores 4N (bases aleatorias en la unión, p. ej. NEXTflex) sí, y si nadie
las recorta, `bowtie -v 1` no alinea un read que trae 4 bases inventadas. Por
eso hay dos mediciones más, y van PRIMERO porque son las que el pipeline puede
arreglar:

  1. Extremos: para los núcleos más abundantes, cuántas variantes distintas de
     los 4 nt de cada punta aparecen. Un isomiR real tiene una variante que
     domina; una base aleatoria, ninguna.
  2. Con --indice: alinea una muestra con bowtie -v 1 (lo mismo que yasma)
     cuatro veces, sacando o no 4 nt de cada extremo. Si recortar dispara la
     fracción alineada, el diagnóstico está hecho y dice qué recortar.

Muestrea las primeras --n lecturas de CADA librería del ledger (recortadas.tsv),
no solo de la primera: un BioProject puede mezclar cosas, y ya mezcló kits. Lee
las PRE-TRIMMED también, que en el ledger apuntan al fastq sin recortar.
"""
import argparse
import collections
import csv
import gzip
import os
import pathlib
import subprocess
import sys
import tempfile

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

# Extremos. K nt de cada punta; se miran los TOP_NUCLEOS núcleos más abundantes
# que tengan al menos MIN_NUCLEO lecturas. Si en promedio la variante más común
# de la punta no llega a DOMINANTE_MAX, la punta es aleatoria: en un isomiR real
# la secuencia templada domina con holgura.
K = 4
TOP_NUCLEOS = 20
MIN_NUCLEO = 20
DOMINANTE_MAX = 0.30
# Recortes que prueba bowtie (5', 3'), y cuánto tiene que mejorar el mejor para
# decir que el problema son los extremos: al menos MEJORA_MIN puntos y MEJORA_X
# veces lo que alinea sin recortar.
RECORTES = [(0, 0), (0, K), (K, 0), (K, K)]
MEJORA_MIN = 10.0
MEJORA_X = 5.0


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


def muestrear(rutas, n, n_bowtie=0):
    """Cuenta secuencias de las primeras n lecturas de cada librería, y guarda
    las primeras n_bowtie de cada una para la prueba de alineamiento.

    El tope se lleva con un entero: sumar el Counter en cada lectura era
    cuadrático, y con 200 000 lecturas casi todas distintas no terminaba."""
    total = collections.Counter()
    por_lib = {}
    para_bowtie = []
    for run, p in rutas:
        c = collections.Counter()
        leidas = 0
        abrir = gzip.open if p.suffix == ".gz" else open
        with abrir(p, "rt") as h:
            for i, linea in enumerate(h):
                if i % 4 == 1:
                    q = linea.strip()
                    c[q] += 1
                    if leidas < n_bowtie:
                        para_bowtie.append(q)
                    leidas += 1
                    if leidas >= n:
                        break
        por_lib[run] = c
        total.update(c)
    return total, por_lib, para_bowtie


def extremos(total, lado):
    """Para los núcleos más abundantes, qué fracción se lleva la variante más
    común de los K nt de la punta `lado` ('5' o '3'). Devuelve
    (dominante_promedio, variantes_promedio, núcleos_mirados)."""
    nucleos = collections.defaultdict(collections.Counter)
    for q, v in total.items():
        if len(q) < K + 15:
            continue
        nucleo, punta = (q[K:], q[:K]) if lado == "5" else (q[:-K], q[-K:])
        nucleos[nucleo][punta] += v
    top = sorted(nucleos.values(), key=lambda c: -sum(c.values()))
    top = [c for c in top if sum(c.values()) >= MIN_NUCLEO][:TOP_NUCLEOS]
    if not top:
        return None, None, 0
    peso = sum(sum(c.values()) for c in top)
    dom = sum(max(c.values()) for c in top) / peso
    var = sum(len(c) for c in top) / len(top)
    return dom, var, len(top)


def prueba_bowtie(seqs, indice, cores):
    """% alineado con bowtie -v 1 (como yasma) para cada recorte de RECORTES.
    Cuenta las líneas de salida con -k 1 en vez de leer el resumen de stderr,
    que cambia de redacción entre versiones de bowtie."""
    res = {}
    with tempfile.TemporaryDirectory() as d:
        fq = pathlib.Path(d) / "muestra.fq"
        with fq.open("w") as h:
            for i, q in enumerate(seqs):
                h.write(f"@m{i}\n{q}\n+\n{'I' * len(q)}\n")
        for t5, t3 in RECORTES:
            r = subprocess.run(["bowtie", "-q", "-v", "1", "-k", "1", "-p", str(cores),
                                "--trim5", str(t5), "--trim3", str(t3), str(indice), str(fq)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"bowtie falló ({r.returncode}): {r.stderr.strip()[-300:]}")
            alineadas = sum(1 for l in r.stdout.splitlines() if l.strip())
            res[(t5, t3)] = 100 * alineadas / max(len(seqs), 1)
    return res


def veredicto_extremos(bt, dom3, dom5):
    """(VEREDICTO, porqué) si el problema son los extremos, o None."""
    if bt:
        base = bt[(0, 0)]
        tope = max(bt.values())
        # El recorte MAS CHICO que llega casi al maximo, no el maximo a secas: si
        # el 5' no tiene bases aleatorias, recortar las dos puntas igual alinea
        # un poco mas (un read mas corto pega por azar mas facil), y elegir el
        # maximo diria "5' y 3'" cuando es solo el 3'.
        (t5, t3), mejor = next(
            (k, v) for k, v in sorted(bt.items(), key=lambda kv: sum(kv[0]))
            if v >= tope - max(2.0, 0.05 * tope))
        if mejor - base >= MEJORA_MIN and mejor >= MEJORA_X * max(base, 0.1):
            cuales = " y ".join(x for x, t in (("5'", t5), ("3'", t3)) if t)
            return ("EXTREMOS SIN RECORTAR",
                    f"sacando {K} nt del {cuales} alinea {mejor:.1f}% contra {base:.1f}% "
                    f"sin sacarlos. El kit deja bases aleatorias pegadas al inserto "
                    f"(adaptadores 4N) y el recorte no las saca: no es el genoma ni el "
                    f"organismo, es el recorte. Hay que recortar {K} nt del {cuales} y "
                    f"re-alinear")
        return None
    aleat = [x for x, d in (("3'", dom3), ("5'", dom5)) if d is not None and d < DOMINANTE_MAX]
    if aleat:
        return ("EXTREMOS ALEATORIOS",
                f"los {K} nt del {' y '.join(aleat)} varían sin una variante que domine, "
                f"como bases aleatorias de un adaptador 4N. Sin índice no lo puedo "
                f"confirmar alineando: corré con --indice")
    return None


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
    ap.add_argument("--indice", type=pathlib.Path,
                    help="prefijo del índice de bowtie: activa la prueba de recortes")
    ap.add_argument("--n-bowtie", type=int, default=20_000,
                    help="lecturas por librería para la prueba de bowtie")
    ap.add_argument("--cores", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)

    reino = reino_de(a.org)
    rutas = librerias(a.proyecto_dir)
    total, por_lib, para_bt = muestrear(rutas, a.n, a.n_bowtie if a.indice else 0)
    n = sum(total.values())
    print(f"{a.proyecto_dir.name}: {len(rutas)} librerías, {n:,} lecturas "
          f"muestreadas ({a.n:,} por librería), reino {reino}\n")

    dom3, var3, n3 = extremos(total, "3")
    dom5, var5, n5 = extremos(total, "5")
    print(f"extremos ({K} nt de cada punta, en los {TOP_NUCLEOS} núcleos más abundantes):")
    for lado, dom, var, nn in (("3'", dom3, var3, n3), ("5'", dom5, var5, n5)):
        if dom is None:
            print(f"  {lado}  sin núcleos con {MIN_NUCLEO}+ lecturas")
        else:
            marca = "  <- ALEATORIO" if dom < DOMINANTE_MAX else ""
            print(f"  {lado}  la variante más común se lleva {100 * dom:5.1f}%  "
                  f"({var:.0f} variantes por núcleo, {nn} núcleos){marca}")
    print()

    bt = None
    if a.indice:
        print(f"bowtie -v 1 sobre {len(para_bt):,} lecturas, recortando o no {K} nt:")
        bt = prueba_bowtie(para_bt, a.indice, a.cores)
        for (t5, t3), pct in bt.items():
            print(f"  --trim5 {t5} --trim3 {t3}   {pct:6.2f}% alineado")
        print()

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

    ext = veredicto_extremos(bt, dom3, dom5)
    if ext:
        print(f"\n>>> {ext[0]}: {ext[1]}")
    ver, porque = veredicto(reino, suma)
    print(f"{'' if ext else chr(10)}>>> {ver}: {porque}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
