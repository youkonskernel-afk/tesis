#!/usr/bin/env python3
"""¿De quién son los reads de un proyecto que casi no alinea?

    ./scripts/de_quien.py <proyecto_dir> <org> [--n 200000]
        [--indice <prefijo> [--registrar data/extremos.tsv]]

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

POR CORRIDA, Y A LA TABLA

La prueba de bowtie se lee también corrida por corrida, con las mismas cuatro
llamadas: cada lectura lleva su número en el nombre y se sabe de qué librería
vino. Un porcentaje agregado de 82% no dice si las 13 corridas son 4N o si 11
lo son y 2 no — y cortar 4 nt de una corrida que no los trae acorta insertos de
verdad sin que nada falle. Fue lo de gadmo PRJNA328800 con el adaptador.

Con --registrar escribe el resultado en data/extremos.tsv (scripts/extremos.py),
que es de donde trim.sh lo toma. Una fila por BioProject cuando sus corridas
medidas son TODAS las del manifiesto y piden lo mismo; si no, una por corrida.
La muestra tiene que venir recortada SIN cortar extremos —es lo que hace
`trim.sh muestra`—, porque lo que se registra es el corte total, medido desde
el adaptador. Por eso solo lo usa `align.sh prueba` y no `verificar`, que mira
el recorte de verdad.
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
sys.path.insert(0, str(RAIZ / "scripts"))
import extremos as ext_tabla  # noqa: E402

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
# Por debajo de esto una corrida no decide nada por su cuenta: 1% de 50 lecturas
# es media lectura.
MIN_LECTURAS_RUN = 200


def reino_de(org, organismos=RAIZ / "data" / "organismos.tsv"):
    lineas = (l for l in organismos.read_text().splitlines() if not l.startswith("#"))
    for fila in csv.DictReader(lineas, delimiter="\t"):
        if fila["org"] == org:
            return fila["reino"]
    raise SystemExit(f"{org} no está en {organismos}")


def librerias(pdir):
    """Los ficheros que entraron al alineamiento, según el ledger del recorte.
    -> ([(run, ruta)], {run: bioproject})"""
    led = pdir / "recortadas.tsv"
    if not led.exists():
        raise SystemExit(f"no hay {led}: ¿se recortó este proyecto?")
    filas = list(csv.DictReader(led.open(), delimiter="\t"))
    rutas = [(f["run"], pdir / f["fichero"]) for f in filas]
    faltan = [r for r, p in rutas if not p.exists()]
    if faltan:
        raise SystemExit(f"faltan en disco: {', '.join(faltan)}")
    return rutas, {f["run"]: f["bioproject"] for f in filas}


def muestrear(rutas, n, n_bowtie=0):
    """Cuenta secuencias de las primeras n lecturas de cada librería, y guarda
    las primeras n_bowtie de cada una para la prueba de alineamiento.

    El tope se lleva con un entero: sumar el Counter en cada lectura era
    cuadrático, y con 200 000 lecturas casi todas distintas no terminaba."""
    total = collections.Counter()
    por_lib = {}
    para_bowtie = []     # (run, secuencia)
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
                        para_bowtie.append((run, q))
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


def prueba_bowtie(muestra, indice, cores):
    """% alineado con bowtie -v 1 (como yasma) para cada recorte de RECORTES.
    muestra: [(run, secuencia)]. -> (agregado, {run: {recorte: %}}, {run: n})

    Cuenta los nombres que salen con -k 1 en vez de leer el resumen de stderr,
    que cambia de redacción entre versiones de bowtie. El nombre es m<i>, y el
    i dice de qué corrida es la lectura: con eso las MISMAS cuatro llamadas dan
    el resultado de cada corrida."""
    res, por_run = {}, collections.defaultdict(dict)
    n_run = collections.Counter(run for run, _ in muestra)
    with tempfile.TemporaryDirectory() as d:
        fq = pathlib.Path(d) / "muestra.fq"
        with fq.open("w") as h:
            for i, (_, q) in enumerate(muestra):
                h.write(f"@m{i}\n{q}\n+\n{'I' * len(q)}\n")
        for t5, t3 in RECORTES:
            r = subprocess.run(["bowtie", "-q", "-v", "1", "-k", "1", "-p", str(cores),
                                "--trim5", str(t5), "--trim3", str(t3), str(indice), str(fq)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"bowtie falló ({r.returncode}): {r.stderr.strip()[-300:]}")
            vistas = set()
            for l in r.stdout.splitlines():
                nombre = l.split("\t", 1)[0].strip()
                if nombre.startswith("m") and nombre[1:].isdigit():
                    vistas.add(int(nombre[1:]))
            res[(t5, t3)] = 100 * len(vistas) / max(len(muestra), 1)
            cuenta = collections.Counter(muestra[i][0] for i in vistas if i < len(muestra))
            for run, n in n_run.items():
                por_run[run][(t5, t3)] = 100 * cuenta[run] / n
    return res, dict(por_run), dict(n_run)


def elegir_recorte(bt):
    """((5', 3'), % con ese recorte) si cortar arregla el alineamiento, o None.

    El recorte MAS CHICO que llega casi al maximo, no el maximo a secas: si el
    5' no tiene bases aleatorias, recortar las dos puntas igual alinea un poco
    mas (un read mas corto pega por azar mas facil), y elegir el maximo diria
    "5' y 3'" cuando es solo el 3'."""
    base = bt[(0, 0)]
    tope = max(bt.values())
    (t5, t3), mejor = next(
        (k, v) for k, v in sorted(bt.items(), key=lambda kv: sum(kv[0]))
        if v >= tope - max(2.0, 0.05 * tope))
    if mejor - base >= MEJORA_MIN and mejor >= MEJORA_X * max(base, 0.1):
        return (t5, t3), mejor
    return None


def veredicto_extremos(bt, dom3, dom5):
    """(VEREDICTO, porqué) si el problema son los extremos, o None."""
    if bt:
        base = bt[(0, 0)]
        elegido = elegir_recorte(bt)
        if elegido:
            (t5, t3), mejor = elegido
            cuales = " y ".join(x for x, t in (("5'", t5), ("3'", t3)) if t)
            return ("EXTREMOS SIN RECORTAR",
                    f"sacando {K} nt del {cuales} alinea {mejor:.1f}% contra {base:.1f}% "
                    f"sin sacarlos. El kit deja bases aleatorias pegadas al inserto "
                    f"(adaptadores 4N) y el recorte no las saca: no es el genoma ni el "
                    f"organismo, es el recorte. Hay que cortar {K} nt del {cuales} antes "
                    f"de alinear: `align.sh prueba <org>/<rol>` lo registra en "
                    f"data/extremos.tsv y `trim.sh correr` rehace lo recortado sin cortar")
        return None
    aleat = [x for x, d in (("3'", dom3), ("5'", dom5)) if d is not None and d < DOMINANTE_MAX]
    if aleat:
        return ("EXTREMOS ALEATORIOS",
                f"los {K} nt del {' y '.join(aleat)} varían sin una variante que domine, "
                f"como bases aleatorias de un adaptador 4N. Sin índice no lo puedo "
                f"confirmar alineando: corré con --indice")
    return None


def decisiones(por_run, n_run):
    """{run: (5', 3')} con lo que pide cada corrida; None si tiene muy pocas
    lecturas para decir algo."""
    return {run: (None if n_run.get(run, 0) < MIN_LECTURAS_RUN
                  else (elegir_recorte(bt) or ((0, 0), None))[0])
            for run, bt in por_run.items()}


def filas_extremos(org, por_run, n_run, proy_de, manifest, existentes, fuente):
    """Lo que hay que escribir en data/extremos.tsv. -> (nuevas, borrar)

    Por BioProject: si se midieron TODAS sus corridas del manifiesto, lo viejo
    se borra y va una fila `-` cuando todas piden lo mismo (ninguna si no piden
    cortar), o una por corrida si no. Si se midieron algunas, solo se escriben
    las que difieren de lo que la tabla ya les aplica — incluido un 0/0
    explícito cuando una fila `-` las cortaría."""
    dec = decisiones(por_run, n_run)
    del_manifest = collections.defaultdict(set)
    with open(manifest) as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["org"] == org:
                del_manifest[r["bioproject"]].add(r["run"])
    por_proy = collections.defaultdict(list)
    for run in por_run:
        por_proy[proy_de.get(run, "-")].append(run)

    def pct(runs, recorte):
        n = sum(n_run[r] for r in runs)
        return sum(por_run[r][recorte] * n_run[r] for r in runs) / max(n, 1), n

    nuevas, borrar = [], set()
    for proy, runs in sorted(por_proy.items()):
        medidas = [r for r in runs if dec[r] is not None]
        todas = del_manifest.get(proy, set())
        if todas and todas <= set(medidas):
            borrar.add((org, proy))
            pedidos = {dec[r] for r in todas}
            if pedidos == {(0, 0)}:
                continue
            if len(pedidos) == 1:
                (c5, c3), = pedidos
                sin, n = pct(sorted(todas), (0, 0))
                con, _ = pct(sorted(todas), (c5, c3))
                nuevas.append(ext_tabla.fila(org, proy, "-", c5, c3, sin, con, n, fuente))
                continue
            elegibles = [r for r in sorted(todas) if dec[r] != (0, 0)]
        else:
            elegibles = [r for r in sorted(medidas)
                         if dec[r] != ext_tabla.requerido(existentes, org, proy, r)]
        for r in elegibles:
            c5, c3 = dec[r]
            nuevas.append(ext_tabla.fila(org, proy, r, c5, c3, por_run[r][(0, 0)],
                                         por_run[r][(c5, c3)], n_run[r], fuente))
    return nuevas, borrar


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
    ap.add_argument("--registrar", type=pathlib.Path,
                    help="con --indice: escribe el corte de cada proyecto/corrida en esta "
                         "tabla (data/extremos.tsv). Solo sobre una muestra recortada SIN "
                         "cortar extremos, que es lo que deja `trim.sh muestra`")
    ap.add_argument("--manifest", type=pathlib.Path, default=ext_tabla.MANIFEST,
                    help="para saber si se midieron todas las corridas de un BioProject")
    a = ap.parse_args(argv)
    if a.registrar and not a.indice:
        ap.error("--registrar necesita --indice: sin alinear no hay medición que registrar")

    reino = reino_de(a.org)
    rutas, proy_de = librerias(a.proyecto_dir)
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
        bt, por_run, n_run = prueba_bowtie(para_bt, a.indice, a.cores)
        for (t5, t3), pct in bt.items():
            print(f"  --trim5 {t5} --trim3 {t3}   {pct:6.2f}% alineado")
        print()
        dec = decisiones(por_run, n_run)
        print(f"por corrida, % alineado con cada corte (5'/3'):")
        print(f"  {'CORRIDA':14s} {'LECT':>7s} " + " ".join(f"{f'{x}/{y}':>7s}" for x, y in RECORTES)
              + "   PIDE")
        for run in por_run:
            d = dec[run]
            pide = ("pocas lecturas" if d is None else
                    "nada" if d == (0, 0) else f"cortar {d[0]}/{d[1]}")
            print(f"  {run:14s} {n_run[run]:7d} "
                  + " ".join(f"{por_run[run][r]:7.2f}" for r in RECORTES) + f"   {pide}")
        pedidos = {d for d in dec.values() if d is not None}
        if len(pedidos) > 1:
            print("  OJO: las corridas no piden todas lo mismo — el proyecto mezcla kits")
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

    if a.registrar:
        existentes = ext_tabla.leer(a.registrar)
        nuevas, borrar = filas_extremos(
            a.org, por_run, n_run, proy_de, a.manifest, existentes,
            f"de_quien.py --indice sobre {a.proyecto_dir.name}")
        sacadas = [f for f in existentes if (f["org"], f["bioproject"]) in borrar]
        ext_tabla.registrar(nuevas, borrar, tabla=a.registrar)
        print(f"\n>>> {a.registrar.name}: {len(nuevas)} fila(s) escrita(s)"
              + (f", {len(sacadas)} vieja(s) reemplazada(s)" if sacadas else ""))
        for f in nuevas:
            print(f"    {f['org']} {f['bioproject']} {f['run']}: cortar "
                  f"{f['cortar_5p']}/{f['cortar_3p']}  ({f['alin_sin_pct']}% -> "
                  f"{f['alin_con_pct']}% alineado, {f['lecturas']} lecturas)")
        if nuevas or sacadas:
            print("    trim.sh correr rehace lo recortado con otro corte; empujá la tabla "
                  "(§5d de 20_alinear) o se pierde con la VM")
        else:
            print("    la tabla ya dice lo que se midió: no hay nada que cambiar")
    return 0


if __name__ == "__main__":
    sys.exit(main())
