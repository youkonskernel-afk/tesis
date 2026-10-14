#!/usr/bin/env python3
"""BAM sintetico para comparar yasma tradeoff original contra el parchado.

    tests/bam_sintetico.py <salida.bam> [--libs 4] [--seed 1] [--escala 1]
                           [--largos chrA:400000,chrB:150000,...]

Tiene lo que ejercita los caminos de tradeoff que el parche toca: pilas tipo
miRNA (con estrella), clusters de siRNA de 24 nt en las dos hebras, phasiRNA,
tRF de 33 nt, reads pegados a los dos extremos de cada cromosoma, un cromosoma
mas corto que la ventana de cobertura (la rama "too short"), un locus en una
sola libreria (la mediana entre librerias lo baja), un empate de largos 21/22
en la misma posicion (el desempate de sizeClass depende del orden), ruido de
fondo y reads sin mapear. La profundidad de cada libreria es distinta, para que
el rpb importe.
"""
import argparse
import random

import pysam

COMP = str.maketrans("ACGT", "TGCA")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("salida")
    ap.add_argument("--libs", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--escala", type=float, default=1.0)
    ap.add_argument("--largos", default="chrA:400000,chrB:150000,chrC:180,chrD:5000,chrE:3000")
    a = ap.parse_args()
    rnd = random.Random(a.seed)

    crs = [(c.split(":")[0], int(c.split(":")[1])) for c in a.largos.split(",")]
    gen = {c: "".join(rnd.choice("ACGT") for _ in range(n)) for c, n in crs}
    libs = [f"SRR{1000 + i}" for i in range(a.libs)]
    hdr = {"HD": {"VN": "1.6", "SO": "unsorted"},
           "SQ": [{"SN": c, "LN": n} for c, n in crs],
           "RG": [{"ID": l, "SM": l} for l in libs]}
    reads = []

    def lectura(c, ini, largo, hebra, lib):
        n = dict(crs)[c]
        if ini < 0 or ini + largo > n:
            return
        reads.append((c, ini, largo, hebra, lib))

    def peso_lib(i):
        # profundidad distinta por libreria, para que el rpb y la mediana importen
        return [1.0, 0.4, 2.3, 0.05, 1.7, 0.9, 3.1, 0.2][i % 8]

    def pila(c, ini, largo, hebra, total, var3=True):
        for i, l in enumerate(libs):
            k = int(total * peso_lib(i) * a.escala * rnd.uniform(0.6, 1.4))
            for _ in range(k):
                ll = largo + (rnd.choice([0, 0, 0, 1, -1]) if var3 else 0)
                lectura(c, ini, ll, hebra, l)

    # miRNA: brazo 5p muy abundante, estrella poco
    for c, ini, h in [("chrA", 10000, "+"), ("chrA", 52000, "-"), ("chrA", 230000, "+"),
                      ("chrB", 40000, "-"), ("chrD", 2000, "+")]:
        pila(c, ini, 22, h, 800)
        pila(c, ini + 50, 21, h, 40)
    # miRNA en una sola libreria (la mediana lo baja)
    for _ in range(int(3000 * a.escala)):
        lectura("chrA", 300000, 22, "+", libs[0])
    # cluster de siRNA de 24 nt, dos hebras
    for c, ini, ln in [("chrA", 80000, 1500), ("chrA", 150000, 900), ("chrB", 100000, 2000)]:
        for _ in range(int(4000 * a.escala)):
            lectura(c, ini + rnd.randrange(ln), rnd.choice([24, 24, 24, 23]),
                     rnd.choice("+-"), rnd.choice(libs))
    # phasiRNA de 21 nt
    for f in range(12):
        pila("chrA", 120000 + 21 * f, 21, "+", 60, var3=False)
    # tRF de 33 nt, muchas copias
    for c, ini in [("chrA", 180000), ("chrA", 181500), ("chrB", 5000)]:
        pila(c, ini, 33, "-", 300, var3=False)
    # bordes de cromosoma
    pila("chrA", 0, 22, "+", 200)
    pila("chrA", 400000 - 22, 22, "-", 200)
    pila("chrB", 150000 - 30, 30, "+", 100)
    pila("chrC", 10, 22, "+", 50)
    pila("chrD", 0, 24, "-", 120)
    # empate de tamaños: 21 y 22 con la misma profundidad en la misma posicion
    for i in range(200):
        for l in libs[:2]:
            lectura("chrA", 260000, 21 if i % 2 else 22, "+", l)
    # Empate de claves de tamaño que solo el ORDEN de llegada desempata, en un
    # cromosoma sin ruido. En la posicion c llegan primero los de 23 nt (empiezan
    # en c-11) y despues los de 21 (empiezan en c-10); los de 22 estan en c+1.
    # Con 21 = 23 y 22 por debajo de la mitad, sizeClass no llama un solo largo y
    # (21,22) empata con (22,23): gana la clave que entro primero al Counter, y
    # eso depende de que el 23 haya llegado antes que el 21. Si el parche
    # reordenara los largos de una posicion, este locus cambiaria de sizecall.
    if "chrE" in dict(crs):
        c = 1500
        # en todas las librerias: con una sola, la mediana entre librerias da 0
        for l in libs:
            for _ in range(int(300 * a.escala)):
                lectura("chrE", c - 11, 23, "+", l)
                lectura("chrE", c - 10, 21, "+", l)
            for _ in range(int(400 * a.escala)):
                lectura("chrE", c + 1 - 11, 22, "+", l)
    # ruido de fondo
    for _ in range(int(30000 * a.escala)):
        c, n = rnd.choice(crs)
        if n < 60 or c == "chrE":
            continue
        lectura(c, rnd.randrange(n - 55), rnd.randrange(15, 51), rnd.choice("+-"),
                rnd.choice(libs))

    rnd.shuffle(reads)
    tmp = a.salida + ".unsorted.bam"
    ids = {c: i for i, (c, _) in enumerate(crs)}
    with pysam.AlignmentFile(tmp, "wb", header=hdr) as out:
        for j, (c, ini, largo, hebra, lib) in enumerate(reads):
            s = pysam.AlignedSegment()
            s.query_name = f"r{j}"
            s.reference_id = ids[c]
            s.reference_start = ini
            s.cigartuples = [(0, largo)]
            s.flag = 16 if hebra == "-" else 0
            s.mapping_quality = 255
            s.query_sequence = gen[c][ini:ini + largo]
            s.query_qualities = pysam.qualitystring_to_array("I" * largo)
            s.set_tag("RG", lib)
            s.set_tag("XY", rnd.choice("UUUPR"))
            out.write(s)
        for j in range(int(500 * a.escala)):
            s = pysam.AlignedSegment()
            s.query_name = f"u{j}"
            s.flag = 4
            s.reference_id = -1
            s.reference_start = -1
            s.query_sequence = "ACGTACGTACGTACGTACGTAC"
            s.query_qualities = pysam.qualitystring_to_array("I" * 22)
            s.set_tag("RG", rnd.choice(libs))
            s.set_tag("XY", "N")
            out.write(s)
    pysam.sort("-o", a.salida, tmp)
    pysam.index(a.salida)
    print(f"{a.salida}: {len(reads)} reads, {len(libs)} librerias")


if __name__ == "__main__":
    main()
