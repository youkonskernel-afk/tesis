#!/usr/bin/env python3
"""El parche de memoria de yasma tradeoff: que guarde disperso y devuelva IGUAL.

`tradeoff.py` de v1.1.1 reserva 8 B por base en `pos_size_d` y 8*C*N B por base
en `pos_d` (N librerias): `-9` en los cuatro proyectos que se probaron en Colab.
El parche (scripts/tradeoff_parche.py + scripts/tradeoff_disperso.py) guarda
solo donde hay reads. Lo unico que importa es que YASMA lea lo mismo.

Lo que este banco protege:
 1. Conteos: cada rebanada que YASMA pide da el mismo arreglo uint32 que el
    denso, con la semantica de rebanada de numpy (negativos, fuera de rango,
    vacia); la cobertura sale bit a bit igual a las cinco lineas del original,
    con 1-3 condiciones, 1-7 librerias, rpb con ceros y con desborde de uint32;
    y la compactacion y los bloques forzados no cambian nada.
 2. Largos: la misma lista por posicion, con el MISMO ORDEN de primera
    aparicion. sizeClass desempata con Counter.most_common(), que respeta el
    orden de insercion: reordenar cambia sizecalls empatados (medido: un locus
    pasa de 22_23 a 21_22, ver tests/bam_sintetico.py, chrE).
 3. El parcheador: aplica las doce lineas, el resultado corre y da lo mismo que
    el original (sobre una version minima de get_kernel_coverage con las
    mismas lineas), es idempotente, y se NIEGA si una linea falta o se repite.
    Detecta un modulo desactualizado.
 4. tradeoff_control.py: compara lo que tiene que comparar e ignora solo la
    hora y las rutas.
 5. Si hay un yasma instalado (YASMA_PYTHON=<python con yasma>), el original
    contra el parchado sobre BAMs sinteticos, byte a byte.
"""
import collections
import contextlib
import io
import importlib
import os
import pathlib
import random
import subprocess
import sys
import tempfile

import numpy as np

RAIZ = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "scripts"))
import tradeoff_control  # noqa: E402
import tradeoff_disperso as td  # noqa: E402
import tradeoff_parche  # noqa: E402

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


def raises(f, exc):
    try:
        f()
    except exc:
        return True
    except Exception:  # noqa: BLE001
        return False
    return False


def original_cobertura(denso, rpbs):
    """Las cinco lineas de tradeoff.py v1.1.1, literales."""
    coverage = np.sum(denso, axis=3, dtype='uint32')
    coverage = np.multiply(coverage, rpbs)
    coverage = np.median(coverage, axis=(2))
    coverage = coverage.astype('float32')
    coverage = np.mean(coverage, axis=(1))
    return coverage


def lecturas(rnd, L, C, N, n):
    """Reads amontonados en pocas posiciones (como sRNA), mas ruido y bordes."""
    pilas = [rnd.randrange(L) for _ in range(8)]
    out = []
    for _ in range(n):
        r = rnd.random()
        if r < 0.6:
            p = rnd.choice(pilas)
        elif r < 0.65:
            p = rnd.choice([0, L - 1])
        else:
            p = rnd.randrange(L)
        out.append((p, rnd.randrange(C), rnd.randrange(N), rnd.randrange(2)))
    return out


print("== 1. Conteos: lo mismo que el arreglo denso")
rnd = random.Random(11)
casos_mal = []
for C, N in [(1, 1), (1, 2), (1, 4), (1, 7), (2, 3), (3, 2)]:
    for compactar, bloque in [(1 << 22, 64 << 20), (5, 64), (97, 300)]:
        td.COMPACTAR, td.BLOQUE_B = compactar, bloque
        L = 3000
        denso = np.zeros((L, C, N, 2), dtype='uint32')
        cs = td.Conteos(L, C, N)
        for p, c, l, s in lecturas(rnd, L, C, N, 4000):
            denso[p, c, l, s] += 1
            cs.sumar(p, c, l, s)
        rpbs = np.array([[rnd.choice([0, 1, 37, 900, 2**31]) for _ in range(N)]
                         for _ in range(C)], dtype='uint32')
        tag = f"C={C} N={N} compactar={compactar} bloque={bloque}"
        ref = original_cobertura(denso, rpbs)
        got = cs.cobertura(rpbs)
        if got.dtype != ref.dtype or got.tobytes() != ref.tobytes():
            casos_mal.append(f"cobertura {tag}")
        if cs.por_posicion().tobytes() != np.sum(denso, axis=(1, 2, 3), dtype='uint32').tobytes():
            casos_mal.append(f"por_posicion {tag}")
        for sl in [slice(0, 10), slice(-5, None), slice(None, None), slice(L - 3, L + 10),
                   slice(10, 5), slice(-L - 5, 3), slice(2990, 3000), slice(1500, 1501),
                   slice(0, 0), slice(5, None), slice(None, 7)]:
            for clave in (sl, (sl,), (sl, Ellipsis)):
                a, b = cs[clave], denso[clave]
                if a.shape != b.shape or a.dtype != b.dtype or not np.array_equal(a, b):
                    casos_mal.append(f"rebanada {clave!r} {tag}")
chk("rebanadas, cobertura y suma por posición iguales al denso en los 18 casos",
    not casos_mal, "; ".join(casos_mal[:4]))
td.COMPACTAR, td.BLOQUE_B = 1 << 22, 64 << 20

# Lo que el desborde de uint32 hace en el original tiene que pasar igual aca:
# cuenta 3 por rpb 2^31 da 2^32 + 2^31, que en uint32 es 2^31.
cs = td.Conteos(300, 1, 1)
denso = np.zeros((300, 1, 1, 2), dtype='uint32')
for _ in range(3):
    cs.sumar(150, 0, 0, 0)
    denso[150, 0, 0, 0] += 1
rp = np.array([[2**31]], dtype='uint32')
chk("el desborde de uint32 da lo mismo que el original",
    cs.cobertura(rp).tobytes() == original_cobertura(denso, rp).tobytes()
    and cs.cobertura(rp)[150] == np.float32(2**31))
chk("sin reads la cobertura es ceros float32 de largo L",
    td.Conteos(50, 2, 3).cobertura(np.ones((2, 3), 'uint32')).tobytes()
    == np.zeros(50, np.float32).tobytes())
cs = td.Conteos(100, 1, 2)
chk("sumar fuera del cromosoma revienta, como el denso",
    raises(lambda: cs.sumar(100, 0, 0, 0), IndexError))
cs.sumar(-1, 0, 1, 1)
chk("y un índice negativo cuenta desde el final, como el denso",
    cs[99:100][0, 0, 1, 1] == 1 and cs.por_posicion()[99] == 1)
chk("un índice que no es rebanada no se atiende a ciegas",
    raises(lambda: cs[5], TypeError) and raises(lambda: cs[0:10:2], TypeError))

# Al pasar de cromosoma, el buffer del anterior se compacta: si no, con cientos
# de scaffolds cada uno retiene hasta COMPACTAR lecturas sin comprimir.
a, b = td.Conteos(100, 1, 1), td.Conteos(100, 1, 1)
for _ in range(3):
    a.sumar(10, 0, 0, 0)
pend = len(a._buf)
b.sumar(5, 0, 0, 0)
chk("al empezar otro cromosoma se compacta el anterior",
    pend == 3 and len(a._buf) == 0 and a.por_posicion()[10] == 3, (pend, len(a._buf)))
la, lb = td.Largos(100), td.Largos(100)
la.agregar(3, 22)
lb.agregar(4, 21)
chk("lo mismo con los largos", len(la._buf) == 0 and la[3] == [22])

print("== 2. Largos: la misma lista, en el mismo orden de primera aparición")
orden_mal = []
for compactar in (1 << 22, 3, 50):
    td.COMPACTAR = compactar
    L = 400
    ref = [None] * L
    ls = td.Largos(L)
    for _ in range(5000):
        p = rnd.choice([0, L - 1, 17, 17, 18, rnd.randrange(L)])
        largo = rnd.choice([15, 21, 22, 22, 23, 24, 33, 50])
        try:
            ref[p].append(largo)
        except AttributeError:
            ref[p] = [largo]
        ls.agregar(p, largo)
    for p in range(-L, L):
        a, b = ls[p], ref[p]
        # Counter conserva el orden de insercion: comparar sus items compara
        # las cuentas Y el orden de primera aparicion, que es lo que sizeClass ve.
        if (a is None) != (b is None) or (b is not None and
                list(collections.Counter(a).items()) != list(collections.Counter(b).items())):
            orden_mal.append((compactar, p))
    # y la secuencia concatenada de posiciones consecutivas, como en una ventana
    for ini in (0, 10, 390):
        flat_a = [x for p in range(ini, ini + 10) for x in (ls[p] or [])]
        flat_b = [x for p in range(ini, ini + 10) for x in (ref[p] or [])]
        if list(collections.Counter(flat_a).items()) != list(collections.Counter(flat_b).items()):
            orden_mal.append((compactar, f"ventana {ini}"))
td.COMPACTAR = 1 << 22
chk("cuentas y orden de primera aparición iguales en cada posición y ventana",
    not orden_mal, orden_mal[:5])

ls = td.Largos(50)
for x in (23, 21, 22, 21, 23):
    ls.agregar(7, x)
chk("[23, 21, 22, 21, 23] sale con el 23 primero y el 21 segundo",
    ls[7] == [23, 23, 21, 21, 22], ls[7])
chk("sin reads devuelve None, como la lista", ls[8] is None)
chk("los largos son int de Python (con numpy, repr diria np.int64(22))",
    all(type(x) is int for x in ls[7]))
chk("fuera de rango revienta igual que la lista",
    raises(lambda: ls[50], IndexError) and raises(lambda: ls[-51], IndexError)
    and raises(lambda: ls.agregar(50, 22), IndexError))
chk("una escritura directa que el parche no esperaba revienta",
    raises(lambda: ls.__setitem__(3, [22]), TypeError))

print("== 3. el parcheador, sobre una get_kernel_coverage mínima con las mismas líneas")
# Las doce lineas que el parche reemplaza, copiadas de tradeoff.py v1.1.1 con su
# indentacion de tabs, en un programa que hace con ellas lo mismo que YASMA.
FUENTE = '''\
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view


def kernel(chromosomes, annotation_conditions, max_lib_counts, reads, rpbs):
\tpos_size_d = dict()
\tpos_d      = dict()
\tfor chrom, chrom_length in chromosomes:
\t\tpos_size_d[chrom] = [None] * chrom_length
\t\tpos_d[chrom] = np.zeros(shape=(chrom_length, len(annotation_conditions), max_lib_counts, 2), dtype='uint32')
\tfor chrom, pos, condition_i, lib_i, strand_i, length in reads:
\t\tpos_d[chrom][pos, condition_i, lib_i, strand_i] += 1
\t\ttry:
\t\t\tpos_size_d[chrom][pos].append(length)
\t\texcept AttributeError:
\t\t\tpos_size_d[chrom][pos] = [length]
\tcov_d = {}
\tout = {}
\tfor chrom, chrom_length in chromosomes:
\t\tcoverage = np.sum(pos_d[chrom], axis=3, dtype='uint32')
\t\tif coverage.shape[0] < 250:
\t\t\tcoverage = np.zeros(shape=(coverage.shape[0]), dtype='uint32')
\t\telse:
\t\t\tcoverage = np.multiply(coverage, rpbs)
\t\t\tcoverage = np.median(coverage, axis=(2))
\t\t\tcoverage = coverage.astype('float32')
\t\t\tcoverage = np.mean(coverage, axis=(1))
\t\tcov_d[chrom] = coverage
\t\tsummed = np.sum(pos_d[chrom], axis=(1,2,3), dtype='uint32')
\t\tout[chrom] = (coverage, summed)
\treturn out, pos_d, pos_size_d
'''
# Las lineas originales que tienen que DESAPARECER (lista propia del banco, no
# la del parche: si alguien borra un reemplazo de CAMBIOS, esto lo ve).
DEBEN_IRSE = [
    "pos_size_d[chrom] = [None] * chrom_length",
    "pos_d[chrom] = np.zeros(shape=(chrom_length, len(annotation_conditions), max_lib_counts, 2), dtype='uint32')",
    "pos_d[chrom][pos, condition_i, lib_i, strand_i] += 1",
    "pos_size_d[chrom][pos].append(length)",
    "coverage = np.sum(pos_d[chrom], axis=3, dtype='uint32')",
    "coverage = np.multiply(coverage, rpbs)",
    "coverage = np.median(coverage, axis=(2))",
    "coverage = coverage.astype('float32')",
    "coverage = np.mean(coverage, axis=(1))",
    "cov_d[chrom] = coverage",
    "summed = np.sum(pos_d[chrom], axis=(1,2,3), dtype='uint32')",
]
PARCHE = [sys.executable, str(RAIZ / "scripts" / "tradeoff_parche.py")]


def paquete(d, texto):
    pkg = d / "yasma"
    pkg.mkdir(parents=True, exist_ok=True)
    (pkg / "__init__.py").write_text("")
    f = pkg / "tradeoff.py"
    f.write_text(texto)
    return f


def correr_parche(f, *extra):
    r = subprocess.run(PARCHE + ["--ruta", str(f), *extra], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def cargar(d):
    for m in [m for m in sys.modules if m == "yasma" or m.startswith("yasma.")]:
        del sys.modules[m]
    sys.path.insert(0, str(d))
    try:
        return importlib.import_module("yasma.tradeoff")
    finally:
        sys.path.remove(str(d))


with tempfile.TemporaryDirectory() as d:
    d = pathlib.Path(d)
    f = paquete(d / "orig", FUENTE)
    fp = paquete(d / "parch", FUENTE)
    rc, txt = correr_parche(fp, "--verificar")
    chk("--verificar sin parche sale 1", rc == 1, txt)
    rc, txt = correr_parche(fp)
    chk("aplica", rc == 0 and "aplicado" in txt, txt)
    nuevo = fp.read_text()
    lineas = [l.strip() for l in nuevo.splitlines()]
    chk("las once líneas densas ya no están", all(x not in lineas for x in DEBEN_IRSE),
        [x for x in DEBEN_IRSE if x in lineas])
    chk("y el módulo disperso quedó al lado", (fp.parent / tradeoff_parche.MODULO).read_text()
        == (RAIZ / "scripts" / "tradeoff_disperso.py").read_text())
    chk("deja el respaldo con el original", (fp.parent / "tradeoff.py.sin_parche").read_text() == FUENTE)
    rc, txt = correr_parche(fp)
    chk("es idempotente", rc == 0 and "ya aplicado" in txt and fp.read_text() == nuevo, txt)
    rc, txt = correr_parche(fp, "--verificar")
    chk("--verificar con parche sale 0", rc == 0, txt)

    # El parchado tiene que correr y dar lo mismo que el original.
    mo, mp = cargar(d / "orig"), cargar(d / "parch")
    rnd = random.Random(5)
    crs = [("c1", 2000), ("c2", 120), ("c3", 900)]
    reads = []
    for _ in range(6000):
        c, L = rnd.choice(crs)
        reads.append((c, rnd.choice([0, L - 1, 50, 51, rnd.randrange(L)]), rnd.randrange(2),
                      rnd.randrange(3), rnd.randrange(2), rnd.choice([21, 22, 22, 24])))
    rpbs = np.array([[40, 0, 913], [7, 2**31, 1]], dtype="uint32")
    oo, po, so = mo.kernel(crs, ["a", "b"], 3, reads, rpbs)
    op, pp, sp = mp.kernel(crs, ["a", "b"], 3, reads, rpbs)
    igual = all(oo[c][0].dtype == op[c][0].dtype and oo[c][0].tobytes() == op[c][0].tobytes()
                and oo[c][1].tobytes() == op[c][1].tobytes() for c, _ in crs)
    chk("el parchado da la misma cobertura y suma por posición (incl. el cromosoma corto)", igual)
    igual = all(np.array_equal(po[c][a:b, ...], pp[c][a:b, ...])
                for c, L in crs for a, b in [(0, L), (40, 60), (L - 5, L + 5)])
    chk("y las mismas rebanadas de pos_d", igual)
    igual = all((so[c][p] is None) == (sp[c][p] is None) and
                list(collections.Counter(so[c][p] or []).items())
                == list(collections.Counter(sp[c][p] or []).items())
                for c, L in crs for p in range(L))
    chk("y los mismos largos por posición, en el mismo orden", igual)

    # Un modulo viejo al lado del tradeoff.py parchado: --verificar lo ve.
    (fp.parent / tradeoff_parche.MODULO).write_text("# viejo\n")
    rc, txt = correr_parche(fp, "--verificar")
    chk("con el módulo desactualizado, --verificar sale 1", rc == 1 and "no coincide" in txt, txt)
    rc, txt = correr_parche(fp)
    chk("y aplicar lo repone", rc == 0 and (fp.parent / tradeoff_parche.MODULO).read_text()
        == (RAIZ / "scripts" / "tradeoff_disperso.py").read_text(), txt)
    (fp.parent / tradeoff_parche.MODULO).unlink()
    rc, txt = correr_parche(fp, "--verificar")
    chk("sin el módulo, --verificar sale 1", rc == 1, txt)

    # Que se NIEGUE si el fuente cambio.
    for nombre, texto in [
        ("falta una línea", FUENTE.replace("\t\tcov_d[chrom] = coverage\n", "")),
        ("una línea aparece dos veces", FUENTE.replace(
            "\t\tout[chrom] = (coverage, summed)\n",
            "\t\tsummed = np.sum(pos_d[chrom], axis=(1,2,3), dtype='uint32')\n"
            "\t\tout[chrom] = (coverage, summed)\n")),
        ("cambió la reserva de pos_d", FUENTE.replace("max_lib_counts, 2), dtype='uint32')",
                                                       "max_lib_counts, 2), dtype='uint16')")),
    ]:
        g = paquete(d / nombre.replace(" ", "_"), texto)
        rc, txt = correr_parche(g)
        chk(f"se niega si {nombre}", rc == 2 and "esperaba 1 vez" in txt
            and g.read_text() == texto and not (g.parent / tradeoff_parche.MODULO).exists(), txt)

print("== 4. tradeoff_control: compara todo salvo la hora y las rutas")
with tempfile.TemporaryDirectory() as d:
    raiz = pathlib.Path(d)
    for nombre in ("original", "parchado"):
        t = raiz / nombre / "p" / "tradeoff"
        t.mkdir(parents=True)
        (t / "loci.txt").write_text("Locus\tName\nchr1:1-50\tlocus_1\n")
        (t / "log.txt").write_text(f"Start time: {nombre}\nruta {raiz / nombre}/p/x.bam\n"
                                   f" elapsed: 0:00:0{len(nombre)}\nfin\n")
        (t / "params.json").write_text(f'{{"alignment_file": "{raiz / nombre}/p/x.bam"}}\n')
    ig, dif, faltan = tradeoff_control.comparar(raiz)
    chk("iguales salvo hora y rutas -> iguales", sorted(ig) == ["loci.txt", "log.txt", "params.json"]
        and not dif and not faltan, (ig, dif, faltan))
    (raiz / "parchado" / "p" / "tradeoff" / "loci.txt").write_text("Locus\tName\nchr1:1-51\tlocus_1\n")
    (raiz / "parchado" / "p" / "tradeoff" / "log.txt").write_text(
        f"Start time: x\nruta {raiz}/parchado/p/x.bam\n elapsed: 9\nOTRA COSA\n")
    (raiz / "parchado" / "p" / "tradeoff" / "reads.txt").write_text("x\n")
    ig, dif, faltan = tradeoff_control.comparar(raiz)
    chk("un loci.txt distinto se marca", "loci.txt" in dif, dif)
    chk("un log.txt distinto en algo más que la hora se marca", "log.txt" in dif, dif)
    chk("un fichero que está en una sola se marca", faltan == ["reads.txt"], faltan)

print("== 4b. tradeoff_control: el recorte del BAM y el aviso de que sigue vivo")
# Sobre los 32 M reads de sclsc/duplicado eran dos tradeoff completos, ~1 h sin
# imprimir nada: en Colab se leyo como una celda colgada.
try:
    import pysam  # noqa: F401
    HAY_PYSAM = True
except ImportError:
    HAY_PYSAM = False
if not HAY_PYSAM:
    print("  (salteado el recorte: sin pysam)")
else:
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        bam = d / "x.bam"
        subprocess.run([sys.executable, str(RAIZ / "tests" / "bam_sintetico.py"), str(bam)],
                       check=True, capture_output=True)
        with pysam.AlignmentFile(str(bam)) as h:
            por_ref = {r: h.count(r) for r in h.references}
            refs, rgs = list(h.references), [g["ID"] for g in h.header.to_dict()["RG"]]
        n = tradeoff_control.recortar_bam(bam, d / "c.bam", 1000)
        with pysam.AlignmentFile(str(d / "c.bam")) as h:
            vistos = {r.reference_name for r in h.fetch()}
            chk("con 1000 entra el primer cromosoma ENTERO y nada más",
                n == por_ref[refs[0]] and vistos == {refs[0]}, (n, vistos))
            chk("la cabecera va entera: secuencias y @RG",
                list(h.references) == refs and [g["ID"] for g in h.header.to_dict()["RG"]] == rgs)
        chk("y queda indexado", (d / "c.bam.bai").exists())
        n2 = tradeoff_control.recortar_bam(bam, d / "c2.bam", 10**9)
        chk("con un tope enorme, todos los mapeados", n2 == sum(por_ref.values()), (n2, por_ref))

with tempfile.TemporaryDirectory() as d:
    sal = pathlib.Path(d) / "salida.txt"
    with open(sal, "w") as h:
        pr = subprocess.Popen([sys.executable, "-c",
                               "import sys, time\n"
                               "print('encoding reads ... 43%', flush=True)\n"
                               "time.sleep(1.5)\nsys.exit(3)"], stdout=h)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        st, ru = tradeoff_control.esperar(pr, sal, cada=0.3)
    chk("espera hasta que termina y devuelve su estado", os.waitstatus_to_exitcode(st) == 3, st)
    chk("mientras tanto avisa que sigue, con la última línea de yasma",
        "min  encoding reads ... 43%" in buf.getvalue(), buf.getvalue())
    sal.write_text("a\rencoding 99%\x1b[1A\x1b[2K\n\n")
    chk("la última línea sin los códigos de terminal",
        tradeoff_control.ultima_linea(sal) == "encoding 99%", tradeoff_control.ultima_linea(sal))

print("== 5. original contra parchado con yasma de verdad")
py = os.environ.get("YASMA_PYTHON")
if not py:
    print("  (salteado: definí YASMA_PYTHON=<python con yasma, pysam y ViennaRNA> para correrlo)")
else:
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        gen = RAIZ / "tests" / "bam_sintetico.py"
        ctl = RAIZ / "scripts" / "tradeoff_control.py"
        for nombre, args_bam, extra, env in [
            ("4 librerías", ["--libs", "4"], [], {}),
            ("3 librerías, compactación forzada", ["--libs", "3", "--seed", "2"], [],
             {"TESIS_DISPERSO_COMPACTAR": "977", "TESIS_DISPERSO_BLOQUE_B": "4096"}),
            ("condiciones desiguales", ["--libs", "4", "--seed", "3"],
             ["--", "-c", "SRR1000:WT", "SRR1001:WT", "SRR1002:WT", "SRR1003:mut"], {}),
        ]:
            bam = d / f"{len(nombre)}.bam"
            subprocess.run([py, str(gen), str(bam), *args_bam], check=True, capture_output=True)
            r = subprocess.run([py, str(ctl), str(bam), *extra], capture_output=True, text=True,
                               env=dict(os.environ, **env))
            chk(f"{nombre}: iguales byte a byte", r.returncode == 0 and "IGUALES" in r.stdout,
                (r.stdout + r.stderr)[-400:])

print()
print("todo bien" if not FALLAS else f"{FALLAS} chequeo(s) MAL")
sys.exit(1 if FALLAS else 0)
