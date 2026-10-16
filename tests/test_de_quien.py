#!/usr/bin/env python3
"""scripts/de_quien.py sobre librerías de mentira.

Paso con rhirr_duplicado: 99.9% sin alinear con recorte y genoma verificados.
La pregunta es si los reads son de una planta (el huésped) o si el pipeline
está mal, y la respuesta tiene que salir de los datos, no de un umbral que mire
otra cosa.

  1. Un hongo con miRNAs de planta en cantidad → HUÉSPED VEGETAL.
  2. Un hongo sin ellos → SIN PLANTA.
  3. Una planta con sus miRNAs → ESPERADO, no un falso huésped.
  4. Lee TODAS las librerías del ledger, PRE-TRIMMED incluidas, no solo la primera.
  5. Un isomiR (3' distinto) cuenta igual: el prefijo es 5'.
"""
import contextlib
import re
import gzip
import io
import pathlib
import sys
import tempfile

RAIZ = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "scripts"))
import de_quien  # noqa: E402

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


MIR166 = "TCGGACCAGGCTTCATTCCCC"
MIR156 = "TGACAGAAGAGAGTGAGCAC"
HONGO = "ACGTTGCAAGGTCAATCGAAT"


def fq(seqs):
    return "".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(seqs))


def proyecto(tmp, libs):
    """libs: (run, fichero relativo, secuencias)."""
    pdir = tmp / "x_dup"
    lin = ["run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc"]
    for run, rel, seqs in libs:
        p = pdir / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(p, "wt") as h:
            h.write(fq(seqs))
        lin.append(f"{run}\tPRJ\t{rel}\t-\t-\tPRE-TRIMMED\tT")
    (pdir / "recortadas.tsv").write_text("\n".join(lin) + "\n")
    return pdir


def correr(pdir, org, n=1000):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = de_quien.main([str(pdir), org, "--n", str(n)])
    return rc, buf.getvalue()


with tempfile.TemporaryDirectory() as d:
    print("== 1. hongo con miRNAs de planta: huésped vegetal")
    tmp = pathlib.Path(d) / "a"
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz", [MIR166] * 30 + [HONGO] * 70)])
    rc, out = correr(pdir, "rhirr")
    chk("sale 0", rc == 0, rc)
    chk("dice HUÉSPED VEGETAL", ">>> HUÉSPED VEGETAL" in out, out[-300:])
    chk("miR166 al 30%", re.search(r"miR166\s+30\.00%", out), out)
    chk("y marca la secuencia en el top", f"{MIR166}  miR166" in out, out)

    print("== 2. hongo sin miRNAs de planta")
    tmp = pathlib.Path(d) / "b"
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz", [HONGO] * 100)])
    rc, out = correr(pdir, "rhirr")
    chk("dice SIN PLANTA", ">>> SIN PLANTA" in out, out[-300:])
    chk("y no HUÉSPED", "HUÉSPED" not in out, out[-300:])

    print("== 3. una planta con sus miRNAs no es un huésped")
    tmp = pathlib.Path(d) / "c"
    pdir = proyecto(tmp, [("P1", "trim/P1.t.fq.gz", [MIR166] * 50 + [MIR156] * 50)])
    rc, out = correr(pdir, "prupe")
    chk("dice ESPERADO", ">>> ESPERADO" in out, out[-300:])
    chk("y no HUÉSPED", "HUÉSPED" not in out, out[-300:])

    print("== 4. todas las librerías del ledger, PRE-TRIMMED incluida")
    tmp = pathlib.Path(d) / "d"
    # La primera es puro hongo; la planta está en la segunda, que es PRE-TRIMMED
    # (fastq sin recortar). Mirando solo la primera saldría SIN PLANTA.
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz", [HONGO] * 100),
                          ("R2", "untrimmed/R2.fastq.gz", [MIR166] * 100)])
    rc, out = correr(pdir, "rhirr")
    chk("dice HUÉSPED VEGETAL", ">>> HUÉSPED VEGETAL" in out, out[-300:])
    chk("miR166 al 50% del total", re.search(r"miR166\s+50\.00%", out), out)
    chk("muestra cada librería", "R1" in out and "R2" in out, out)
    chk("R2 al 100%, R1 al 0%",
        re.search(r"R2\s+100\.00%", out) and re.search(r"R1\s+0\.00%", out), out)

    print("== 5. el tope de lecturas es por librería")
    tmp = pathlib.Path(d) / "e"
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz", [HONGO] * 500),
                          ("R2", "trim/R2.t.fq.gz", [MIR166] * 500)])
    rc, out = correr(pdir, "rhirr", n=100)
    chk("200 muestreadas (100 + 100)", "200 lecturas" in out, out[:200])

    print("== 6. un isomiR con otro 3' cuenta igual")
    tmp = pathlib.Path(d) / "f"
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz",
                           [MIR166[:-2]] * 10 + [MIR166 + "A"] * 10 + [HONGO] * 80)])
    rc, out = correr(pdir, "rhirr")
    chk("miR166 al 20%", re.search(r"miR166\s+20\.00%", out), out)

    print("== 7. si falta un fichero del ledger, lo dice")
    tmp = pathlib.Path(d) / "g"
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz", [HONGO] * 10)])
    (pdir / "trim" / "R1.t.fq.gz").unlink()
    try:
        correr(pdir, "rhirr")
        e = None
    except SystemExit as x:
        e = x
    chk("corta nombrando la corrida", e is not None and "R1" in str(e.code), repr(e))

print("== 8. la celda §5b del notebook llama al script con el proyecto de §2")
import json  # noqa: E402
NB = RAIZ / "notebooks" / "20_alinear.ipynb"
CELDA = [("".join(c["source"])) for c in json.loads(NB.read_text())["cells"]
         if c["cell_type"] == "code" and "'de_quien.py'" in "".join(c["source"])]
chk("hay una sola celda que lo llama", len(CELDA) == 1, len(CELDA))
llamadas = []
ns = {"PROYECTO": "rhirr/duplicado", "PROY_DIR": pathlib.Path("/content/proyectos"),
      "CLON": RAIZ, "correr": lambda *a, **k: llamadas.append(a) or 0}
exec(CELDA[0], ns)
chk("con el dir del proyecto y el org",
    llamadas == [("de_quien.py", "/content/proyectos/rhirr_duplicado", "rhirr")], llamadas)
for falta in ("PROYECTO", "correr"):
    try:
        exec(CELDA[0], {k: v for k, v in ns.items() if k != falta})
        e = None
    except Exception as x:  # noqa: BLE001
        e = x
    chk(f"sin {falta}: dice cuál falta, no NameError",
        isinstance(e, RuntimeError) and falta in str(e), repr(e))
with tempfile.TemporaryDirectory() as d:
    try:
        exec(CELDA[0], dict(ns, CLON=pathlib.Path(d)))
        e = None
    except Exception as x:  # noqa: BLE001
        e = x
    chk("con un clon viejo manda a re-clonar", isinstance(e, RuntimeError)
        and "clonar" in str(e), repr(e))

print("== 9. extremos: un inserto con 4 nt aleatorios en el 3' (lo de rhirr_duplicado)")
import itertools  # noqa: E402
import os  # noqa: E402
import random  # noqa: E402
rnd = random.Random(7)
NUCLEO_R = "TGCTGAGATTAAGCCCGTGTTCTAAGATTTGT"   # el de verdad, del top 15 de rhirr
OTRO = "GATTCGAAGGTCAATCGAATCCGTAGCATGCA"


def cuatro():
    return "".join(rnd.choice("ACGT") for _ in range(4))


def con_4n(nucleo, n, cinco=False):
    return [(cuatro() if cinco else "") + nucleo + cuatro() for _ in range(n)]


with tempfile.TemporaryDirectory() as d:
    tmp = pathlib.Path(d) / "h"
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz", con_4n(NUCLEO_R, 600) + con_4n(OTRO, 400))])
    rc, out = correr(pdir, "rhirr")
    chk("el 3' sale ALEATORIO", re.search(r"3'\s+la variante más común se lleva\s+\d+\.\d%.*ALEATORIO", out), out)
    chk("el 5' no", not re.search(r"5'.*ALEATORIO", out), out)
    chk("sin índice: EXTREMOS ALEATORIOS, que nombra el 3'",
        ">>> EXTREMOS ALEATORIOS: los 4 nt del 3' varían" in out, out[-500:])
    chk("y pide --indice para confirmarlo", "--indice" in out, out[-400:])
    chk("el de planta sigue saliendo debajo", ">>> SIN PLANTA" in out, out[-300:])

    print("== 10. una librería normal, con isomiRs, no es ALEATORIO")
    # miR-like de 22 nt: la forma templada domina, y hay isomiRs 3' de mismo largo.
    MIR = "TGAGGTAGTAGGTTGTATAGTT"
    seqs = ([MIR] * 800 + [MIR[:-1] + "A"] * 100 + [MIR[:-2] + "TT"] * 50
            + [OTRO] * 300 + [OTRO[:-1] + "G"] * 30)
    tmp = pathlib.Path(d) / "i"
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz", seqs)])
    rc, out = correr(pdir, "rhirr")
    chk("ningún extremo ALEATORIO", "ALEATORIO" not in out, out)
    chk("ni veredicto de extremos", "EXTREMOS" not in out, out[-400:])

    print("== 11. con --indice: bowtie decide, y dice QUÉ recortar")
    # bowtie falso: alinea mucho solo si se le recorta lo que hay que recortar.
    # BUEN_3/BUEN_5 dicen qué recortes lo arreglan; lee el fq real y cuenta.
    bindir = pathlib.Path(d) / "bin"
    bindir.mkdir()
    (bindir / "bowtie").write_text(r"""#!/usr/bin/env bash
echo "$*" >> "$LOG_BT"
t5=0; t3=0; prev=""
for a in "$@"; do
  [[ $prev == --trim5 ]] && t5=$a; [[ $prev == --trim3 ]] && t3=$a; prev=$a
done
fq="${@: -1}"; n=$(( $(wc -l < "$fq") / 4 ))
ok=1
[[ ${BUEN_3:-0} -gt 0 && $t3 -lt $BUEN_3 ]] && ok=0
[[ ${BUEN_5:-0} -gt 0 && $t5 -lt $BUEN_5 ]] && ok=0
[[ ${BUEN_3:-0} -eq 0 && ${BUEN_5:-0} -eq 0 ]] && ok=0
if [[ $ok == 1 ]]; then k=$(( n * 85 / 100 + (t5 + t3) * n / 200 )); else k=$(( n / 500 + (t5 + t3) * n / 400 )); fi
for ((i=0; i<k; i++)); do echo "m$i	+	chr1	1	ACGT"; done
""")
    (bindir / "bowtie").chmod(0o755)
    viejo_path = os.environ["PATH"]
    os.environ["PATH"] = f"{bindir}:{viejo_path}"
    os.environ["LOG_BT"] = str(pathlib.Path(d) / "bt.log")
    tmp = pathlib.Path(d) / "j"
    pdir = proyecto(tmp, [("R1", "trim/R1.t.fq.gz", con_4n(NUCLEO_R, 600) + con_4n(OTRO, 400)),
                          ("R2", "trim/R2.t.fq.gz", con_4n(OTRO, 1000))])

    def con_indice(**env):
        for k in ("BUEN_3", "BUEN_5"):
            os.environ.pop(k, None)
        os.environ.update({k: str(v) for k, v in env.items()})
        pathlib.Path(os.environ["LOG_BT"]).write_text("")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            de_quien.main([str(pdir), "rhirr", "--n", "1000", "--n-bowtie", "300",
                           "--indice", "/x/GCF_X.1", "--cores", "2"])
        return buf.getvalue(), pathlib.Path(os.environ["LOG_BT"]).read_text()

    try:
        out, log = con_indice(BUEN_3=4)
        chk("prueba los cuatro recortes", len(log.splitlines()) == 4, log)
        chk("con -v 1 y -k 1, como yasma", all("-v 1" in l and "-k 1" in l for l in log.splitlines()), log)
        chk("contra el índice que se le pasa", all("/x/GCF_X.1" in l for l in log.splitlines()), log)
        chk("sobre 300 por librería (600)", "bowtie -v 1 sobre 600 lecturas" in out, out)
        chk("dice EXTREMOS SIN RECORTAR", ">>> EXTREMOS SIN RECORTAR" in out, out[-600:])
        chk("del 3' solo, aunque las dos puntas alinean un poco más",
            "sacando 4 nt del 3' alinea" in out, out[-600:])
        chk("con los dos porcentajes", re.search(r"alinea \d+\.\d% contra \d+\.\d%", out), out[-600:])

        out, log = con_indice(BUEN_3=4, BUEN_5=4)
        chk("con 4N en las dos puntas, dice las dos", "sacando 4 nt del 5' y 3' alinea" in out, out[-600:])

        out, log = con_indice()
        chk("si recortar no cambia nada, no culpa a los extremos",
            "EXTREMOS" not in out.split("bowtie -v 1")[1], out[-600:])
        chk("aunque la firma de extremos diga ALEATORIO", "ALEATORIO" in out, out)
    finally:
        os.environ["PATH"] = viejo_path

print("== 12. por corrida, y lo que va a data/extremos.tsv")
# Un bowtie falso que alinea de verdad: corta lo que le piden y busca la
# secuencia en un "genoma". Así cada corrida se comporta según lo que trae, y el
# reparto por corrida sale de los nombres, que es lo que se prueba.
import csv  # noqa: E402
import extremos  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    d = pathlib.Path(d)
    genoma = "".join(rnd.choice("ACGT") for _ in range(20000))
    (d / "genoma.txt").write_text(genoma)
    bindir = d / "bin"
    bindir.mkdir()
    (bindir / "bowtie").write_text(f"""#!{sys.executable}
import sys
a = sys.argv[1:]
t5 = int(a[a.index('--trim5') + 1]); t3 = int(a[a.index('--trim3') + 1])
g = open({str(d / "genoma.txt")!r}).read()
l = open(a[-1]).read().split('\\n')
for i in range(0, len(l) - 1, 4):
    q = l[i + 1][t5:len(l[i + 1]) - t3]
    if len(q) >= 15 and q in g:
        print(l[i][1:] + '\\t+\\tchr1\\t1\\t' + q)
""")
    (bindir / "bowtie").chmod(0o755)

    def nucleos(n):
        return [genoma[p:p + 22] for p in (rnd.randrange(0, len(genoma) - 22) for _ in range(n))]

    def lib4n(n):
        return [cuatro() + q + cuatro() for q in nucleos(n)]

    def manifest(runs_por_proy):
        f = d / "man.tsv"
        with f.open("w") as h:
            h.write("org\trun\tbioproject\trol\n")
            for proy, runs in runs_por_proy.items():
                for r in runs:
                    h.write(f"rhirr\t{r}\t{proy}\tduplicado\n")
        return f

    def proy_libs(nombre, libs):
        """libs: (run, bioproject, secuencias)"""
        pdir = d / nombre
        lin = ["run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc"]
        for run, proy, seqs in libs:
            q = pdir / "trim" / f"{run}.t.fq.gz"
            q.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(q, "wt") as h:
                h.write(fq(seqs))
            lin.append(f"{run}\t{proy}\ttrim/{run}.t.fq.gz\t1\t1\t100.0\tT")
        (pdir / "recortadas.tsv").write_text("\n".join(lin) + "\n")
        return pdir

    def prueba(pdir, man, tabla, n=400):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            de_quien.main([str(pdir), "rhirr", "--n", str(n), "--n-bowtie", str(n),
                           "--indice", "/x/idx", "--cores", "1",
                           "--registrar", str(tabla), "--manifest", str(man)])
        return buf.getvalue()

    viejo_path = os.environ["PATH"]
    os.environ["PATH"] = f"{bindir}:{viejo_path}"
    try:
        tabla = d / "extremos.tsv"
        man = manifest({"PRJ_R": ["R1", "R2", "R3"], "PRJ_O": ["O1"]})

        print("  -- 12a. las tres corridas del BioProject piden 4/4: UNA fila de proyecto")
        pdir = proy_libs("a", [("R1", "PRJ_R", lib4n(300)), ("R2", "PRJ_R", lib4n(300)),
                               ("R3", "PRJ_R", lib4n(300))])
        out = prueba(pdir, man, tabla)
        filas = extremos.leer(tabla)
        chk("una sola fila", len(filas) == 1, filas)
        f = filas[0] if filas else {}
        chk("del proyecto (run -)", (f.get("bioproject"), f.get("run")) == ("PRJ_R", "-"), f)
        chk("con 4 y 4", (f.get("cortar_5p"), f.get("cortar_3p")) == (4, 4), f)
        chk("y la evidencia: casi nada sin cortar, todo cortando",
            float(f.get("alin_sin_pct", 99)) < 5 and float(f.get("alin_con_pct", 0)) > 95, f)
        chk("sobre las 900 lecturas", str(f.get("lecturas")) == "900", f)
        chk("la tabla por corrida nombra las tres",
            all(re.search(rf"{r}\s+300 .*cortar 4/4", out) for r in ("R1", "R2", "R3")), out)
        chk("y lo que escribió", "rhirr PRJ_R -: cortar 4/4" in out, out[-500:])
        chk("trim.sh lo ve para cada corrida",
            extremos.requerido(filas, "rhirr", "PRJ_R", "R2") == (4, 4), filas)
        chk("la cabecera de la tabla queda",
            tabla.read_text().startswith("# Cuantos nt sacar"), tabla.read_text()[:80])

        print("  -- 12b. una corrida que no es 4N: filas por corrida, y avisa la mezcla")
        pdir = proy_libs("b", [("R1", "PRJ_R", lib4n(300)), ("R2", "PRJ_R", lib4n(300)),
                               ("R3", "PRJ_R", nucleos(300))])
        out = prueba(pdir, man, tabla)
        filas = extremos.leer(tabla)
        claves = sorted((f["run"], f["cortar_5p"], f["cortar_3p"]) for f in filas)
        chk("R1 y R2 con 4/4, R3 sin fila, y la de proyecto se fue",
            claves == [("R1", 4, 4), ("R2", 4, 4)], claves)
        chk("R3 no se corta", extremos.requerido(filas, "rhirr", "PRJ_R", "R3") == (0, 0), filas)
        chk("R3 dice nada", re.search(r"R3\s+300 .*nada", out), out)
        chk("avisa que mezcla kits", "mezcla kits" in out, out)

        print("  -- 12c. medir algunas corridas no borra lo de las otras")
        extremos.registrar([extremos.fila("rhirr", "PRJ_R", "-", 4, 4, 0.1, 90, 900, "t")],
                           {("rhirr", "PRJ_R")}, tabla=tabla)
        pdir = proy_libs("c", [("R3", "PRJ_R", nucleos(300))])
        out = prueba(pdir, man, tabla)
        filas = extremos.leer(tabla)
        chk("la fila del proyecto sigue",
            any(f["run"] == "-" and f["cortar_5p"] == 4 for f in filas), filas)
        chk("y R3 lleva un 0/0 explícito, porque la del proyecto la cortaría",
            extremos.requerido(filas, "rhirr", "PRJ_R", "R3") == (0, 0), filas)
        chk("R1 sigue cortándose", extremos.requerido(filas, "rhirr", "PRJ_R", "R1") == (4, 4), filas)
        pdir = proy_libs("c2", [("R1", "PRJ_R", lib4n(300))])
        n_antes = len(extremos.leer(tabla))
        out = prueba(pdir, man, tabla)
        chk("medir una que ya pide lo que dice la tabla no agrega nada",
            len(extremos.leer(tabla)) == n_antes and "no hay nada que cambiar" in out,
            extremos.leer(tabla))

        print("  -- 12d. todo el BioProject sin 4N: las filas viejas se van")
        pdir = proy_libs("d", [(r, "PRJ_R", nucleos(300)) for r in ("R1", "R2", "R3")])
        out = prueba(pdir, man, tabla)
        filas = extremos.leer(tabla)
        chk("no queda ninguna de PRJ_R", not [f for f in filas if f["bioproject"] == "PRJ_R"], filas)
        chk("y lo dice", "reemplazada" in out, out[-400:])

        print("  -- 12e. las filas de otro proyecto no se tocan")
        extremos.registrar([extremos.fila("rhirr", "PRJ_O", "-", 0, 4, 1, 80, 300, "t")],
                           tabla=tabla)
        pdir = proy_libs("e", [(r, "PRJ_R", lib4n(300)) for r in ("R1", "R2", "R3")])
        prueba(pdir, man, tabla)
        filas = extremos.leer(tabla)
        chk("PRJ_O sigue", extremos.requerido(filas, "rhirr", "PRJ_O", "O1") == (0, 4), filas)

        print("  -- 12f. una corrida con pocas lecturas no decide")
        tabla.unlink()
        pdir = proy_libs("f", [("R1", "PRJ_R", lib4n(300)), ("R2", "PRJ_R", lib4n(300)),
                               ("R3", "PRJ_R", lib4n(50))])
        out = prueba(pdir, man, tabla)
        filas = extremos.leer(tabla)
        chk("R3 dice pocas lecturas", re.search(r"R3\s+50 .*pocas lecturas", out), out)
        chk("y como no se midieron todas, van por corrida",
            sorted(f["run"] for f in filas) == ["R1", "R2"], filas)

        print("  -- 12g. --registrar sin --indice no registra nada")
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                de_quien.main([str(pdir), "rhirr", "--registrar", str(tabla)])
            e = None
        except SystemExit as x:
            e = x
        chk("corta", e is not None and e.code != 0, repr(e))
    finally:
        os.environ["PATH"] = viejo_path

print("== 13. la tabla se niega a un corte que no es de un kit")
with tempfile.TemporaryDirectory() as d:
    t = pathlib.Path(d) / "x.tsv"
    t.write_text("\t".join(extremos.COLS) + "\nrhirr\tP\t-\t4\t40\t0\t0\t0\tt\tF\n")
    try:
        extremos.leer(t)
        e = None
    except ValueError as x:
        e = x
    chk("40 nt no", e is not None and "fuera de" in str(e), repr(e))
    t.write_text("\t".join(extremos.COLS) + "\nrhirr\tP\t-\tcuatro\t4\t0\t0\t0\tt\tF\n")
    try:
        extremos.leer(t)
        e = None
    except ValueError as x:
        e = x
    chk("ni un texto", e is not None and "no entero" in str(e), repr(e))

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
