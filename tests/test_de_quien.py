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

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
