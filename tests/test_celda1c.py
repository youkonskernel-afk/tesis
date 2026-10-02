#!/usr/bin/env python3
"""§1c y §2 de 20_alinear con varias máquinas: la TANDA y los claims.

Se corren las celdas DE VERDAD con exec, contra scripts/reparto.py real y un
Drive de mentira en un tmpdir. Lo único que se toca del fuente es la línea
`TANDA = None`, que es donde la persona pega la tanda que le dio reparto.py.

Hasta acá ningún banco las cubría, y no se sabía porque mutar.py no copiaba
data/ y todo salía "detectado". Lo que tienen que garantizar:

  1. Lo que ya está en Drive no se vuelve a proponer (horas tiradas).
  2. Lo que otra máquina tiene tomado tampoco.
  3. §2 toma el claim ANTES de gastar horas, y se niega si otra lo tiene.
"""
import contextlib
import io
import json
import pathlib
import sys
import tempfile

RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB = RAIZ / "notebooks" / "20_alinear.ipynb"
sys.path.insert(0, str(RAIZ / "scripts"))
import reparto  # noqa: E402

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


def celda(marca):
    hall = ["".join(c["source"]) for c in json.loads(NB.read_text())["cells"]
            if c["cell_type"] == "code" and marca in "".join(c["source"])]
    assert len(hall) == 1, f"{marca!r} aparece {len(hall)} veces en {NB.name}"
    return hall[0]


C1C = celda("import reparto")
C2 = celda("# El de la TANDA si")
assert "TANDA = None\n" in C1C, "§1c ya no tiene la línea TANDA = None"


def correr_1c(drive, tanda):
    ns = {"CLON": RAIZ, "DRIVE": drive}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(C1C.replace("TANDA = None\n", f"TANDA = {tanda!r}\n", 1), ns)
    return buf.getvalue(), ns


def correr_2(ns, tmp):
    genomes = tmp / "genomes"
    genomes.mkdir(exist_ok=True)
    llamadas = []
    import shutil  # lo trae la celda de configuracion
    ns.setdefault("SIGUIENTE", None)
    ns.setdefault("PRUEBA", None)
    ns.update(GENOMES=genomes, shutil=shutil,
              correr=lambda *a, **k: llamadas.append(a) or 0)
    buf, err = io.StringIO(), None
    try:
        with contextlib.redirect_stdout(buf):
            exec(C2, ns)
    except Exception as e:  # noqa: BLE001
        err = e
    return buf.getvalue(), err, llamadas


with tempfile.TemporaryDirectory() as d:
    tmp = pathlib.Path(d)
    drive = tmp / "drive"
    for org in ("aa", "bb", "cc"):
        (drive / "70_genomas" / org).mkdir(parents=True)
    (drive / "10_bam" / "aa").mkdir(parents=True)
    (drive / "10_bam" / "aa" / "primario.bam").write_bytes(b"x")   # ya hecho
    claims = drive / "90_claims"
    reparto.tomar(claims, "bb/primario", "otra-maquina")             # tomado

    print("== 1. §1c saltea lo hecho y lo tomado por otra")
    out, ns = correr_1c(drive, ["aa/primario", "bb/primario", "cc/primario"])
    chk("lo hecho dice YA EN DRIVE", "aa/primario" in out and "YA EN DRIVE" in out, out)
    chk("lo tomado dice por quién", "tomado por otra-maquina" in out, out)
    chk("PROYECTO es el primero libre", ns.get("PROYECTO") == "cc/primario",
        ns.get("PROYECTO"))

    print("== 2. §2 toma el claim antes de empezar")
    out2, err, ll = correr_2(ns, tmp)
    chk("arranca", err is None, repr(err))
    chk("con el proyecto de la tanda", "proyecto: cc/primario" in out2, out2)
    chk("y el claim queda a su nombre",
        reparto.tomados(claims).get("cc/primario", ("",))[0] == ns["YO"],
        reparto.tomados(claims))

    print("== 3. §2 se niega si otra máquina lo tomó en el medio")
    # Entre §1c y §2 pueden pasar minutos: otra maquina lo pudo tomar.
    reparto.tomar(claims, "aa/duplicado", "otra-maquina")
    ns["PROYECTO"] = "aa/duplicado"
    out3, err3, ll3 = correr_2(ns, tmp)
    chk("revienta", isinstance(err3, RuntimeError), repr(err3))
    chk("diciendo que lo tiene otra", "otra máquina" in str(err3), str(err3))
    chk("sin preparar el genoma", ll3 == [], ll3)

    print("== 4. sin TANDA, §1c no propone nada y §2 no toma claims")
    out4, ns4 = correr_1c(drive, None)
    chk("PROYECTO queda en None", ns4.get("PROYECTO") is None, ns4.get("PROYECTO"))
    ns4["SIGUIENTE"] = "cc/duplicado"
    _, err4, _ = correr_2(ns4, tmp)
    chk("§2 usa SIGUIENTE", err4 is None and ns4.get("PROYECTO") == "cc/duplicado",
        (repr(err4), ns4.get("PROYECTO")))
    chk("y no deja claim", "cc/duplicado" not in reparto.tomados(claims),
        reparto.tomados(claims))

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
