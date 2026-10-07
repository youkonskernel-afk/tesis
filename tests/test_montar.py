#!/usr/bin/env python3
"""La celda de montar Drive y los guardias de 20_alinear, sin Colab.

Se corre la celda DE VERDAD —la fuente del notebook, con exec— contra un
`google.colab.drive` falso que falla cuando se le pide. Paso en serio: un
`drive.mount` pelado dio 'mount failed' sin decir que hacer, la persona salto a
§1, y §1 revento con `NameError: CLON` — un error que apunta a la celda
equivocada.

  1. Si el primer montaje falla, se desmonta y se reintenta forzando.
  2. Si el segundo tambien, el error dice QUE hacer (cuenta, casillas), no un
     traceback de google.colab.
  3. Montar con otra cuenta se distingue de no montar.
  4. §1 y §2 sin las celdas de arriba dicen "Ejecutar anteriores", no NameError.
  5. La celda de montaje es la misma en todos los notebooks.
"""
import io
import json
import pathlib
import sys
import tempfile
import types
from contextlib import redirect_stdout

RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB = RAIZ / "notebooks" / "20_alinear.ipynb"
sys.path.insert(0, str(RAIZ / "scripts"))
import validate_notebooks  # noqa: E402

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


def celda(marca, nb=NB):
    hall = ["".join(c["source"]) for c in json.loads(nb.read_text())["cells"]
            if c["cell_type"] == "code" and marca in "".join(c["source"])]
    assert len(hall) == 1, f"{marca!r} aparece {len(hall)} veces en {nb.name}"
    return hall[0]


MONTAR = celda("drive.mount(")


def montar(fallas, tesis_existe=True):
    """Corre la celda con un drive falso que falla las primeras `fallas` veces.
    Solo se cambia la ruta de Drive, a un tmpdir: el resto es la fuente tal cual."""
    llamadas = []

    def mount(punto, force_remount=False, timeout_ms=None):
        llamadas.append(("mount", force_remount))
        if len([x for x in llamadas if x[0] == "mount"]) <= fallas:
            raise ValueError("mount failed")

    def flush_and_unmount():
        llamadas.append(("flush", None))

    colab = types.ModuleType("google.colab")
    colab.drive = types.SimpleNamespace(mount=mount, flush_and_unmount=flush_and_unmount)
    google = types.ModuleType("google")
    google.colab = colab
    viejos = {k: sys.modules.get(k) for k in ("google", "google.colab")}
    sys.modules.update({"google": google, "google.colab": colab})

    with tempfile.TemporaryDirectory() as d:
        tesis = pathlib.Path(d) / "MyDrive" / "tesis"
        if tesis_existe:
            tesis.mkdir(parents=True)
        fuente = MONTAR.replace("/content/drive/MyDrive/tesis", str(tesis))
        buf, err = io.StringIO(), None
        try:
            with redirect_stdout(buf):
                exec(fuente, {})
        except Exception as e:  # noqa: BLE001
            err = e
        finally:
            for k, v in viejos.items():
                if v is None:
                    sys.modules.pop(k, None)
                else:
                    sys.modules[k] = v
    return llamadas, buf.getvalue(), err


print("== 1. monta a la primera")
ll, out, err = montar(0)
chk("sin error", err is None, repr(err))
chk("un solo mount, sin forzar", ll == [("mount", False)], ll)
chk("dice Drive OK", "Drive OK" in out, out)

print("== 2. 'mount failed' una vez: desmonta y reintenta forzando")
ll, out, err = montar(1)
chk("sin error", err is None, repr(err))
chk("desmontó antes de reintentar", ("flush", None) in ll, ll)
chk("el segundo fuerza", ll[-1] == ("mount", True), ll)
chk("avisa que reintentó", "reintento" in out, out)
chk("y terminó montado", "Drive OK" in out, out)

print("== 3. falla las dos veces: dice qué hacer")
ll, out, err = montar(2)
chk("lanza RuntimeError", isinstance(err, RuntimeError), repr(err))
msg = str(err)
chk("nombra la cuenta", "seb.ugazm@gmail.com" in msg, msg)
chk("habla de las casillas", "casillas" in msg, msg)
chk("manda a borrar el entorno", "borrar entorno" in msg, msg)
chk("sin el traceback de google.colab encadenado",
    err is not None and err.__suppress_context__, repr(err))
chk("intentó dos veces, no más", [x for x in ll if x[0] == "mount"] ==
    [("mount", False), ("mount", True)], ll)

print("== 4. monta, pero con otra cuenta: no ve tesis/")
ll, out, err = montar(0, tesis_existe=False)
chk("lanza AssertionError", isinstance(err, AssertionError), repr(err))
chk("dice que es la cuenta", "otra cuenta" in str(err), str(err))

print("== 5. §1 y §2 sin las celdas de arriba: Ejecutar anteriores, no NameError")
for nombre, marca in (("§1", "B_FQGZ = 22"), ("§2", "# El de la TANDA si"),
                      ("§8", "cola.correr_cola(")):
    try:
        exec(celda(marca), {})
        e = None
    except Exception as x:  # noqa: BLE001
        e = x
    chk(f"{nombre}: RuntimeError", isinstance(e, RuntimeError), repr(e))
    chk(f"{nombre}: dice Ejecutar anteriores", "Ejecutar anteriores" in str(e), str(e))
    chk(f"{nombre}: nombra lo que falta", "CLON" in str(e) or "correr" in str(e), str(e))

# §8 es la excepcion: parado ahi, "Ejecutar anteriores" corre tambien §1 a §7, y
# §3-§4 recortan y alinean un proyecto suelto durante horas antes de la cola.
# Paso dos veces: Preambulo corrido, Herramientas y Configuracion no. Ahora §8
# las corre sola, con la fuente del notebook del clon; el Preambulo no puede
# (montar pide autorizar en una ventana), asi que ese lo pide.
print("== 5b. §8 sin Preámbulo: lo pide, y NO Ejecutar anteriores")
COLA = celda("cola.correr_cola(")
try:
    exec(COLA, {})
    e = None
except Exception as x:  # noqa: BLE001
    e = x
chk("§8 sin nada: RuntimeError", isinstance(e, RuntimeError), repr(e))
chk("§8 sin nada: pide el Preámbulo", "Preámbulo" in str(e), str(e))
chk("§8 sin nada: advierte no usar Ejecutar anteriores",
    'NO uses "Ejecutar anteriores"' in str(e), str(e))

print("== 5c. §8 con Preámbulo: corre Herramientas y Configuración y sigue a la cola")
HERR = "YASMA_REF = 'x'\nORDEN.append('Herramientas')\n"
CONF = ("import os\nENV = dict(os.environ)\nORDEN.append('Configuración')\n"
        "GENOMES = PROY_DIR = BAM_DIR = CLON\n")
COLA_FALSA = ("LLAMADAS = []\n"
              "def correr_cola(**kw):\n    LLAMADAS.append(kw)\n    return {'ok': True}\n")


def clon_falso(d, celdas):
    (d / "notebooks").mkdir(parents=True)
    (d / "scripts").mkdir()
    (d / "notebooks" / "20_alinear.ipynb").write_text(json.dumps(
        {"cells": [{"cell_type": "code", "source": [c]} for c in celdas]}))
    (d / "scripts" / "cola.py").write_text(COLA_FALSA)
    (d / "scripts" / "reparto.py").write_text("")


def correr_cola_con(d, ns):
    viejos = {k: sys.modules.pop(k, None) for k in ("cola", "reparto")}
    path = list(sys.path)
    buf, err = io.StringIO(), None
    try:
        with redirect_stdout(buf):
            exec(COLA, ns)
    except Exception as x:  # noqa: BLE001
        err = x
    finally:
        sys.path[:] = path
        for k, v in viejos.items():
            sys.modules.pop(k, None)
            if v is not None:
                sys.modules[k] = v
    return buf.getvalue(), err


with tempfile.TemporaryDirectory() as d:
    d = pathlib.Path(d)
    # Configuracion ANTES que Herramientas en el notebook: el orden lo pone §8.
    clon_falso(d / "a", [CONF, "print('§1')", HERR])
    ns = {"CLON": d / "a", "DRIVE": d / "a", "ORDEN": []}
    out, err = correr_cola_con(d / "a", ns)
    chk("§8 con Preámbulo: sin error", err is None, repr(err))
    chk("corre Herramientas y después Configuración",
        ns["ORDEN"] == ["Herramientas", "Configuración"], ns["ORDEN"])
    chk("y no otras celdas (§1)", "§1" not in out, out)
    chk("avisa que las corrió", "la corro" in out, out)
    llam = ns.get("cola") and ns["cola"].LLAMADAS
    chk("llega a la cola con el ENV de Configuración",
        bool(llam) and llam[0]["env"] is ns.get("ENV"), llam)

    # Ya corridas: no las repite (Herramientas instala; no hace falta otra vuelta).
    ns2 = dict(ns, ORDEN=[])
    out, err = correr_cola_con(d / "a", ns2)
    chk("ya corridas: no las repite", err is None and ns2["ORDEN"] == [], (err, ns2["ORDEN"]))

    # Un notebook donde no esta la celda: error claro, no NameError mas abajo.
    clon_falso(d / "b", [HERR])
    ns3 = {"CLON": d / "b", "DRIVE": d / "b", "ORDEN": []}
    out, err = correr_cola_con(d / "b", ns3)
    chk("sin la celda de Configuración: RuntimeError que la nombra",
        isinstance(err, RuntimeError) and "Configuración" in str(err), repr(err))

# Y en el notebook DE VERDAD cada marca encuentra exactamente una celda: si no,
# en Colab §8 cortaria en vez de correrlas.
# Mismo filtro que §8: la propia celda de la cola lleva las marcas como texto.
_codigo = ["".join(c["source"]) for c in json.loads(NB.read_text())["cells"]
           if c["cell_type"] == "code"]
reales = [[c for c in _codigo if m in c and "correr_cola(" not in c]
          for m in ("YASMA_REF = ", "ENV = dict(os.environ)")]
chk("las marcas de §8 encuentran UNA celda cada una en el notebook real",
    [len(r) for r in reales] == [1, 1], [len(r) for r in reales])
chk("y son Herramientas y Configuración",
    all(reales) and "apt-get" in reales[0][0] and "PROY_DIR" in reales[1][0])

print("== 6. la celda de montaje es la misma en todos los notebooks")
todas = {nb.name: celda("drive.mount(", nb)
         for nb in sorted((RAIZ / "notebooks").glob("*.ipynb"))
         if "drive.mount(" in nb.read_text()}
chk("cinco notebooks montan", len(todas) == 5, sorted(todas))
chk("con la misma celda", len(set(todas.values())) == 1)

# Y validate_notebooks lo detecta si una deriva.
with tempfile.TemporaryDirectory() as d:
    d = pathlib.Path(d)
    for i, src in enumerate((MONTAR, "from google.colab import drive\ndrive.mount('/content/drive')\n")):
        (d / f"n{i}.ipynb").write_text(json.dumps(
            {"cells": [{"cell_type": "code", "source": [src]}]}))
    viejo, validate_notebooks.NB_DIR = validate_notebooks.NB_DIR, d
    try:
        fallos = []
        validate_notebooks.preambulos(fallos)
    finally:
        validate_notebooks.NB_DIR = viejo
chk("validate_notebooks grita si deriva", any("montaje" in f for f in fallos), fallos)

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
