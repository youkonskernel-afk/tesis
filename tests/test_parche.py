#!/usr/bin/env python3
"""El parche a yasma/nativealign.py: que corrija el bug y que se niegue si no.

El bug: `bowtie_generator` hace el Popen sin mirar `mmap`, asi que en la etapa
`over` levanta un bowtie por libreria cuya salida no lee y que no espera. Cada
uno se queda vivo con el indice del genoma en RAM. Con 12 librerias y 670 Mb eso
es un OOM (`gadmo_duplicado`, `Killed` al 96.5%); con 95 librerias y 1.05 Gb
(`galga_duplicado`) no hay maquina donde entre.

Lo que este banco protege:
 1. Que despues del parche el Popen NO se ejecute en 'over' — que es todo el
    punto— y que el fichero siga siendo Python valido.
 2. Que el parche no toque las otras dos etapas: si 'unique' o 'multi' dejaran
    de levantar bowtie, no habria alineamiento y el BAM saldria vacio.
 3. Que se NIEGUE cuando el fuente cambio. Un parche que aplica a ciegas sobre
    una version distinta es peor que no tenerlo: no falla ruidosamente.
 4. Que se niegue si desaparece `if mmap != 'over':`, que es la linea que
    vuelve seguro saltear el Popen. Si en 'over' se esperara al proceso,
    saltearlo SI cambiaria el comportamiento.
"""
import ast
import pathlib
import subprocess
import sys
import tempfile

RAIZ = pathlib.Path(__file__).resolve().parent.parent
PARCHE = RAIZ / "scripts" / "yasma_parche.py"

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


# Copiado literal de nativealign.py de YASMA v1.1.1, lineas 270-345 y 435-437,
# con tabs como en el original. Es la unica forma de que el banco afirme sobre
# el fuente real y no sobre una maqueta comoda.
FUENTE = '''\
from subprocess import Popen, PIPE

LEVANTADOS = []


class _P:
\tpass


def Popen(*a, **k):
\tLEVANTADOS.append(k.get('stage'))
\treturn _P()


def bowtie_generator(lib, mmap, errf=None):

\t\tbowtie_call = ['bowtie']

\t\tsuff = '.fq'

\t\tif mmap == 'unique':
\t\t\tbowtie_call += ['-v', '1', '-m', '1', '--best', '--strata']

\t\telif mmap == 'multi':
\t\t\tbowtie_call += ['-v', '1', '-m', '50', '-a', '--best', '--strata']


\t\tprint(f"stage: {mmap}", file=errf)

\t\tif ".gz" in lib.suffixes:

\t\t\tgzip = Popen(stage=[mmap, 'gzip'])

\t\t\tbowtie_call.append("-")
\t\t\tp = Popen(stage=[mmap, 'bowtie'])

\t\telse:
\t\t\tbowtie_call.append(str(lib))
\t\t\tp = Popen(stage=[mmap, 'bowtie'])



\t\tif mmap == 'over':
\t\t\treturn 'leido del .maxN, sin tocar p.stdout'

\t\tif mmap != 'over':
\t\t\tp.wait_llamado = True

\t\treturn 'de bowtie'
'''

ARRANQUE = '''
import sys

class L:
\tsuffixes = ['.t', '.fq', '.gz']


for etapa in ['unique', 'multi', 'over']:
\tbowtie_generator(L(), etapa, errf=sys.stderr)
print(sorted(x for x in LEVANTADOS if x))
'''


def correr(*args):
    r = subprocess.run([sys.executable, str(PARCHE), *args],
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def levantados(f):
    """Qué procesos arranca el generador, corriendo el fichero de verdad."""
    r = subprocess.run([sys.executable, "-c", f.read_text() + ARRANQUE],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None, r.stderr
    return ast.literal_eval(r.stdout.strip()), ""


with tempfile.TemporaryDirectory() as d:
    f = pathlib.Path(d) / "nativealign.py"

    print("== 1. sin parche, la etapa 'over' levanta un bowtie que nadie lee")
    f.write_text(FUENTE)
    antes, err = levantados(f)
    chk("el fixture corre", antes is not None, err)
    chk("y 'over' arranca bowtie (el bug)", ['over', 'bowtie'] in (antes or []), antes)

    print("== 2. --verificar no toca nada y lo dice")
    rc, out = correr("--verificar", "--ruta", str(f))
    chk("exit 1", rc == 1, rc)
    chk("dice que no está", "NO aplicado" in out, out)
    chk("y no escribió", f.read_text() == FUENTE)

    print("== 3. aplicado: 'over' ya no levanta nada")
    rc, out = correr("--ruta", str(f))
    chk("exit 0", rc == 0, out)
    despues, err = levantados(f)
    chk("sigue siendo Python válido", despues is not None, err)
    chk("'over' no arranca bowtie", ['over', 'bowtie'] not in (despues or []), despues)
    chk("ni el gzip que lo alimenta", ['over', 'gzip'] not in (despues or []), despues)

    print("== 4. y las otras dos etapas siguen igual")
    # Si 'unique' o 'multi' dejaran de levantar bowtie no habria alineamiento y
    # el BAM saldria vacio — el mismo fallo silencioso que el adaptador
    # equivocado, un paso mas arriba.
    for etapa in ("unique", "multi"):
        chk(f"{etapa} sigue levantando bowtie",
            [etapa, 'bowtie'] in (despues or []), despues)
    chk("nada más cambió",
        despues == [x for x in (antes or []) if x[0] != 'over'], despues)

    print("== 5. es idempotente y deja respaldo")
    copia = f.read_text()
    rc, out = correr("--ruta", str(f))
    chk("exit 0 la segunda vez", rc == 0, out)
    chk("y no cambió nada", f.read_text() == copia)
    chk("hay respaldo sin parche",
        (pathlib.Path(d) / "nativealign.py.sin_parche").read_text() == FUENTE)
    rc, out = correr("--verificar", "--ruta", str(f))
    chk("--verificar ahora sale 0", rc == 0, out)

    print("== 6. si el fuente cambió, se niega en vez de adivinar")
    g = pathlib.Path(d) / "otra.py"
    g.write_text(FUENTE.replace('if ".gz" in lib.suffixes:', 'if lib.esta_comprimido():'))
    rc, out = correr("--ruta", str(g))
    chk("exit != 0", rc != 0, rc)
    chk("dice qué pasó", "ancla" in out and "encontré 0" in out, out)
    chk("y no escribió", "PARCHE_TESIS_OVER" not in g.read_text())

    print("== 7. si el ancla aparece dos veces, tampoco")
    h = pathlib.Path(d) / "dos.py"
    h.write_text(FUENTE + '\n\t\tif ".gz" in lib.suffixes:\n\t\t\tpass\n')
    rc, out = correr("--ruta", str(h))
    chk("exit != 0", rc != 0, rc)
    chk("dice cuántas encontró", "encontré 2" in out, out)

    print("== 8. si ya no vale que en 'over' no se espere, se niega")
    # Saltear el Popen es seguro SOLO porque abajo dice `if mmap != 'over'`. Si
    # esa linea cambia, el parche deja de ser equivalente y hay que releerlo.
    i = pathlib.Path(d) / "espera.py"
    i.write_text(FUENTE.replace("if mmap != 'over':", "if True:"))
    rc, out = correr("--ruta", str(i))
    chk("exit != 0", rc != 0, rc)
    chk("nombra la suposición", "mmap != 'over'" in out, out)

    print("== 9. un fichero que no existe no es un parche aplicado")
    rc, out = correr("--verificar", "--ruta", str(pathlib.Path(d) / "no_existe.py"))
    chk("exit != 0", rc != 0, rc)
    chk("y lo dice", "no existe" in out, out)

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    raise SystemExit(1)
print("TODO OK")
