#!/usr/bin/env python3
"""Parche a `yasma/nativealign.py`: no levantar un bowtie por libreria en `over`.

EL BUG, EN EL FUENTE DE v1.1.1

`bowtie_generator(lib, mmap)` arma `bowtie_call` en tres ramas —`unique`,
`multi` y `over`— pero solo las dos primeras le agregan flags. Despues, sin
mirar `mmap`, hace el `Popen` (linea 324-336 de nativealign.py). O sea que en
`over` tambien arranca un bowtie: sin `-v`, sin `-m`, sin `-S`, sobre la
libreria entera.

Esa salida NO SE USA. La rama `over` lee `<RG>.max50.fq` del disco y construye
los `AlignedSegment` a mano; `p.stdout` no se toca nunca. Y como abajo dice
`if mmap != 'over': p.wait()`, el proceso tampoco se espera: se queda vivo,
bloqueado escribiendo a un pipe que nadie lee, con el indice de bowtie entero
en memoria. Uno por libreria, todos a la vez.

POR QUE IMPORTA ACA

`gadmo_duplicado` son 12 librerias contra un genoma de ~670 Mb: 12 indices
residentes encima de `unique_d` (8 bytes por base del genoma, medido) es mas
memoria que la que tiene la VM. Murio con `Killed` al 96.5%, en la etapa `over`,
que es justo cuando ya hay once orfanos vivos. `sclsc_duplicado` habia pasado
porque son 2 librerias y 39 Mb.
`galga_duplicado` son 95 librerias contra 1.05 Gb: no hay maquina donde entre.

QUE CAMBIA EN LA SALIDA: NADA. El proceso que se deja de levantar es el que no
se lee. El parche lo unico que hace es no llamar a Popen cuando mmap == 'over'.

COMO SE APLICA

    ./scripts/yasma_parche.py              aplica (idempotente)
    ./scripts/yasma_parche.py --verificar  sale 0 si esta, 1 si no
    ./scripts/yasma_parche.py --ruta <f>   sobre un fichero concreto (bancos)

`align.sh correr` lo exige antes de alinear: sin el, los proyectos grandes
mueren a las horas y el mensaje no dice por que.
"""
import argparse
import importlib.util
import pathlib
import sys

MARCA = "PARCHE_TESIS_OVER"

# El ancla es el Popen que hay que saltear. Tiene que aparecer UNA sola vez:
# si aparece cero, la version cambio y el parche ya no corresponde; si aparece
# mas de una, no se cual es.
ANCLA = '\t\tif ".gz" in lib.suffixes:\n'

NUEVO = (
    "\t\tif mmap == 'over':\n"
    f"\t\t\t# {MARCA}: en 'over' la salida de bowtie NO se lee (se lee el .maxN\n"
    "\t\t\t# del disco) y abajo no se hace p.wait(), asi que el proceso queda\n"
    "\t\t\t# vivo con el indice entero en RAM. Uno por libreria = OOM.\n"
    "\t\t\tp = None\n"
    '\t\telif ".gz" in lib.suffixes:\n'
)

# La condicion que vuelve seguro el parche: si en 'over' se esperara al proceso,
# saltearlo cambiaria el comportamiento. Esta linea es la que dice que no.
GUARDIA = "if mmap != 'over':"


def ruta_instalada():
    spec = importlib.util.find_spec("yasma")
    if spec is None or not spec.submodule_search_locations:
        return None
    return pathlib.Path(list(spec.submodule_search_locations)[0]) / "nativealign.py"


def estado(f):
    """-> ('aplicado'|'pendiente'|<motivo de que no se pueda>, texto)"""
    try:
        texto = f.read_text()
    except OSError as e:
        return f"no pude leer {f}: {e}", ""
    if MARCA in texto:
        return "aplicado", texto
    n = texto.count(ANCLA)
    if n != 1:
        return (f"esperaba 1 vez el ancla del Popen y encontré {n}: "
                "nativealign.py cambió, revisá el parche contra la versión nueva"), texto
    if GUARDIA not in texto:
        return ("no encontré `if mmap != 'over':` antes de p.wait(): el parche "
                "asume que en 'over' el proceso no se espera, y eso ya no vale"), texto
    return "pendiente", texto


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verificar", action="store_true",
                    help="no toca nada; sale 0 si el parche está aplicado")
    ap.add_argument("--ruta", help="nativealign.py a parchar (por defecto, el instalado)")
    args = ap.parse_args()

    f = pathlib.Path(args.ruta) if args.ruta else ruta_instalada()
    if f is None:
        print("ERROR: no encontré el paquete yasma instalado (ver docs/yasma.md)", file=sys.stderr)
        return 2
    if not f.is_file():
        print(f"ERROR: no existe {f}", file=sys.stderr)
        return 2

    est, texto = estado(f)

    if est == "aplicado":
        print(f"parche ya aplicado en {f}")
        return 0
    if args.verificar:
        print(f"parche NO aplicado en {f}"
              + ("" if est == "pendiente" else f"\n  {est}"), file=sys.stderr)
        print("  aplicalo con: ./scripts/yasma_parche.py", file=sys.stderr)
        return 1
    if est != "pendiente":
        print(f"ERROR: {est}", file=sys.stderr)
        return 2

    respaldo = f.with_suffix(".py.sin_parche")
    if not respaldo.exists():
        respaldo.write_text(texto)
    f.write_text(texto.replace(ANCLA, NUEVO, 1))

    est2, _ = estado(f)
    if est2 != "aplicado":
        f.write_text(texto)
        print(f"ERROR: el parche no quedó ({est2}); dejé el fichero como estaba", file=sys.stderr)
        return 2
    print(f"parche aplicado en {f}\n  respaldo: {respaldo}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
