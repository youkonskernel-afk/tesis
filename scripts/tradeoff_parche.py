#!/usr/bin/env python3
"""Parche a `yasma/tradeoff.py`: la misma anotacion sin arreglos densos por base.

EL PROBLEMA, EN EL FUENTE DE v1.1.1

`get_kernel_coverage()` reserva, por cada cromosoma y antes de leer un solo read:

    pos_size_d[chrom] = [None] * chrom_length                    # 8 B por base
    pos_d[chrom] = np.zeros((chrom_length, C, N, 2), 'uint32')   # 8*C*N B por base

con N librerias. Y `numpy_method()` guarda en `cov_d` la cobertura de todo el
genoma (4 B por base) sin leerla nunca despues. Medido en Colab: `-9` (SIGKILL,
sin memoria) en gadmo_duplicado, gadmo_primario, phypa_duplicado y
prupe_primario. Las cuentas, con los genomas y librerias del manifiesto:

    gadmo_duplicado    670 Mb x 12 lib    64 GB de pos_d + 10.7 fijos
    prupe_primario     228 Mb x  8 lib    15 GB          +  3.6
    galga_duplicado   1053 Mb x 95 lib   800 GB          + 16.8

`np.zeros` reserva perezoso, asi que lo residente es lo que se toca; pero cada
read toca una pagina distinta y con decenas de millones de reads alcanza.

QUE CAMBIA EN LA SALIDA: NADA. Cambia donde se guardan los conteos, no que se
calcula. `tradeoff_disperso.py` devuelve, para cada rebanada que YASMA pide, el
mismo arreglo uint32 que devolveria el denso, y para cada posicion la misma
lista de largos en el mismo orden; la cobertura sale de las mismas cinco
operaciones numpy, fila por fila. Se verifico corriendo el original y el
parchado sobre BAMs sinteticos y comparando byte a byte todo lo que escribe
(tests/test_tradeoff_parche.py). `cov_d` deja de guardarse: no se lee.

QUE NO ARREGLA

- `claim_d`: un dict por base DENTRO de las regiones. Es memoria proporcional a
  lo anotado (unos % del genoma), no al genoma.
- `ker_d`: el kernel de cobertura, 4 B por base del genoma, float32. Hace falta
  entero para buscar regiones despues de fijar el umbral global.
- El tiempo: los bucles de Python por base siguen igual.

COMO SE APLICA

    ./scripts/tradeoff_parche.py              aplica (idempotente)
    ./scripts/tradeoff_parche.py --verificar  sale 0 si esta, 1 si no
    ./scripts/tradeoff_parche.py --ruta <f>   sobre un tradeoff.py concreto (bancos)

Copia `tradeoff_disperso.py` al paquete como `_tesis_disperso.py` y reemplaza
once lineas de `tradeoff.py`. Se niega si alguna no aparece exactamente una vez:
un parche que aplica a ciegas sobre otra version es peor que no tenerlo.
`matriz.py anotar` lo exige antes de llamar a tradeoff.
"""
import argparse
import importlib.util
import pathlib
import sys

MARCA = "PARCHE_TESIS_MEMORIA"
MODULO = "_tesis_disperso.py"
FUENTE_MODULO = pathlib.Path(__file__).resolve().parent / "tradeoff_disperso.py"

# (linea original sin indentacion, reemplazo). Cada original tiene que aparecer
# UNA vez en tradeoff.py; se reemplaza conservando su indentacion.
CAMBIOS = [
    ("import numpy as np",
     "import numpy as np\nfrom . import _tesis_disperso as _td  # " + MARCA),
    ("pos_size_d[chrom] = [None] * chrom_length",
     "pos_size_d[chrom] = _td.Largos(chrom_length)  # " + MARCA),
    ("pos_d[chrom] = np.zeros(shape=(chrom_length, len(annotation_conditions), "
     "max_lib_counts, 2), dtype='uint32')",
     "pos_d[chrom] = _td.Conteos(chrom_length, len(annotation_conditions), "
     "max_lib_counts)  # " + MARCA),
    ("pos_d[chrom][pos, condition_i, lib_i, strand_i] += 1",
     "pos_d[chrom].sumar(pos, condition_i, lib_i, strand_i)  # " + MARCA),
    ("pos_size_d[chrom][pos].append(length)",
     "pos_size_d[chrom].agregar(pos, length)  # " + MARCA),
    # coverage.shape[0] sigue valiendo: Conteos tiene .shape = (L, C, N, 2)
    ("coverage = np.sum(pos_d[chrom], axis=3, dtype='uint32')",
     "coverage = pos_d[chrom]  # " + MARCA),
    ("coverage = np.multiply(coverage, rpbs)",
     "coverage = coverage.cobertura(rpbs)  # " + MARCA + ": suma, rpb, mediana, f32 y media"),
    ("coverage = np.median(coverage, axis=(2))", "pass  # " + MARCA),
    ("coverage = coverage.astype('float32')", "pass  # " + MARCA),
    ("coverage = np.mean(coverage, axis=(1))", "pass  # " + MARCA),
    ("cov_d[chrom] = coverage",
     "cov_d[chrom] = None  # " + MARCA + ": no se lee nunca y son 4 B por base"),
    ("summed = np.sum(pos_d[chrom], axis=(1,2,3), dtype='uint32')",
     "summed = pos_d[chrom].por_posicion()  # " + MARCA),
]


def ruta_instalada():
    spec = importlib.util.find_spec("yasma")
    if spec is None or not spec.submodule_search_locations:
        return None
    return pathlib.Path(list(spec.submodule_search_locations)[0]) / "tradeoff.py"


def _lineas_con(texto, original):
    return [i for i, l in enumerate(texto.split("\n")) if l.strip() == original]


def estado(f):
    """-> ('aplicado'|'pendiente'|<motivo de que no se pueda>, texto)"""
    try:
        texto = f.read_text()
    except OSError as e:
        return f"no pude leer {f}: {e}", ""
    if MARCA in texto:
        mod = f.parent / MODULO
        if not mod.is_file():
            return f"tradeoff.py está parchado pero falta {mod.name}", texto
        if mod.read_text() != FUENTE_MODULO.read_text():
            return (f"{mod.name} no coincide con scripts/tradeoff_disperso.py: "
                    "volvé a aplicar el parche"), texto
        return "aplicado", texto
    for original, _ in CAMBIOS:
        n = len(_lineas_con(texto, original))
        if n != 1:
            return (f"esperaba 1 vez `{original}` y encontré {n}: tradeoff.py "
                    "cambió, revisá el parche contra la versión nueva"), texto
    return "pendiente", texto


def parchar(texto):
    for original, nuevo in CAMBIOS:
        lineas = texto.split("\n")
        (i,) = _lineas_con(texto, original)
        sangria = lineas[i][: len(lineas[i]) - len(lineas[i].lstrip())]
        lineas[i] = "\n".join(sangria + n for n in nuevo.split("\n"))
        texto = "\n".join(lineas)
    return texto


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verificar", action="store_true",
                    help="no toca nada; sale 0 si el parche está aplicado y al día")
    ap.add_argument("--ruta", help="tradeoff.py a parchar (por defecto, el instalado)")
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
        print(f"parche de memoria ya aplicado en {f}")
        return 0
    if args.verificar:
        print(f"parche de memoria NO aplicado en {f}"
              + ("" if est == "pendiente" else f"\n  {est}"), file=sys.stderr)
        print("  aplicalo con: ./scripts/tradeoff_parche.py", file=sys.stderr)
        return 1
    mod = f.parent / MODULO
    if MARCA in texto:
        # Parchado pero con el modulo viejo o ausente: solo hay que copiarlo.
        mod.write_text(FUENTE_MODULO.read_text())
    else:
        if est != "pendiente":
            print(f"ERROR: {est}", file=sys.stderr)
            return 2
        respaldo = f.with_suffix(".py.sin_parche")
        if not respaldo.exists():
            respaldo.write_text(texto)
        mod.write_text(FUENTE_MODULO.read_text())
        f.write_text(parchar(texto))

    est2, _ = estado(f)
    if est2 != "aplicado":
        if MARCA not in texto:
            f.write_text(texto)
        print(f"ERROR: el parche no quedó ({est2}); dejé tradeoff.py como estaba",
              file=sys.stderr)
        return 2
    print(f"parche de memoria aplicado en {f}\n  módulo: {mod}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
