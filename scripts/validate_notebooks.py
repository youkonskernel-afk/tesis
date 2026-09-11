#!/usr/bin/env python3
"""
Valida los notebooks del repo: JSON bien formado, estructura de celdas, y que
cada celda de codigo parsee como Python.

Reemplaza a la validacion que hacia gen_colab_notebook.py al generar, y cubre
todos los notebooks en vez de solo el generado.

    ./scripts/validate_notebooks.py
"""
import ast
import json
import pathlib
import sys

RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB_DIR = RAIZ / "notebooks"


def valida(path):
    fallos = []
    try:
        nb = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        return [f"JSON invalido: {e}"]

    if nb.get("nbformat") != 4:
        fallos.append(f"nbformat={nb.get('nbformat')}, se espera 4")
    celdas = nb.get("cells")
    if not isinstance(celdas, list) or not celdas:
        return fallos + ["sin celdas"]

    n_code = 0
    for i, c in enumerate(celdas):
        tipo = c.get("cell_type")
        if tipo not in ("code", "markdown"):
            fallos.append(f"celda {i}: cell_type '{tipo}'")
            continue
        if not isinstance(c.get("source"), list):
            fallos.append(f"celda {i}: 'source' tiene que ser lista de lineas")
            continue
        if tipo != "code":
            continue
        n_code += 1
        if "outputs" not in c or "execution_count" not in c:
            fallos.append(f"celda {i}: falta 'outputs' o 'execution_count'")
        lineas = c["source"]
        sin_salto = [j for j, ln in enumerate(lineas[:-1]) if not ln.endswith("\n")]
        if sin_salto:
            fallos.append(
                f"celda {i}: {len(sin_salto)} lineas de 'source' sin \\n final "
                f"(la primera, la {sin_salto[0]}). Jupyter las concatena sin "
                f"separador y el codigo queda pegado en una sola linea.")
        # nbformat concatena 'source' SIN separador: cada linea tiene que
        # terminar en \n. Unir con "\n" aca enmascararia justamente el bug de
        # lineas sin salto, que en Colab llega como codigo pegado en una linea.
        src = "".join(c["source"])
        # Las lineas de shell (! y %) no son Python: se neutralizan para el
        # parse, que es lo unico que se quiere verificar.
        limpio = "\n".join(
            "pass" if ln.lstrip().startswith(("!", "%")) else ln
            for ln in src.split("\n")
        )
        try:
            ast.parse(limpio)
        except SyntaxError as e:
            fallos.append(f"celda {i}: sintaxis Python, linea {e.lineno}: {e.msg}")

    if n_code == 0:
        fallos.append("no tiene celdas de codigo")
    return fallos


def main():
    nbs = sorted(NB_DIR.glob("*.ipynb"))
    if not nbs:
        print(f"no hay notebooks en {NB_DIR}", file=sys.stderr)
        return 1
    malos = 0
    for nb in nbs:
        f = valida(nb)
        if f:
            malos += 1
            print(f"[MAL] {nb.name}")
            for x in f:
                print(f"       {x}")
        else:
            n = len(json.loads(nb.read_text())["cells"])
            print(f"[OK ] {nb.name}  ({n} celdas)")
    print()
    print(f"{len(nbs)} notebooks, {malos} con problemas")
    return 1 if malos else 0


if __name__ == "__main__":
    sys.exit(main())
