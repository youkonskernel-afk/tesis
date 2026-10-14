#!/usr/bin/env python3
"""Corre `yasma tradeoff` original y parchado sobre el MISMO BAM y compara todo.

    ./scripts/tradeoff_control.py <bam> [--trabajo DIR] [-- <args extra de tradeoff>]

Es la prueba de que el parche de memoria (scripts/tradeoff_parche.py) no cambia
la anotacion: arma dos copias del paquete yasma instalado —una con el
tradeoff.py original, sacado del respaldo `.py.sin_parche` si el instalado ya
esta parchado, y otra parchada— corre las dos con los parametros de
`matriz.TRADEOFF` y compara byte a byte cada fichero que escriben en tradeoff/.
Solo se normalizan las lineas con la hora y la duracion en log.txt y las rutas
de trabajo en log.txt y params.json.

Imprime la memoria pico de cada corrida, que es el otro dato que importa.
El original tiene que entrar en RAM: correlo sobre un proyecto chico
(sclsc/duplicado son 39 Mb y 2 librerias). Sale 0 si son iguales, 1 si no,
2 si no se pudo correr.
"""
import argparse
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

AQUI = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
import tradeoff_parche  # noqa: E402

# Lineas de log.txt que cambian de una corrida a otra aunque la anotacion sea igual.
VOLATILES = ("Start time", "elapsed", "Run completed")


def tradeoff_args():
    """Los parametros con los que anota el proyecto. Sale de matriz.py y no de
    una copia: si la lista cambia alla, el control mide la que se usa."""
    import matriz  # noqa: PLC0415
    return list(matriz.TRADEOFF)


def paquete_instalado():
    import importlib.util  # noqa: PLC0415
    spec = importlib.util.find_spec("yasma")
    if spec is None or not spec.submodule_search_locations:
        raise SystemExit("no encontré el paquete yasma instalado")
    return pathlib.Path(list(spec.submodule_search_locations)[0])


def armar_copias(pkg, raiz):
    """-> {'original': raiz/original, 'parchado': raiz/parchado}, cada una con un
    paquete yasma adentro."""
    texto = (pkg / "tradeoff.py").read_text()
    if tradeoff_parche.MARCA in texto:
        respaldo = pkg / "tradeoff.py.sin_parche"
        if not respaldo.is_file():
            raise SystemExit(f"el tradeoff.py instalado está parchado y no está {respaldo.name}: "
                             "no tengo el original para comparar")
        original = respaldo.read_text()
    else:
        original = texto
    out = {}
    for nombre in ("original", "parchado"):
        dst = raiz / nombre / "yasma"
        shutil.copytree(pkg, dst, ignore=shutil.ignore_patterns("__pycache__"))
        (dst / "tradeoff.py").write_text(original)
        (dst / tradeoff_parche.MODULO).unlink(missing_ok=True)
        out[nombre] = raiz / nombre
    f = out["parchado"] / "yasma" / "tradeoff.py"
    f.write_text(tradeoff_parche.parchar(original))
    (f.parent / tradeoff_parche.MODULO).write_text(tradeoff_parche.FUENTE_MODULO.read_text())
    est, _ = tradeoff_parche.estado(f)
    if est != "aplicado":
        raise SystemExit(f"no pude parchar la copia: {est}")
    return out


def correr(pypath, bam, d, extra):
    """yasma tradeoff con el paquete de `pypath`. -> (rc, segundos, pico en bytes)"""
    d.mkdir(parents=True)
    # El BAM se copia: tradeoff escribe <bam>.depth.txt al lado, y no quiero que
    # la segunda corrida lea el de la primera.
    b = d / "x.bam"
    shutil.copy2(bam, b)
    for suf in (".bai", ".csi"):
        idx = pathlib.Path(str(bam) + suf)
        if idx.exists():
            shutil.copy2(idx, pathlib.Path(str(b) + suf))
    env = dict(os.environ, PYTHONPATH=str(pypath) + os.pathsep + os.environ.get("PYTHONPATH", ""))
    cmd = [sys.executable, "-m", "yasma", "tradeoff", "-a", str(b), "-o", str(d),
           *tradeoff_args(), *extra]
    t = time.time()
    with open(d / "salida.txt", "w") as h:
        p = subprocess.Popen(cmd, cwd=d, env=env, stdout=h, stderr=subprocess.STDOUT)
        _, st, ru = os.wait4(p.pid, 0)
    rc = os.waitstatus_to_exitcode(st)
    return rc, time.time() - t, ru.ru_maxrss * 1024


def normalizar(texto, raiz, nombre):
    texto = texto.replace(str(raiz / nombre), "<TRABAJO>")
    return [l for l in texto.splitlines() if not any(v in l for v in VOLATILES)]


def comparar(raiz):
    """-> (iguales, distintos, faltan) con los nombres de fichero."""
    a, b = raiz / "original" / "p" / "tradeoff", raiz / "parchado" / "p" / "tradeoff"
    fa = {f.name for f in a.iterdir()} if a.is_dir() else set()
    fb = {f.name for f in b.iterdir()} if b.is_dir() else set()
    iguales, distintos = [], []
    for n in sorted(fa & fb):
        if n in ("log.txt", "params.json"):
            ok = (normalizar((a / n).read_text(errors="replace"), raiz, "original")
                  == normalizar((b / n).read_text(errors="replace"), raiz, "parchado"))
        else:
            ok = (a / n).read_bytes() == (b / n).read_bytes()
        (iguales if ok else distintos).append(n)
    return iguales, distintos, sorted(fa ^ fb)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("--trabajo", help="directorio de trabajo (por defecto, uno temporal)")
    ap.add_argument("extra", nargs="*", help="argumentos extra para tradeoff, después de --")
    a = ap.parse_args()
    bam = pathlib.Path(a.bam).resolve()
    if not bam.is_file():
        print(f"ERROR: no existe {bam}", file=sys.stderr)
        return 2

    tmp = None
    if a.trabajo:
        raiz = pathlib.Path(a.trabajo).resolve()
        if raiz.exists():
            shutil.rmtree(raiz)
        raiz.mkdir(parents=True)
    else:
        tmp = tempfile.TemporaryDirectory()
        raiz = pathlib.Path(tmp.name)
    try:
        copias = armar_copias(paquete_instalado(), raiz)
        res = {}
        for nombre, pypath in copias.items():
            print(f"== {nombre}: yasma tradeoff sobre {bam.name}", flush=True)
            res[nombre] = correr(pypath, bam, raiz / nombre / "p", a.extra)
            rc, s, pico = res[nombre]
            print(f"   rc={rc}  {s/60:.1f} min  memoria pico {pico/1e9:.2f} GB", flush=True)
            if rc != 0:
                print(f"   salió con {rc}; últimas líneas:", file=sys.stderr)
                print("\n".join((raiz / nombre / "p" / "salida.txt").read_text(
                    errors="replace").splitlines()[-15:]), file=sys.stderr)
                return 2
        iguales, distintos, faltan = comparar(raiz)
        print(f"\n{len(iguales)} ficheros iguales: {', '.join(iguales)}")
        for n in distintos:
            print(f"DIFIERE  {n}")
        for n in faltan:
            print(f"FALTA en una de las dos  {n}")
        po, pp = res["original"][2], res["parchado"][2]
        print(f"\nmemoria pico: original {po/1e9:.2f} GB, parchado {pp/1e9:.2f} GB "
              f"({po/max(pp, 1):.1f}x)")
        if distintos or faltan or not iguales:
            print("\n>>> EL PARCHE CAMBIA LA SALIDA: no lo uses hasta entender por qué")
            return 1
        print("\n>>> IGUALES: el parche no cambia la anotación de este BAM")
        return 0
    finally:
        if tmp:
            tmp.cleanup()


if __name__ == "__main__":
    sys.exit(main())
