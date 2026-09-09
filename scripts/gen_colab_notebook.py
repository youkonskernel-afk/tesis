#!/usr/bin/env python3
"""
Genera notebooks/descarga_genomas.ipynb a partir de data/genomas.tsv.

El notebook corre en Google Colab, que sí tiene salida a NCBI y monta Drive
nativamente: el FASTA va de NCBI a Mi unidad/tesis/70_genomas/ sin pasar por
el disco local. Ver docs/colab.md.

Se genera en vez de escribirse a mano para que la tabla de candidatos no se
desincronice de data/genomas.tsv.

    ./scripts/gen_colab_notebook.py            # escribe el .ipynb
    ./scripts/gen_colab_notebook.py --check    # falla si esta desactualizado
"""
import argparse
import ast
import json
import pathlib
import sys

RAIZ = pathlib.Path(__file__).resolve().parent.parent
SPEC = RAIZ / "data" / "genomas.tsv"
SALIDA = RAIZ / "notebooks" / "descarga_genomas.ipynb"
DRIVE_DIR = "/content/drive/MyDrive/tesis/70_genomas"


def lee_spec():
    filas = []
    with SPEC.open() as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if f[0] == "org":
                continue
            # org especie fuente assembly accession estado confianza nota
            filas.append({
                "org": f[0], "especie": f[1], "assembly": f[3],
                "accession": f[4], "estado": f[5], "confianza": f[6],
            })
    return filas


def md(src):
    return {"cell_type": "markdown", "metadata": {}, "source": src.strip().split("\n")}


def code(src):
    ast.parse(src)  # no generar un notebook que no parsea
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src.strip().split("\n")}


def construye():
    filas = lee_spec()
    pendientes = [f for f in filas if f["estado"] == "candidato"]
    tabla = json.dumps(
        [{k: f[k] for k in ("org", "especie", "assembly", "accession", "confianza")}
         for f in pendientes],
        indent=2, ensure_ascii=False,
    )

    celdas = [
        md(f"""
# Descarga de ensamblados → Google Drive

Baja los genomas de referencia desde NCBI y los escribe directo en
`Mi unidad/tesis/70_genomas/`. No pasan por tu disco local.

**Por qué acá y no en la sesión de Claude:** esa sesión tiene bloqueada la
salida a NCBI y Ensembl por política del proxy. Colab no.

**Orden:** montar Drive → verificar candidatos → descargar → pegar el
resultado de vuelta en `data/genomas.tsv`.

Los {len(pendientes)} ensamblados de abajo están como `candidato`: son
propuestas **sin comprobar**. La celda de verificación los contrasta contra
NCBI antes de que se baje nada.
"""),
        md("## 1. Montar Drive"),
        code("""
from google.colab import drive
drive.mount('/content/drive')

import os, hashlib, json, urllib.request, zipfile, io, pathlib

DESTINO = pathlib.Path(%r)
DESTINO.mkdir(parents=True, exist_ok=True)
print('destino:', DESTINO)
print('existe:', DESTINO.exists())
""" % DRIVE_DIR),
        md("""
## 2. Candidatos a verificar

Snapshot de `data/genomas.tsv`. Si lo cambiaste en el repo, volvé a generar el
notebook con `./scripts/gen_colab_notebook.py`.
"""),
        code(f"""
CANDIDATOS = {tabla}

for c in CANDIDATOS:
    print(f"{{c['org']:8}} {{c['confianza']:6}} {{c['accession'] or '(sin candidato)':20}} {{c['assembly']}}")
"""),
        md("""
## 3. Verificar contra NCBI

Para cada organismo pregunta dos cosas: si el accession propuesto existe, y
cuál es el ensamblado de **referencia vigente** de la especie. La segunda
importa más que la primera — un accession puede existir y no ser el que
corresponde.
"""),
        code('''
import urllib.parse

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

def get(url):
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            return json.load(r)
    except Exception as e:
        return {"_error": str(e)}

def verifica(c):
    out = {"org": c["org"], "candidato": c["accession"], "coincide": None}

    if c["accession"]:
        d = get(f"{API}/genome/accession/{c['accession']}/dataset_report")
        reps = d.get("reports") or []
        if reps:
            r = reps[0]
            out["ncbi_nombre"] = r.get("assembly_info", {}).get("assembly_name")
            out["ncbi_organismo"] = r.get("organism", {}).get("organism_name")
        else:
            out["ncbi_nombre"] = None
            out["error"] = "accession no encontrado"

    esp = urllib.parse.quote(c["especie"])
    d = get(f"{API}/genome/taxon/{esp}/dataset_report"
            "?filters.reference_only=true&page_size=3")
    reps = d.get("reports") or []
    if reps:
        r = reps[0]
        out["referencia_vigente"] = r.get("accession")
        out["referencia_nombre"] = r.get("assembly_info", {}).get("assembly_name")
        out["nivel"] = r.get("assembly_info", {}).get("assembly_level")
        out["coincide"] = (out["referencia_vigente"] == c["accession"])
    else:
        out["referencia_vigente"] = None
        out["nota"] = "sin ensamblado de referencia; buscar por cepa"
    return out

RESULTADOS = [verifica(c) for c in CANDIDATOS]

for r in RESULTADOS:
    marca = "OK " if r.get("coincide") else "REVISAR"
    print(f"\\n[{marca}] {r['org']}")
    print(f"   candidato : {r.get('candidato') or '(ninguno)'}  {r.get('ncbi_nombre') or ''}")
    print(f"   vigente   : {r.get('referencia_vigente')}  {r.get('referencia_nombre') or ''}"
          f"  nivel={r.get('nivel')}")
    if r.get("nota"):
        print(f"   ! {r['nota']}")
'''),
        md("""
## 4. Elegir qué bajar

**Mirá la salida de arriba antes de correr esto.** Por defecto se baja el
ensamblado **de referencia vigente**, no el candidato — si difieren, gana NCBI.
Para forzar otro accession, editá `A_BAJAR` a mano.
"""),
        code('''
A_BAJAR = {r["org"]: r.get("referencia_vigente")
           for r in RESULTADOS if r.get("referencia_vigente")}

# Editá acá si querés forzar alguno, por ejemplo:
# A_BAJAR["cloro"] = "GCA_XXXXXXXXX.1"

faltan = [r["org"] for r in RESULTADOS if not r.get("referencia_vigente")]
if faltan:
    print("SIN RESOLVER (hay que buscarlos por cepa):", faltan)
print()
for o, a in sorted(A_BAJAR.items()):
    print(f"{o:8} -> {a}")
'''),
        md("""
## 5. Descargar a Drive

Escribe `<org>/<accession>.fna.gz` en `70_genomas/` y calcula el `sha256`.
Si el archivo ya está, lo saltea, así que se puede re-ejecutar sin problema.
"""),
        code('''
def baja(org, acc):
    destino = DESTINO / org
    destino.mkdir(parents=True, exist_ok=True)
    final = destino / f"{acc}.fna.gz"
    if final.exists() and final.stat().st_size > 0:
        print(f"== {org}: ya está ({final.stat().st_size/1e6:.0f} MB)")
        return None

    url = (f"{API}/genome/accession/{acc}/download"
           "?include_annotation_type=GENOME_FASTA")
    print(f"== {org}: bajando {acc} ...")
    with urllib.request.urlopen(url, timeout=1800) as r:
        blob = r.read()

    import gzip
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        nombres = [n for n in z.namelist()
                   if n.endswith("_genomic.fna") and f"/{acc}/" in n]
        if not nombres:
            print(f"   ERROR: no encontré el FASTA en el zip de {acc}")
            print("  ", z.namelist()[:10])
            return None
        crudo = z.read(nombres[0])

    # Se escribe comprimido: son cientos de MB y Drive cobra por byte.
    tmp = final.with_suffix(".gz.parcial")
    with gzip.open(tmp, "wb", compresslevel=6) as fh:
        fh.write(crudo)
    tmp.rename(final)

    h = hashlib.sha256()
    with final.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    sha = h.hexdigest()
    (destino / f"{acc}.fna.gz.sha256").write_text(f"{sha}  {acc}.fna.gz\\n")
    print(f"   {final.stat().st_size/1e6:.0f} MB  sha256={sha[:16]}...")
    return sha

SHAS = {}
for org, acc in sorted(A_BAJAR.items()):
    s = baja(org, acc)
    if s:
        SHAS[org] = (acc, s)
'''),
        md("""
## 6. Cerrar el círculo con git

Imprime las líneas para pegar en `data/genomas.tsv` y en
`data/genomas.sha256`. **Ese paso es el que deja constancia versionada de qué
genoma se usó** — el checksum guardado solo al lado del FASTA en Drive no
prueba nada, porque quien reemplace el genoma reemplaza el checksum con él.
"""),
        code('''
import datetime
hoy = datetime.date.today().isoformat()

print("--- pegar en data/genomas.sha256 ---")
print("org\\taccession\\tassembly\\tsha256\\tfecha_utc")
for org, (acc, sha) in sorted(SHAS.items()):
    nombre = next((r.get("referencia_nombre") for r in RESULTADOS
                   if r["org"] == org), "")
    print(f"{org}\\t{acc}\\t{nombre}\\t{sha}\\t{hoy}")

print()
print("--- en data/genomas.tsv, poner estado=verificado en:", sorted(SHAS), "---")
'''),
        md(f"""
## Lo que esto NO resuelve

Sirve para los genomas, que son unos pocos GB. **No sirve para los BAMs**: son
del orden de 340 GB y las sesiones de Colab son efímeras y con límite de
tiempo. Ese volumen sigue yendo por `rclone` desde tu máquina, con
`./scripts/drive_push.sh bam <org> --go`.

Tampoco sirve para el alineamiento: 30-40 h de bowtie no entran en una sesión
de Colab.
"""),
    ]

    return {
        "cells": celdas,
        "metadata": {
            "colab": {"provenance": [], "toc_visible": True},
            "kernelspec": {"display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="falla si el notebook esta desactualizado")
    a = ap.parse_args()

    nb = construye()
    txt = json.dumps(nb, indent=1, ensure_ascii=False) + "\n"

    if a.check:
        if not SALIDA.exists():
            print(f"falta {SALIDA}; corré {sys.argv[0]}", file=sys.stderr)
            return 1
        if SALIDA.read_text() != txt:
            print(f"{SALIDA} está desactualizado respecto de {SPEC.name}; "
                  f"corré {sys.argv[0]}", file=sys.stderr)
            return 1
        print("notebook al día")
        return 0

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(txt)
    n_code = sum(1 for c in nb["cells"] if c["cell_type"] == "code")
    print(f"{SALIDA}: {len(nb['cells'])} celdas ({n_code} de código)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
