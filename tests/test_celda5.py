#!/usr/bin/env python3
"""§5 de 20_alinear: lo medido que va a data/calibracion.tsv.

Se corre la celda DE VERDAD con exec, sobre un proyecto de mentira en un
tmpdir. Paso en serio con gadmo_duplicado: 6 de sus 12 corridas son
PRE-TRIMMED, que en recortadas.tsv dicen '-' en READS_OUT, y §5 contaba los
reads de ahi. Calibro con 16 M reads en vez de 55: 360 s/M en vez de ~105 y
B_BAM 44 en vez de ~13. El cronograma de los 4 proyectos grandes se habria
triplicado por un conteo, no por una medicion.
"""
import contextlib
import csv
import io
import json
import pathlib
import shutil
import sys
import tempfile

RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB = RAIZ / "notebooks" / "20_alinear.ipynb"
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


FUENTE = celda("medido en")


def proyecto(tmp, filas_ledger, libs, bam_bytes, fq_bytes):
    """filas_ledger: (run, fichero, reads_out|'-'); libs: (run, total)."""
    pdir = tmp / "proyectos" / "aa_dup"
    (pdir / "align").mkdir(parents=True)
    (pdir / "trim").mkdir()
    (pdir / "untrimmed").mkdir()
    lin = ["run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc"]
    for run, fich, out in filas_ledger:
        (pdir / fich).write_bytes(b"x" * fq_bytes)
        pct = "PRE-TRIMMED" if out == "-" else "50.0"
        lin.append(f"{run}\tPRJ\t{fich}\t{'-' if out == '-' else 2 * out}\t{out}\t{pct}\tT")
    (pdir / "recortadas.tsv").write_text("\n".join(lin) + "\n")
    st = ["project\tlibrary\tumap\tmmap_wg\tmmap_nw\txmap_nw\txmap_ma\txmap_nv\txmap_fr"]
    for run, tot in libs:
        # se reparte el total entre las 7 columnas: tiene que sumarlas todas
        st.append(f"aa_dup\t{run}\t{tot - 6}\t1\t1\t1\t1\t1\t1")
    (pdir / "align" / "library_stats.txt").write_text("\n".join(st) + "\n")
    (pdir / "align" / "alignment.bam").write_bytes(b"b" * bam_bytes)
    (pdir / "alineado.tsv").write_text("proyecto\nx\n")
    (tmp / "bams" / "aa").mkdir(parents=True)
    return pdir


def correr_celda(tmp, t_alin=1000.0):
    calib = tmp / "calibracion.tsv"
    ns = {
        "PROYECTO": "aa/dup", "correr": lambda *a, **k: 0,
        "DRIVE": tmp / "drive", "BAM_DIR": tmp / "bams", "PROY_DIR": tmp / "proyectos",
        "TANDA": None, "shutil": shutil, "csv": csv,
        "bases_de": lambda org: (100_000_000,), "proy": {("aa", "dup"): (99,)},
        "B_FQGZ": 22, "B_BAM": 16, "s_por_m": lambda mb: 60.0,
        "CALIB": calib, "_acc": {"aa": "GCF_X.1"},
        "T_ALIN": t_alin, "INDICE_EN_RELOJ": False,
    }
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(FUENTE, ns)
    filas = list(csv.DictReader(
        (l for l in calib.read_text().splitlines() if not l.startswith("#")),
        delimiter="\t")) if calib.exists() else []
    return buf.getvalue(), filas


with tempfile.TemporaryDirectory() as d:
    print("== 1. un proyecto mixto: las PRE-TRIMMED cuentan")
    # 2 recortadas (READS_OUT 1000 c/u) y 2 PRE-TRIMMED ('-'): el BAM tiene 4
    # librerias y 10 000 reads. Con READS_OUT habrian sido 2000.
    tmp = pathlib.Path(d) / "a"
    proyecto(tmp,
             [("SRR1", "trim/SRR1.t.fq.gz", 1000), ("SRR2", "trim/SRR2.t.fq.gz", 1000),
              ("SRR3", "untrimmed/SRR3.fastq.gz", "-"), ("SRR4", "untrimmed/SRR4.fastq.gz", "-")],
             [("SRR1", 1000), ("SRR2", 1000), ("SRR3", 4000), ("SRR4", 4000)],
             bam_bytes=150_000, fq_bytes=50_000)
    out, filas = correr_celda(tmp, t_alin=1000.0)
    r = filas[0] if filas else {}
    chk("escribe una fila", len(filas) == 1, out)
    chk("reads = lo que entró al BAM (10 000)", r.get("reads_trim") == "10000", r)
    chk("librerías = las 4 del BAM, no las 2 .t.fq.gz", r.get("librerias") == "4", r)
    chk("s/M sobre todos los reads (1000 s / 0.01 M)", r.get("s_por_m") == "100000.0", r)
    chk("B_BAM sobre todos los reads (150 000 / 10 000)", r.get("b_bam") == "15.0", r)
    chk("el disco incluye el fastq de las PRE-TRIMMED (4 x 50 000)",
        r.get("b_fqgz") == "20.0", r)
    chk("y lo imprime así", "0 M reads" in out and "B_BAM   15.0" in out, out)

    print("== 2. sin library_stats no inventa una calibración")
    tmp = pathlib.Path(d) / "b"
    pdir = proyecto(tmp, [("SRR1", "trim/SRR1.t.fq.gz", 1000)], [("SRR1", 1000)],
                    bam_bytes=10, fq_bytes=10)
    (pdir / "align" / "library_stats.txt").unlink()
    out, filas = correr_celda(tmp)
    chk("no escribe fila", filas == [], filas)
    chk("y dice de dónde faltan", "library_stats" in out, out)

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
