#!/usr/bin/env python3
"""La celda §1 de 20_alinear.ipynb: que entra en el disco de esta VM.

Es la celda que existe para que un proyecto que no entra se sepa EN UN SEGUNDO y
no a las seis horas. Un proyecto grande —galga_duplicado son ~1500 M reads tras
el recorte— necesita recortado + 2 x BAM al mismo tiempo, porque pysam.sort
escribe el BAM ordenado antes de borrar el sin ordenar. En una VM de Colab Free
eso no entra, y `yasma align` no avisa: se queda sin disco a mitad de camino.

Se afirma sobre las LISTAS que la celda arma (ENTRAN / NO_ENTRAN), que es lo que
la persona usa para elegir el proyecto, no sobre el texto que imprime al lado.
"""
import io
import json
import pathlib
import shutil
import sys
import tempfile
from contextlib import redirect_stdout

RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB = RAIZ / "notebooks" / "20_alinear.ipynb"
CELDA = 9

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


HDR = ("org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\t"
       "avg_len\tstrategy\tlayout\tsource")


def corrida(org, run, proy, rol, reads):
    return f"{org}\t{run}\t{proy}\t{rol}\tapl\t{reads}\t{reads*50}\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC"


def correr(tmp, corridas, libre_gb, retenciones):
    clon = tmp / "clon"
    (clon / "data").mkdir(parents=True, exist_ok=True)
    (clon / "data" / "srr_manifest.tsv").write_text("\n".join([HDR] + corridas) + "\n")
    (clon / "data" / "adaptadores.tsv").write_text(
        "# tabla\n"
        "org\tbioproject\tfamilia\tsecuencia\tadapt_pct\tinserto_modal\tretencion_est\tveredicto\tfecha\n"
        + "".join(f"{o}\t{p}\tRA3\tACGT\t100\t22\t{r}\tPARECE sRNA-seq\t2026-09-23\n"
                 for (o, p), r in retenciones.items()))

    # Se parchea shutil.disk_usage de verdad en vez de inyectar un falso en el
    # namespace: la celda hace `import shutil`, asi que un falso inyectado lo
    # pisaria el import y el chequeo correria contra el disco de esta maquina.
    class _U:
        free = int(libre_gb * 1e9)

    real = shutil.disk_usage
    shutil.disk_usage = lambda _: _U()
    try:
        src = "".join(json.loads(NB.read_text())["cells"][CELDA]["source"])
        ns = {"CLON": clon}
        buf = io.StringIO()
        with redirect_stdout(buf):
            exec(src, ns)
    finally:
        shutil.disk_usage = real
    return buf.getvalue(), ns


with tempfile.TemporaryDirectory() as d:
    tmp = pathlib.Path(d)

    print("== 1. uno chico entra y uno enorme no")
    # 10 M reads -> pico ~1.2 GB ; 1500 M -> ~170 GB
    out, ns = correr(
        tmp / "a",
        [corrida("sclsc", "SRR1", "PRJ_S", "duplicado", 10_000_000),
         corrida("galga", "SRR2", "PRJ_G", "duplicado", 1_500_000_000)],
        libre_gb=78,
        retenciones={("sclsc", "PRJ_S"): 100, ("galga", "PRJ_G"): 100})
    chk("el chico entra", "sclsc/duplicado" in ns["ENTRAN"], ns["ENTRAN"])
    chk("el grande NO", "galga/duplicado" in ns["NO_ENTRAN"], ns["NO_ENTRAN"])
    chk("y lo dice", "NO ENTRA en esta VM" in out)
    chk("y dice qué hacer", "Colab Pro" in out or "máquina local" in out, out)

    print("== 2. el orden es del más chico al más grande")
    # Empezar por el mas chico es lo que hace que una sesion que se muere haya
    # dejado algo hecho: la unidad reanudable es el proyecto entero.
    lineas = [l for l in out.splitlines() if l.startswith(("sclsc", "galga"))]
    chk("el chico va primero", lineas and lineas[0].startswith("sclsc"), lineas)

    print("== 3. el mismo proyecto entra o no según el disco")
    _, ns_poco = correr(
        tmp / "b", [corrida("sclsc", "SRR1", "PRJ_S", "duplicado", 10_000_000)],
        libre_gb=5, retenciones={("sclsc", "PRJ_S"): 100})
    chk("con 5 GB no entra", ns_poco["NO_ENTRAN"] == ["sclsc/duplicado"], ns_poco["NO_ENTRAN"])
    _, ns_mucho = correr(
        tmp / "c", [corrida("sclsc", "SRR1", "PRJ_S", "duplicado", 10_000_000)],
        libre_gb=78, retenciones={("sclsc", "PRJ_S"): 100})
    chk("con 78 GB sí", ns_mucho["ENTRAN"] == ["sclsc/duplicado"], ns_mucho["ENTRAN"])

    print("== 4. cuenta con lo que SOBREVIVE al recorte, no con los reads crudos")
    # galga PRJEB12164 tiene 95% de adaptador y retiene 52%: usar el crudo
    # sobreestimaria el disco casi al doble y marcaria como imposible algo que
    # entra.
    # Los 1200 M estan elegidos para que la retencion decida con las constantes
    # MEDIDAS (22 B/read en fq.gz, 16 en BAM = 54 B/read de pico). Eran 400 M
    # cuando B_BAM era la estimacion de 45, y al corregirla a lo medido el
    # escenario dejo de discriminar: los dos entraban.
    _, r100 = correr(
        tmp / "d", [corrida("xx", "SRR1", "PRJ_X", "primario", 1_200_000_000)],
        libre_gb=50, retenciones={("xx", "PRJ_X"): 100})
    _, r50 = correr(
        tmp / "e", [corrida("xx", "SRR1", "PRJ_X", "primario", 1_200_000_000)],
        libre_gb=50, retenciones={("xx", "PRJ_X"): 50})
    chk("con 100% de retención no entra", r100["NO_ENTRAN"] == ["xx/primario"])
    chk("con 50% sí", r50["ENTRAN"] == ["xx/primario"])

    print("== 5. el primario y el duplicado se miden por separado")
    # Son proyectos YASMA distintos y se alinean por separado: sumarlos daria un
    # pico que nunca ocurre.
    _, ns5 = correr(
        tmp / "f",
        [corrida("aa", "SRR1", "PRJ_A", "primario", 200_000_000),
         corrida("aa", "SRR2", "PRJ_B", "duplicado", 200_000_000)],
        libre_gb=20,
        retenciones={("aa", "PRJ_A"): 100, ("aa", "PRJ_B"): 100})
    chk("los dos entran por separado",
        sorted(ns5["ENTRAN"]) == ["aa/duplicado", "aa/primario"], ns5["ENTRAN"])

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
