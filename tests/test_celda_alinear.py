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
def celda(marca):
    """La celda por lo que dice, no por su indice: agregar una celda arriba no
    tiene por que romper el banco de otra — y ya paso al insertar §1b."""
    hall = [i for i, c in enumerate(json.loads(NB.read_text())["cells"])
            if c["cell_type"] == "code" and marca in "".join(c["source"])]
    assert len(hall) == 1, f"{marca!r} aparece {len(hall)} veces en {NB.name}"
    return hall[0]


CELDA = celda("B_FQGZ = 22")
CELDA_1B = [CELDA, celda("PRUEBA = None")]   # §1b usa el namespace de §1

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


def correr(tmp, corridas, libre_gb, retenciones, en_drive=(),
           genomas=None, calib=None, ram_gb=64, celda=CELDA, ns_extra=None):
    """Ejecuta la celda con un clon y un Drive de mentira.

    `genomas` es {org: Mb}: se escribe el .bases que align.sh cachea al lado del
    FASTA, asi que la celda lo lee como EXACTO. Un org que no este ahi no tiene
    genoma medible, que es el caso de una VM recien arrancada.
    """
    clon = tmp / "clon"
    # Drive es lo unico que persiste entre sesiones de Colab: la VM arranca
    # con /content vacio, asi que "que falta" no se puede leer del disco local.
    drive = tmp / "drive"
    genomes = tmp / "genomes"
    for p in en_drive:
        org, rol = p.split("/")
        (drive / "10_bam" / org).mkdir(parents=True, exist_ok=True)
        (drive / "10_bam" / org / f"{rol}.bam").write_bytes(b"BAM")
    (clon / "data").mkdir(parents=True, exist_ok=True)
    (clon / "data" / "srr_manifest.tsv").write_text("\n".join([HDR] + corridas) + "\n")
    (clon / "data" / "adaptadores.tsv").write_text(
        "# tabla\n"
        "org\tbioproject\tfamilia\tsecuencia\tadapt_pct\tinserto_modal\tretencion_est\tveredicto\tfecha\n"
        + "".join(f"{o}\t{p}\tRA3\tACGT\t100\t22\t{r}\tPARECE sRNA-seq\t2026-09-23\n"
                 for (o, p), r in retenciones.items()))

    orgs = sorted({c.split("\t")[0] for c in corridas})
    (clon / "data" / "genomas.sha256").write_text(
        "org\taccession\tassembly\tsha256\tfecha_utc\n"
        + "".join(f"{o}\tGCF_{o}.1\tAsm\tdead\t2026-09-23\n" for o in orgs))
    for o, mb in (genomas or {}).items():
        (genomes / o).mkdir(parents=True, exist_ok=True)
        (genomes / o / f"GCF_{o}.1.fna.bases").write_text(str(int(mb * 1e6)))

    if calib is not None:
        (clon / "data" / "calibracion.tsv").write_text(
            "proyecto\torg\trol\taccession\tgenoma_mb\tlibrerias\treads_trim\t"
            "seg_reloj\ts_por_m\tb_bam\tb_fqgz\tfecha_utc\n"
            + "".join(f"p{i}\to\tr\tA\t{mb}\t2\t1000\t100\t{spm}\t14\t22\t2026-09-23T00:00:00Z\n"
                      for i, (mb, spm) in enumerate(calib)))

    # Se parchea shutil.disk_usage de verdad en vez de inyectar un falso en el
    # namespace: la celda hace `import shutil`, asi que un falso inyectado lo
    # pisaria el import y el chequeo correria contra el disco de esta maquina.
    class _U:
        free = int(libre_gb * 1e9)

    # Idem /proc/meminfo: la celda lo abre por nombre.
    meminfo = tmp / "meminfo"
    meminfo.write_text(f"MemTotal:  99999999 kB\nMemAvailable: {int(ram_gb * 2**30 / 1024)} kB\n")
    real_open, real_du = open, shutil.disk_usage

    def _open(f, *a, **k):
        return real_open(meminfo if f == "/proc/meminfo" else f, *a, **k)

    shutil.disk_usage = lambda _: _U()
    try:
        celdas = celda if isinstance(celda, (list, tuple)) else [celda]
        fuentes = ["".join(json.loads(NB.read_text())["cells"][c]["source"])
                   for c in celdas]
        ns = {"CLON": clon, "DRIVE": drive, "GENOMES": genomes, "open": _open}
        ns.update(ns_extra or {})
        buf = io.StringIO()
        with redirect_stdout(buf):
            for fuente in fuentes:
                exec(fuente, ns)
    finally:
        shutil.disk_usage = real_du
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
    chk("y lo dice", "NO ENTRA: disco" in out, out)
    chk("y dice qué hacer", "máquina local" in out, out)

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

    print("== 6. un proyecto que ya tiene BAM en Drive no se rehace")
    # La VM es efimera: align.sh estado ve /content vacio y diria que falta.
    # Se perdieron 31 min re-alineando sclsc_duplicado por esto.
    _, ns6 = correr(
        tmp / "g",
        [corrida("aa", "SRR1", "PRJ_A", "duplicado", 10_000_000),
         corrida("aa", "SRR2", "PRJ_B", "primario", 20_000_000)],
        libre_gb=100,
        retenciones={("aa", "PRJ_A"): 100, ("aa", "PRJ_B"): 100},
        en_drive=["aa/duplicado"])
    chk("el hecho queda aparte", ns6["HECHOS"] == ["aa/duplicado"], ns6["HECHOS"])
    chk("y no cuenta como pendiente", "aa/duplicado" not in ns6["ENTRAN"], ns6["ENTRAN"])
    chk("SIGUIENTE saltea el hecho", ns6["SIGUIENTE"] == "aa/primario", ns6["SIGUIENTE"])

    print("== 7. sin nada en Drive, SIGUIENTE es el mas chico")
    # Es el mismo orden de la tabla: una sesion que se muere pierde un proyecto,
    # asi que conviene que sea el chico.
    _, ns7 = correr(
        tmp / "h",
        [corrida("aa", "SRR1", "PRJ_A", "duplicado", 10_000_000),
         corrida("aa", "SRR2", "PRJ_B", "primario", 20_000_000)],
        libre_gb=100,
        retenciones={("aa", "PRJ_A"): 100, ("aa", "PRJ_B"): 100})
    chk("el mas chico primero", ns7["SIGUIENTE"] == "aa/duplicado", ns7["SIGUIENTE"])

    print("== 8. todo hecho: SIGUIENTE queda en None y §2 corta")
    out8, ns8 = correr(
        tmp / "i", [corrida("aa", "SRR1", "PRJ_A", "duplicado", 10_000_000)],
        libre_gb=100, retenciones={("aa", "PRJ_A"): 100},
        en_drive=["aa/duplicado"])
    chk("SIGUIENTE es None", ns8["SIGUIENTE"] is None, ns8["SIGUIENTE"])
    chk("y lo dice", "No queda ninguno pendiente" in out8, out8)

    print("== 9. uno que no entra tampoco puede ser SIGUIENTE")
    # Marcarlo como el que sigue mandaria a gastar horas en algo que se queda
    # sin disco a la mitad.
    _, ns9 = correr(
        tmp / "j",
        [corrida("aa", "SRR1", "PRJ_A", "duplicado", 2_000_000_000),
         corrida("aa", "SRR2", "PRJ_B", "primario", 10_000_000)],
        libre_gb=50,
        retenciones={("aa", "PRJ_A"): 100, ("aa", "PRJ_B"): 100})
    chk("el grande no es SIGUIENTE", ns9["SIGUIENTE"] == "aa/primario", ns9["SIGUIENTE"])
    chk("y queda en NO_ENTRAN", ns9["NO_ENTRAN"] == ["aa/duplicado"], ns9["NO_ENTRAN"])

    print("== 10. la RAM decide igual que el disco, y por separado")
    # unique_d es un entero de Python por base del genoma: 8 B medidos, y se
    # arma entero antes del primer read. No depende de los reads NI de las
    # librerias. Un proyecto chico contra un genoma grande entra en disco y
    # muere por memoria a las horas — que es como murio gadmo_duplicado.
    out10, ns10 = correr(
        tmp / "j", [corrida("gadmo", "SRR1", "PRJ_D", "primario", 10_000_000)],
        libre_gb=200, ram_gb=4, retenciones={("gadmo", "PRJ_D"): 100},
        genomas={"gadmo": 670})
    chk("no entra por RAM", ns10["NO_ENTRAN"] == ["gadmo/primario"], ns10["NO_ENTRAN"])
    chk("y dice que es la RAM", "NO ENTRA: RAM" in out10, out10)
    chk("no el disco", "NO ENTRA: disco" not in out10, out10)
    chk("ni es SIGUIENTE", ns10["SIGUIENTE"] is None, ns10["SIGUIENTE"])

    out10b, ns10b = correr(
        tmp / "k", [corrida("gadmo", "SRR1", "PRJ_D", "primario", 10_000_000)],
        libre_gb=200, ram_gb=32, retenciones={("gadmo", "PRJ_D"): 100},
        genomas={"gadmo": 670})
    chk("con 32 GB de RAM sí entra", ns10b["ENTRAN"] == ["gadmo/primario"], ns10b["ENTRAN"])

    print("== 11. un genoma que no se pudo medir NO se declara imposible")
    # Mandar a la maquina local un proyecto que a lo mejor entra cuesta tanto
    # como lo contrario. Sin genoma medible se avisa y se deja pasar.
    out11, ns11 = correr(
        tmp / "l", [corrida("gadmo", "SRR1", "PRJ_D", "primario", 10_000_000)],
        libre_gb=200, ram_gb=1, retenciones={("gadmo", "PRJ_D"): 100})
    chk("entra igual", ns11["ENTRAN"] == ["gadmo/primario"], ns11["ENTRAN"])
    chk("y el genoma sale como ?", " ?" in out11, out11)

    print("== 12. con un solo punto el s/M es plano; con dos, ajustado")
    # El cronograma de ~98 h salia de UN proyecto medido contra un genoma de
    # 39 Mb, y galga es 1.05 Gb. Si el s/M escala con el genoma, el numero se va
    # -- y los 4 proyectos mas caros son el 54% de los reads del set.
    corr12 = [corrida("galga", "SRR1", "PRJ_G", "primario", 100_000_000)]
    ret12 = {("galga", "PRJ_G"): 100}
    _, uno = correr(tmp / "m", corr12, libre_gb=200, retenciones=ret12,
                    genomas={"galga": 1000}, calib=[(39, 58)])
    chk("un punto -> sin ajuste", uno["AJUSTE"] is None, uno["AJUSTE"])
    chk("y usa el plano", uno["s_por_m"](1000) == uno["S_POR_M_PLANO"])

    out12, dos = correr(tmp / "n", corr12, libre_gb=200, retenciones=ret12,
                        genomas={"galga": 1000}, calib=[(39, 58), (639, 108)])
    chk("dos puntos -> recta", dos["AJUSTE"] is not None, dos["AJUSTE"])
    # 58 s/M a 39 Mb y 108 a 639: la pendiente es 50/600 = 0.0833 s/M por Mb,
    # asi que a 1000 Mb da 58 + 0.0833*961 = 138.
    chk("interpola donde se midió", round(dos["s_por_m"](39)) == 58, dos["s_por_m"](39))
    chk("y extrapola al genoma grande", round(dos["s_por_m"](1000)) == 138,
        dos["s_por_m"](1000))
    chk("lo dice", "AJUSTADO" in out12, out12)
    chk("y muestra los dos estimados", "plano" in out12, out12)

    print("== 13. si el ajustado se va del plano, avisa")
    # Dos puntos no son un modelo. Cuando la recta manda el cronograma lejos del
    # numero plano, lo que hay que hacer es medir un tercero, no creerle.
    out13, _ = correr(tmp / "o", corr12, libre_gb=200, retenciones=ret12,
                      genomas={"galga": 1000}, calib=[(39, 58), (639, 300)])
    chk("avisa de la extrapolación", "extrapolación" in out13, out13)
    out13b, _ = correr(tmp / "p", corr12, libre_gb=200, retenciones=ret12,
                       genomas={"galga": 1000}, calib=[(39, 58), (639, 59)])
    chk("y no avisa cuando coinciden", "extrapolación" not in out13b, out13b)

    print("== 14. §1b propone el que mide más genoma nuevo por hora")
    # No el mas barato ni el mas grande: el que mas reduce la incertidumbre por
    # unidad de tiempo. Y solo entre los que ENTRAN -- proponer algo que se
    # queda sin disco a la mitad es peor que no proponer nada.
    corr14 = [corrida("sclsc", "SRR1", "PRJ_S", "duplicado", 30_000_000),
              corrida("gadmo", "SRR2", "PRJ_D", "primario", 56_000_000),
              corrida("cloro", "SRR3", "PRJ_C", "primario", 40_000_000)]
    ret14 = {("sclsc", "PRJ_S"): 100, ("gadmo", "PRJ_D"): 100, ("cloro", "PRJ_C"): 100}
    gen14 = {"sclsc": 39, "gadmo": 670, "cloro": 71}
    out14, ns14 = correr(tmp / "q", corr14, libre_gb=200, retenciones=ret14,
                         genomas=gen14, calib=[(39, 58)], celda=CELDA_1B,
                         ns_extra=None)
    chk("elige el de genoma lejano", ns14["PRUEBA"] == "gadmo/primario", ns14["PRUEBA"])
    chk("no el más barato", ns14["SIGUIENTE"] == "sclsc/duplicado", ns14["SIGUIENTE"])
    chk("y explica la diferencia", "mide más genoma nuevo" in out14, out14)

    print("== 15. §1b no propone uno que no entra")
    out15, ns15 = correr(tmp / "r", corr14, libre_gb=200, ram_gb=4,
                         retenciones=ret14, genomas=gen14, calib=[(39, 58)],
                         celda=CELDA_1B)
    chk("gadmo queda afuera por RAM", "gadmo/primario" in ns15["NO_ENTRAN"], ns15["NO_ENTRAN"])
    chk("y no es la PRUEBA", ns15["PRUEBA"] != "gadmo/primario", ns15["PRUEBA"])
    chk("pero propone alguno", ns15["PRUEBA"] in ns15["ENTRAN"], ns15["PRUEBA"])

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
