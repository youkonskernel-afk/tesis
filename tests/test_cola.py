#!/usr/bin/env python3
"""scripts/cola.py contra un align.sh y un trim.sh falsos, sobre un Drive de mentira.

La cola corre sola durante horas sin que nadie mire cada paso, asi que lo que
tiene que garantizar se prueba aca y no en Colab:

  1. Un recorte que no verifica NO se alinea (gadmo_duplicado se alineo con 6
     librerias vacias y llego a Drive como bueno), y la cola sigue.
  2. Un BAM que no verifica no sube a Drive.
  3. Lo hecho (BAM en Drive) no se rehace; lo que no entra en RAM se saltea;
     lo tomado por otra maquina tambien.
  4. No arranca un proyecto que no termina antes del tope de horas.
  5. Lo medido queda en Drive al lado del BAM, por si la sesion muere.
"""
import contextlib
import io
import os
import pathlib
import sys
import tempfile

RAIZ = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "scripts"))
import cola  # noqa: E402
import reparto  # noqa: E402

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


ALIGN = r'''#!/usr/bin/env bash
echo "align.sh $*" >> "$LOG"
case "$1" in
  genoma)
    acc=$(awk -F'\t' -v o="$2" '$1==o{print $2}' "$GENOMAS_LEDGER")
    mkdir -p "$GENOMES_DIR/$2"
    printf 'c1\t5000000\n' > "$GENOMES_DIR/$2/$acc.fna.fai"
    echo 60 > "$GENOMES_DIR/$2/$acc.indice_s" ;;
  correr)
    org=${2%/*}; rol=${2#*/}; d="$PROY_DIR/${org}_${rol}"
    mkdir -p "$d/align" "$BAM_DIR/$org"
    { printf 'project\tlibrary\tumap\tmmap_wg\tmmap_nw\txmap_nw\txmap_ma\txmap_nv\txmap_fr\n'
      awk -F'\t' -v o="$org" -v r="$rol" 'NR>1 && $1==o && $4==r {print o"_"r"\t"$2"\t994\t1\t1\t1\t1\t1\t1"}' "$MANIFEST"
    } > "$d/align/library_stats.txt"
    head -c 30000 /dev/zero > "$d/align/alignment.bam"
    cp "$d/align/alignment.bam" "$BAM_DIR/$org/$rol.bam"
    printf 'proyecto\tx\n%s_%s\t1\n' "$org" "$rol" > "$d/alineado.tsv" ;;
  verificar) [[ " ${ALIN_FALLA:-} " == *" $2 "* ]] && exit 1; exit 0 ;;
  ledger) touch "$2" ;;
esac
'''

TRIM = r'''#!/usr/bin/env bash
echo "trim.sh $*" >> "$LOG"
case "$1" in
  correr)
    org=${2%/*}; rol=${2#*/}; d="$PROY_DIR/${org}_${rol}"
    mkdir -p "$d/trim"
    { printf 'run\tbioproject\tfichero\treads_in\treads_out\tretencion_pct\tfecha_utc\n'
      awk -F'\t' -v o="$org" -v r="$rol" 'NR>1 && $1==o && $4==r {print $2"\tPRJ\ttrim/"$2".t.fq.gz\t2000\t1000\t50.0\tT"}' "$MANIFEST"
    } > "$d/recortadas.tsv"
    for run in $(awk -F'\t' -v o="$org" -v r="$rol" 'NR>1 && $1==o && $4==r {print $2}' "$MANIFEST"); do
      head -c 20000 /dev/zero > "$d/trim/$run.t.fq.gz"
    done ;;
  verificar) [[ " ${VERIF_FALLA:-} " == *" $2 "* ]] && exit 1; exit 0 ;;
esac
'''

# org run bioproject rol read_count. Del mas chico al mas grande, en reads:
# sclsc/primario < rhirr/duplicado < sclsc/duplicado < galga/duplicado.
MANIFIESTO = [
    ("sclsc", "S1", "PRJS", "primario", 500_000), ("sclsc", "S2", "PRJS", "primario", 500_000),
    ("rhirr", "R1", "PRJR", "duplicado", 2_000_000),
    ("sclsc", "S3", "PRJT", "duplicado", 3_000_000),
    ("galga", "G1", "PRJG", "duplicado", 4_000_000),
    ("prupe", "P1", "PRJP", "primario", 100_000),
]


def escenario(tmp):
    clon = tmp / "clon"
    (clon / "data").mkdir(parents=True)
    (clon / "scripts").mkdir()
    for nom, txt in (("align.sh", ALIGN), ("trim.sh", TRIM)):
        f = clon / "scripts" / nom
        f.write_text(txt)
        f.chmod(0o755)
    (clon / "data" / "srr_manifest.tsv").write_text(
        "org\trun\tbioproject\trol\tread_count\n"
        + "".join(f"{o}\t{r}\t{b}\t{ro}\t{n}\n" for o, r, b, ro, n in MANIFIESTO))
    (clon / "data" / "adaptadores.tsv").write_text(
        "org\tbioproject\trun\tretencion_est\n"
        + "".join(f"{o}\t{b}\t-\t100\n" for o, _, b, _, _ in MANIFIESTO))
    (clon / "data" / "genomas.sha256").write_text(
        "org\taccession\tassembly\tsha256\tfecha_utc\n"
        + "".join(f"{o}\tGCF_{o.upper()}.1\tx\tdead\t2026\n"
                  for o in ("sclsc", "rhirr", "galga", "prupe")))
    drive = tmp / "drive"
    (drive / "10_bam" / "prupe").mkdir(parents=True)
    (drive / "10_bam" / "prupe" / "primario.bam").write_bytes(b"x")   # ya hecho
    for o in ("sclsc", "rhirr", "galga"):
        (drive / "70_genomas" / o).mkdir(parents=True)
    dirs = {k: tmp / k for k in ("genomes", "proyectos", "bams")}
    for v in dirs.values():
        v.mkdir()
    env = dict(os.environ, LOG=str(tmp / "log"), MANIFEST=str(clon / "data" / "srr_manifest.tsv"),
               GENOMAS_LEDGER=str(clon / "data" / "genomas.sha256"),
               GENOMES_DIR=str(dirs["genomes"]), PROY_DIR=str(dirs["proyectos"]),
               BAM_DIR=str(dirs["bams"]))
    (tmp / "log").write_text("")
    return clon, drive, env, dirs


def correr(tmp, clon, drive, env, dirs, **kw):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = cola.correr_cola(clon=clon, drive=drive, env=env, genomes=dirs["genomes"],
                               proy_dir=dirs["proyectos"], bam_dir=dirs["bams"],
                               yo="esta-vm", ram_b=8e9, disco_libre=500e9, **kw)
    return res, buf.getvalue(), (tmp / "log").read_text()


with tempfile.TemporaryDirectory() as d:
    tmp = pathlib.Path(d) / "a"
    clon, drive, env, dirs = escenario(tmp)
    env.update(VERIF_FALLA="rhirr/duplicado", ALIN_FALLA="sclsc/duplicado")
    res, out, log = correr(tmp, clon, drive, env, dirs)

    print("== 1. el que anda se hace entero y sube a Drive")
    chk("sclsc/primario hecho", res["hechos"] == ["sclsc/primario"], res)
    chk("su BAM en Drive", (drive / "10_bam" / "sclsc" / "primario.bam").exists())
    chk("y su alineado.tsv al lado",
        (drive / "10_bam" / "sclsc" / "primario.alineado.tsv").exists())
    cal = drive / "10_bam" / "sclsc" / "primario.calibracion.tsv"
    chk("y la calibración al lado del BAM", cal.exists())
    filas = cal.read_text().splitlines() if cal.exists() else []
    chk("con los reads del BAM (2 libs x 1000)", len(filas) == 2 and "\t2000\t" in filas[1], filas)
    chk("también en data/calibracion.tsv del clon",
        "sclsc_primario" in (clon / "data" / "calibracion.tsv").read_text())

    print("== 2. un recorte que no verifica NO se alinea")
    chk("rhirr/duplicado falló", any(p == "rhirr/duplicado" for p, _ in res["fallados"]),
        res["fallados"])
    chk("dice por qué", "VACIA" in out and "3b" in out, out[-800:])
    chk("y no se llamó a align.sh correr", "align.sh correr rhirr/duplicado" not in log, log)
    chk("ni subió nada", not (drive / "10_bam" / "rhirr").exists())

    print("== 3. un BAM que no verifica no sube")
    chk("sclsc/duplicado falló", any(p == "sclsc/duplicado" for p, _ in res["fallados"]))
    chk("y su BAM no está en Drive", not (drive / "10_bam" / "sclsc" / "duplicado.bam").exists())

    print("== 4. lo hecho no se rehace; lo que no entra en RAM se saltea")
    chk("prupe no aparece en el log", "prupe" not in log, log)
    chk("galga saltado por RAM", ("galga/duplicado", "RAM") in res["saltados"], res["saltados"])
    chk("y nunca se tocó", "galga" not in log, log)

    print("== 5. el orden es del más chico al más grande")
    orden = [ln.split()[-1] for ln in log.splitlines() if ln.startswith("trim.sh correr")]
    chk("sclsc/primario, rhirr/duplicado, sclsc/duplicado",
        orden == ["sclsc/primario", "rhirr/duplicado", "sclsc/duplicado"], orden)

    print("== 6. limpia el disco y suelta los claims")
    chk("sin directorios de proyecto", list(dirs["proyectos"].iterdir()) == [],
        list(dirs["proyectos"].iterdir()))
    chk("sin BAM locales", not any(dirs["bams"].rglob("*.bam")))
    chk("sin genomas de los que ya no siguen", list(dirs["genomes"].iterdir()) == [],
        list(dirs["genomes"].iterdir()))
    chk("ningún claim queda tomado", reparto.tomados(drive / "90_claims") == {},
        reparto.tomados(drive / "90_claims"))

    print("== 7. relanzar sigue desde lo que falta")
    env.update(VERIF_FALLA="", ALIN_FALLA="")
    (tmp / "log").write_text("")
    res2, out2, log2 = correr(tmp, clon, drive, env, dirs)
    chk("no rehace sclsc/primario", "trim.sh correr sclsc/primario" not in log2, log2)
    chk("y hace los dos que fallaron",
        sorted(res2["hechos"]) == ["rhirr/duplicado", "sclsc/duplicado"], res2)

    print("== 8. no arranca uno que no termina antes del tope de horas")
    tmp = pathlib.Path(d) / "b"
    clon, drive, env, dirs = escenario(tmp)
    res3, out3, log3 = correr(tmp, clon, drive, env, dirs, max_horas=0.0001)
    chk("se corta", res3["cortado"], res3)
    chk("sin empezar ninguno", "trim.sh correr" not in log3, log3)
    chk("y dice que hay que relanzar", "relanz" in out3, out3[-400:])

    print("== 9. lo que otra máquina tiene tomado se saltea")
    tmp = pathlib.Path(d) / "c"
    clon, drive, env, dirs = escenario(tmp)
    reparto.tomar(drive / "90_claims", "sclsc/primario", "otra-vm")
    res4, out4, log4 = correr(tmp, clon, drive, env, dirs, solo=["sclsc/primario"])
    chk("saltado", any(p == "sclsc/primario" for p, _ in res4["saltados"]), res4)
    # Que diga QUIEN lo tiene: el resumen de la cola es lo que se lee al volver.
    chk("diciendo quién lo tiene", ("sclsc/primario", "tomado por otra-vm") in res4["saltados"],
        res4["saltados"])
    chk("sin tocarlo", "sclsc/primario" not in log4, log4)

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
