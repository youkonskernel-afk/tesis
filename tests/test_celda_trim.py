#!/usr/bin/env python3
"""Las celdas §3b y §3c de 20_alinear.ipynb: qué corridas re-medir y rehacer.

Existen por gadmo/duplicado: 12 corridas de un mismo BioProject, 6 recortadas
con la secuencia equivocada porque `perfil --proyectos` mide UNA corrida y le
tocó una de las buenas. El guardia de §3 cortó bien, pero el paso siguiente
quedaba en transcribir seis accessions a mano desde una tabla formateada — que
es exactamente donde se cuela un error que después no falla ruidosamente.

Se afirma sobre la lista MALAS que la celda arma y sobre los argumentos con los
que llama a `rehacer`, no sobre el texto que imprime al lado.
"""
import io
import json
import pathlib
import tempfile
from contextlib import redirect_stdout

RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB = RAIZ / "notebooks" / "20_alinear.ipynb"
CELDA_3B = 15
CELDA_3C = 17

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


# La tabla real de `trim.sh verificar`, con el ancho real. Copiada de la corrida
# de gadmo/duplicado: 6 VACIA y 6 ok.
TABLA = """PROYECTO           BIOPROJECT     CORRIDA        READS_IN  READS_OUT      MEDIDA  ESPERA  VEREDICTO
gadmo_duplicado    PRJNA328800    SRR3884830      4213889    2072357       49.2%     51%  ok
gadmo_duplicado    PRJNA328800    SRR3884837      8322376      48445        0.6%     51%  VACIA — casi seguro el adaptador equivocado
gadmo_duplicado    PRJNA328800    SRR3884824      6511402    2491213       38.3%     51%  ok
gadmo_duplicado    PRJNA328800    SRR3884834      6821307     115103        1.7%     51%  VACIA — casi seguro el adaptador equivocado
gadmo_duplicado    PRJNA328800    SRR3884866      7108218      51125        0.7%     51%  VACIA — casi seguro el adaptador equivocado
gadmo_duplicado    PRJNA328800    SRR3884814      5143154    2981804       58.0%     51%  ok
gadmo_duplicado    PRJNA328800    SRR3884815      6695993    2558498       38.2%     51%  ok
gadmo_duplicado    PRJNA328800    SRR3884828      6066301    3081639       50.8%     51%  ok
gadmo_duplicado    PRJNA328800    SRR3884833      5078178      59783        1.2%     51%  VACIA — casi seguro el adaptador equivocado
gadmo_duplicado    PRJNA328800    SRR3884832      5450651    2843846       52.2%     51%  ok
gadmo_duplicado    PRJNA328800    SRR3884835      6350751     124714        2.0%     51%  VACIA — casi seguro el adaptador equivocado
gadmo_duplicado    PRJNA328800    SRR3884836      5283583      51000        1.0%     51%  VACIA — casi seguro el adaptador equivocado

12 corridas verificadas, 6 fuera de lo esperado
"""

MALAS_REALES = ["SRR3884837", "SRR3884834", "SRR3884866",
                "SRR3884833", "SRR3884835", "SRR3884836"]

TODO_OK = """PROYECTO           BIOPROJECT     CORRIDA        READS_IN  READS_OUT      MEDIDA  ESPERA  VEREDICTO
sclsc_duplicado    PRJNA1135930   SRR31851668    12367469   12248728       99.0%     97%  ok
sclsc_duplicado    PRJNA1135930   SRR31851669    20014299   19897362       99.4%     97%  ok

2 corridas verificadas, ninguna desviada
"""

# Una DESVIADA (no VACIA) tambien tiene que entrar: es el otro veredicto que
# verificar emite, y el que aparece cuando el adaptador es casi correcto.
DESVIADA = """PROYECTO           BIOPROJECT     CORRIDA        READS_IN  READS_OUT      MEDIDA  ESPERA  VEREDICTO
aa_primario        PRJ_A          SRR_A1          1000000     300000       30.0%     95%  DESVIADA — 65 puntos por debajo
aa_primario        PRJ_A          SRR_A2          1000000     940000       94.0%     95%  ok
"""


def fuente(i):
    return "".join(json.loads(NB.read_text())["cells"][i]["source"])


def correr_3b(salida_verificar, rc):
    """Ejecuta §3b con un `correr` falso. Devuelve (ns, stdout, llamadas)."""
    llamadas = []

    def correr(script, *args, env_extra=None, capturar=False):
        llamadas.append((script, args))
        if capturar:
            return (rc, salida_verificar)
        return 0

    ns = {"correr": correr, "PROYECTO": "gadmo/duplicado"}
    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            exec(fuente(CELDA_3B), ns)
    except AssertionError as e:
        ns["_error"] = str(e)
    return ns, buf.getvalue(), llamadas


def correr_3c(filas, malas, tabla):
    """Ejecuta §3c con un `correr` falso y una data/adaptadores.tsv de mentira."""
    llamadas = []

    def correr(script, *args, env_extra=None, capturar=False):
        llamadas.append((script, args))
        return 0

    with tempfile.TemporaryDirectory() as d:
        clon = pathlib.Path(d)
        (clon / "data").mkdir()
        (clon / "data" / "adaptadores.tsv").write_text(tabla)
        src = fuente(CELDA_3C).replace("FILAS = '''\n'''", f"FILAS = '''\n{filas}\n'''")
        ns = {"correr": correr, "PROYECTO": "gadmo/duplicado",
              "CLON": clon, "MALAS": malas}
        buf = io.StringIO()
        try:
            with redirect_stdout(buf):
                exec(src, ns)
        except AssertionError as e:
            ns["_error"] = str(e)
        ns["_tabla"] = (clon / "data" / "adaptadores.tsv").read_text()
    return ns, buf.getvalue(), llamadas


print("== 1. MALAS sale de la tabla, no de transcribirla")
ns, out, llam = correr_3b(TABLA, 1)
chk("las 6 VACIA, en orden", ns.get("MALAS") == MALAS_REALES, ns.get("MALAS"))
chk("y ninguna de las ok", not set(ns.get("MALAS", [])) & {"SRR3884830", "SRR3884814"})

print("== 2. y después perfila TODAS, una por una")
chk("llama a perfil --corridas",
    ("fetch_runs.sh", ("perfil", "--corridas", "gadmo/duplicado", "--tsv")) in llam, llam)

print("== 3. si verificar pasa, no re-mide nada")
ns2, out2, llam2 = correr_3b(TODO_OK, 0)
chk("MALAS vacía", ns2.get("MALAS") == [], ns2.get("MALAS"))
chk("no llama a perfil", not any(s == "fetch_runs.sh" for s, _ in llam2), llam2)
chk("y lo dice", "no hay nada que re-medir" in out2, out2)

print("== 4. DESVIADA cuenta igual que VACIA")
# Son los dos veredictos que verificar emite. Mirar solo VACIA dejaria pasar el
# caso en que el adaptador es casi correcto, que es el mas dificil de ver.
ns3, _, _ = correr_3b(DESVIADA, 1)
chk("agarra la DESVIADA", ns3.get("MALAS") == ["SRR_A1"], ns3.get("MALAS"))

print("== 5. si falla pero no se puede leer la tabla, corta en vez de seguir")
# Una tabla que cambie de forma no puede terminar en "0 corridas para rehacer"
# y un rehacer vacio que parece exito.
ns4, _, llam4 = correr_3b("algo salio mal y no hay tabla\n", 1)
chk("MALAS vacía", ns4.get("MALAS") == [], ns4.get("MALAS"))
chk("y revienta", "_error" in ns4, ns4.get("_error"))
chk("sin llamar a perfil", not any(s == "fetch_runs.sh" for s, _ in llam4), llam4)

print("== 6. §3c pega las filas y rehace SOLO las malas")
CAB = "org\tbioproject\trun\tfamilia\tsecuencia\tadapt_pct\tinserto_modal\tretencion_est\tveredicto\tfecha_utc\n"
BASE = CAB + "gadmo\tPRJNA328800\t-\tRA3\tTGGAATTCTCGGGTGCCAAGG\t84\t10\t51\tPARECE sRNA-seq\t2026-09-22\n"
NUEVA = "gadmo\tPRJNA328800\tSRR3884837\tIllumina_universal\tAGATCGGAAGAGCACACGTCT\t91\t22\t77\tPARECE sRNA-seq\t2026-09-24"
ns5, out5, llam5 = correr_3c(NUEVA, MALAS_REALES, BASE)
chk("la fila nueva queda en la tabla", NUEVA in ns5["_tabla"])
chk("y la del proyecto sigue", "PRJNA328800\t-\tRA3" in ns5["_tabla"])
chk("rehacer con las 6 corridas",
    ("trim.sh", tuple(["rehacer", "gadmo/duplicado"] + MALAS_REALES)) in llam5, llam5)

print("== 7. §3c no duplica una fila que ya está")
ns6, _, _ = correr_3c(NUEVA, MALAS_REALES, BASE + NUEVA + "\n")
chk("una sola vez", ns6["_tabla"].count(NUEVA) == 1, ns6["_tabla"].count(NUEVA))

print("== 8. sin FILAS no rehace nada")
# Rehacer sin haber corregido la tabla vuelve a recortar con la MISMA secuencia
# equivocada: el mismo resultado, otra vez, y pareciendo que se hizo algo.
ns7, out7, llam7 = correr_3c("", MALAS_REALES, BASE)
chk("no llama a rehacer", not any(s == "trim.sh" for s, _ in llam7), llam7)
chk("y lo dice", "FILAS está vacío" in out7, out7)

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    raise SystemExit(1)
print("TODO OK")
