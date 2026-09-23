#!/usr/bin/env python3
"""
Chequea que lo que AFIRMAN los docs coincida con lo que hay en el repo.

Por que existe: en este proyecto ya paso cinco veces que un limite o un dato
viejo de la herramienta se viera identico a un hecho — la paginacion cortada en
20 de 78 ensamblados, la cepa que no se leia del biosample, el 0% de adaptador,
el ledger que no espejaba el manifiesto, y el total del manifiesto que coincidia
por casualidad. Los docs tienen el mismo problema: la tabla del dataset en
CLAUDE.md decia que el duplicado de sclsc era PRJNA985401 cuando la spec ya
decia PRJNA1135930, y el README mandaba a correr cinco comandos de los cuales
cuatro no existen. Nadie recalcula una afirmacion escrita a mano.

    ./scripts/check_docs.py

Lo que se puede comprobar sin ambiguedad falla con exit 1. Los conteos en prosa
se imprimen para comparar a ojo, porque reescribir una frase no deberia romper
el chequeo.
"""
import pathlib
import re
import sys

RAIZ = pathlib.Path(__file__).resolve().parent.parent
DOCS = ["README.md", "CLAUDE.md"]

# Ficheros citados en los docs que a proposito NO estan en el repo.
# Cada entrada es deuda: el dia que el fichero llegue, se borra la linea y el
# chequeo vuelve a exigirlo.
CITAS_SIN_FICHERO = {
    "gen_manifest.sh": "script de R1, reemplazado por scripts/fetch_runs.sh; "
                       "se cita como historia, no como algo que se corra",
    "config.sh": "vive en el main local, todavia sin subir",
    "orchestrate.sh": "idem config.sh",
    "verify.sh": "idem config.sh",
    "check_env.sh": "idem config.sh",
    "environment.yml": "idem config.sh",
    "environment_pip.txt": "idem config.sh",
}

# Ficheros que el pipeline PRODUCE al correr y que, por la regla de ubicacion del
# proyecto, no van a estar nunca en git. No son deuda: no hay nada que esperar.
# Se listan aparte para que CITAS_SIN_FICHERO siga queriendo decir "esto tendria
# que llegar algun dia".
CITAS_DE_SALIDA = {
    "recortadas.tsv": "ledger que escribe scripts/trim.sh en trim/<org>_<rol>/",
    "log.txt": "el log de cutadapt que deja yasma trim",
    "inputs.json": "lo escribe trim.sh y lo pisa yasma; es estado de corrida",
    "loci.gff3": "salida de yasma tradeoff",
    "library_stats.txt": "conteos por read group que escribe yasma align",
    "alineado.tsv": "contra que genoma alineo cada proyecto; lo escribe align.sh",
}

# Ficheros de una dependencia externa, citados para decir DONDE se midio algo.
# Tampoco son deuda: no van a estar nunca en este repo. Llevan el repo al lado
# para que se sepa contra que se verifican.
CITAS_EXTERNAS = {
    "nativealign.py": "NateyJay/YASMA@v1.1.1 — el `yasma align` real",
    "align.py": "NateyJay/YASMA@v1.1.1 — el wrapper de ShortStack, comentado",
    "__init__.py": "NateyJay/YASMA@v1.1.1 — es donde align.py esta comentado",
    "trim.py": "NateyJay/YASMA@v1.1.1",
    "generics.py": "NateyJay/YASMA@v1.1.1",
}

ORGS = ["rhirr", "sclsc", "cloro", "phypa", "prupe", "maldo", "gadmo", "galga",
        "maggi"]


def sin_comentarios(path):
    """organismos.tsv y excluidas.tsv arrancan con comentarios '#', asi que la
    primera linea del fichero no es el encabezado."""
    return [ln for ln in path.read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")]


def tsv(path, con_comentarios=True):
    lineas = sin_comentarios(path) if con_comentarios else path.read_text().splitlines()
    cab = lineas[0].split("\t")
    return [dict(zip(cab, ln.split("\t"))) for ln in lineas[1:]]


def citas_rotas(fallos):
    """Un fichero citado que no existe. El README mandaba a correr ./verify.sh,
    ./check_env.sh y ./gen_manifest.sh, ninguno de los tres en el repo."""
    # sha256 va ANTES que sh, o la alternancia corta 'genomas.sha256' en
    # 'genomas.sh' y el chequeo inventa un fichero que nadie cito. El \b final
    # evita el mismo error con cualquier extension que sea prefijo de otra.
    pat = re.compile(r"[A-Za-z0-9_/.-]+\.(?:sha256|sh|py|yml|txt|ipynb|tsv)\b")
    vistos = {}
    for d in DOCS:
        for m in pat.finditer((RAIZ / d).read_text()):
            vistos.setdefault(m.group(0).lstrip("./"), set()).add(d)

    for cita, donde in sorted(vistos.items()):
        base = cita.rsplit("/", 1)[-1]
        existe = any((RAIZ / p / cita).exists() or (RAIZ / p / base).exists()
                     for p in ("", "scripts", "data", "docs", "notebooks", "tests"))
        if base in CITAS_DE_SALIDA or base in CITAS_EXTERNAS:
            continue
        if existe:
            if base in CITAS_SIN_FICHERO:
                fallos.append(
                    f"{cita} ya existe pero sigue en CITAS_SIN_FICHERO de este "
                    f"script — borra esa linea y el chequeo vuelve a exigirlo")
        elif base not in CITAS_SIN_FICHERO:
            fallos.append(f"{', '.join(sorted(donde))} cita {cita}, que no existe")


def tabla_dataset(fallos):
    """La tabla de CLAUDE.md contra data/organismos.tsv.

    Es el bug que disparo este script: la tabla decia PRJNA985401 para el
    duplicado de sclsc y la prosa 130 lineas mas abajo decia PRJNA1135930.
    Quien lee la tabla no llega a la prosa."""
    spec = {}
    for r in tsv(RAIZ / "data" / "organismos.tsv"):
        spec.setdefault((r["org"], r["rol"]), []).append(r["bioproject"])

    texto = (RAIZ / "CLAUDE.md").read_text()
    for org in ORGS:
        m = re.search(rf"^\|\s*{org}\s*\|(.*)$", texto, re.M)
        if not m:
            fallos.append(f"CLAUDE.md no tiene fila de tabla para {org}")
            continue
        cols = [c.strip() for c in m.group(1).split("|")]
        if len(cols) < 4:
            fallos.append(f"CLAUDE.md: la fila de {org} tiene {len(cols)} columnas")
            continue
        for i, rol in ((2, "primario"), (3, "duplicado")):
            # El primario de maggi son dos BioProjects: la tabla los escribe
            # 'PRJNA154615 + PRJNA232734', igual que dos filas en la spec.
            dice = {a.strip() for a in cols[i].split("+") if a.strip()}
            esta = set(spec.get((org, rol), []))
            if dice != esta:
                fallos.append(
                    f"CLAUDE.md: {org} {rol} dice {sorted(dice) or ['nada']}, "
                    f"organismos.tsv dice {sorted(esta) or ['nada']}")


def data_consistente(fallos, hechos):
    """El manifiesto contra la spec, las exclusiones y el ledger de md5."""
    man = tsv(RAIZ / "data" / "srr_manifest.tsv", con_comentarios=False)
    corridas = {r["run"] for r in man}
    proy_man = {r["bioproject"] for r in man}

    spec = tsv(RAIZ / "data" / "organismos.tsv")
    proy_spec = {r["bioproject"] for r in spec}
    huerfanos = sorted(proy_spec - proy_man)

    excl = [ln.split("\t")[0] for ln in
            sin_comentarios(RAIZ / "data" / "excluidas.tsv")[1:]]
    adentro = [r for r in excl if r in corridas]

    led = tsv(RAIZ / "data" / "sra_md5.tsv", con_comentarios=False)
    en_led = {r["run"] for r in led}

    hechos += [
        f"organismos.tsv : {len(spec)} filas, {len(proy_spec)} BioProjects",
        f"manifiesto     : {len(man)} corridas, {len(proy_man)} BioProjects",
        f"ledger md5     : {len(led)} filas "
        f"({sum(1 for r in led if r['formato'] == 'sralite')} sralite)",
        f"excluidas      : {len(excl)}",
        f"genomas        : {len(tsv(RAIZ / 'data' / 'genomas.sha256', False))} sha256",
    ]

    if adentro:
        fallos.append(f"el manifiesto tiene corridas excluidas adentro: {adentro}")
    if huerfanos:
        # Es el estado de hoy: sclsc gano PRJNA1135930 en la spec y el
        # manifiesto todavia no se regenero contra la ENA.
        fallos.append(
            f"BioProjects de la spec sin ninguna corrida en el manifiesto: "
            f"{huerfanos} — regenerar con fetch_runs.sh manifest (necesita red)")
    if corridas - en_led:
        fallos.append(f"en el manifiesto y no en el ledger: {sorted(corridas - en_led)}")
    if en_led - corridas:
        fallos.append(f"en el ledger y no en el manifiesto: {sorted(en_led - corridas)}")

    gen = tsv(RAIZ / "data" / "genomas.tsv")
    con_sha = {r["org"] for r in tsv(RAIZ / "data" / "genomas.sha256", False)}
    sin_sha = sorted({r["org"] for r in gen} - con_sha)
    if sin_sha:
        fallos.append(f"organismos en genomas.tsv sin sha256: {sin_sha}")


def numeros_en_prosa():
    """Los conteos escritos a mano, para comparar contra los hechos de arriba.

    No falla: reescribir una frase no deberia romper el chequeo, y un chequeo
    que falla por algo que no importa se termina ignorando."""
    pat = re.compile(r"\b(\d{2,4})\b(?=[^.\n]{0,40}\b"
                     r"(?:corridas?|md5|proyectos?|filas|\.sra)\b)")
    filas = []
    for d in DOCS:
        for i, ln in enumerate((RAIZ / d).read_text().splitlines(), 1):
            for m in pat.finditer(ln):
                filas.append((d, i, m.group(1), ln.strip()[:88]))
    return filas


def main():
    fallos, hechos = [], []
    citas_rotas(fallos)
    tabla_dataset(fallos)
    data_consistente(fallos, hechos)

    print("== hechos, recalculados de data/")
    for h in hechos:
        print(f"   {h}")

    print()
    print("== conteos escritos a mano en los docs (comparar con lo de arriba)")
    for d, i, n, ln in numeros_en_prosa():
        print(f"   {d}:{i}  {n:>4}  {ln}")

    print()
    if fallos:
        print(f"== {len(fallos)} problema(s)")
        for f in fallos:
            print(f"   MAL  {f}")
        return 1
    print("== los docs coinciden con data/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
