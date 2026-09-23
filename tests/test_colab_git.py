#!/usr/bin/env python3
"""scripts/colab_git.py contra un repo git DE VERDAD, sin red.

Se levanta un repo bare en un tmpdir y se clona contra el, asi que el push se
prueba empujando en serio en vez de mirar los argumentos que se le pasan a un
git falso. Eso importa: las tres reglas que este modulo tiene que garantizar se
afirman sobre el CONTENIDO del commit que llego al remoto, no sobre el texto que
se imprime al lado.

  1. Nunca `git add -A`. Una VM de Colab tiene Drive montado con ~190 GB de
     .sra; un add -A es exactamente donde se cuela lo que no va al repo.
  2. El token no aparece en ningun mensaje. git mete la URL —con el token
     adentro— en sus errores.
  3. Un push rechazado falla fuerte. Si no, el resultado queda en Drive y no en
     git, que es el estado que este modulo existe para evitar.
"""
import os
import pathlib
import subprocess
import sys
import tempfile

RAIZ = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "scripts"))
import colab_git  # noqa: E402

FALLAS = 0


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


def git(d, *a):
    return subprocess.run(["git", "-C", str(d), *a], capture_output=True, text=True)


def escenario(tmp):
    """Un bare como 'GitHub' y un clon de trabajo, los dos de verdad."""
    bare = tmp / "remoto.git"
    subprocess.run(["git", "init", "-q", "--bare", "-b", "main", str(bare)], check=True)

    semilla = tmp / "semilla"
    semilla.mkdir()
    subprocess.run(["git", "init", "-q", "-b", "main", str(semilla)], check=True)
    for k, v in (("user.name", "t"), ("user.email", "t@t")):
        git(semilla, "config", k, v)
    (semilla / "data").mkdir()
    (semilla / "data" / "sra_md5.tsv").write_text("org\trun\n")
    (semilla / "scripts").mkdir()
    (semilla / "scripts" / "x.sh").write_text("#!/bin/sh\n")
    git(semilla, "add", "-A")
    git(semilla, "commit", "-qm", "inicial")
    git(semilla, "push", "-q", str(bare), "main")

    clon = tmp / "clon"
    subprocess.run(["git", "clone", "-q", str(bare), str(clon)], check=True)
    for k, v in (("user.name", "c"), ("user.email", "c@c")):
        git(clon, "config", k, v)
    return bare, clon


def en_remoto(bare, ruta, rama="main"):
    r = git(bare, "show", f"{rama}:{ruta}")
    return r.stdout if r.returncode == 0 else None


with tempfile.TemporaryDirectory() as td:
    tmp = pathlib.Path(td)
    os.environ.pop("GITHUB_TOKEN", None)

    print("== 1. revisar=True muestra el diff y NO empuja")
    bare, clon = escenario(tmp / "a")
    (clon / "data" / "sra_md5.tsv").write_text("org\trun\naa\tSRR1\n")
    s = colab_git.empujar(clon, ["data/sra_md5.tsv"], "prueba", revisar=True)
    chk("lo dice", "no empuje nada" in s, s)
    chk("y muestra qué iría", "sra_md5.tsv" in s, s)
    chk("el remoto sigue igual", en_remoto(bare, "data/sra_md5.tsv") == "org\trun\n")
    chk("y el índice queda limpio",
        git(clon, "diff", "--cached", "--name-only").stdout.strip() == "")

    print("== 2. revisar=False empuja de verdad")
    s = colab_git.empujar(clon, ["data/sra_md5.tsv"], "agrega SRR1", revisar=False)
    chk("lo dice", "empujado a main" in s, s)
    chk("y el remoto lo tiene", en_remoto(bare, "data/sra_md5.tsv") == "org\trun\naa\tSRR1\n")

    print("== 3. NUNCA add -A: lo que no se pide no se commitea")
    # La regla que mas importa: una VM de Colab tiene Drive montado al lado.
    bare, clon = escenario(tmp / "b")
    (clon / "data" / "sra_md5.tsv").write_text("org\trun\nbb\tSRR2\n")
    (clon / "scripts" / "x.sh").write_text("#!/bin/sh\nrm -rf /\n")   # sucio
    (clon / "enorme.sra").write_bytes(b"x" * 1024)                    # data suelta
    colab_git.empujar(clon, ["data/sra_md5.tsv"], "solo el ledger", revisar=False)
    chk("el fichero pedido llegó",
        en_remoto(bare, "data/sra_md5.tsv") == "org\trun\nbb\tSRR2\n")
    chk("el script sucio NO",  en_remoto(bare, "scripts/x.sh") == "#!/bin/sh\n")
    chk("y el .sra tampoco",   en_remoto(bare, "enorme.sra") is None)

    print("== 4. rutas que no se pueden empujar")
    for ruta, por in [("scripts/x.sh", "fuera de data/"),
                      ("data/../scripts/x.sh", "con .."),
                      ("data/noexiste.tsv", "que no existe")]:
        try:
            colab_git.empujar(clon, [ruta], "no", revisar=False)
            chk(f"rechaza una ruta {por}", False, "no lanzó")
        except ValueError as e:
            chk(f"rechaza una ruta {por}", True)
            chk(f"  y dice cuál ({por})", ruta.split("/")[-1] in str(e), str(e))

    print("== 5. sin cambios no commitea")
    s = colab_git.empujar(clon, ["data/sra_md5.tsv"], "nada", revisar=False)
    chk("lo dice", "nada que commitear" in s, s)
    antes = git(bare, "rev-list", "--count", "main").stdout.strip()
    colab_git.empujar(clon, ["data/sra_md5.tsv"], "nada", revisar=False)
    chk("y no deja un commit vacío",
        git(bare, "rev-list", "--count", "main").stdout.strip() == antes)

    print("== 6. si el remoto se movió, rebasa y empuja igual")
    # Dos sesiones de Colab escribiendo el mismo ledger. Sin el reintento, el
    # push se rechaza y el resultado queda en Drive y no en git.
    bare, clon = escenario(tmp / "c")
    otro = tmp / "c" / "otro"
    subprocess.run(["git", "clone", "-q", str(bare), str(otro)], check=True)
    for k, v in (("user.name", "o"), ("user.email", "o@o")):
        git(otro, "config", k, v)
    (otro / "data" / "otro.tsv").write_text("de la otra sesión\n")
    git(otro, "add", "data/otro.tsv"); git(otro, "commit", "-qm", "otra sesión")
    git(otro, "push", "-q", "origin", "main")

    (clon / "data" / "sra_md5.tsv").write_text("org\trun\ncc\tSRR3\n")
    s = colab_git.empujar(clon, ["data/sra_md5.tsv"], "desde la primera", revisar=False)
    chk("lo dice", "rebasado" in s, s)
    chk("el remoto tiene lo mío",
        en_remoto(bare, "data/sra_md5.tsv") == "org\trun\ncc\tSRR3\n")
    chk("y NO perdió lo del otro",
        en_remoto(bare, "data/otro.tsv") == "de la otra sesión\n")

    print("== 7. un push que no se puede resolver revienta")
    bare2 = tmp / "c" / "no_existe.git"
    (clon / "data" / "sra_md5.tsv").write_text("org\trun\ndd\tSRR4\n")
    git(clon, "remote", "set-url", "origin", str(bare2))
    try:
        colab_git.empujar(clon, ["data/sra_md5.tsv"], "al vacío", revisar=False)
        chk("lanza", False, "no lanzó")
    except RuntimeError as e:
        chk("lanza", True)
        chk("y avisa que el commit ya está hecho", "commit esta HECHO" in str(e), str(e))

    print("== 8. el token no se filtra, y no se empuja a GitHub por error")
    TOK = "ghp_UnTokenDeVerdadQueNoDebeAparecer"
    os.environ["GITHUB_TOKEN"] = TOK
    chk("token() lo lee del entorno", colab_git.token() == TOK)
    chk("sin_token lo tapa", colab_git.sin_token(f"url https://x:{TOK}@gh", TOK)
        == "url https://x:***@gh")
    # git mete la URL —con el token adentro— en sus errores. localhost:1 falla
    # al instante y sin red, que es justo lo que hace falta para verlo.
    r = colab_git._git(clon, "fetch", f"https://x-access-token:{TOK}@localhost:1/x.git",
                       tok=TOK)
    chk("git de verdad falló", r.returncode != 0)
    chk("y su stderr no trae el token", TOK not in r.stderr, r.stderr[:200])

    # Con un origin local, un GITHUB_TOKEN en el entorno NO puede hacer que se
    # empuje al repo de verdad. Con la URL cableada, sí podía.
    bare, clon8 = escenario(tmp / "d")
    (clon8 / "data" / "sra_md5.tsv").write_text("org\trun\nee\tSRR5\n")
    s = colab_git.empujar(clon8, ["data/sra_md5.tsv"], "con token y origin local",
                          revisar=False)
    chk("empuja al origin local igual", "empujado a main" in s, s)
    chk("y el bare local lo tiene",
        en_remoto(bare, "data/sra_md5.tsv") == "org\trun\nee\tSRR5\n")
    chk("_destino elige origin, no github", colab_git._destino(clon8, TOK) == "origin")
    os.environ.pop("GITHUB_TOKEN", None)

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    sys.exit(1)
print("TODO OK")
