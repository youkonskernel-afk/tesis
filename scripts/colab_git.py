#!/usr/bin/env python3
"""Traer el repo a Colab y empujar de vuelta lo que hay que versionar.

Hasta ahora Colab solo LEIA de GitHub. El preambulo de los notebooks clona con
`fetch --depth 1 && reset --hard` y lo dice explicito: "el clon es un CACHE del
repo, no un espacio de trabajo". Lo que Colab produce y tiene que quedar
versionado —las filas nuevas de data/srr_manifest.tsv, data/sra_md5.tsv,
data/genomas.sha256— se imprimia para COPIAR Y PEGAR a mano. Este modulo cierra
ese circulo.

Vive en scripts/ y no adentro de una celda por dos motivos: el preambulo de los
notebooks ya declara que "los notebooks llaman a los scripts del repo en vez de
reimplementarlos", y asi se puede probar sin Colab y sin red
(tests/test_colab_git.py levanta un repo bare de verdad y empuja contra el).

TRES REGLAS QUE NO SON NEGOCIABLES

1. NUNCA `git add -A`. La regla del proyecto es que la data no entra al repo, y
   una VM de Colab tiene Drive montado en /content/drive con ~190 GB de .sra.
   .gitignore cubre lo obvio, pero un `add -A` es exactamente donde se cuela lo
   que no. Solo se agregan las rutas que se pasan, y se rechaza cualquiera que
   no este bajo data/.

2. EL TOKEN NUNCA SE IMPRIME. git mete la URL en sus mensajes de error, y esa
   URL lleva el token adentro. Todo lo que se muestre pasa por sin_token().

3. UN PUSH RECHAZADO FALLA FUERTE. Dos sesiones de Colab que escriben el mismo
   ledger divergen. Un push rechazado que nadie mira deja el resultado en Drive
   y no en git, que es justo el estado que este modulo existe para evitar. Se
   reintenta UNA vez rebasando sobre lo que haya, y si vuelve a fallar, revienta.
"""
import os
import pathlib
import shutil
import subprocess

REPO = "youkonskernel-afk/tesis"
URL_ANON = f"https://github.com/{REPO}.git"

AYUDA_TOKEN = (
    "No hay GITHUB_TOKEN en los Secrets de Colab.",
    "Para empujar hace falta un PAT CON permiso de escritura:",
    "  1. github.com -> Settings -> Developer settings -> Personal access",
    "     tokens -> Fine-grained tokens -> Generate new token",
    "  2. Repository access: solo " + REPO,
    "  3. Permissions -> Repository permissions -> Contents: Read and write",
    "  4. Guardalo como GITHUB_TOKEN en el panel de Secrets de Colab (la llave",
    "     a la izquierda) y habilita el acceso para este notebook.",
    "El token NO va al repo ni a Drive: vive en los Secrets de la sesion.",
)


def sin_token(txt, secreto):
    """git incluye la URL en sus mensajes de error, y esa URL lleva el token."""
    return txt.replace(secreto, "***") if secreto else txt


def token():
    """El PAT, de los Secrets de Colab o del entorno. None si no hay."""
    try:
        from google.colab import userdata  # noqa: PLC0415
        t = userdata.get("GITHUB_TOKEN")
        if t:
            return t.strip()
    except Exception:
        pass
    t = os.environ.get("GITHUB_TOKEN", "").strip()
    return t or None


def _git(clon, *args, tok=None):
    r = subprocess.run(["git", "-C", str(clon), *args],
                       capture_output=True, text=True)
    r.stdout = sin_token(r.stdout, tok)
    r.stderr = sin_token(r.stderr, tok)
    return r


# EL CLON NO ESTA ACA, Y ES A PROPOSITO
# Este modulo vive ADENTRO del repo, asi que no se puede importar antes de
# clonarlo: el clon es el bootstrap y tiene que estar en la celda. Duplicarlo
# aca dejaria dos implementaciones de lo mismo —la trampa que este proyecto ya
# se comio con las tres rutas de los .sra— y la del modulo no correria nunca.
# La celda de clon esta en los 5 notebooks y es identica en todos;
# scripts/validate_notebooks.py falla si una deriva, y tests/test_clon.py la
# ejercita con un git falso.


def rama_actual(clon):
    r = _git(clon, "rev-parse", "--abbrev-ref", "HEAD")
    n = r.stdout.strip()
    # Un clon --depth 1 puede quedar en HEAD desacoplado tras el reset.
    if r.returncode != 0 or not n or n == "HEAD":
        r = _git(clon, "rev-parse", "--abbrev-ref", "origin/HEAD")
        n = r.stdout.strip().replace("origin/", "")
    return n or "main"


# --- empujar ------------------------------------------------------------------

def _validar(clon, rutas):
    """Solo ficheros de texto bajo data/. Ver la regla 1 de la cabecera."""
    limpias, malas = [], []
    for r in rutas:
        p = pathlib.PurePosixPath(str(r).replace("\\", "/"))
        if p.is_absolute():
            try:
                p = pathlib.PurePosixPath(
                    os.path.relpath(str(p), str(pathlib.Path(clon).resolve())))
            except ValueError:
                malas.append((str(r), "esta fuera del clon"))
                continue
        s = str(p)
        if ".." in p.parts:
            malas.append((s, "sube de directorio"))
        elif p.parts[:1] != ("data",):
            malas.append((s, "no esta bajo data/"))
        elif not (pathlib.Path(clon) / s).is_file():
            malas.append((s, "no existe en el clon"))
        else:
            limpias.append(s)
    if malas:
        raise ValueError(
            "rutas que no se pueden empujar:\n"
            + "\n".join(f"  {s}: {por}" for s, por in malas)
            + "\nSolo ficheros de texto bajo data/. La data pesada va a Drive:"
              " ver la regla de ubicacion en CLAUDE.md.")
    return limpias


def _destino(clon, tok):
    """A donde empujar: `origin`, o la URL con el token si origin ES GitHub.

    Mirar el origin y no cablear la URL no es cosmetico. Con la URL cableada,
    cualquier clon —el de un banco, el de una prueba a mano— empujaria al repo
    de verdad apenas hubiera un GITHUB_TOKEN en el entorno, sin importar contra
    que remoto se clono. Lo encontro su propio banco al primer escenario.
    """
    url = _git(clon, "remote", "get-url", "origin").stdout.strip()
    if tok and "github.com" in url:
        return f"https://x-access-token:{tok}@github.com/{REPO}.git"
    return "origin"


def _deshallow(clon, tok):
    # Un clon --depth 1 no tiene base de merge, asi que el rebase del reintento
    # no podria correr. El repo pesa menos de 1 MB: traerlo entero sale gratis.
    if (pathlib.Path(clon) / ".git" / "shallow").exists():
        _git(clon, "fetch", "--unshallow", "origin", tok=tok)


def empujar(clon, rutas, mensaje, rama=None, revisar=True, autor=None):
    """Commitea SOLO `rutas` (bajo data/) y empuja.

    revisar=True (el default) imprime el diff y NO empuja, igual que el dry-run
    de drive_push.sh: mover cosas no puede ser lo que pasa si te equivocas de
    celda. Devuelve una linea para imprimir.
    """
    clon = pathlib.Path(clon)
    limpias = _validar(clon, rutas)
    rama = rama or rama_actual(clon)
    tok = token()

    # `add --` y las rutas explicitas: nunca -A. Ver la regla 1.
    r = _git(clon, "add", "--", *limpias, tok=tok)
    if r.returncode != 0:
        raise RuntimeError("git add fallo: " + r.stderr)

    hay = _git(clon, "diff", "--cached", "--name-only", tok=tok).stdout.strip()
    if not hay:
        return "nada que commitear: git ya tiene lo que hay en Drive"

    diff = _git(clon, "diff", "--cached", "--stat", tok=tok).stdout.strip()
    if revisar:
        _git(clon, "reset", tok=tok)     # dejar el indice como estaba
        return ("REVISAR_PRIMERO esta en True, asi que no empuje nada.\n"
                f"Esto es lo que iria a la rama {rama}:\n{diff}\n"
                "Si esta bien, pone REVISAR_PRIMERO = False y corre de nuevo.")

    if autor:
        _git(clon, "config", "user.name", autor[0], tok=tok)
        _git(clon, "config", "user.email", autor[1], tok=tok)
    r = _git(clon, "commit", "-m", mensaje, tok=tok)
    if r.returncode != 0:
        raise RuntimeError("git commit fallo: " + r.stderr + r.stdout)

    destino = _destino(clon, tok)
    r = _git(clon, "push", destino, f"HEAD:{rama}", tok=tok)
    if r.returncode == 0:
        return f"empujado a {rama}:\n{diff}"

    # Rechazado: casi siempre otra sesion de Colab movio la rama. Se reintenta
    # UNA vez rebasando encima, y si vuelve a fallar revienta (regla 3).
    _deshallow(clon, tok)
    if _git(clon, "fetch", destino, rama, tok=tok).returncode == 0:
        _git(clon, "rebase", "FETCH_HEAD", tok=tok)
        r2 = _git(clon, "push", destino, f"HEAD:{rama}", tok=tok)
        if r2.returncode == 0:
            return f"empujado a {rama} (rebasado sobre lo que habia):\n{diff}"
        r = r2

    raise RuntimeError(
        "el push fue rechazado dos veces. El commit esta HECHO en el clon, que "
        "es efimero: si cerras la sesion se pierde.\n" + r.stderr + r.stdout
        + ("\n" + "\n".join(AYUDA_TOKEN) if not tok else ""))
