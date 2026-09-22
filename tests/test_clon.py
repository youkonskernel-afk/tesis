#!/usr/bin/env python3
"""Banco de pruebas de la celda de clon, con un git falso en el PATH.
Misma tecnica que ya se uso para el filtro de RNA y el flujo de prefetch."""
import json, os, pathlib, shutil, subprocess, sys, tempfile

import pathlib
RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB = str(RAIZ / 'notebooks' / '10_descarga_runs.ipynb')
nb = json.load(open(NB, encoding='utf-8'))
CODE = next(''.join(c['source']) for c in nb['cells']
            if c['cell_type'] == 'code' and 'URL_ANON' in ''.join(c['source']))

FALLAS = 0
def chk(n, ok, extra=''):
    global FALLAS
    print(('  ok   ' if ok else '  MAL  ') + n + (f'  [{extra}]' if not ok and extra else ''))
    if not ok: FALLAS += 1

GIT = r'''#!/usr/bin/env python3
import os, pathlib, sys
a = sys.argv[1:]
if a[0] == '-C': a = a[2:]             # git -C <dir> ...
modo = os.environ.get('GIT_FAKE', 'ok')
log = pathlib.Path(os.environ['GIT_LOG'])
log.write_text(log.read_text() + ' '.join(a) + '\n')
if a[0] == 'clone':
    if modo == 'clone_falla_anon' and 'x-access-token' not in a[-2]:
        sys.stderr.write('remote: Repository not found.\n'); sys.exit(128)
    if modo == 'clone_falla_todo':
        sys.stderr.write('fatal: could not read from ' + a[-2] + '\n'); sys.exit(128)
    d = pathlib.Path(a[-1]); (d / 'data').mkdir(parents=True, exist_ok=True)
    (d / 'data' / 'sra_md5.tsv').write_text('limpio\n')
    sys.exit(0)
if a[0] == 'fetch':
    sys.exit(1 if modo == 'fetch_falla' else 0)
if a[0] == 'reset':
    # reset --hard: restaura el archivo versionado, pise lo que pise
    for p in pathlib.Path.cwd().rglob('sra_md5.tsv'): p.write_text('limpio\n')
    sys.exit(0)
if a[0] == 'log':
    print('abc1234 commit falso'); sys.exit(0)
sys.exit(0)
'''

def correr(modo, clon_previo=None, token=None):
    tmp = pathlib.Path(tempfile.mkdtemp())
    (tmp / 'bin').mkdir()
    g = tmp / 'bin' / 'git'; g.write_text(GIT); g.chmod(0o755)
    clon = tmp / 'tesis'
    if clon_previo is not None:
        (clon / 'data').mkdir(parents=True)
        (clon / 'data' / 'sra_md5.tsv').write_text(clon_previo)
    logf = tmp / 'git.log'; logf.write_text('')
    env = dict(os.environ, PATH=f"{tmp/'bin'}:{os.environ['PATH']}",
               GIT_FAKE=modo, GIT_LOG=str(logf))
    salida, err = [], None
    g_ns = {'CLON': clon, 'print': lambda *a: salida.append(' '.join(map(str, a)))}
    viejo = os.environ.copy(); os.environ.update(env)
    cwd = os.getcwd(); os.chdir(tmp)
    if token is not None:
        mod = type(sys)('google.colab')
        mod.userdata = type(sys)('userdata'); mod.userdata.get = lambda k: token
        sys.modules['google.colab'] = mod
        sys.modules.setdefault('google', type(sys)('google')).colab = mod
    try:
        exec(CODE, g_ns)
    except Exception as e:
        err = str(e)
    finally:
        os.chdir(cwd); os.environ.clear(); os.environ.update(viejo)
        sys.modules.pop('google.colab', None); sys.modules.pop('google', None)
    return '\n'.join(salida), err, logf.read_text(), clon

print('== 1. no hay clon -> clona de cero')
out, err, log, clon = correr('ok')
chk('sin excepcion', err is None, err)
chk('clona', 'clone' in log)
chk('reporta el clon anonimo', 'anonimo' in out, out)
chk('imprime el commit', 'abc1234' in out, out)

print('== 2. clon limpio -> se actualiza, no re-clona')
out, err, log, clon = correr('ok', clon_previo='limpio\n')
chk('sin excepcion', err is None, err)
chk('usa fetch', 'fetch --depth 1 origin HEAD' in log)
chk('usa reset --hard', 'reset --hard FETCH_HEAD' in log)
chk('no re-clona', 'clone' not in log)
chk('dice actualizado', 'actualizado' in out, out)

print('== 3. REGRESION: clon con un archivo versionado modificado')
out, err, log, clon = correr('ok', clon_previo='SUCIO: 414 filas viejas\n')
chk('sin excepcion', err is None, err)
chk('se actualiza igual', 'actualizado' in out, out)
chk('el archivo quedo restaurado',
    (clon / 'data' / 'sra_md5.tsv').read_text() == 'limpio\n',
    (clon / 'data' / 'sra_md5.tsv').read_text())

print('== 4. no se puede actualizar -> tira el clon y re-clona, no sigue viejo')
out, err, log, clon = correr('fetch_falla', clon_previo='SUCIO\n')
chk('sin excepcion', err is None, err)
chk('re-clona', 'clone' in log)
chk('el clon quedo limpio',
    (clon / 'data' / 'sra_md5.tsv').read_text() == 'limpio\n')

print('== 5. repo privado, sin token -> levanta excepcion, no sigue')
out, err, log, clon = correr('clone_falla_todo')
chk('levanta excepcion', err is not None)
chk('el mensaje explica las dos salidas', err and 'GITHUB_TOKEN' in err, err)

print('== 6. repo privado con token -> clona y NO filtra el token')
TOK = 'ghp_TOKENSUPERSECRETO123'
out, err, log, clon = correr('clone_falla_anon', token=TOK)
chk('sin excepcion', err is None, err)
chk('clona con token', 'con token' in out, out)
chk('el token NO aparece en la salida', TOK not in out)
chk('resetea el remoto sin token', 'remote set-url origin https://github.com/' in log)

print('== 7. el token tampoco se filtra cuando el clon con token falla')
out, err, log, clon = correr('clone_falla_todo', token=TOK)
chk('levanta excepcion', err is not None)
chk('el token NO aparece en el error', err and TOK not in err, err)

print()
print('TODO OK' if FALLAS == 0 else f'{FALLAS} fallas')
sys.exit(1 if FALLAS else 0)
