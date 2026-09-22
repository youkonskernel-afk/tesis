#!/usr/bin/env python3
"""
La celda §1 de 10_descarga_runs.ipynb, ejercitada sin red ni Colab.

Por que tiene banco propio: esa celda estuvo mal tres veces seguidas.
  1. Reusaba la copia de Drive SIEMPRE, asi que cambiar organismos.tsv o
     excluidas.tsv no tenia ningun efecto visible.
  2. Despues detectaba la copia vieja pero solo la imprimia, y la celda seguia:
     las tres siguientes corrian sobre el manifiesto viejo y §4 terminaba
     ofreciendo la fila excluida para commitear.
  3. Despues cortaba con RuntimeError, que frena el dano pero deja a la persona
     moviendo un flag a mano para algo que la celda ya sabe hacer.

Ahora regenera sola cuando el chequeo falla, y solo falla si la ENA sigue sin
coincidir con la spec despues de consultarla — que ya no es una copia vieja.

    tests/test_celda1.py
"""
import io
import json
import pathlib
import sys
import tempfile
from contextlib import redirect_stdout

RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB = RAIZ / 'notebooks' / '10_descarga_runs.ipynb'
CELDA = 7

FALLAS = 0


def ok(m):
    print(f'  ok   {m}')


def mal(m):
    global FALLAS
    FALLAS += 1
    print(f'  MAL  {m}')


def fuente_celda():
    nb = json.loads(NB.read_text())
    return ''.join(nb['cells'][CELDA]['source'])


HDR_MAN = ('org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\t'
           'avg_len\tstrategy\tlayout\tsource')


def fila_man(org, run, proy):
    return f'{org}\t{run}\t{proy}\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC'


def escenario(tmp, manifiesto_drive, salida_ena, regenerar=False):
    """Arma un clon falso y corre la celda. Devuelve (stdout, excepcion o None).

    `manifiesto_drive`: filas (sin cabecera) de la copia de Drive, o None si no
    hay copia. `salida_ena`: lo que el fetch_runs.sh falso escribe cuando se lo
    llama, o None para que falle."""
    clon = tmp / 'clon'
    (clon / 'data').mkdir(parents=True, exist_ok=True)
    (clon / 'scripts').mkdir(parents=True, exist_ok=True)

    # spec: 2 proyectos, y una corrida excluida
    (clon / 'data' / 'organismos.tsv').write_text(
        '# comentario largo que no es el encabezado\n'
        '# otra linea de comentario\n'
        'org\tespecie\treino\tclado\trol\tbioproject\truns\tspots_M\tstrategy\tassembly\tnota\tset_modelo\n'
        'aa\tE uno\tFungi\tAsco\tprimario\tPRJ_A\t2\t10\tmiRNA-Seq\tA\tn\taplicacion\n'
        'bb\tE dos\tFungi\tAsco\tduplicado\tPRJ_B\t1\t5\tmiRNA-Seq\tA\tn\taplicacion\n')
    (clon / 'data' / 'excluidas.tsv').write_text(
        '# corridas que no entran\n'
        'run\tmotivo\tfecha_utc\n'
        'SRR_MALA\tes mRNA\t2026-09-22\n')

    # fetch_runs.sh falso: escribe el manifiesto que le digamos
    sh = clon / 'scripts' / 'fetch_runs.sh'
    if salida_ena is None:
        sh.write_text('#!/usr/bin/env bash\necho "la ENA no respondio" >&2\nexit 1\n')
    else:
        cuerpo = '\\n'.join([HDR_MAN] + salida_ena)
        sh.write_text('#!/usr/bin/env bash\n'
                      f'printf \'{cuerpo}\\n\' > data/srr_manifest.tsv\n'
                      'echo "consulte la ENA"\n')
    sh.chmod(0o755)

    man_drive = tmp / 'drive' / 'srr_manifest.tsv'
    man_drive.parent.mkdir(parents=True, exist_ok=True)
    if manifiesto_drive is not None:
        man_drive.write_text('\n'.join([HDR_MAN] + manifiesto_drive) + '\n')
    elif man_drive.exists():
        man_drive.unlink()

    ns = {'CLON': clon, 'MANIFIESTO_DRIVE': man_drive, 'REGENERAR': regenerar}
    src = fuente_celda()
    # REGENERAR se pasa por el namespace, asi que se saca la asignacion literal
    src = '\n'.join(ln for ln in src.split('\n')
                    if not ln.startswith('REGENERAR = '))
    buf = io.StringIO()
    exc = None
    try:
        with redirect_stdout(buf):
            exec(src, ns)
    except Exception as e:      # noqa: BLE001 — cualquier fallo se reporta
        exc = e
    return buf.getvalue(), exc


AL_DIA = [fila_man('aa', 'SRR_A1', 'PRJ_A'), fila_man('bb', 'SRR_B1', 'PRJ_B')]
# la copia vieja: tiene la excluida adentro y le falta PRJ_B
VIEJA = [fila_man('aa', 'SRR_A1', 'PRJ_A'), fila_man('aa', 'SRR_MALA', 'PRJ_A')]


def main():
    with tempfile.TemporaryDirectory() as d:
        tmp = pathlib.Path(d)

        print('== 1. copia al dia: la reusa y no consulta la ENA')
        out, exc = escenario(tmp, AL_DIA, AL_DIA)
        (ok if exc is None else mal)(f'sin excepcion ({exc})')
        (ok if 'lo reuso.' in out else mal)('dice que la reusa')
        if 'consulte la ENA' in out:
            mal('NO consulta la ENA')
        else:
            ok('NO consulta la ENA')
        if 'coincide con organismos.tsv' in out:
            ok('y confirma que coincide')
        else:
            mal('y confirma que coincide')

        print('== 2. copia vieja: REGENERA SOLA, sin pedir que se mueva un flag')
        out, exc = escenario(tmp, VIEJA, AL_DIA)
        (ok if exc is None else mal)(f'no corta ({exc})')
        for frag, etiq in (
            ('hay corridas excluidas adentro: SRR_MALA', 'detecta la excluida'),
            ('BioProjects de la spec sin ninguna corrida: PRJ_B', 'y el proyecto que falta'),
            ('La regenero sola', 'lo dice y regenera'),
            ('consulte la ENA', 'consulta la ENA de verdad'),
            ('coincide con organismos.tsv', 'y termina coincidiendo'),
        ):
            (ok if frag in out else mal)(etiq)
        if 'REGENERAR = True' in out:
            mal('no manda a mover el flag a mano')
        else:
            ok('no manda a mover el flag a mano')

        print('== 3. sin copia en Drive: la genera')
        out, exc = escenario(tmp, None, AL_DIA)
        (ok if exc is None else mal)(f'sin excepcion ({exc})')
        (ok if 'no hay manifiesto en Drive' in out else mal)('lo dice')
        (ok if 'consulte la ENA' in out else mal)('consulta la ENA')

        print('== 4. REGENERAR = True fuerza la consulta aunque la copia sirva')
        out, exc = escenario(tmp, AL_DIA, AL_DIA, regenerar=True)
        (ok if exc is None else mal)(f'sin excepcion ({exc})')
        (ok if 'consulto la ENA igual' in out else mal)('lo dice')
        (ok if 'consulte la ENA' in out else mal)('consulta la ENA')

        print('== 5. la ENA sigue sin traer un proyecto: AHORA si corta')
        # regenerar no alcanza porque la ENA no devuelve PRJ_B
        out, exc = escenario(tmp, VIEJA, [fila_man('aa', 'SRR_A1', 'PRJ_A')])
        if isinstance(exc, RuntimeError) and 'no tienen ninguna corrida' in str(exc):
            ok('corta con el motivo correcto')
        else:
            mal(f'corta con el motivo correcto (dio {exc!r})')
        # Lo que importa es que NO mande a mover el flag: regenerar ya se
        # intento y no alcanzo. (Buscar 'copia vieja' aca daba falso positivo:
        # el mensaje dice 'No es una copia vieja', que es justo lo correcto.)
        if 'REGENERAR' in str(exc):
            mal('y no manda a poner REGENERAR = True')
        else:
            ok('y no manda a poner REGENERAR = True')
        if 'fetch_runs.sh buscar' in str(exc):
            ok('sino a buscar un reemplazo, que es lo que corresponde')
        else:
            mal('sino a buscar un reemplazo, que es lo que corresponde')

        print('== 6. la ENA sigue devolviendo la excluida: corta distinto')
        out, exc = escenario(tmp, VIEJA, AL_DIA + [fila_man('aa', 'SRR_MALA', 'PRJ_A')])
        if isinstance(exc, RuntimeError) and 'corridas excluidas' in str(exc):
            ok('corta señalando excluidas.tsv')
        else:
            mal(f'corta señalando excluidas.tsv (dio {exc!r})')

        print('== 7. si la ENA falla, no pisa la copia buena de Drive')
        out, exc = escenario(tmp, VIEJA, None)
        if isinstance(exc, RuntimeError) and 'manifest fallo' in str(exc):
            ok('corta antes de copiar')
        else:
            mal(f'corta antes de copiar (dio {exc!r})')
        drive = (tmp / 'drive' / 'srr_manifest.tsv').read_text()
        if 'SRR_MALA' in drive:
            ok('la copia de Drive quedo intacta')
        else:
            mal('la copia de Drive quedo intacta')

    print()
    if FALLAS:
        print(f'{FALLAS} fallas')
        return 1
    print('TODO OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
