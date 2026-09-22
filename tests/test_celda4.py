#!/usr/bin/env python3
"""
La celda §4 de 10_descarga_runs.ipynb, ejercitada sin red ni Colab.

Es la celda cuya salida se copia a git, o sea donde un error hace dano de
verdad: ya ofrecio TRES veces la fila de una corrida excluida para commitear,
que la habria devuelto al manifiesto sin que nadie lo notara.

Dos cosas que este banco fija:
  1. El guardia: con una corrida excluida adentro del manifiesto no imprime
     NADA, corta.
  2. Que compare contra lo que tiene GIT, no contra el working tree. §1 pisa
     CLON/data/srr_manifest.tsv con la copia de Drive, asi que leer el fichero
     del clon da el manifiesto de Drive y el diff sale vacio siempre. Por eso
     el escenario arma un repo git real y commitea la version vieja.

    tests/test_celda4.py
"""
import io
import json
import pathlib
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout

RAIZ = pathlib.Path(__file__).resolve().parent.parent
NB = RAIZ / 'notebooks' / '10_descarga_runs.ipynb'
CELDA = 17

FALLAS = 0


def ok(m):
    print(f'  ok   {m}')


def mal(m):
    global FALLAS
    FALLAS += 1
    print(f'  MAL  {m}')


HDR_MAN = ('org\trun\tbioproject\trol\tset_modelo\tread_count\tbase_count\t'
           'avg_len\tstrategy\tlayout\tsource')
HDR_LED = 'org\trun\tmd5\tformato\tfecha_utc'


def fila_man(run, proy='PRJ_A', org='aa'):
    return f'{org}\t{run}\t{proy}\tprimario\tapl\t100\t100\t50\tmiRNA-Seq\tSINGLE\tTRANSCRIPTOMIC'


def fila_led(run, org='aa'):
    return f'{org}\t{run}\t{"a" * 32}\tsra\t2026-09-22'


def git(clon, *args):
    subprocess.run(['git', '-C', str(clon)] + list(args),
                   capture_output=True, text=True, check=True)


def escenario(tmp, man_git, led_git, man_uso, led_drive):
    """man_git/led_git son lo COMMITEADO; man_uso/led_drive lo que hay en Drive."""
    clon = tmp / 'clon'
    if clon.exists():
        subprocess.run(['rm', '-rf', str(clon)], check=True)
    (clon / 'data').mkdir(parents=True)

    (clon / 'data' / 'excluidas.tsv').write_text(
        '# corridas que no entran\n'
        'run\tmotivo\tfecha_utc\n'
        'SRR_MALA\tes mRNA\t2026-09-22\n')
    (clon / 'data' / 'srr_manifest.tsv').write_text('\n'.join([HDR_MAN] + man_git) + '\n')
    (clon / 'data' / 'sra_md5.tsv').write_text('\n'.join([HDR_LED] + led_git) + '\n')

    git(clon, 'init', '-q')
    git(clon, 'config', 'user.email', 't@t')
    git(clon, 'config', 'user.name', 't')
    git(clon, 'add', '-A')
    git(clon, 'commit', '-q', '-m', 'estado commiteado')

    # §1 pisa el fichero del clon con la copia de Drive: el working tree deja de
    # decir que hay en git, y por eso la celda tiene que usar `git show`.
    (clon / 'data' / 'srr_manifest.tsv').write_text('\n'.join([HDR_MAN] + man_uso) + '\n')

    ledger = tmp / 'ledger_drive.tsv'
    ledger.write_text('\n'.join([HDR_LED] + led_drive) + '\n')

    nb = json.loads(NB.read_text())
    src = ''.join(nb['cells'][CELDA]['source'])
    ns = {'CLON': clon, 'LEDGER_DRIVE': ledger}
    buf, exc = io.StringIO(), None
    try:
        with redirect_stdout(buf):
            exec(src, ns)
    except Exception as e:      # noqa: BLE001 — cualquier fallo se reporta
        exc = e
    return buf.getvalue(), exc


def main():
    with tempfile.TemporaryDirectory() as d:
        tmp = pathlib.Path(d)

        print('== 1. el caso de sclsc: Drive tiene 2 corridas que git no')
        out, exc = escenario(
            tmp,
            man_git=[fila_man('SRR_A1')],
            led_git=[fila_led('SRR_A1')],
            man_uso=[fila_man('SRR_A1'), fila_man('SRR_B1', 'PRJ_B', 'bb'),
                     fila_man('SRR_B2', 'PRJ_B', 'bb')],
            led_drive=[fila_led('SRR_A1'), fila_led('SRR_B1', 'bb'),
                       fila_led('SRR_B2', 'bb')])
        (ok if exc is None else mal)(f'sin excepcion ({exc})')
        for frag, etiq in (
            ('manifiesto en uso : 3 corridas', 'cuenta el manifiesto en uso'),
            ('manifiesto en git : 1 corridas', 'y el de git, que es el viejo'),
            ('AGREGÁ estas filas a data/srr_manifest.tsv', 'pide las filas del manifiesto'),
            ('AGREGÁ estas filas a data/sra_md5.tsv', 'y las del ledger'),
        ):
            (ok if frag in out else mal)(etiq)
        # las 2 filas nuevas del manifiesto, y NINGUNA de las viejas
        bloque_man = out.split('AGREGÁ estas filas a data/srr_manifest.tsv ---')[1]
        bloque_man = bloque_man.split('---')[0]
        (ok if 'SRR_B1' in bloque_man and 'SRR_B2' in bloque_man else mal)(
            'las dos corridas nuevas')
        (ok if 'SRR_A1' not in bloque_man else mal)(
            'y NO repite la que git ya tiene')

        print('== 2. compara contra git, no contra el working tree')
        # Si leyera el fichero del clon (que §1 ya piso con la copia de Drive),
        # el diff daria vacio y no pediria nada. Es el bug que esto fija.
        (ok if 'nada que commitear' not in out else mal)(
            'no dice "nada que commitear" teniendo 2 filas de diferencia')

        print('== 3. todo sincronizado: no pide nada')
        out, exc = escenario(
            tmp,
            man_git=[fila_man('SRR_A1')], led_git=[fila_led('SRR_A1')],
            man_uso=[fila_man('SRR_A1')], led_drive=[fila_led('SRR_A1')])
        (ok if exc is None else mal)(f'sin excepcion ({exc})')
        (ok if 'nada que commitear' in out else mal)('lo dice')
        (ok if 'AGREGÁ' not in out else mal)('y no pide ninguna fila')

        print('== 4. una corrida excluida adentro: CORTA sin imprimir filas')
        out, exc = escenario(
            tmp,
            man_git=[fila_man('SRR_A1')], led_git=[fila_led('SRR_A1')],
            man_uso=[fila_man('SRR_A1'), fila_man('SRR_MALA')],
            led_drive=[fila_led('SRR_A1'), fila_led('SRR_MALA')])
        if isinstance(exc, RuntimeError) and 'corridas excluidas adentro' in str(exc):
            ok('corta con RuntimeError')
        else:
            mal(f'corta con RuntimeError (dio {exc!r})')
        if 'AGREGÁ' in out or 'SRR_MALA' in out:
            mal('y NO llega a imprimir ninguna fila')
        else:
            ok('y NO llega a imprimir ninguna fila')

        print('== 5. una fila del ledger que ya no esta en el manifiesto')
        # El caso de SRR23277331: excluida del manifiesto, pero su md5 sigue en
        # el ledger de Drive. No se commitea, se reporta aparte.
        out, exc = escenario(
            tmp,
            man_git=[fila_man('SRR_A1')], led_git=[fila_led('SRR_A1')],
            man_uso=[fila_man('SRR_A1')],
            led_drive=[fila_led('SRR_A1'), fila_led('SRR_VIEJA')])
        (ok if exc is None else mal)(f'sin excepcion ({exc})')
        (ok if 'NO commitear' in out else mal)('la reporta como no commiteable')
        (ok if 'Corré §3' in out else mal)('y dice como sacarla')
        bloque_led = out.split('AGREGÁ estas filas a data/sra_md5.tsv')
        if len(bloque_led) > 1 and 'SRR_VIEJA' in bloque_led[1]:
            mal('y NO la ofrece para commitear')
        else:
            ok('y NO la ofrece para commitear')

        print('== 6. una fila que git tiene y el manifiesto nuevo ya no')
        out, exc = escenario(
            tmp,
            man_git=[fila_man('SRR_A1'), fila_man('SRR_FUERA')],
            led_git=[fila_led('SRR_A1')],
            man_uso=[fila_man('SRR_A1')], led_drive=[fila_led('SRR_A1')])
        (ok if exc is None else mal)(f'sin excepcion ({exc})')
        (ok if 'SACÁ estas de data/srr_manifest.tsv' in out else mal)(
            'pide sacarla, en vez de ignorarla')
        (ok if 'SRR_FUERA' in out else mal)('y la nombra')

    print()
    if FALLAS:
        print(f'{FALLAS} fallas')
        return 1
    print('TODO OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
