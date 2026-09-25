#!/usr/bin/env python3
"""Repartir los 18 proyectos entre varias máquinas, sin que dos hagan el mismo.

POR QUÉ NO SE PARALELIZA DENTRO DE UNA MÁQUINA

`yasma align` ya corre `bowtie -p <cores>`: en una VM de Colab Free con 2 vCPU,
un solo proyecto ya las satura. Lanzar dos a la vez no divide el tiempo, lo
reparte — y encima **suma las dos `unique_d`**, que son 8 bytes por base del
genoma cada una. Dos proyectos de `gadmo` a la vez piden ~12.4 GB de los ~11.4
que tiene la VM: los dos mueren, y más tarde que uno solo.

Y un proyecto tampoco se puede partir en varias llamadas: `nativealign.py`
acumula `unique_d` sobre **todas** las librerías antes de que la etapa `multi`
la use para pesar, así que partirlo no cambia la contabilidad sino a qué locus
se asigna cada read multimapeado.

O sea que la unidad paralelizable es **el proyecto entero en otra máquina**.

EL TECHO ES 4x, Y NO LO MUEVE AGREGAR MÁQUINAS

Son ~101 h de trabajo total, pero `galga_duplicado` solo son ~25 h y no se
puede partir. El makespan nunca baja de ahí:

    1 worker  101 h      3 workers  34 h      5 workers  25 h
    2 workers  51 h      4 workers  25 h      8 workers  25 h

A partir de 4 máquinas, la quinta no hace nada. Eso se dice acá y no se
descubre habiendo conseguido cinco.

EL REPARTO ES LPT (largest processing time first)

Asignar de a uno del más caro al menos caro, cada uno a la máquina menos
cargada **que pueda hospedarlo**. Al revés —del más chico al más grande— todas
terminan los baratos y una queda sola con `galga_duplicado` al final.

La restricción que importa es la RAM: `galga` pide ~9.8 GB y no entra en Colab
Free. Una máquina que no puede con un proyecto no lo recibe, y si ninguna
puede, se dice cuál y por qué en vez de repartirlo igual.

DENTRO de cada máquina el orden vuelve a ser del más chico al más grande: no
cambia el makespan —es el mismo total— y hace que una sesión que se muere haya
dejado algo terminado.

Uso:

    ./scripts/reparto.py --maquinas 3
    ./scripts/reparto.py --maquinas colab:12:220 colab:12:220 local:64:1600
    ./scripts/reparto.py --maquinas 3 --hechos sclsc/duplicado

`<nombre>:<RAM GB>:<disco GB>`, con la RAM **disponible** (la que reporta
`MemAvailable`, no la nominal). Con un número suelto se asumen máquinas
iguales a las de Colab Free: 11.4 GB de RAM y 220 de disco, medidos.
"""
import argparse
import collections
import csv
import datetime
import pathlib
import sys

# Las mismas constantes medidas que usa §1 de 20_alinear.ipynb.
B_FQGZ = 22    # bytes por read en .t.fq.gz   (medido 21.8)
B_BAM = 16     # bytes por read en BAM        (medido 14.1)
B_BASE = 10    # bytes de RAM por base        (8 medidos de unique_d + ~1.5 del indice)
S_POR_M_PLANO = 60
MARGEN_DISCO = 5e9
MARGEN_RAM = 2e9

# RAM DISPONIBLE y disco libre, medidos en la VM real de Colab Free. La RAM va
# como la reporta MemAvailable de /proc/meminfo y NO como la nominal: una VM
# que se anuncia con 12 GB deja ~11.4 utilizables, y galga cae justo en esa
# diferencia — con 12 entra y con 11.4 no. El numero que hay que poner acá es
# el que imprime §1 del notebook o `align.sh plan`.
COLAB_FREE = (11.4, 220.0)


# ---------------------------------------------------------------- el modelo
def ajuste_s_por_m(calib):
    """Recta s/M = a + b*genoma_Mb, o None si no hay con que ajustarla.

    Con UN punto no se ajusta nada: el s/M plano sale de sclsc_duplicado, un
    genoma de 39 Mb, y galga es 1.05 Gb. Extrapolar de un punto a 27x es
    justamente lo que esto existe para no hacer.
    """
    pts = [(g, s) for g, s in calib]
    if len({g for g, _ in pts}) < 2:
        return None
    n = len(pts)
    sx = sum(g for g, _ in pts); sy = sum(s for _, s in pts)
    sxx = sum(g * g for g, _ in pts); sxy = sum(g * s for g, s in pts)
    den = n * sxx - sx * sx
    if not den:
        return None
    b = (n * sxy - sx * sy) / den
    return ((sy - b * sx) / n, b)


def s_por_m(mb, aj):
    return S_POR_M_PLANO if aj is None else max(1.0, aj[0] + aj[1] * mb)


def leer_calibracion(clon):
    pts = []
    f = pathlib.Path(clon) / 'data' / 'calibracion.tsv'
    if f.is_file():
        with open(f) as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t'):
                try:
                    pts.append((float(r['genoma_mb']), float(r['s_por_m'])))
                except (KeyError, ValueError, TypeError):
                    pass
    return pts


def costos(clon, genomas_mb, calib=None):
    """-> {(org, rol): dict} con runs, trim, horas, ram_b, pico_b, genoma_mb.

    `genomas_mb` es {org: Mb}. Un org que no este ahi queda con genoma 0, que
    quien llama tiene que tratar como "no se pudo medir" y NO como "no entra".
    """
    clon = pathlib.Path(clon)
    aj = ajuste_s_por_m(calib if calib is not None else leer_calibracion(clon))

    ret = {}
    with open(clon / 'data' / 'adaptadores.tsv') as f:
        for r in csv.DictReader((l for l in f if not l.startswith('#')), delimiter='\t'):
            ret[(r['org'], r['bioproject'])] = float(r['retencion_est'])

    proy = collections.defaultdict(lambda: [0, 0.0])
    with open(clon / 'data' / 'srr_manifest.tsv') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            k = (r['org'], r['rol'])
            proy[k][0] += 1
            proy[k][1] += int(r['read_count']) * ret.get((r['org'], r['bioproject']), 80) / 100

    out = {}
    for (org, rol), (n, kept) in proy.items():
        mb = float(genomas_mb.get(org, 0))
        out[(org, rol)] = dict(
            runs=n, trim=kept, genoma_mb=mb,
            horas=kept / 1e6 * s_por_m(mb, aj) / 3600,
            ram_b=mb * 1e6 * B_BASE,
            pico_b=kept * B_FQGZ + 2 * kept * B_BAM)
    return out, aj


# ---------------------------------------------------------------- el reparto
def repartir(trabajos, maquinas):
    """LPT con restricciones. Devuelve (asignacion, sin_lugar).

    `trabajos`: [(nombre, horas, ram_b, pico_b)]
    `maquinas`: [(nombre, ram_b, disco_b)]
    `asignacion`: {maquina: [(nombre, horas), ...]} ya ordenado de chico a
    grande, que es el orden en que conviene correrlos.
    `sin_lugar`: [(nombre, motivo)] — ninguna maquina puede con ellos.
    """
    carga = {m[0]: 0.0 for m in maquinas}
    asign = {m[0]: [] for m in maquinas}
    sin_lugar = []
    for nom, h, ram, pico in sorted(trabajos, key=lambda t: -t[1]):
        # Un genoma que no se pudo medir (ram == 0) NO se declara imposible:
        # mandar a la maquina grande algo que entra en cualquiera cuesta tanto
        # como lo contrario.
        aptas = [m for m in maquinas
                 if (ram == 0 or ram + MARGEN_RAM <= m[1])
                 and pico + MARGEN_DISCO <= m[2]]
        if not aptas:
            falta = 'RAM' if not any(ram == 0 or ram + MARGEN_RAM <= m[1] for m in maquinas) else 'disco'
            sin_lugar.append((nom, falta))
            continue
        elegida = min(aptas, key=lambda m: (carga[m[0]], m[0]))[0]
        carga[elegida] += h
        asign[elegida].append((nom, h))
    # De chico a grande DENTRO de cada maquina: no cambia el makespan y deja
    # algo terminado si la sesion se muere.
    for m in asign:
        asign[m].sort(key=lambda t: t[1])
    return asign, sin_lugar


def makespan(asign):
    return max((sum(h for _, h in v) for v in asign.values()), default=0.0)


# ---------------------------------------------------------------- los claims
# Coordinacion entre maquinas via Drive, que es lo unico que las tres ven. Un
# claim es un fichero con quien lo tomo y cuando; vence para que una sesion que
# se murio no deje el proyecto bloqueado para siempre.
#
# NO es exclusion mutua de verdad: el FUSE de Drive no da atomicidad. Lo que
# hace es releer lo que escribio y comprobar que siga siendo suyo, que atrapa
# el caso comun —dos maquinas arrancando con minutos de diferencia— y deja
# afuera el de dos escrituras en el mismo instante. El costo de una colision
# son horas perdidas, no datos corruptos: el BAM se escribe local y se copia a
# Drive al final, asi que lo peor es que el segundo pise al primero con el
# mismo contenido.
TTL_H = 12


def _f_claim(dir_claims, proyecto):
    return pathlib.Path(dir_claims) / (proyecto.replace('/', '_') + '.claim')


def _edad_h(f, ahora=None):
    ahora = ahora or datetime.datetime.now(datetime.timezone.utc)
    try:
        cuando = datetime.datetime.strptime(
            f.read_text().strip().split('\t')[1], '%Y-%m-%dT%H:%M:%SZ'
        ).replace(tzinfo=datetime.timezone.utc)
    except (OSError, IndexError, ValueError):
        return None
    return (ahora - cuando).total_seconds() / 3600


def tomar(dir_claims, proyecto, quien, ttl_h=TTL_H, ahora=None):
    """True si el proyecto queda tomado por `quien`."""
    d = pathlib.Path(dir_claims)
    d.mkdir(parents=True, exist_ok=True)
    f = _f_claim(d, proyecto)
    if f.exists():
        try:
            duenio = f.read_text().strip().split('\t')[0]
        except OSError:
            return False
        edad = _edad_h(f, ahora)
        if duenio == quien:
            pass                      # renovarlo es idempotente
        elif edad is not None and edad < ttl_h:
            return False              # de otro y vigente
        # vencido: se puede robar, y se avisa arriba
    ahora = ahora or datetime.datetime.now(datetime.timezone.utc)
    f.write_text(f"{quien}\t{ahora.strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
    try:
        return f.read_text().strip().split('\t')[0] == quien
    except OSError:
        return False


def soltar(dir_claims, proyecto, quien):
    f = _f_claim(dir_claims, proyecto)
    try:
        if f.exists() and f.read_text().strip().split('\t')[0] == quien:
            f.unlink()
            return True
    except OSError:
        pass
    return False


def tomados(dir_claims, ttl_h=TTL_H, ahora=None):
    """{proyecto: (quien, edad_h)} de los claims vigentes."""
    d = pathlib.Path(dir_claims)
    out = {}
    if not d.is_dir():
        return out
    for f in sorted(d.glob('*.claim')):
        try:
            quien = f.read_text().strip().split('\t')[0]
        except OSError:
            continue
        edad = _edad_h(f, ahora)
        if edad is not None and edad < ttl_h:
            org, _, rol = f.stem.partition('_')
            out[f'{org}/{rol}'] = (quien, edad)
    return out


# ---------------------------------------------------------------- CLI
GENOMAS_MB = {'sclsc': 39, 'cloro': 71, 'rhirr': 147, 'prupe': 228, 'phypa': 472,
              'maggi': 588, 'maldo': 650, 'gadmo': 670, 'galga': 1053}


def _maquina(spec, i):
    if ':' in spec:
        nom, ram, disco = (spec.split(':') + ['', ''])[:3]
        return (nom, float(ram) * 2**30, float(disco) * 1e9)
    return (f'w{i}', COLAB_FREE[0] * 2**30, COLAB_FREE[1] * 1e9)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--maquinas', nargs='+', default=['3'],
                    help='un número (máquinas tipo Colab Free) o <nombre>:<RAM disponible GB>:<disco libre GB>')
    ap.add_argument('--hechos', nargs='*', default=[],
                    help='proyectos ya terminados, como org/rol')
    ap.add_argument('--clon', default=str(pathlib.Path(__file__).resolve().parent.parent))
    args = ap.parse_args()

    if len(args.maquinas) == 1 and args.maquinas[0].isdigit():
        maquinas = [_maquina('', i) for i in range(1, int(args.maquinas[0]) + 1)]
    else:
        maquinas = [_maquina(s, i) for i, s in enumerate(args.maquinas, 1)]

    c, aj = costos(args.clon, GENOMAS_MB)
    hechos = {h.replace('_', '/') for h in args.hechos}
    trabajos = [(f'{o}/{r}', v['horas'], v['ram_b'], v['pico_b'])
                for (o, r), v in c.items() if f'{o}/{r}' not in hechos]

    tot = sum(t[1] for t in trabajos)
    mayor = max(trabajos, key=lambda t: t[1]) if trabajos else ('-', 0, 0, 0)
    print(f's/M read: ' + ('plano ' + str(S_POR_M_PLANO) if aj is None
                           else f'{aj[0]:.0f} + {aj[1]:.3f} x genoma_Mb (ajustado)'))
    print(f'{len(trabajos)} proyectos, {tot:.0f} h de trabajo total.')
    print(f'El más largo es {mayor[0]} con {mayor[1]:.1f} h, y NO se puede partir:')
    print(f'  ningún reparto baja de ahí. Con {len(maquinas)} máquinas la cota es '
          f'{max(mayor[1], tot/len(maquinas)):.1f} h.\n')

    asign, sin_lugar = repartir(trabajos, maquinas)
    ms = makespan(asign)
    for nom, ram, disco in maquinas:
        js = asign[nom]
        h = sum(x for _, x in js)
        print(f'--- {nom}  ({ram/2**30:.0f} GB RAM, {disco/1e9:.0f} GB disco)  '
              f'{h:.1f} h en {len(js)} proyecto(s)')
        print("    TANDA = [" + ', '.join(repr(n) for n, _ in js) + "]")
        for n, x in js:
            print(f'      {n:<20}{x:>6.1f} h')
    # OJO CON EL RESUMEN. Lo que NO se pudo asignar sigue habiendo que correrlo,
    # asi que dividir el trabajo total por el makespan de lo asignado da un
    # speedup inventado — daba 4.5x mientras los dos galga, que son el 33% de
    # los reads del set, no estaban asignados a ninguna maquina. Aca el
    # speedup se calcula SOLO sobre lo repartido, y lo de afuera se dice
    # aparte con sus horas.
    h_sin = sum(t[1] for t in trabajos if t[0] in {n for n, _ in sin_lugar})
    h_asig = tot - h_sin
    if sin_lugar:
        print('\nNINGUNA máquina puede con estos, y siguen habiendo que correrlos:')
        for n, falta in sin_lugar:
            hn = next(t[1] for t in trabajos if t[0] == n)
            print(f'  {n:<20}{hn:>6.1f} h   le falta {falta}')
        print(f'  Son {h_sin:.0f} h que este reparto NO cubre. Van a una máquina más')
        print('  grande (Colab Pro high-RAM o la local); volvé a correr esto con')
        print('  esa máquina en --maquinas para que entren en el plan.')

    print(f'\nmakespan de lo repartido: {ms:.1f} h'
          f'  ({h_asig/ms:.1f}x contra las {h_asig:.0f} h en serie)')
    if sin_lugar:
        print(f'De punta a punta son al menos {max(ms, h_sin):.1f} h: '
              f'las {h_sin:.0f} h de arriba van en paralelo SOLO si conseguís esa máquina.')

    sobran = sum(1 for m in maquinas if not asign[m[0]])
    asignados = [(n, h) for v in asign.values() for n, h in v]
    if sobran:
        print(f'{sobran} máquina(s) sin trabajo: sobran para este set.')
    elif asignados and len(maquinas) > 1:
        # El piso lo marca el mas largo DE LO ASIGNADO, no el del set entero.
        mayor_asig = max(asignados, key=lambda t: t[1])
        if ms <= mayor_asig[1] * 1.01:
            print(f'Ya estás en el piso de lo repartido: lo marca {mayor_asig[0]} '
                  f'sola ({mayor_asig[1]:.1f} h). Otra máquina igual no ayuda.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
