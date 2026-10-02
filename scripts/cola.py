#!/usr/bin/env python3
"""Alinear en cola todos los proyectos que faltan, uno detras de otro.

Hace lo mismo que §2 a §6 de 20_alinear.ipynb, proyecto por proyecto, sin que
haya que volver a tocar el notebook entre uno y otro:

    genoma -> recorte -> verificar recorte -> alinear -> verificar BAM
           -> BAM y registros a Drive -> calibracion -> (push a git) -> limpiar

POR QUE ES REANUDABLE Y NO UNA SOLA CORRIDA LARGA
Lo que falta son ~110 h y Colab Free corta la sesion a las ~12. Cada proyecto
sube su BAM a Drive apenas termina, y la cola arranca mirando Drive: relanzarla
sigue desde el primero que no esta. Una sesion que se muere pierde el proyecto
en curso, no la cola. Y antes de empezar uno, si no alcanza a terminar dentro
de --max-horas, se detiene: arrancar algo de 8 h a la hora 10 es tirarlo.

LO QUE NO HACE A CIEGAS
- Si `trim.sh verificar` no pasa (una corrida VACIA o DESVIADA), ese proyecto
  NO se alinea: es el caso de gadmo_duplicado, que se alineo con 6 librerias
  vacias y llego a Drive como bueno. Se reporta y la cola sigue con el
  siguiente; ese se arregla a mano con §3b/§3c del notebook.
- Si `align.sh verificar` no pasa, el BAM no sube.
- Lo que no entra en la RAM de esta maquina (galga en Colab Free) se saltea y
  se nombra, en vez de morir con `Killed` a las seis horas.
- Lo que otra maquina tiene tomado (claims en Drive/90_claims) se saltea.

Desde el notebook (§8):   cola.correr_cola(clon=CLON, drive=DRIVE, env=ENV, ...)
Desde una maquina local:  ./scripts/cola.py --drive <ruta a tesis/> [--max-horas H]
                          (con las mismas variables de entorno que align.sh)
"""
import argparse
import csv
import datetime
import os
import pathlib
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import reparto  # noqa: E402

COLS_CALIB = ['proyecto', 'org', 'rol', 'accession', 'genoma_mb', 'librerias',
              'reads_trim', 'seg_reloj', 's_por_m', 'b_bam', 'b_fqgz', 'fecha_utc',
              'indice_en_reloj']


# ------------------------------------------------------------------ utilidades
def correr(clon, env, script, *args, capturar=False):
    """Corre un script del repo mostrando su salida en vivo. -> (rc, salida)."""
    p = subprocess.Popen([str(pathlib.Path(clon) / 'scripts' / script), *args],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, env=env, bufsize=1)
    lineas = []
    for ln in p.stdout:
        print(ln, end='')
        if capturar:
            lineas.append(ln)
    return p.wait(), ''.join(lineas)


def ram_disponible_b():
    """MemAvailable, no MemTotal: galga cae justo en esa diferencia."""
    try:
        for ln in open('/proc/meminfo'):
            if ln.startswith('MemAvailable:'):
                return int(ln.split()[1]) * 1024
    except OSError:
        pass
    return None


def accesiones(clon):
    out = {}
    led = pathlib.Path(clon) / 'data' / 'genomas.sha256'
    for ln in led.read_text().splitlines():
        c = ln.split('\t')
        if len(c) >= 2 and not ln.startswith('#') and c[0] != 'org':
            out[c[0]] = c[1]
    return out


def bases_genoma(genomes, org, acc):
    """Bases del ensamblado, del .fai o del cache que deja align.sh. 0 = no se sabe."""
    d = pathlib.Path(genomes) / org
    try:
        fai = d / f'{acc}.fna.fai'
        if fai.exists():
            return sum(int(l.split('\t')[1]) for l in fai.read_text().splitlines() if l)
        cache = d / f'{acc}.fna.bases'
        if cache.exists():
            return int(cache.read_text().strip())
    except (OSError, ValueError, IndexError):
        pass
    return 0


# ---------------------------------------------------------- la calibracion
def registrar_calibracion(pdir, org, rol, acc, gen_mb, t_alin, indice_en_reloj,
                          calib, nlib_def=0, b_fqgz_ref=None, b_bam_ref=None,
                          spm_pred=None):
    """Mide el proyecto recien alineado y escribe su fila en `calib`.

    UNA implementacion, la usan §5 del notebook y la cola. Devuelve la fila
    escrita, o None con el motivo impreso.

    Los reads son los que ENTRARON al alineamiento, de align/library_stats.txt
    —la misma fuente que `align.sh verificar`—, no los READS_OUT de
    recortadas.tsv: una PRE-TRIMMED no pasa por cutadapt y ahi dice '-'.
    gadmo_duplicado tiene 6 asi y se calibro con 16 M reads en vez de 55: 360 s/M
    en vez de 105, y B_BAM 44 en vez de 13. Por el mismo motivo las librerias se
    cuentan del BAM y el disco del recortado sale de las rutas del ledger, que
    para una PRE-TRIMMED es su fastq original y no un .t.fq.gz.
    """
    pdir, calib = pathlib.Path(pdir), pathlib.Path(calib)
    stats = pdir / 'align' / 'library_stats.txt'
    libs = [l.split('\t') for l in
            (stats.read_text().splitlines()[1:] if stats.exists() else [])
            if l.strip()]
    reads = sum(sum(int(x) for x in r[2:9]) for r in libs if len(r) >= 9)
    rec = pdir / 'recortadas.tsv'
    rutas = [l.split('\t')[2] for l in
             (rec.read_text().splitlines()[1:] if rec.exists() else []) if l.strip()]
    fq = sum((pdir / r).stat().st_size for r in rutas if (pdir / r).is_file())
    bam_f = pdir / 'align' / 'alignment.bam'
    bam = bam_f.stat().st_size if bam_f.exists() else 0
    nlib = len(libs) or nlib_def

    if not reads:
        print('\nsin reads en align/library_stats.txt: no puedo calibrar este proyecto')
        return None
    if t_alin is None:
        print('\nel alineamiento no corrio en esta sesion: no escribo calibracion '
              '(el reloj seria el de otro proyecto)')
        return None
    if not gen_mb:
        print('\nno pude medir el genoma: no escribo calibracion '
              '(la fila sin genoma no sirve para ajustar la recta)')
        return None

    spm = t_alin / (reads / 1e6)
    print(f'\nmedido en {org}_{rol}: {reads/1e6:.0f} M reads, genoma {gen_mb:.0f} Mb')
    print(f'  B_FQGZ  {fq/reads:.1f} B/read' + (f'     (§1 usa {b_fqgz_ref})' if b_fqgz_ref else ''))
    print(f'  B_BAM   {bam/reads:.1f} B/read' + (f'     (§1 usa {b_bam_ref})' if b_bam_ref else ''))
    print(f'  s/M     {spm:.0f} s/M read' + (f'   (§1 predijo {spm_pred:.0f})' if spm_pred else ''))

    # indice_en_reloj = si: el reloj incluye bowtie-build y el s/M no sirve
    # para la recta (gadmo_primario lo tuvo: 156 s/M, inflado). §1 lo saltea.
    en_reloj = 'si' if indice_en_reloj else 'no'
    nueva = dict(zip(COLS_CALIB, [
        f'{org}_{rol}', org, rol, acc or '?', f'{gen_mb:.0f}', str(nlib), str(reads),
        f'{t_alin:.0f}', f'{spm:.1f}', f'{bam/reads:.1f}', f'{fq/reads:.1f}',
        datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        en_reloj]))
    if en_reloj == 'si':
        print('  OJO: el reloj incluye la construccion del indice; la fila queda '
              'marcada y §1 no la usa para la recta')

    calib.parent.mkdir(parents=True, exist_ok=True)
    comentarios, filas = [], []
    if calib.exists():
        lineas = calib.read_text().splitlines()
        comentarios = [l for l in lineas if l.startswith('#')]
        # Por nombre de columna y no por posicion: una fila escrita con el
        # formato viejo (sin indice_en_reloj) no puede quedar corrida.
        for r in csv.DictReader((l for l in lineas if l and not l.startswith('#')),
                                delimiter='\t'):
            if r.get('proyecto') and r['proyecto'] != f'{org}_{rol}':
                filas.append({c: (r.get(c) or '?') for c in COLS_CALIB})
    filas.append(nueva)
    filas.sort(key=lambda r: r['proyecto'])
    calib.write_text('\n'.join(comentarios + ['\t'.join(COLS_CALIB)]
                               + ['\t'.join(r[c] for c in COLS_CALIB) for r in filas]) + '\n')
    print(f'  -> {calib.name}: {len(filas)} proyecto(s) calibrado(s)')
    return nueva


# ------------------------------------------------------------------ la cola
def pendientes(clon, drive, solo=()):
    """[(org, rol, costo)] que faltan en Drive, de menor a mayor."""
    cst, _aj = reparto.costos(clon, reparto.GENOMAS_MB)
    hechos = {k for k in cst
              if (pathlib.Path(drive) / '10_bam' / k[0] / f'{k[1]}.bam').exists()}
    quiero = {tuple(s.split('/')) for s in solo} if solo else None
    out = [(o, r, c) for (o, r), c in cst.items()
           if (o, r) not in hechos and (quiero is None or (o, r) in quiero)]
    return sorted(out, key=lambda t: t[2]['horas']), len(hechos)


def _bytes_por_read(clon):
    """El mayor medido, con las constantes de reparto como piso (igual que §1)."""
    bfq, bbam = reparto.B_FQGZ, reparto.B_BAM
    calib = pathlib.Path(clon) / 'data' / 'calibracion.tsv'
    if calib.exists():
        for r in csv.DictReader((l for l in calib.read_text().splitlines()
                                 if l and not l.startswith('#')), delimiter='\t'):
            try:
                bfq = max(bfq, float(r.get('b_fqgz') or 0))
                bbam = max(bbam, float(r.get('b_bam') or 0))
            except ValueError:
                pass
    return bfq, bbam


def _empujar(clon, mensaje):
    try:
        import colab_git  # noqa: PLC0415
        rutas = [p for p in ('data/alineamientos.tsv', 'data/calibracion.tsv',
                             'data/adaptadores.tsv')
                 if (pathlib.Path(clon) / p).exists()]
        print(colab_git.empujar(clon, rutas, mensaje, revisar=False))
        return True
    except Exception as e:  # noqa: BLE001
        # Lo que no llega a git queda en Drive (alineado.tsv y calibracion.tsv
        # al lado del BAM): se recupera de ahi, como ya se hizo una vez.
        print(f'   no pude empujar a git ({e.__class__.__name__}): las filas quedan '
              f'en Drive al lado del BAM')
        return False


def correr_cola(clon, drive, env, genomes, proy_dir, bam_dir, max_horas=11.0,
                empujar=False, solo=(), yo=None, ram_b=None, disco_libre=None,
                presupuesto_gb='20', reloj=time.monotonic):
    """Corre la cola. -> dict con hechos, fallados, saltados y si se corto por tiempo."""
    clon, drive = pathlib.Path(clon), pathlib.Path(drive)
    genomes, proy_dir, bam_dir = map(pathlib.Path, (genomes, proy_dir, bam_dir))
    yo = yo or os.uname().nodename
    claims = drive / '90_claims'
    acc = accesiones(clon)
    ram_b = ram_disponible_b() if ram_b is None else ram_b
    bfq, bbam = _bytes_por_read(clon)
    t0 = reloj()
    res = dict(hechos=[], fallados=[], saltados=[], cortado=False)

    cola, n_hechos = pendientes(clon, drive, solo)
    print(f'== cola: {len(cola)} pendiente(s), {n_hechos} ya en Drive. '
          f'Tope de esta sesion: {max_horas:.1f} h.')
    for o, r, c in cola:
        print(f'   {o}/{r:<10} {c["runs"]:>3} corridas  ~{c["horas"]:.1f} h  '
              f'RAM ~{c["ram_b"]/1e9:.1f} GB')

    for i, (org, rol, c) in enumerate(cola):
        p = f'{org}/{rol}'
        # --- lo que no se puede o no se debe empezar -------------------------
        if ram_b is not None and c['ram_b'] and c['ram_b'] + reparto.MARGEN_RAM > ram_b:
            print(f'\n== {p}: NO ENTRA EN RAM (pide ~{c["ram_b"]/1e9:.1f} GB, hay '
                  f'{ram_b/1e9:.1f}). Va a otra maquina.')
            res['saltados'].append((p, 'RAM'))
            continue
        libre = shutil.disk_usage(proy_dir).free if disco_libre is None else disco_libre
        pico = c['trim'] * (bfq + 2 * bbam)
        if pico + reparto.MARGEN_DISCO > libre:
            print(f'\n== {p}: NO ENTRA EN DISCO (pico ~{pico/1e9:.0f} GB, hay '
                  f'{libre/1e9:.0f}).')
            res['saltados'].append((p, 'disco'))
            continue
        tom = reparto.tomados(claims)
        if p in tom and tom[p][0] != yo:
            print(f'\n== {p}: lo tiene {tom[p][0]}; sigo con el proximo.')
            res['saltados'].append((p, f'tomado por {tom[p][0]}'))
            continue
        usado_h = (reloj() - t0) / 3600
        if usado_h + c['horas'] > max_horas:
            print(f'\n== {p}: ~{c["horas"]:.1f} h no entran en lo que queda de la sesion '
                  f'({max_horas - usado_h:.1f} h). Corto aca: relanzá la cola en una '
                  f'sesion nueva y sigue desde este.')
            res['cortado'] = True
            break
        if not reparto.tomar(claims, p, yo):
            res['saltados'].append((p, 'claim'))
            continue

        print(f'\n{"=" * 70}\n== {p}   ({i + 1} de {len(cola)}, ~{c["horas"]:.1f} h)\n{"=" * 70}')
        ok, motivo = _un_proyecto(clon, drive, env, genomes, proy_dir, bam_dir,
                                  org, rol, acc.get(org), c, presupuesto_gb, reloj)
        if ok:
            res['hechos'].append(p)
            if empujar:
                _empujar(clon, f'Alineado {org}_{rol} desde la cola')
        else:
            res['fallados'].append((p, motivo))
            print(f'\n!! {p}: {motivo}. Sigo con el proximo.')
        # El claim se suelta despues de subir el BAM (o de fallar): si no, queda
        # una ventana en la que nadie lo tiene y nadie lo hizo.
        reparto.soltar(claims, p, yo)
        _limpiar(proy_dir, bam_dir, genomes, org, rol,
                 sigue_org=any(o == org for o, _, _ in cola[i + 1:]))

    print(f'\n{"=" * 70}\n== resumen de la cola')
    print(f'   hechos   : {", ".join(res["hechos"]) or "-"}')
    for p, m in res['fallados']:
        print(f'   FALLO    : {p} — {m}')
    for p, m in res['saltados']:
        print(f'   saltado  : {p} — {m}')
    if res['cortado']:
        print('   cortado por tiempo: relanzá la cola en una sesion nueva.')
    return res


def _un_proyecto(clon, drive, env, genomes, proy_dir, bam_dir, org, rol, acc, c,
                 presupuesto_gb, reloj):
    p = f'{org}/{rol}'
    # Genoma: el .gz se copia de Drive al disco de la VM (align.sh descomprime y
    # bowtie indexa al lado: eso no se escribe al FUSE de Drive).
    src = drive / '70_genomas' / org
    if src.exists():
        shutil.copytree(src, genomes / org, dirs_exist_ok=True)
    rc, _ = correr(clon, env, 'align.sh', 'genoma', org)
    if rc != 0:
        return False, 'el genoma no se pudo preparar (sha256 o descompresion)'

    rc, _ = correr(clon, dict(env, PRESUPUESTO_GB=str(presupuesto_gb)),
                   'trim.sh', 'correr', p)
    if rc != 0:
        return False, 'el recorte fallo'
    rc, _ = correr(clon, env, 'trim.sh', 'verificar', p)
    if rc != 0:
        return False, ('el recorte no verifica (alguna corrida VACIA o DESVIADA): '
                       'no lo alineo. Arreglalo con §3b/§3c del notebook')

    # El indice esta completo solo si existe .indice_s (lo escribe align.sh
    # genoma DESPUES de que bowtie-build sale bien).
    en_reloj = not (genomes / org / f'{acc}.indice_s').exists() if acc else True
    t = reloj()
    rc, _ = correr(clon, env, 'align.sh', 'correr', p)
    t_alin = reloj() - t
    if rc != 0:
        return False, ('se quedo sin memoria (Killed)' if rc == 137
                       else f'el alineamiento fallo (rc={rc})')
    print(f'\n   alineamiento: {t_alin/60:.0f} min de reloj')

    rc, _ = correr(clon, env, 'align.sh', 'verificar', p)
    if rc != 0:
        return False, 'el BAM no paso la verificacion: no lo subo a Drive'

    destino = drive / '10_bam' / org
    destino.mkdir(parents=True, exist_ok=True)
    for suf in ('', '.bai'):
        f = bam_dir / org / f'{rol}.bam{suf}'
        if f.exists():
            shutil.copy2(f, destino / f.name)
            print(f'   {f.name} -> {destino}  ({f.stat().st_size/1e9:.1f} GB)')
    pdir = proy_dir / f'{org}_{rol}'
    led = pdir / 'alineado.tsv'
    if led.exists():
        shutil.copy2(led, destino / f'{rol}.alineado.tsv')

    gen_mb = bases_genoma(genomes, org, acc) / 1e6 if acc else 0
    calib = clon / 'data' / 'calibracion.tsv'
    fila = registrar_calibracion(pdir, org, rol, acc, gen_mb, t_alin, en_reloj, calib,
                                 nlib_def=c['runs'])
    if fila:
        # Copia en Drive al lado del BAM: si la sesion muere antes del push,
        # la fila se recupera de ahi y no de lo que se haya impreso.
        (destino / f'{rol}.calibracion.tsv').write_text(
            '\t'.join(COLS_CALIB) + '\n' + '\t'.join(fila[k] for k in COLS_CALIB) + '\n')
    correr(clon, env, 'align.sh', 'ledger', str(clon / 'data' / 'alineamientos.tsv'))
    return True, ''


def _limpiar(proy_dir, bam_dir, genomes, org, rol, sigue_org):
    """Libera el disco de la VM: lo que importa ya esta en Drive."""
    shutil.rmtree(proy_dir / f'{org}_{rol}', ignore_errors=True)
    for suf in ('', '.bai'):
        (bam_dir / org / f'{rol}.bam{suf}').unlink(missing_ok=True)
    if not sigue_org:
        shutil.rmtree(genomes / org, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--drive', required=True, help='ruta a Mi unidad/tesis (o su copia)')
    ap.add_argument('--max-horas', type=float, default=11.0)
    ap.add_argument('--empujar', action='store_true')
    ap.add_argument('--solo', nargs='*', default=[])
    a = ap.parse_args()
    raiz = pathlib.Path(__file__).resolve().parent.parent
    env = dict(os.environ)
    res = correr_cola(clon=raiz, drive=a.drive, env=env,
                      genomes=env.get('GENOMES_DIR', raiz / 'genomes'),
                      proy_dir=env.get('PROY_DIR', raiz / 'proyectos'),
                      bam_dir=env.get('BAM_DIR', raiz / 'bams'),
                      max_horas=a.max_horas, empujar=a.empujar, solo=a.solo)
    return 1 if res['fallados'] else 0


if __name__ == '__main__':
    sys.exit(main())
