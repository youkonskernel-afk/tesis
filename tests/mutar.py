#!/usr/bin/env python3
"""
Rompe el codigo a proposito y verifica que algun banco se ponga rojo.

Un banco que pasa no prueba nada: prueba algo el dia que se rompe lo que cubre
y el banco se queja. En este proyecto ya aparecieron dos huecos que solo se
vieron asi — el banco de `perfil` verificaba el mensaje `>>> YA RECORTADA` pero
no la variable que decide el exit code, y `estado` no tenia banco.

    tests/mutar.py              # todas las mutaciones
    tests/mutar.py adaptador    # solo las que matcheen ese texto

Cada mutacion que salga HUECO es un chequeo que falta. Una que diga PATRON NO
ENCONTRADO es una mutacion que quedo vieja: el codigo cambio y hay que
reescribirla o sacarla (una mutacion que no se aplica tampoco prueba nada).
"""
import pathlib
import shutil
import subprocess
import sys
import tempfile

RAIZ = pathlib.Path(__file__).resolve().parent.parent

# (nombre, fichero, [(viejo, nuevo), ...]). Varios pares = mutacion combinada,
# para el caso en que dos rutas del codigo se cubren entre si.
MUTACIONES = [
    ("filtro: acepta GENOMIC", "scripts/fetch_runs.sh",
     [('[[ "$src" == "$FUENTE_OK" ]] || return 1', ': # roto')]),
    ("exclusiones: no se aplican", "scripts/fetch_runs.sh",
     [('if grep -qx -- "$run" "$excl"; then', 'if false; then')]),
    ("ledger: no poda lo que no esta en el manifiesto", "scripts/fetch_runs.sh",
     [('NR > 1 && !($2 in en) { print $2 }', 'NR > 1 && 0 { print $2 }')]),
    ("prefetch: salida a /dev/null otra vez", "scripts/fetch_runs.sh",
     [('eco_log', 'true #')]),
    # Las dos rutas al .sralite se tapan entre si, asi que hay que matar ambas.
    ("sralite: ni candidatos ni find", "scripts/fetch_runs.sh",
     [('"$dir/$org/$run.sralite" "$dir/$run/$run.sralite" "$dir/$run.sralite"',
       '"$dir/$run.NADA"'),
      ("-o -name '*.sralite'", "-o -name '*.NADA'")]),
    ("registrar_md5: sin columna formato", "scripts/fetch_runs.sh",
     [('fmt="${4:-sra}"', 'fmt=""')]),
    ("perfil: YA RECORTADA vuelve a NO PARECE", "scripts/fetch_runs.sh",
     [('        ver = "YA RECORTADA"', '        ver = "NO PARECE sRNA-seq"')]),
    ("perfil: el adaptador no se separa del nombre", "scripts/fetch_runs.sh",
     [('split(campo[i], par, ":"); A[i] = par[1]; NOM[i] = par[2]',
       'A[i] = campo[i]; NOM[i] = campo[i]')]),
    ("perfil: el adaptador se cae del RESUMEN", "scripts/fetch_runs.sh",
     [('ver, ad_top > "/dev/stderr"', 'ver > "/dev/stderr"')]),
    ("estado: ruta_sra pierde una disposicion", "scripts/fetch_runs.sh",
     [('for c in "$dir/$org/$run.sra" "$dir/$run/$run.sra" "$dir/$run.sra"; do',
       'for c in "$dir/$org/$run.sra"; do')]),
    ("estado: acepta ficheros vacios", "scripts/fetch_runs.sh",
     [('for c in "$dir/$org/$run.sra" "$dir/$run/$run.sra" "$dir/$run.sra"; do\n'
       '    [[ -s "$c" ]]',
       'for c in "$dir/$org/$run.sra" "$dir/$run/$run.sra" "$dir/$run.sra"; do\n'
       '    [[ -e "$c" ]]')]),
    ("estado: sin manifiesto sigue adelante", "scripts/fetch_runs.sh",
     [('[[ -f "$MANIFEST" ]] || die "no existe $MANIFEST — corré: $0 manifest"\n'
       '  echo "manifiesto: $MANIFEST"', 'echo "manifiesto: $MANIFEST"')]),
    ("cepas: vuelve a una sola pagina de 20", "scripts/fetch_genomes.sh",
     [('page_size=100${tok:+&page_token=$tok}', 'page_size=20')]),
    ("cepas: la cepa solo del campo organism", "scripts/fetch_genomes.sh",
     [('def bioattr', 'def _bioattr_sin_usar')]),
]


def correr(mut_dir):
    r = subprocess.run(['./tests/run_all.sh'], cwd=mut_dir,
                       capture_output=True, text=True)
    for ln in r.stdout.splitlines():
        if 'chequeos:' in ln:
            # 'N bancos, M chequeos: X ok, Y fallas' -> Y
            return int(ln.rsplit(' ok, ', 1)[1].split()[0]), ln
    return None, '(run_all.sh no imprimio resumen)'


def main():
    filtro = sys.argv[1] if len(sys.argv) > 1 else ''
    base = tempfile.mkdtemp(prefix='mutar.')
    try:
        shutil.copytree(RAIZ / 'tests', pathlib.Path(base) / 'tests')
        huecos = viejas = 0
        elegidas = [m for m in MUTACIONES if filtro.lower() in m[0].lower()]
        if not elegidas:
            print(f"ninguna mutacion matchea '{filtro}'", file=sys.stderr)
            return 1

        for nombre, fich, pares in elegidas:
            shutil.rmtree(pathlib.Path(base) / 'scripts', ignore_errors=True)
            shutil.copytree(RAIZ / 'scripts', pathlib.Path(base) / 'scripts')
            p = pathlib.Path(base) / fich
            texto = p.read_text()
            falta = next((v for v, _ in pares if v not in texto), None)
            if falta is not None:
                print(f"[VIEJA ] {nombre}")
                print(f"          el patron ya no esta en {fich}: {falta[:60]}")
                viejas += 1
                continue
            for viejo, nuevo in pares:
                texto = texto.replace(viejo, nuevo, 1)
            p.write_text(texto)

            # Comparar contra un numero, no buscar '0 fallas' como substring:
            # '20 fallas' lo contiene y la primera version de esto lo reporto
            # como hueco. Mismo error que el resto del proyecto viene cazando.
            fallas, linea = correr(pathlib.Path(base))
            if fallas is None:
                print(f"[ERROR ] {nombre}: {linea}")
                huecos += 1
            elif fallas == 0:
                print(f"[HUECO ] {nombre}  -> ningun banco se quejo")
                huecos += 1
            else:
                print(f"[OK    ] {nombre}  -> {fallas} chequeo(s) en rojo")

        print()
        print(f"{len(elegidas)} mutaciones: {len(elegidas) - huecos - viejas} "
              f"detectadas, {huecos} huecos, {viejas} patrones viejos")
        return 1 if (huecos or viejas) else 0
    finally:
        shutil.rmtree(base, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())
