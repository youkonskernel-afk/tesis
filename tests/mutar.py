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
     [('A[i] = substr(campo[i], 1, j - 1); NOM[i] = substr(campo[i], j + 1)',
       'A[i] = campo[i]; NOM[i] = campo[i]')]),
    ("perfil: el 5p deja de avisar que no sirve", "scripts/fetch_runs.sh",
     [('if (ad_top ~ /^5p:/) {', 'if (0) {')]),
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
    ("trim: recorta sin fila en la tabla", "scripts/trim.sh",
     [('if [[ -z "$sec" ]]; then', 'if false; then')]),
    ("trim: acepta un adaptador 5p", "scripts/trim.sh",
     [("5p:*|'??:'*) no_recortable+=", "xx:*) no_recortable+=")]),
    ("trim: pierde las rutas relativas de inputs.json", "scripts/trim.sh",
     [('libs.append(f"untrimmed/{nombre}")', 'libs.append(str(dir_proy / "untrimmed" / nombre))')]),
    ("trim: no le pasa -t a fasterq-dump", "scripts/trim.sh",
     [('-t "$TMP_FASTERQ" ', '')]),
    ("trim: la ventana deja de ir explicita", "scripts/trim.sh",
     [('--min_length 15 --max_length 50 </dev/null', '</dev/null')]),
    ("trim: vuelve a adivinar el nombre de salida", "scripts/trim.sh",
     [("print(q.name.split('.')[0])", "print(q.name.split('.tfq')[0])")]),
    # §1 del notebook. Estuvo mal tres veces, asi que se muta igual que el
    # codigo de scripts/.
    ("celda1: vuelve a reusar la copia siempre",
     "notebooks/10_descarga_runs.ipynb",
     [('if any(revisar()):', 'if False:')]),
    ("celda1: no chequea las exclusiones",
     "notebooks/10_descarga_runs.ipynb",
     [('mal = [r for r in _excluidas() if r in corridas]', 'mal = []')]),
    ("celda1: no chequea los proyectos de la spec",
     "notebooks/10_descarga_runs.ipynb",
     [('faltan = sorted(_proyectos_spec() - proy_man)', 'faltan = []')]),
    ("celda1: pisa la copia de Drive aunque la ENA falle",
     "notebooks/10_descarga_runs.ipynb",
     [("raise RuntimeError('fetch_runs.sh manifest fallo; no piso la copia de Drive')",
       'pass')]),
    ("celda4: lee el working tree en vez de git",
     "notebooks/10_descarga_runs.ipynb",
     [("'git', '-C', str(CLON), 'show', f'HEAD:{ruta_rel}'",
       "'cat', str(CLON / ruta_rel)")]),
    ("celda4: deja de pedir las filas del manifiesto",
     "notebooks/10_descarga_runs.ipynb",
     [('man_faltan = [l for l in man_uso[1:] if l not in man_git[1:]]',
       'man_faltan = []')]),
    ("celda4: ofrece filas del ledger que no estan en el manifiesto",
     "notebooks/10_descarga_runs.ipynb",
     # OJO: en un .ipynb el codigo va JSON-escapado, asi que un \\t del codigo
     # Python aparece como \\\\t en el fichero. mutar.py muta el texto crudo.
     [(r"if l not in led_git[1:] and l.split('\\t')[1] in en_manifiesto]",
       "if l not in led_git[1:]]")]),
    ("celda4: el guardia de excluidas deja de cortar",
     "notebooks/10_descarga_runs.ipynb",
     [('_mal = [r for r in _excl if r in en_manifiesto]', '_mal = []')]),
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
        shutil.copytree(RAIZ / 'notebooks', pathlib.Path(base) / 'notebooks')
        huecos = viejas = 0
        elegidas = [m for m in MUTACIONES if filtro.lower() in m[0].lower()]
        if not elegidas:
            print(f"ninguna mutacion matchea '{filtro}'", file=sys.stderr)
            return 1

        for nombre, fich, pares in elegidas:
            # Se restauran los dos arboles mutables antes de cada mutacion:
            # las de celda1 tocan notebooks/, no scripts/.
            for sub in ('scripts', 'notebooks'):
                shutil.rmtree(pathlib.Path(base) / sub, ignore_errors=True)
                shutil.copytree(RAIZ / sub, pathlib.Path(base) / sub)
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
