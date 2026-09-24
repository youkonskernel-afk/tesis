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
     [('ver, ad_top,\n             100*dentro/total > "/dev/stderr"',
       '100*dentro/total > "/dev/stderr"')]),
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
     [("run = q.name.split('.')[0]", "run = q.name.split('.tfq')[0]")]),
    # yasma trim v1.1.1 PISA trimmed_libraries con lo de su propia llamada. Sin
    # el ledger acumulativo, la tanda 2 desmiente a la 1 y todo se re-recorta.
    ("trim: la tanda N vuelve a pisar el registro de la N-1", "scripts/trim.sh",
     [('previas, orden = {}, []\nif led.is_file():',
       'previas, orden = {}, []\nif False:')]),
    # El duplicado es la validacion independiente: comparte directorio de
    # proyecto con el primario y deja de serlo.
    # gadmo PRJNA328800: 6 de 12 corridas con otro kit. Una fila por proyecto
    # no puede expresar eso y el fallo es el silencioso de --trimmed-only.
    ("trim: la fila de la corrida deja de ganarle a la del proyecto",
     "scripts/trim.sh",
     [('$3 == r && r != "" { print $c; hecho = 1; exit }',
       '$3 == r && 0        { print $c; hecho = 1; exit }')]),
    ("trim: sin el flag, exit imprime tambien la fila del proyecto",
     "scripts/trim.sh",
     [('END { if (!hecho && gen != "") print gen }',
       'END { if (gen != "") print gen }')]),
    ("trim: validar vuelve a deduplicar por proyecto y no mira las corridas",
     "scripts/trim.sh",
     [('      [[ " $vistos " == *" $etiqueta|"* ]] && continue\n'
       '      vistos="$vistos $etiqueta|"',
       '      [[ " $vistos " == *" $org/$proy "* ]] && continue\n'
       '      vistos="$vistos $org/$proy"')]),
    ("trim: rehacer no saca la corrida del ledger", "scripts/trim.sh",
     [("    led.write_text('\\n'.join([cab] + quedan) + '\\n')", "    pass")]),
    ("trim: rehacer borra la salida de una PRE-TRIMMED", "scripts/trim.sh",
     [("    if q.is_file() and q.parent.name == 'trim':", "    if q.is_file():")]),
    ("trim: rehacer no limpia inputs.json", "scripts/trim.sh",
     [("    datos['trimmed_libraries'] = [r for r in prev\n"
       "                                  if not any(run in r for run in runs)]",
       "    pass")]),
    ("perfil: --corridas mide una sola, como --proyectos",
     "scripts/fetch_runs.sh",
     [("""    'NR>1 && $1==o && $4==r {print $1"\\t"$3"\\t"$4"\\t"$9"\\t"$2}' "$MANIFEST")""",
       """    'NR>1 && $1==o && $4==r && !seen++ {print $1"\\t"$3"\\t"$4"\\t"$9"\\t"$2}' "$MANIFEST")""")]),
    ("perfil: la columna run sale siempre en '-'", "scripts/fetch_runs.sh",
     [('clave = (pc == "1") ? $5 : "-"', 'clave = "-"')]),
    ("trim: el primario y el duplicado vuelven a un solo directorio",
     "scripts/trim.sh",
     [('dir="$PROY_DIR/${org}_${rol}"', 'dir="$PROY_DIR/${org}"')] * 4),
    ("align: el primario y el duplicado comparten BAM", "scripts/align.sh",
     [('dir="$PROY_DIR/${org}_${rol}"', 'dir="$PROY_DIR/${org}"')] * 5),
    # maggi_primario cruza dos BioProjects con 62% y 87% de retencion esperada.
    ("trim: el ledger pone el bioproject del primero de la tanda",
     "scripts/trim.sh",
     [("de_proyecto.get(run, '-')", "list(de_proyecto.values() or ['-'])[0]")] * 2),
    # Recortar con la secuencia equivocada no da error: --trimmed-only descarta
    # todo lo que no matchea y deja un .t.fq.gz casi vacio.
    ("trim: verificar deja pasar una retencion vacia", "scripts/trim.sh",
     [('if (m < 5)       print "VACIA', 'if (0)           print "VACIA')]),
    # Para una PRE-TRIMMED el fastq sin recortar ES la salida que YASMA anoto.
    ("trim: borra el fastq de una PRE-TRIMMED", "scripts/trim.sh",
     [('        [[ "$sec" == "PRE-TRIMMED" ]] && continue\n        rm -f',
       '        rm -f')]),
    # --- align.sh ---
    # El ledger de sha256 es el unico registro de contra que se alineo.
    ("align: alinea sin comprobar el sha256 del genoma", "scripts/align.sh",
     [('[[ "$real" == "$esp" ]] || die "el sha256', '[[ 1 ]] || die "el sha256')]),
    # ic.check() hace relative_to(output_directory) sin protegerlo: un genoma de
    # afuera tira ValueError. Y el symlink tiene que ser AL DIRECTORIO, o
    # bowtie-build corre una vez por proyecto en vez de una por organismo.
    ("align: el genoma vuelve a ir de afuera del -o", "scripts/align.sh",
     [('-g "$dir/genome/$(basename "$fna")"', '-g "$fna"')]),
    ("align: el symlink del genoma apunta al fichero, no al directorio",
     "scripts/align.sh",
     [('ln -sfn "$GENOMES_DIR/$org" "$dir/genome"',
       'mkdir -p "$dir/genome" && ln -sfn "$fna" "$dir/genome/$(basename "$fna")"')]),
    # Alinear contra el genoma equivocado no falla: yasma sale con 0 y 0% alineado.
    ("align: verificar deja pasar una fraccion alineada muy baja",
     "scripts/align.sh",
     [('if (al < 10)       v = "MUY BAJA', 'if (0)             v = "MUY BAJA')]),
    # 17.2% colocado con 67.1% sin alineamiento valido pasaba como "ok": el
    # unico umbral miraba lo colocado, y 17.2 > 10. Medido en sclsc_duplicado.
    ("align: verificar no mira lo que no alinea en ninguna parte",
     "scripts/align.sh",
     [('else if (sa > 50)', 'else if (0)      ')]),
    # Y que los dos diagnosticos no se confundan: el de -m es un genoma
    # repetitivo (los tRF de danre), el de SIN_AL son reads de otro genoma.
    ("align: SIN_AL y >m se calculan del mismo conteo", "scripts/align.sh",
     [('sa = 100*n/tot', 'sa = 100*h/tot')]),
    ("align: no reporta una corrida que no llego al BAM", "scripts/align.sh",
     [('    while read -r run; do', '    while false; do')]),
    ("align: alinea un proyecto a medio recortar", "scripts/align.sh",
     [('    if [[ "$n_rec" -lt "$n_man" ]]; then', '    if false; then')]),
    ("align: el -m 50 deja de ir explicito", "scripts/align.sh",
     [('--max_multi "$MAX_MULTI" ', '')]),
    ("align: el BAM no llega a donde drive_push lo busca", "scripts/align.sh",
     [('    enlazar_bam "$org" "$rol" "$dir"', '    :')]),
    ("align: no queda registrado contra que genoma se alineo", "scripts/align.sh",
     [('    registrar "$dir" "$org" "$rol" "$acc" "$esp"', '    :')]),
    # El parche de yasma: sin el, la etapa `over` levanta un bowtie por
    # libreria que nadie lee ni espera y el proyecto muere por OOM a las horas.
    ("align: alinea aunque falte el parche de yasma", "scripts/align.sh",
     [('"$ROOT/scripts/yasma_parche.py" --verificar >/dev/null 2>&1 \\\n    ||',
       'true \\\n    ||')]),
    ("align: la RAM no se dice antes de alinear", "scripts/align.sh",
     [('    n_bases=$(bases_de "$fna"); ram_pide=$(ram_gb_de "$n_bases"); ram_hay=$(ram_libre_gb)',
       '    n_bases=0; ram_pide=0; ram_hay=""; if false; then')]),
    ("align: plan pierde la columna de RAM", "scripts/align.sh",
     [("%10d %10d %7s  %s\\n' \"${org}_${rol}\"", "%10d %10d %.0s  %s\\n' \"${org}_${rol}\"")]),
    # --- yasma_parche.py ---
    # Un parche que aplica a ciegas sobre una version que cambio es peor que no
    # tenerlo: no falla ruidosamente, deja el fuente en cualquier estado.
    ("parche: aplica aunque el ancla no este una sola vez", "scripts/yasma_parche.py",
     [('    if n != 1:', '    if False:')]),
    ("parche: no comprueba que en 'over' no se espere al proceso", "scripts/yasma_parche.py",
     [('    if GUARDIA not in texto:', '    if False:')]),
    ("parche: saltea el Popen en todas las etapas", "scripts/yasma_parche.py",
     [("\"\\t\\tif mmap == 'over':\\n\"", "\"\\t\\tif True:\\n\"")]),
    ("parche: --verificar sale 0 sin estar aplicado", "scripts/yasma_parche.py",
     [('    if args.verificar:', '    if False:')]),
    # --- colab_git.py: empujar a GitHub desde Colab ---
    # Una VM de Colab tiene Drive montado al lado con ~190 GB de .sra. Un
    # `add -A` desde ahi es exactamente donde se cuela lo que no va al repo.
    ("colab_git: vuelve a `git add -A`", "scripts/colab_git.py",
     [('"add", "--", *limpias', '"add", "-A"')]),
    ("colab_git: deja pasar rutas fuera de data/", "scripts/colab_git.py",
     [('elif p.parts[:1] != ("data",):', 'elif False:')]),
    # git mete la URL —con el token adentro— en sus mensajes de error.
    ("colab_git: el token se filtra", "scripts/colab_git.py",
     [('return txt.replace(secreto, "***") if secreto else txt', 'return txt')]),
    # Un push rechazado que nadie mira deja el resultado en Drive y no en git,
    # que es el estado que este modulo existe para evitar.
    ("colab_git: un push rechazado no revienta", "scripts/colab_git.py",
     [('    raise RuntimeError(\n        "el push fue rechazado dos veces.',
       '    return (\n        "el push fue rechazado dos veces.')]),
    ("colab_git: revisar=True empuja igual", "scripts/colab_git.py",
     [('    if revisar:', '    if False:')]),
    # Con la URL cableada, cualquier clon empujaria al repo de verdad apenas
    # hubiera un GITHUB_TOKEN en el entorno.
    ("colab_git: empuja a github aunque el origin sea otro", "scripts/colab_git.py",
     [('    if tok and "github.com" in url:', '    if tok:')]),
    # La celda de clon esta copiada en todos los notebooks: si deriva, alguien
    # arregla el bug en uno y los otros siguen rotos.
    ("notebooks: la celda de clon puede derivar", "scripts/validate_notebooks.py",
     [('    if len(vistas) > 1:', '    if False:')]),
    # §1 de 20_alinear: un proyecto que no entra en disco tiene que saberse en
    # un segundo, no a las seis horas.
    ("alinear: un proyecto que no entra se da por bueno",
     "notebooks/20_alinear.ipynb",
     [('    entra = pico + MARGEN < libre', '    entra = True')]),
    ("alinear: mide los reads crudos y no lo que sobrevive al recorte",
     "notebooks/20_alinear.ipynb",
     [("int(r['read_count']) * _ret.get((r['org'], r['bioproject']), 80) / 100",
       "int(r['read_count'])")]),
    # B_BAM estuvo en 45 sin haberse medido nunca: el pico de BAM salia 3x. Un
    # numero inventado en una constante no falla ruidosamente, solo manda
    # proyectos que entran a la maquina local.
    # Se perdieron 31 min re-alineando un proyecto que ya estaba en Drive: la
    # VM es efimera, asi que "que falta" no se puede leer del disco local.
    ("alinear: no mira si el proyecto ya esta en Drive",
     "notebooks/20_alinear.ipynb",
     [("hecho = (DRIVE / '10_bam' / org / f'{rol}.bam').exists()",
       "hecho = False")]),
    ("alinear: SIGUIENTE puede caer en uno que no entra",
     "notebooks/20_alinear.ipynb",
     [('    elif entra:', '    elif True:')]),
    # §3b/§3c: el paso siguiente a que el recorte corte. Transcribir seis
    # accessions a mano desde una tabla formateada es donde se cuela el error.
    ("celda3b: MALAS deja de salir de la tabla", "notebooks/20_alinear.ipynb",
     [("        MALAS.append(_c[2])", "        pass")]),
    ("celda3b: una DESVIADA deja de contar", "notebooks/20_alinear.ipynb",
     [("if len(_c) > 3 and ('VACIA' in _ln or 'DESVIADA' in _ln):",
       "if len(_c) > 3 and ('VACIA' in _ln):")]),
    ("celda3b: lee mal la columna y agarra el bioproject",
     "notebooks/20_alinear.ipynb",
     [("        MALAS.append(_c[2])", "        MALAS.append(_c[1])")]),
    ("celda3b: sigue aunque no haya podido leer la tabla",
     "notebooks/20_alinear.ipynb",
     [("    assert MALAS, ('verificar falló pero no pude leer qué corridas",
       "    assert True, ('verificar falló pero no pude leer qué corridas")]),
    ("celda3c: rehace sin haber corregido la tabla", "notebooks/20_alinear.ipynb",
     [("if FILAS:", "if True:")]),
    ("alinear: B_BAM vuelve a la estimacion sin medir",
     "notebooks/20_alinear.ipynb",
     [('B_BAM    = 16', 'B_BAM    = 45')]),
    # §1 del notebook. Estuvo mal tres veces, asi que se muta igual que el
    # codigo de scripts/.
    ("celda1: vuelve a reusar la copia siempre",
     "notebooks/10_descarga_runs.ipynb",
     [('if any(visto):', 'if False:')]),
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
    ("perfil --tsv: emite el prefijo en vez de la secuencia completa",
     "scripts/fetch_runs.sh",
     [('else if (fam in SEC)       sec = SEC[fam]', 'else                       sec = fam')]),
    ("perfil --tsv: una ya recortada deja de ir PRE-TRIMMED",
     "scripts/fetch_runs.sh",
     [('sec = "PRE-TRIMMED"; ret = "100"', 'ret = "100"')]),
    ("perfil: reporta familia sobre un punado de reads", "scripts/fetch_runs.sh",
     [('hay_ad = (con >= 0.2 * total)', 'hay_ad = (con > 0)')]),
    ("perfil: la retencion vuelve a ser el adapt_pct", "scripts/fetch_runs.sh",
     [('100*dentro/total > "/dev/stderr"', '100*con/total > "/dev/stderr"')]),
    ("perfil: la PRE-TRIMMED retiene 0 en vez de 100", "scripts/fetch_runs.sh",
     [('sec = "PRE-TRIMMED"; ret = "100"', 'sec = "PRE-TRIMMED"')]),
    ("perfil: deja de avisar que el inserto modal esta fuera", "scripts/fetch_runs.sh",
     [('if (modal + 0 < vmin || modal + 0 > vmax) {', 'if (0) {')]),
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
     # Es el caso de SRR23277331: excluida del manifiesto, pero su md5 sigue en
     # el ledger de Drive. Empujarla desharia la exclusion sin que nadie lo vea.
     # OJO: en un .ipynb el codigo va JSON-escapado, asi que un \t del codigo
     # Python aparece como \\t en el fichero. mutar.py muta el texto crudo.
     [(r"led_bueno = [l for l in lineas[1:] if l.split('\\t')[1] in en_manifiesto]",
       r"led_bueno = list(lineas[1:])")]),
    ("celda4: el guardia de excluidas deja de cortar",
     "notebooks/10_descarga_runs.ipynb",
     [('_mal = [r for r in _excl if r in en_manifiesto]', '_mal = []')]),
    # Las rutas locales: tres scripts tenian tres respuestas distintas.
    ("rutas: fetch_runs vuelve a su ruta propia", "scripts/fetch_runs.sh",
     [('CACHE="${SRA_CACHE:-$(ruta_local sra)}"', 'CACHE="/home/dev/sra_cache"')]),
    ("rutas: trim.sh vuelve a su ruta propia", "scripts/trim.sh",
     [('SRA_DEST="${SRA_DEST:-$(ruta_local sra)}"',
       'SRA_DEST="$HOME/tesis_data/80_sra"')]),
    ("rutas: ruta_local ignora LOCAL_ROOT", "scripts/_drive_lib.sh",
     [('raiz="${LOCAL_ROOT:-$(cd "${BASH_SOURCE[0]%/*}/.." && pwd)}"',
       'raiz="$(cd "${BASH_SOURCE[0]%/*}/.." && pwd)"')]),
    # purge borra la copia local, que para lo producido aca puede ser la unica.
    ("purge: borra sin confirmar contra Drive", "scripts/drive_pull.sh",
     [('    sra|genomas) ;;   # Drive es la fuente; nada que confirmar',
       '    sra|genomas|bam|yasma|qc|features|modelos|figuras) ;;')]),
    ("purge: deja de exigir un ORG", "scripts/drive_pull.sh",
     [('[[ -n "$ORG" ]] || { echo "purge necesita un ORG explícito" >&2; exit 2; }',
       ':')]),
    ("drive_check: deja de mirar el scope", "scripts/drive_check.sh",
     [('  drive.file)', '  drive.file_NUNCA)')]),
    ("drive_check: deja de mirar root_folder_id", "scripts/drive_check.sh",
     [('elif [[ -n "$DRIVE_ROOT" ]]; then', 'elif false; then')]),
    ("drive_check: una raiz vacia pasa desapercibida", "scripts/drive_check.sh",
     [('if [[ "$n" -eq 0 ]]; then', 'if false; then')]),
    # `Error 401: invalid_client` muere en el navegador y no vuelve a rclone,
    # que igual ofrece guardar el remoto. Queda uno sin token.
    ("drive_check: un client_id malformado pasa", "scripts/drive_check.sh",
     [("      'propio' if c.endswith('.apps.googleusercontent.com') else 'malo')",
       "      'propio')")]),
    ("drive_check: un remoto sin autorizar pasa", "scripts/drive_check.sh",
     [('if [[ "${tok_len:-0}" -eq 0 ]]; then', 'if false; then')]),
    # DRIVE_ROOT='' es una config valida; `:-` la pisaria con el default.
    ("DRIVE_ROOT vacio se pisa con el default", "scripts/drive_pull.sh",
     [('DRIVE_ROOT="${DRIVE_ROOT-tesis}"', 'DRIVE_ROOT="${DRIVE_ROOT:-tesis}"')]),
    ("push: vuelve a --size-only", "scripts/drive_push.sh",
     [('--checksum --progress --transfers 4 --drive-chunk-size 64M',
       '--size-only --progress --transfers 4 --drive-chunk-size 64M')]),
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
