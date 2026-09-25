#!/usr/bin/env python3
"""El reparto entre máquinas y los claims que evitan que dos hagan lo mismo.

Lo que este banco protege, por orden de gravedad:

 1. Que un proyecto que NO entra en ninguna máquina salga como tal y no se
    reparta igual. Son los dos `galga` —el 33% de los reads del set— y su modo
    de fallo es un `Killed` a las horas, no un error.
 2. Que el resumen NO cuente en el speedup lo que no pudo asignar. La primera
    version decia 4.5x mientras las 34 h de `galga` no estaban en ninguna
    maquina: el numero que la persona lee para decidir cuantas conseguir.
 3. Que el reparto sea LPT y no cualquier cosa. Del mas chico al mas grande,
    todas terminan los baratos y una queda sola con `galga_duplicado` 25 h al
    final — el makespan se va al doble sin que nada falle.
 4. Que un genoma que no se pudo medir no se declare imposible.
 5. Que dos maquinas no tomen el mismo proyecto, y que una que se murio no lo
    deje bloqueado para siempre.
"""
import datetime
import importlib.util
import pathlib
import subprocess
import sys
import tempfile

RAIZ = pathlib.Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("reparto", RAIZ / "scripts" / "reparto.py")
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)

FALLAS = 0
GB = 2**30


def chk(n, ok, extra=""):
    global FALLAS
    print(("  ok   " if ok else "  MAL  ") + n + (f"  [{extra}]" if not ok and extra else ""))
    if not ok:
        FALLAS += 1


CHICA = ('chica', 11.4 * GB, 220e9)
GRANDE = ('grande', 64 * GB, 1600e9)


def trabajo(nom, h, ram_gb=0.5, pico_gb=5):
    return (nom, h, ram_gb * 1e9, pico_gb * 1e9)


print("== 1. LPT: el reparto queda cerca de la cota, no en cualquier lado")
# La cota inferior de cualquier reparto es max(trabajo mas largo, total/m).
# LPT garantiza quedar cerca; empezar por los chicos no.
JOBS = [trabajo(f'p{i}', h) for i, h in
        enumerate([25.2, 11.0, 10.0, 8.8, 8.4, 7.4, 7.2, 4.6, 3.8, 3.0,
                   2.8, 2.4, 1.8, 1.4, 1.3, 0.9, 0.6, 0.5])]
tot = sum(j[1] for j in JOBS)
mayor = max(j[1] for j in JOBS)
for m in (2, 3, 4):
    maqs = [(f'w{i}', 64 * GB, 1600e9) for i in range(m)]
    asign, sin = R.repartir(JOBS, maqs)
    ms = R.makespan(asign)
    cota = max(mayor, tot / m)
    chk(f"con {m} máquinas queda a menos de 5% de la cota ({ms:.1f} vs {cota:.1f})",
        ms <= cota * 1.05, f"{ms:.2f} > {cota*1.05:.2f}")
    chk(f"  y reparte los {len(JOBS)} sin perder ninguno",
        sum(len(v) for v in asign.values()) == len(JOBS) and not sin)

print("== 2. el techo lo marca el trabajo más largo, y no lo mueve agregar máquinas")
# galga_duplicado son 25.2 h y no se puede partir: nativealign acumula unique_d
# sobre TODAS las librerias antes de pesar los multimapeados.
ms5 = R.makespan(R.repartir(JOBS, [(f'w{i}', 64 * GB, 1600e9) for i in range(5)])[0])
ms9 = R.makespan(R.repartir(JOBS, [(f'w{i}', 64 * GB, 1600e9) for i in range(9)])[0])
chk("con 5 ya está en el piso", abs(ms5 - mayor) < 0.01, ms5)
chk("y con 9 no baja", abs(ms9 - ms5) < 0.01, (ms5, ms9))

print("== 3. un proyecto que no entra en NINGUNA máquina no se reparte")
gordo = trabajo('galga/duplicado', 25.2, ram_gb=10.53, pico_gb=81)
asign, sin = R.repartir([gordo, trabajo('chico', 1.0)], [CHICA, CHICA[:1] + CHICA[1:]])
chk("queda en sin_lugar", [n for n, _ in sin] == ['galga/duplicado'], sin)
chk("y dice que es la RAM", sin and sin[0][1] == 'RAM', sin)
chk("no lo asignó a nadie",
    all('galga/duplicado' not in [n for n, _ in v] for v in asign.values()), asign)
chk("el chico sí entra", sum(len(v) for v in asign.values()) == 1, asign)

print("== 3b. y sí se reparte si hay una máquina que puede")
asign2, sin2 = R.repartir([gordo, trabajo('chico', 1.0)], [CHICA, GRANDE])
chk("va a la grande", ('galga/duplicado', 25.2) in asign2['grande'], asign2)
chk("y nada queda sin lugar", sin2 == [], sin2)

print("== 4. el disco se distingue de la RAM")
ancho = trabajo('ancho', 5.0, ram_gb=0.5, pico_gb=500)
_, sin3 = R.repartir([ancho], [CHICA])
chk("le falta disco, no RAM", sin3 and sin3[0][1] == 'disco', sin3)

print("== 5. un genoma que no se pudo medir NO se declara imposible")
# Mandar a la maquina grande algo que entra en cualquiera cuesta tanto como lo
# contrario. Sin genoma medible (ram_b == 0) pasa igual.
asign4, sin4 = R.repartir([trabajo('sin_genoma', 3.0, ram_gb=0)], [CHICA])
chk("entra igual", sum(len(v) for v in asign4.values()) == 1, asign4)
chk("y no queda sin lugar", sin4 == [], sin4)

print("== 6. dentro de cada máquina el orden es del más chico al más grande")
# No cambia el makespan y deja algo terminado si la sesion se muere.
asign5, _ = R.repartir(JOBS, [('w', 64 * GB, 1600e9)])
hs = [h for _, h in asign5['w']]
chk("ordenado ascendente", hs == sorted(hs), hs)

print("== 7. los claims: dos máquinas no toman el mismo proyecto")
with tempfile.TemporaryDirectory() as d:
    chk("la primera lo toma", R.tomar(d, 'galga/duplicado', 'colab1'))
    chk("la segunda NO", not R.tomar(d, 'galga/duplicado', 'colab2'))
    chk("y sigue siendo de la primera",
        R.tomados(d).get('galga/duplicado', ('?',))[0] == 'colab1', R.tomados(d))
    chk("renovar el propio es idempotente", R.tomar(d, 'galga/duplicado', 'colab1'))
    chk("otro proyecto sí se puede", R.tomar(d, 'cloro/primario', 'colab2'))
    chk("y quedan los dos", len(R.tomados(d)) == 2, R.tomados(d))

    print("== 8. soltar lo libera, y sólo el dueño puede")
    chk("un ajeno no lo suelta", not R.soltar(d, 'galga/duplicado', 'colab2'))
    chk("sigue tomado", 'galga/duplicado' in R.tomados(d))
    chk("el dueño sí", R.soltar(d, 'galga/duplicado', 'colab1'))
    chk("y queda libre", 'galga/duplicado' not in R.tomados(d), R.tomados(d))

print("== 9. un claim vencido se puede robar: una VM que se murió no bloquea nada")
with tempfile.TemporaryDirectory() as d:
    ahora = datetime.datetime(2026, 9, 25, 12, 0, tzinfo=datetime.timezone.utc)
    viejo = ahora - datetime.timedelta(hours=20)
    R.tomar(d, 'galga/duplicado', 'muerta', ahora=viejo)
    chk("a las 20 h ya no cuenta como vigente",
        R.tomados(d, ttl_h=12, ahora=ahora) == {}, R.tomados(d, ttl_h=12, ahora=ahora))
    chk("y otra lo puede tomar",
        R.tomar(d, 'galga/duplicado', 'viva', ttl_h=12, ahora=ahora))
    chk("dentro del TTL no", not R.tomar(d, 'galga/duplicado', 'tercera',
                                         ttl_h=12, ahora=ahora))

print("== 10. un claim ilegible no bloquea ni miente")
with tempfile.TemporaryDirectory() as d:
    p = pathlib.Path(d); p.mkdir(exist_ok=True)
    (p / 'galga_duplicado.claim').write_text('basura sin fecha\n')
    chk("no aparece como vigente", R.tomados(d) == {}, R.tomados(d))
    chk("y se puede tomar", R.tomar(d, 'galga/duplicado', 'colab1'))

print("== 11. el ajuste de s/M: un punto no alcanza")
chk("un solo genoma -> None", R.ajuste_s_por_m([(39, 58)]) is None)
chk("dos iguales tampoco", R.ajuste_s_por_m([(39, 58), (39, 60)]) is None)
aj = R.ajuste_s_por_m([(39, 58), (639, 108)])
chk("dos distintos sí", aj is not None, aj)
chk("pasa por los dos puntos",
    round(R.s_por_m(39, aj)) == 58 and round(R.s_por_m(639, aj)) == 108,
    (R.s_por_m(39, aj), R.s_por_m(639, aj)))
chk("y nunca devuelve <= 0", R.s_por_m(0, (-100, 0.001)) > 0)

print("== 12. el resumen NO cuenta en el speedup lo que no repartió")
# Decia 4.5x mientras las 34 h de galga no estaban en ninguna maquina. Es el
# numero que la persona lee para decidir cuantas maquinas conseguir.
r = subprocess.run([sys.executable, str(RAIZ / 'scripts' / 'reparto.py'),
                    '--maquinas', '3'], capture_output=True, text=True)
out = r.stdout + r.stderr
chk("sale 0", r.returncode == 0, out[-400:])
chk("nombra los que no entran", 'galga/duplicado' in out and 'le falta RAM' in out, out[-600:])
chk("dice cuántas horas deja afuera", 'NO cubre' in out, out[-600:])
chk("y da el de punta a punta", 'punta a punta' in out, out[-600:])
_sp = [l for l in out.splitlines() if 'makespan de lo repartido' in l]
chk("el speedup es sobre lo repartido", _sp and '3.0x' in _sp[0], _sp)
chk("no infla a 4.5x", not any('4.5x' in l for l in _sp), _sp)

print()
if FALLAS:
    print(f"{FALLAS} fallas")
    raise SystemExit(1)
print("TODO OK")
