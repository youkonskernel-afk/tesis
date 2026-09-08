#!/usr/bin/env python3
"""
Validación cruzada dejando un organismo afuera, entre los organismos con
positivos curados (set_modelo == entrenamiento en data/organismos.tsv).

Pregunta que responde: ¿el modelo transfiere entre taxones? Se entrena en dos
organismos y se predice en el tercero. Si no transfiere de pez a molusco, no va
a transferir de pez a musgo, y conviene saberlo antes de correr los nueve.
Ver docs/positivos.md.

Uso:
    ./scripts/loo_cv.py features.tsv                 # todos los métodos
    ./scripts/loo_cv.py features.tsv --metodo bagging_pu
    ./scripts/loo_cv.py --self-test                  # sin datos, verifica la lógica

Entrada: TSV con columnas
    org        código de organismo (gadmo, galga, maggi)
    locus_id   identificador del locus
    y          1 = positivo conocido (coincide con MirGeneDB), 0 = unlabeled
    <resto>    features numéricas

IMPORTANTE — 'y == 0' significa NO ETIQUETADO, no negativo. Toda la evaluación
de acá está construida sobre ese supuesto: no se reporta precisión ni AUC
contra "negativos", porque no hay negativos conocidos.
"""

import argparse
import sys

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold

SEED = 0
METODOS = ("naive", "elkan_noto", "bagging_pu")

# El self-test baja estos valores para correr en segundos en vez de minutos.
N_ARBOLES = 300
N_BOLSAS = 25


# --------------------------------------------------------------------------
# Métodos
# --------------------------------------------------------------------------
def _rf(seed):
    return RandomForestClassifier(
        n_estimators=N_ARBOLES, min_samples_leaf=2, n_jobs=-1, random_state=seed
    )


def naive(X_tr, s_tr, X_te, seed=SEED):
    """Trata los no etiquetados como negativos. Es lo que NO hay que hacer, y
    va incluido justamente porque es la referencia contra la que se compara."""
    clf = _rf(seed).fit(X_tr, s_tr)
    return clf.predict_proba(X_te)[:, 1]


def elkan_noto(X_tr, s_tr, X_te, seed=SEED):
    """Elkan-Noto: g(x)=P(s=1|x), y P(y=1|x)=g(x)/c con c=P(s=1|y=1) estimado
    fuera de muestra sobre los positivos etiquetados.

    OJO: dividir por la constante c es una transformación monótona, así que
    el RANKING es idéntico al de naive. Elkan-Noto sirve para calibrar la
    probabilidad y elegir umbral, no para reordenar. Si se busca cambiar el
    orden hace falta un método que cambie la pérdida (nnPU) o el muestreo
    (bagging_pu)."""
    n_pos = int(s_tr.sum())
    n_folds = min(5, n_pos) if n_pos >= 2 else 0
    if n_folds < 2:
        c = 1.0
    else:
        oof = np.zeros(len(s_tr))
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
        for tr, te in skf.split(X_tr, s_tr):
            oof[te] = _rf(seed).fit(X_tr[tr], s_tr[tr]).predict_proba(X_tr[te])[:, 1]
        c = float(oof[s_tr == 1].mean())
        c = max(c, 1e-6)
    g = _rf(seed).fit(X_tr, s_tr).predict_proba(X_te)[:, 1]
    # SIN recortar a [0,1]. Si c queda subestimado, g/c pasa de 1 y recortarlo
    # satura el tope del ranking en empates — justo donde miramos (top 5%).
    # Para leer el valor como probabilidad sí conviene recortar; para ordenar,
    # nunca.
    return g / c


def bagging_pu(X_tr, s_tr, X_te, seed=SEED, n_bags=None):
    """Bagging PU (Mordelet & Vert): cada bolsa usa todos los positivos y una
    submuestra de no etiquetados del mismo tamaño, tratada como negativa. Al
    promediar, los no etiquetados que se parecen a positivos son penalizados en
    menos bolsas y suben. A diferencia de Elkan-Noto, esto SÍ cambia el orden."""
    n_bags = N_BOLSAS if n_bags is None else n_bags
    rng = np.random.default_rng(seed)
    idx_pos = np.flatnonzero(s_tr == 1)
    idx_unl = np.flatnonzero(s_tr == 0)
    if len(idx_pos) == 0 or len(idx_unl) == 0:
        return np.zeros(len(X_te))
    n_draw = min(len(idx_pos), len(idx_unl))
    acc = np.zeros(len(X_te))
    for b in range(n_bags):
        sub = rng.choice(idx_unl, size=n_draw, replace=False)
        idx = np.concatenate([idx_pos, sub])
        y = np.concatenate([np.ones(len(idx_pos)), np.zeros(len(sub))])
        acc += _rf(seed + b).fit(X_tr[idx], y).predict_proba(X_te)[:, 1]
    return acc / n_bags


# --------------------------------------------------------------------------
# Métricas — todas válidas sin negativos conocidos
# --------------------------------------------------------------------------
def metricas(y_conocidos, scores, k_frac=0.05):
    n = len(scores)
    n_pos = int(y_conocidos.sum())
    if n_pos == 0:
        return None
    k = max(1, int(round(k_frac * n)))
    top = np.zeros(n, dtype=bool)
    top[np.argsort(-scores, kind="stable")[:k]] = True

    recall_k = float(top[y_conocidos == 1].mean())
    # Enriquecimiento: 1.0 = azar. Es la métrica que dice si transfiere.
    enriquecimiento = recall_k / (k / n)
    # Lee & Liu: recall^2 / P(predicho positivo). No necesita negativos.
    lee_liu = recall_k**2 / (k / n)
    # Percentil de los positivos conocidos; 0 = arriba de todo, 0.5 = azar.
    orden = np.argsort(-scores, kind="stable")
    pctl = np.empty(n)
    pctl[orden] = np.arange(n) / max(n - 1, 1)
    pctl_mediano = float(np.median(pctl[y_conocidos == 1]))

    return {
        "n": n,
        "n_pos": n_pos,
        "prior_etiquetado": n_pos / n,
        f"recall@{int(k_frac*100)}%": recall_k,
        "enriquecimiento": enriquecimiento,
        "lee_liu": lee_liu,
        "pctl_mediano_pos": pctl_mediano,
    }


# --------------------------------------------------------------------------
# LOO por organismo
# --------------------------------------------------------------------------
def normaliza_por_organismo(X, orgs):
    """Estandariza dentro de cada organismo. Sin esto, diferencias de
    profundidad de secuenciación entre organismos —galga tiene 32 M spots y
    maggi ~353 M— se cuelan como señal y el modelo aprende el organismo en vez
    del sRNA."""
    Z = np.empty_like(X, dtype=float)
    for o in np.unique(orgs):
        m = orgs == o
        bloque = X[m]
        mu = bloque.mean(axis=0)
        sd = bloque.std(axis=0)
        sd[sd < 1e-12] = 1.0
        Z[m] = (bloque - mu) / sd
    return Z


def loo_cv(orgs, X, y, metodo="bagging_pu", k_frac=0.05, seed=SEED, normalizar=True):
    if normalizar:
        X = normaliza_por_organismo(X, orgs)
    fn = {"naive": naive, "elkan_noto": elkan_noto, "bagging_pu": bagging_pu}[metodo]
    out = {}
    for held in np.unique(orgs):
        tr = orgs != held
        te = orgs == held
        if y[tr].sum() == 0 or y[te].sum() == 0:
            out[held] = None
            continue
        scores = fn(X[tr], y[tr], X[te], seed=seed)
        out[held] = metricas(y[te], scores, k_frac=k_frac)
    return out


def imprime(res, metodo):
    print(f"\n=== método: {metodo} ===")
    cab = f"{'held-out':<10} {'n':>7} {'n_pos':>6} {'recall@5%':>10} {'enriq':>7} {'lee_liu':>8} {'pctl_med':>9}"
    print(cab)
    print("-" * len(cab))
    enr = []
    for org, m in sorted(res.items()):
        if m is None:
            print(f"{org:<10} {'sin positivos en train o test':>50}")
            continue
        enr.append(m["enriquecimiento"])
        print(
            f"{org:<10} {m['n']:>7} {m['n_pos']:>6} "
            f"{m[[k for k in m if k.startswith('recall')][0]]:>10.3f} "
            f"{m['enriquecimiento']:>7.2f} {m['lee_liu']:>8.3f} "
            f"{m['pctl_mediano_pos']:>9.3f}"
        )
    if enr:
        print(f"\nenriquecimiento medio: {np.mean(enr):.2f}   (1.0 = azar)")
        if np.mean(enr) < 1.5:
            print("=> NO transfiere entre estos organismos. Ver el plan B en")
            print("   docs/positivos.md antes de aplicar a plantas y hongos.")
    return enr


# --------------------------------------------------------------------------
# Self-test: verifica que la validación detecta AMBOS casos
# --------------------------------------------------------------------------
def _sintetico(compartido, n_por_org=600, n_pos=60, n_feat=6, seed=1):
    """compartido=True  -> la señal es la misma feature en los 3 organismos.
    compartido=False -> cada organismo usa una feature distinta, así que lo
    aprendido en dos no sirve en el tercero."""
    rng = np.random.default_rng(seed)
    orgs, Xs, ys = [], [], []
    for i, o in enumerate(["gadmo", "galga", "maggi"]):
        X = rng.normal(size=(n_por_org, n_feat))
        f = 0 if compartido else i
        # los positivos se corren en la feature informativa de ese organismo
        idx = rng.choice(n_por_org, size=n_pos, replace=False)
        X[idx, f] += 3.0
        y = np.zeros(n_por_org, dtype=int)
        y[idx] = 1
        orgs.append(np.repeat(o, n_por_org))
        Xs.append(X)
        ys.append(y)
    return np.concatenate(orgs), np.vstack(Xs), np.concatenate(ys)


def self_test():
    global N_ARBOLES, N_BOLSAS
    N_ARBOLES, N_BOLSAS = 60, 8   # el test verifica la lógica, no la precisión
    fallos = []

    for metodo in METODOS:
        # Caso 1: señal compartida -> tiene que transferir.
        orgs, X, y = _sintetico(compartido=True)
        r = loo_cv(orgs, X, y, metodo=metodo)
        enr = [m["enriquecimiento"] for m in r.values() if m]
        ok1 = len(enr) == 3 and min(enr) > 5.0
        print(f"[{'OK ' if ok1 else 'MAL'}] {metodo}: señal compartida -> "
              f"enriquecimiento min={min(enr):.2f} (se espera >5)")
        if not ok1:
            fallos.append(f"{metodo}: no detecta señal compartida")

        # Caso 2: señal específica de cada organismo -> NO tiene que transferir.
        # Este es el test que importa: si la validación no puede distinguir
        # este caso del anterior, no sirve para nada.
        orgs, X, y = _sintetico(compartido=False)
        r = loo_cv(orgs, X, y, metodo=metodo)
        enr2 = [m["enriquecimiento"] for m in r.values() if m]
        ok2 = len(enr2) == 3 and np.mean(enr2) < 2.0
        print(f"[{'OK ' if ok2 else 'MAL'}] {metodo}: señal por organismo -> "
              f"enriquecimiento medio={np.mean(enr2):.2f} (se espera <2)")
        if not ok2:
            fallos.append(f"{metodo}: no detecta fallo de transferencia")

    # Elkan-Noto es monótono respecto de naive: mismo ranking, misma métrica.
    orgs, X, y = _sintetico(compartido=True)
    a = loo_cv(orgs, X, y, metodo="naive")
    b = loo_cv(orgs, X, y, metodo="elkan_noto")
    iguales = all(
        abs(a[o]["pctl_mediano_pos"] - b[o]["pctl_mediano_pos"]) < 1e-9 for o in a
    )
    print(f"[{'OK ' if iguales else 'MAL'}] elkan_noto y naive dan el mismo "
          f"ranking (esperado: c es una constante)")
    if not iguales:
        fallos.append("elkan_noto deberia empatar a naive en ranking")

    # La normalización por organismo tiene que borrar un corrimiento global.
    orgs, X, y = _sintetico(compartido=True)
    X2 = X.copy()
    X2[orgs == "galga"] += 50.0          # simula profundidad muy distinta
    r = loo_cv(orgs, X2, y, metodo="bagging_pu", normalizar=True)
    enr3 = [m["enriquecimiento"] for m in r.values() if m]
    ok4 = min(enr3) > 5.0
    print(f"[{'OK ' if ok4 else 'MAL'}] normalización por organismo absorbe un "
          f"corrimiento de escala -> min={min(enr3):.2f} (se espera >5)")
    if not ok4:
        fallos.append("la normalizacion por organismo no absorbe el corrimiento")

    print()
    if fallos:
        for f in fallos:
            print("FALLO:", f)
        return 1
    print("todo OK")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("features", nargs="?", help="TSV con org, locus_id, y, features")
    ap.add_argument("--metodo", choices=METODOS + ("todos",), default="todos")
    ap.add_argument("--k-frac", type=float, default=0.05)
    ap.add_argument("--sin-normalizar", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        return self_test()
    if not a.features:
        ap.error("falta el TSV de features (o usá --self-test)")

    import pandas as pd
    df = pd.read_csv(a.features, sep="\t")
    for col in ("org", "locus_id", "y"):
        if col not in df.columns:
            sys.exit(f"error: falta la columna '{col}' en {a.features}")
    feats = [c for c in df.columns if c not in ("org", "locus_id", "y")]
    if not feats:
        sys.exit("error: no hay columnas de features")

    orgs = df["org"].to_numpy()
    X = df[feats].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=int)
    print(f"{len(df)} loci, {len(feats)} features, "
          f"{sorted(set(orgs))}, {int(y.sum())} positivos conocidos")

    for m in (METODOS if a.metodo == "todos" else (a.metodo,)):
        imprime(loo_cv(orgs, X, y, metodo=m, k_frac=a.k_frac,
                       normalizar=not a.sin_normalizar), m)
    return 0


if __name__ == "__main__":
    sys.exit(main())
