"""Lo que `yasma tradeoff` guarda por base del genoma, guardado solo donde hay reads.

Este modulo NO se importa desde el repo: `scripts/tradeoff_parche.py` lo copia
adentro del paquete yasma como `_tesis_disperso.py` y hace que `tradeoff.py` lo
use. Ver el docstring de ese script para el por que.

EL PROBLEMA, EN EL FUENTE DE v1.1.1 (tradeoff.py, get_kernel_coverage)

    pos_size_d[chrom] = [None] * chrom_length                    # 8 B por base
    pos_d[chrom] = np.zeros((chrom_length, C, N, 2), 'uint32')   # 8*C*N B por base

con C condiciones y N librerias. Para gadmo_duplicado (670 Mb, 12 librerias)
son 64 GB de `pos_d` y 5.4 de `pos_size_d`; para galga_duplicado (1.05 Gb, 95
librerias), 800 GB. Encima, `cov_d` guarda la cobertura de todo el genoma
(4 B/base) sin usarla nunca. La VM de Colab tiene ~12 GB: `-9` en los cuatro
proyectos que se probaron.

LO QUE HACE ESTE MODULO

Guarda lo mismo de forma dispersa y lo DEVUELVE IGUAL a lo que leeria el
original:
- `Conteos[a:b, ...]` da el mismo arreglo denso uint32 (n, C, N, 2) que daba
  `pos_d[chrom][a:b, ...]`, con la misma semantica de rebanada (negativos,
  fuera de rango, vacia). Todo el codigo de YASMA que suma rebanadas sigue igual.
- `Largos[p]` da la misma lista de largos que `pos_size_d[chrom][p]`, o None si
  no hay reads, y con el MISMO ORDEN de primera aparicion. Importa:
  `sizeClass.get()` desempata con `Counter.most_common()`, que conserva el orden
  de insercion, asi que reordenar los largos podria cambiar un sizecall empatado.
- `cobertura()` y `por_posicion()` hacen las mismas operaciones numpy que el
  original (suma de hebras en uint32, por rpb en uint32, mediana entre
  librerias, float32, media entre condiciones), fila por fila, pero solo sobre
  las posiciones con reads. Las otras dan 0 en el original y 0 aca.
"""
import array
import os

import numpy as np

# Lecturas pendientes antes de compactar. 4 M * 8 B = 32 MB por buffer.
# Las variables de entorno existen para el banco: con pocos reads, forzar
# compactaciones y bloques chicos es la unica forma de ejercitar la fusion.
COMPACTAR = int(os.environ.get("TESIS_DISPERSO_COMPACTAR", 1 << 22))
# Memoria maxima del bloque denso que se arma para la mediana.
BLOQUE_B = int(os.environ.get("TESIS_DISPERSO_BLOQUE_B", 64 << 20))
# Largo maximo de read que entra en la clave de Largos.
_BITS_LARGO = 16

_activo = {}


def _cambiar_activo(obj):
    """Los reads llegan cromosoma por cromosoma: al pasar al siguiente, el buffer
    del anterior se compacta. Si no, cada cromosoma retendria hasta COMPACTAR
    lecturas sin comprimir y con cientos de scaffolds eso vuelve a ser mucho."""
    clase = type(obj)
    prev = _activo.get(clase)
    if prev is not None and prev is not obj:
        prev._compactar()
    _activo[clase] = obj


def _indice_lista(p, n):
    """La semantica de indexar una lista de largo n."""
    if p >= n or p < -n:
        raise IndexError("list index out of range")
    return p + n if p < 0 else p


def _fusionar(k, cols, reducir):
    """Ordena por clave y junta las repetidas. -> (claves unicas, columnas)."""
    orden = np.argsort(k, kind="stable")
    k = k[orden]
    cols = [c[orden] for c in cols]
    if k.size == 0:
        return k, cols
    idx = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    return k[idx], [r.reduceat(c, idx) for r, c in zip(reducir, cols)]


class Conteos:
    """En lugar de np.zeros((L, C, N, 2), uint32): (clave, cuenta) con
    clave = ((pos*C + cond)*N + lib)*2 + hebra, ordenadas."""

    def __init__(self, largo, n_cond, n_lib):
        self.L, self.C, self.N = int(largo), int(n_cond), int(n_lib)
        self.paso = self.C * self.N * 2
        self.shape = (self.L, self.C, self.N, 2)
        self._buf = array.array("q")
        self._k = np.empty(0, np.int64)
        self._n = np.empty(0, np.int64)

    def sumar(self, pos, cond, lib, hebra):
        """pos_d[chrom][pos, cond, lib, hebra] += 1"""
        if _activo.get(Conteos) is not self:
            _cambiar_activo(self)
        pos = _indice_lista(pos, self.L)
        self._buf.append(((pos * self.C + cond) * self.N + lib) * 2 + hebra)
        if len(self._buf) >= COMPACTAR:
            self._compactar()

    def _compactar(self):
        if not len(self._buf):
            return
        b = np.frombuffer(self._buf, dtype=np.int64)
        k, n = np.unique(b, return_counts=True)
        self._buf = array.array("q")
        if self._k.size:
            k = np.concatenate((self._k, k))
            n = np.concatenate((self._n, n.astype(np.int64)))
            k, (n,) = _fusionar(k, [n], [np.add])
        self._k, self._n = k, n.astype(np.int64)

    def _datos(self):
        self._compactar()
        return self._k, self._n

    def __getitem__(self, clave):
        if isinstance(clave, tuple):
            if not clave or any(c is not Ellipsis for c in clave[1:]):
                raise TypeError(f"indice no soportado por el parche: {clave!r}")
            clave = clave[0]
        if not isinstance(clave, slice):
            raise TypeError(f"indice no soportado por el parche: {clave!r}")
        r = range(self.L)[clave]
        if r.step != 1:
            raise TypeError(f"rebanada con paso no soportada por el parche: {clave!r}")
        n = len(r)
        out = np.zeros((n, self.C, self.N, 2), dtype=np.uint32)
        if n == 0:
            return out
        k, c = self._datos()
        lo, hi = r.start * self.paso, (r.start + n) * self.paso
        i0, i1 = np.searchsorted(k, lo, "left"), np.searchsorted(k, hi, "left")
        out.reshape(-1)[k[i0:i1] - lo] = c[i0:i1].astype(np.uint32)
        return out

    def por_posicion(self):
        """np.sum(pos_d[chrom], axis=(1,2,3), dtype='uint32')"""
        out = np.zeros(self.L, dtype=np.uint32)
        k, c = self._datos()
        if k.size:
            pos = k // self.paso
            idx = np.flatnonzero(np.r_[True, pos[1:] != pos[:-1]])
            out[pos[idx]] = np.add.reduceat(c, idx).astype(np.uint32)
        return out

    def cobertura(self, rpbs):
        """Lo que el original calcula sobre todo el cromosoma con cinco lineas
        (suma de hebras en uint32, por rpb, mediana entre librerias, float32,
        media entre condiciones). Son las mismas cinco lineas, copiadas tal
        cual abajo, pero por bloques de posiciones con reads."""
        out = np.zeros(self.L, dtype=np.float32)
        k, c = self._datos()
        if not k.size:
            # Sin reads el original da ceros float32 de largo L: lo mismo.
            return out
        pos = k // self.paso
        resto = k - pos * self.paso
        ini = np.flatnonzero(np.r_[True, pos[1:] != pos[:-1]])
        upos = pos[ini]
        fin = np.r_[ini[1:], k.size]
        por_bloque = max(1, BLOQUE_B // (self.paso * 4))
        for b0 in range(0, upos.size, por_bloque):
            b1 = min(b0 + por_bloque, upos.size)
            e0, e1 = ini[b0], fin[b1 - 1]
            fila = np.repeat(np.arange(b1 - b0), fin[b0:b1] - ini[b0:b1])
            denso = np.zeros((b1 - b0, self.C, self.N, 2), dtype=np.uint32)
            denso.reshape(-1)[fila * self.paso + resto[e0:e1]] = c[e0:e1].astype(np.uint32)
            coverage = np.sum(denso, axis=3, dtype='uint32')
            coverage = np.multiply(coverage, rpbs)
            coverage = np.median(coverage, axis=(2))
            coverage = coverage.astype('float32')
            coverage = np.mean(coverage, axis=(1))
            out[upos[b0:b1]] = coverage
        return out


class Largos:
    """En lugar de [None] * L con una lista de largos por posicion:
    (clave, cuenta, primera aparicion) con clave = pos << 16 | largo."""

    def __init__(self, largo):
        self.L = int(largo)
        self._buf = array.array("q")
        self._t = 0
        self._k = np.empty(0, np.int64)
        self._n = np.empty(0, np.int64)
        self._p = np.empty(0, np.int64)

    def agregar(self, pos, largo):
        """pos_size_d[chrom][pos].append(largo), o = [largo] si estaba en None"""
        if _activo.get(Largos) is not self:
            _cambiar_activo(self)
        pos = _indice_lista(pos, self.L)
        if not 0 <= largo < (1 << _BITS_LARGO):
            raise ValueError(f"largo de read fuera de rango para el parche: {largo}")
        self._buf.append((pos << _BITS_LARGO) | largo)
        if len(self._buf) >= COMPACTAR:
            self._compactar()

    def _compactar(self):
        if not len(self._buf):
            return
        b = np.frombuffer(self._buf, dtype=np.int64)
        k, p, n = np.unique(b, return_index=True, return_counts=True)
        p = p.astype(np.int64) + self._t
        self._t += len(b)
        self._buf = array.array("q")
        if self._k.size:
            k, (n, p) = _fusionar(np.concatenate((self._k, k)),
                                  [np.concatenate((self._n, n.astype(np.int64))),
                                   np.concatenate((self._p, p))],
                                  [np.add, np.minimum])
        self._k, self._n, self._p = k, n.astype(np.int64), p

    def __getitem__(self, p):
        if isinstance(p, slice):
            raise TypeError("rebanada no soportada por el parche")
        p = _indice_lista(int(p), self.L)
        self._compactar()
        i0 = np.searchsorted(self._k, p << _BITS_LARGO, "left")
        i1 = np.searchsorted(self._k, (p + 1) << _BITS_LARGO, "left")
        if i0 == i1:
            return None
        largos = (self._k[i0:i1] & ((1 << _BITS_LARGO) - 1)).tolist()
        cuentas = self._n[i0:i1].tolist()
        out = []
        for j in np.argsort(self._p[i0:i1], kind="stable").tolist():
            out.extend([largos[j]] * cuentas[j])
        return out

    def __setitem__(self, p, v):
        raise TypeError("Largos se llena con agregar(); el parche no esperaba esta escritura")
