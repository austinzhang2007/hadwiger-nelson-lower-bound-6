"""Independent rational arithmetic for the primary 32-term coordinate field.

No search module, SymPy arithmetic, or floating-point predicate is imported.
The basis is 3^(a/2) 5^(b/2) 11^(c/2) h_minus^d h_plus^e,
where h_minus^2=(25-3*sqrt(33))/8 and h_plus^2=(25+3*sqrt(33))/8.
The multiplication table is reduced directly from these five relations.
"""

from functools import lru_cache
import re

from gmpy2 import mpq


class PrimaryField:
    def __init__(self):
        self.zero = (mpq(0),) * 32
        self.one = (mpq(1),) + self.zero[1:]
        self.table = tuple(tuple(self._basis_product(i, j) for j in range(32))
                           for i in range(32))

    @staticmethod
    def element(values):
        texts = tuple(str(v) for v in values)
        if any(re.fullmatch(r"[+-]?\d+(?:/[+-]?\d+)?", v) is None for v in texts):
            raise ValueError("coefficients must be explicit exact rationals")
        result = tuple(mpq(v) for v in texts)
        if len(result) != 32:
            raise ValueError("expected 32 exact rational coefficients")
        return result

    @staticmethod
    @lru_cache(maxsize=None)
    def _reduce(exponents):
        """Reduce a monomial, eliminating h powers before base radicals."""
        for bit in (4, 3, 2, 1, 0):
            if exponents[bit] < 2:
                continue
            lower = list(exponents)
            lower[bit] -= 2
            if bit < 3:
                square = (3, 5, 11)[bit]
                return tuple((m, c * square) for m, c in
                             PrimaryField._reduce(tuple(lower)))
            out = {}
            for m, c in PrimaryField._reduce(tuple(lower)):
                out[m] = out.get(m, mpq(0)) + c * mpq(25, 8)
            lower[0] += 1
            lower[2] += 1
            coefficient = mpq(3 if bit == 4 else -3, 8)
            for m, c in PrimaryField._reduce(tuple(lower)):
                out[m] = out.get(m, mpq(0)) + c * coefficient
            return tuple((m, c) for m, c in sorted(out.items()) if c)
        mask = sum(e << i for i, e in enumerate(exponents))
        return ((mask, mpq(1)),)

    @classmethod
    def _basis_product(cls, i, j):
        return cls._reduce(tuple(((i >> b) & 1) + ((j >> b) & 1)
                                 for b in range(5)))

    def mul(self, a, b):
        out = [mpq(0) for _ in range(32)]
        right = [(j, y) for j, y in enumerate(b) if y]
        for i, x in enumerate(a):
            if not x:
                continue
            row = self.table[i]
            for j, y in right:
                xy = x * y
                for mask, coefficient in row[j]:
                    out[mask] += xy * coefficient
        return tuple(out)

    @staticmethod
    def sub(a, b):
        return tuple(x - y for x, y in zip(a, b))

    def unit(self, left, right):
        """Prove the homogeneous squared-distance identity coefficientwise."""
        x, y, d = left
        u, v, e = right
        if d == self.zero or e == self.zero:
            raise ValueError("zero projective denominator")
        dx = self.sub(self.mul(x, e), self.mul(u, d))
        dy = self.sub(self.mul(y, e), self.mul(v, d))
        de = self.mul(d, e)
        return all(a + b == c for a, b, c in zip(
            self.mul(dx, dx), self.mul(dy, dy), self.mul(de, de)
        ))
