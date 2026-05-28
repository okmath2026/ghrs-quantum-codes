"""Finite-field search utilities for GHRS quantum-code examples.

This script supports the small parameter table in the main manuscript.  It
constructs GHRS generator matrices over prime fields and computes relative
Hamming distances by rank:

    d(D/C) = min{|S| : dim D_S > dim C_S}.

The script is intended for exact exploratory searches over small fields.  The
characteristic-zero determinant certificates used for infinite families are
handled separately by ``s2_bad_primes.py``.
"""

from itertools import combinations, product
from math import comb


def inv(a, p):
    return pow(a % p, p - 2, p)


def rank_mod(mat, p):
    a = [row[:] for row in mat if any(x % p for x in row)]
    if not a:
        return 0
    m, n = len(a), len(a[0])
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if a[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        scale = inv(a[r][c], p)
        a[r] = [(x * scale) % p for x in a[r]]
        for i in range(m):
            if i != r and a[i][c] % p:
                fac = a[i][c] % p
                a[i] = [(a[i][j] - fac * a[r][j]) % p for j in range(n)]
        r += 1
        if r == m:
            break
    return r


def null_vector_one_dim(mat, p):
    a = [row[:] for row in mat]
    m = len(a)
    n = len(a[0])
    pivots = []
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if a[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        scale = inv(a[r][c], p)
        a[r] = [(x * scale) % p for x in a[r]]
        for i in range(m):
            if i != r and a[i][c] % p:
                fac = a[i][c] % p
                a[i] = [(a[i][j] - fac * a[r][j]) % p for j in range(n)]
        pivots.append(c)
        r += 1
    free = [c for c in range(n) if c not in pivots]
    if len(free) != 1:
        return None
    x = [0] * n
    x[free[0]] = 1
    for row, c in enumerate(pivots):
        x[c] = (-a[row][free[0]]) % p
    return x


def nullspace_basis(mat, p):
    if not mat:
        return []
    a = [row[:] for row in mat]
    m = len(a)
    n = len(a[0])
    pivots = []
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if a[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        scale = inv(a[r][c], p)
        a[r] = [(x * scale) % p for x in a[r]]
        for i in range(m):
            if i != r and a[i][c] % p:
                fac = a[i][c] % p
                a[i] = [(a[i][j] - fac * a[r][j]) % p for j in range(n)]
        pivots.append(c)
        r += 1
    free = [c for c in range(n) if c not in pivots]
    basis = []
    for f in free:
        x = [0] * n
        x[f] = 1
        for row, c in enumerate(pivots):
            x[c] = (-a[row][f]) % p
        basis.append(x)
    return basis


def hasse_monomial(k, a, alpha, p):
    if k < a:
        return 0
    return (comb(k, a) * pow(alpha, k - a, p)) % p


def lambda_vector(p, r, s, alphas):
    n = r * s
    mat = []
    for k in range(n - 1):
        row = []
        for j in range(r):
            for i in range(s):
                row.append(hasse_monomial(k, i, alphas[j], p))
        mat.append(row)
    return null_vector_one_dim(mat, p)


def square_root_table(p):
    roots = {}
    for x in range(p):
        roots.setdefault((x * x) % p, x)
    return roots


def admissible_multiplier(p, lam):
    if any(x % p == 0 for x in lam):
        return None
    roots = square_root_table(p)
    for mu in range(1, p):
        vals = [(mu * x) % p for x in lam]
        if all(v in roots for v in vals):
            return [roots[v] for v in vals], mu
    return None


def generator(p, r, s, alphas, mult, m):
    rows = []
    for k in range(m):
        row = []
        for j in range(r):
            for i in range(s):
                row.append((mult[j * s + i] * hasse_monomial(k, i, alphas[j], p)) % p)
        rows.append(row)
    return rows


def min_distance_by_rank(G, p):
    k = rank_mod(G, p)
    n = len(G[0]) if G else 0
    if k == 0:
        return 0
    cols = range(n)
    for w in range(1, n + 1):
        for supp in combinations(cols, w):
            outside = [c for c in cols if c not in supp]
            restricted = [[row[c] for c in outside] for row in G]
            if rank_mod(restricted, p) < k:
                return w
    return n + 1


def supported_subcode_dim(G, p, supp):
    k = rank_mod(G, p)
    n = len(G[0]) if G else 0
    outside = [c for c in range(n) if c not in supp]
    restricted = [[row[c] for c in outside] for row in G]
    return k - rank_mod(restricted, p)


def relative_distance_by_rank(G_big, G_sub, p):
    n = len(G_big[0]) if G_big else 0
    for w in range(1, n + 1):
        for supp in combinations(range(n), w):
            if supported_subcode_dim(G_big, p, supp) > supported_subcode_dim(G_sub, p, supp):
                return w
    return None


def gram_rank(G, p):
    gram = []
    for row in G:
        gram.append([sum(row[c] * other[c] for c in range(len(row))) % p for other in G])
    return rank_mod(gram, p)


def block_distance_poly(r, s, m):
    return r - ((m - 1) // s)


def search():
    rows = []
    for p in [5, 7, 11, 13]:
        for s in [2, 3]:
            for r in range(2, min(p, 6) + 1):
                n = r * s
                if n > 12:
                    continue
                for alphas in combinations(range(p), r):
                    lam = lambda_vector(p, r, s, list(alphas))
                    if lam is None:
                        continue
                    adm = admissible_multiplier(p, lam)
                    if adm is None:
                        continue
                    mult, mu = adm
                    for t in range(1, n // 2 + 1):
                        Gt = generator(p, r, s, list(alphas), mult, t)
                        Gdual = generator(p, r, s, list(alphas), mult, n - t)
                        d_dual = min_distance_by_rank(Gdual, p)
                        d_rel = relative_distance_by_rank(Gdual, Gt, p)
                        d_c = min_distance_by_rank(Gt, p)
                        c_ea = gram_rank(Gdual, p)
                        rows.append(
                            {
                                "q": p,
                                "r": r,
                                "s": s,
                                "n": n,
                                "t": t,
                                "alphas": alphas,
                                "mu": mu,
                                "d_css": d_dual,
                                "d_rel": d_rel,
                                "bound": (t + s) // s,
                                "d_Ct": d_c,
                                "d_blk_dual": block_distance_poly(r, s, n - t),
                                "c_ea_for_Ct": c_ea,
                            }
                        )
                    break
    return rows


def search_low_ea():
    rows = []
    for p, r, s, alphas in [
        (5, 3, 2, (0, 1, 2)),
    ]:
        n = r * s
        for t in range(2, n):
            best = None
            checked = 0
            for mult in product(range(1, p), repeat=n):
                G = generator(p, r, s, list(alphas), list(mult), t)
                H = nullspace_basis(G, p)
                c = gram_rank(H, p)
                d = min_distance_by_rank(G, p)
                checked += 1
                cand = (c, -d, mult, d)
                if best is None or cand < best:
                    best = cand
                    if c == 0:
                        break
                if checked >= 5000:
                    break
            rows.append(
                {
                    "q": p,
                    "r": r,
                    "s": s,
                    "n": n,
                    "t": t,
                    "alphas": alphas,
                    "best_c": best[0],
                    "d": best[3],
                    "k_ea": 2 * t - n + best[0],
                    "checked": checked,
                    "mult": best[2],
                }
            )
    return rows


if __name__ == "__main__":
    rows = search()
    interesting = [x for x in rows if x["d_css"] > x["bound"]]
    print("total admissible rows:", len(rows))
    print("improved Hamming cases:", len(interesting))
    print()
    for x in rows[:40]:
        print(x)
    print("\nImproved:")
    for x in interesting[:80]:
        print(x)
    print("\nLow-entanglement EA search:")
    for x in search_low_ea():
        print(x)
