"""Exact certificates for the s=2 GHRS Hamming-QMDS rank criterion.

The script works over the rational numbers.  Given integer evaluation points
beta_1,...,beta_r, it verifies all rank conditions from Proposition
``Rank criterion for the s=2 Hamming-QMDS property`` in the main manuscript.

For each pair (A,B) with |A|+|B|=r+1 it records two witnesses:

* det M_{A,B}^{(r)} != 0;
* an (r-1)x(r-1) minor of M_{A,B}^{(r-2)} != 0.

Together with the denominators/numerators of the sums S_j and the evaluation
point differences, these witnesses give the finite exceptional prime set used
in the rational lifting theorem.  No floating-point arithmetic is used.
"""

from fractions import Fraction
from itertools import combinations
import argparse


def primes_of(n):
    n = abs(n)
    out = set()
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.add(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out.add(n)
    return out


def det_fraction(mat):
    a = [[Fraction(x) for x in row] for row in mat]
    n = len(a)
    det = Fraction(1)
    for c in range(n):
        piv = None
        for i in range(c, n):
            if a[i][c]:
                piv = i
                break
        if piv is None:
            return Fraction(0)
        if piv != c:
            a[c], a[piv] = a[piv], a[c]
            det = -det
        pivot = a[c][c]
        det *= pivot
        for i in range(c + 1, n):
            fac = a[i][c] / pivot
            if fac:
                for j in range(c, n):
                    a[i][j] -= fac * a[c][j]
    return det


def rank_fraction(mat):
    a = [[Fraction(x) for x in row] for row in mat]
    if not a:
        return 0
    m, n = len(a), len(a[0])
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if a[i][c]:
                piv = i
                break
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        pivot = a[r][c]
        a[r] = [x / pivot for x in a[r]]
        for i in range(m):
            if i != r and a[i][c]:
                fac = a[i][c]
                a[i] = [a[i][j] - fac * a[r][j] for j in range(n)]
        r += 1
    return r


def first_nonzero_minor(mat, size):
    row_count = len(mat)
    col_count = len(mat[0]) if mat else 0
    for rows in combinations(range(row_count), size):
        for cols in combinations(range(col_count), size):
            sub = [[mat[i][j] for j in cols] for i in rows]
            d = det_fraction(sub)
            if d:
                return d, rows, cols
    return Fraction(0), None, None


def rows_for(alphas, A, B, degree):
    rows = []
    for j in A:
        a = alphas[j]
        rows.append([a**k for k in range(degree + 1)])
    for j in B:
        a = alphas[j]
        rows.append([0 if k == 0 else k * a ** (k - 1) for k in range(degree + 1)])
    return rows


def square_class_rank(rationals):
    primes = []
    rows = []

    def factor(n):
        n = abs(n)
        out = []
        d = 2
        while d * d <= n:
            while n % d == 0:
                out.append(d)
                n //= d
            d += 1 if d == 2 else 2
        if n > 1:
            out.append(n)
        return out

    for x in rationals:
        for p in factor(x.numerator) + factor(x.denominator):
            if p not in primes:
                primes.append(p)
    basis = ["-1"] + primes
    for x in rationals:
        row = [1 if x < 0 else 0]
        num_fac = factor(x.numerator)
        den_fac = factor(x.denominator)
        for p in primes:
            row.append((num_fac.count(p) + den_fac.count(p)) % 2)
        rows.append(row)

    rank = 0
    a = [row[:] for row in rows]
    for col in range(len(basis)):
        piv = None
        for i in range(rank, len(a)):
            if a[i][col]:
                piv = i
                break
        if piv is None:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        for i in range(len(a)):
            if i != rank and a[i][col]:
                a[i] = [x ^ y for x, y in zip(a[i], a[rank])]
        rank += 1
    return rank, basis, rows


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def is_square_mod(x, p):
    if x.numerator % p == 0 or x.denominator % p == 0:
        return False
    val = (x.numerator % p) * pow(x.denominator % p, p - 2, p) % p
    return pow(val, (p - 1) // 2, p) == 1


def first_prime_field_hits(gammas, bad, limit=10, search_bound=100000):
    hits = []
    for p in range(3, search_bound):
        if p in bad or not is_prime(p):
            continue
        if all(is_square_mod(g, p) for g in gammas):
            hits.append(p)
            if len(hits) >= limit:
                return hits
    return hits


def good_minor_primes(alphas):
    r = len(alphas)
    bad = set()
    determinants = []
    indices = range(r)
    for a_size in range(r + 1):
        b_size = r + 1 - a_size
        if b_size > r:
            continue
        for A in combinations(indices, a_size):
            for B in combinations(indices, b_size):
                small = rows_for(alphas, A, B, r - 2)
                big = rows_for(alphas, A, B, r)
                small_rank = rank_fraction(small)
                big_rank = rank_fraction(big)
                if big_rank != small_rank + 2:
                    raise SystemExit(f"rank criterion fails over Q for A={A}, B={B}")
                # Pick one nonzero full-rank minor of the big matrix. Here rows are
                # r+1 and columns are r+1, so the determinant suffices.
                d = det_fraction(big)
                if d == 0:
                    raise SystemExit(f"unexpected zero determinant A={A}, B={B}")
                determinants.append(d)
                bad |= primes_of(d.numerator)
                bad |= primes_of(d.denominator)
                small_d, _, _ = first_nonzero_minor(small, r - 1)
                if small_d == 0:
                    raise SystemExit(f"no small rank witness A={A}, B={B}")
                determinants.append(small_d)
                bad |= primes_of(small_d.numerator)
                bad |= primes_of(small_d.denominator)
    # Distinctness of evaluation points and nonzero S_j denominators.
    for i, a in enumerate(alphas):
        for b in alphas[i + 1 :]:
            bad |= primes_of(a - b)
    for j, a in enumerate(alphas):
        S = sum(Fraction(1, a - b) for b in alphas if b != a)
        if S == 0:
            raise SystemExit(f"S_j vanishes over Q at {a}")
        bad |= primes_of(S.numerator)
        bad |= primes_of(S.denominator)
    gammas = []
    for j, a in enumerate(alphas):
        gammas.append(-2 * sum(Fraction(1, a - b) for b in alphas if b != a))
    return sorted(bad), determinants, gammas


def print_report(alphas, prime_hits=10):
    bad, dets, gammas = good_minor_primes(alphas)
    rank, basis, rows = square_class_rank(gammas)
    hits = first_prime_field_hits(gammas, set(bad), limit=prime_hits)
    print("evaluation set:", alphas)
    print("r:", len(alphas))
    print("rank conditions:", sum(1 for _ in range(len(dets) // 2)))
    print("witness determinants/minors:", len(dets))
    print("bad primes:", bad)
    print("gamma=-2S_j:", [str(g) for g in gammas])
    print("square-class basis:", basis)
    print("square-class vectors:", rows)
    print("square-class rank:", rank)
    print("prime-field density:", f"1/{2**rank}")
    print("first prime-field hits:", hits)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Exact rational certificates for s=2 GHRS Hamming-QMDS rank tests."
    )
    parser.add_argument(
        "alphas",
        nargs="*",
        type=int,
        help="integer evaluation points, e.g. 3 6 12 17 24",
    )
    parser.add_argument("--prime-hits", type=int, default=10)
    args = parser.parse_args()

    examples = [args.alphas] if args.alphas else [[3, 6, 12, 17, 24], [1, -1, 2, -2]]
    for alphas in examples:
        print_report(alphas, prime_hits=args.prime_hits)
