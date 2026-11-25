"""Core algorithm: uniform embeddings for Robinson matrices (SPEC.md).

Pure stdlib; trusts its input; no validation, no try/except/raise.
Each function's docstring states its input/output contract; inline
comments reference the numbered pseudocode lines (L#) of SPEC.md.
All indices are 0-based.
"""

import math
from fractions import Fraction

INF = math.inf


# ---- helpers (SPEC.md §2) ----

def unit(k, j):
    """Length-k int tuple with 1 at index j, 0 elsewhere (0 <= j < k)."""
    b = [0] * k
    b[j] = 1
    return tuple(b)


def zero(k):
    """Length-k int tuple of zeros."""
    return (0,) * k


def vec_add(x, y):
    """Componentwise sum of two equal-length int tuples."""
    return tuple(xi + yi for xi, yi in zip(x, y))


def neg(x):
    """Componentwise negation of an int tuple."""
    return tuple(-xi for xi in x)


def dot(x, y):
    """Inner product x·y; exact (Fraction) or float according to y's type."""
    return sum(xi * yi for xi, yi in zip(x, y))


# ---- primitives (SPEC.md §3) ----

def beta_plus(a, k, u, v):
    """Upper-bound coefficient vector of the edge (u, v); u != v.

    a: Robinson matrix (tuple of tuples of int); k = a[0][0]; u, v: vertices.
    Returns a length-k int tuple b such that pi[v] - pi[u] < b·d is the
    edge's upper-bound inequality, or None if the edge imposes no upper
    bound (u < v and a[u][v] = 0).
    """
    if u < v:                                  # L1
        t = a[u][v]                            # L2
        return unit(k, t - 1) if t >= 1 else None   # L3
    return neg(beta_minus(a, k, v, u))         # L4


def beta_minus(a, k, u, v):
    """Lower-bound coefficient vector of the edge (u, v); u != v.

    Returns a length-k int tuple b such that pi[v] - pi[u] > b·d is the
    edge's lower-bound inequality (zero(k) when the lower bound is 0),
    or None if there is none (u > v and a[v][u] = 0).
    """
    if u < v:                                  # L1
        t = a[u][v]                            # L2
        return unit(k, t) if t <= k - 1 else zero(k)   # L3
    b = beta_plus(a, k, v, u)                  # L4
    return neg(b) if b is not None else None   # L5


def walk_bound(a, k, w):
    """Bound of an upper-bound walk w (every up-step has a > 0).

    Returns the length-k int tuple b = sum of the edge bounds, so
    pi[w[-1]] - pi[w[0]] < b·d. Trusts w to be a UBW (no None edges).
    """
    b = zero(k)                                # L1
    for i in range(1, len(w)):                 # L2
        b = vec_add(b, beta_plus(a, k, w[i - 1], w[i]))   # L3
    return b                                   # L4


def simple(w):
    """True iff the walk w is simple: open with all vertices distinct,
    or closed with len(w) >= 3, w[0] == w[-1] and no other repetition.
    A one-vertex walk and (v, v) are not simple.
    """
    if len(w) < 2:                             # L1
        return False
    if w[0] == w[-1]:                          # L2
        return len(w) >= 3 and len(set(w[:-1])) == len(w) - 1
    return len(set(w)) == len(w)               # L3


def preceq(x, y):
    """Prefix order x ≼ y on equal-length int tuples: True iff every
    prefix sum of x is <= the corresponding prefix sum of y, i.e.
    x·d <= y·d for every threshold vector d.
    """
    sx = sy = 0
    for j in range(len(x)):                    # L1  j = 1 .. k (prefix lengths)
        sx += x[j]
        sy += y[j]
        if sx > sy:                            # L2
            return False
    return True                                # L3


def cell_items(cell):
    """The (walk, bound) pairs of a cell: the values of the k=2 slot
    dict, or the list itself for k != 2. The only way to iterate a cell.
    """
    return cell.values() if isinstance(cell, dict) else cell   # L1


# ---- bound generation (SPEC.md §4) ----

def slot_key(b):
    """Chain slot of a k=2 bound: (|b[0]| + |b[1]|, 1 if b[1] > 0 else 2).
    Bounds in one slot form a chain under ≼, so only the slot minimum
    is kept.
    """
    return (abs(b[0]) + abs(b[1]), 1 if b[1] > 0 else 2)


def minimal_in(b, cell):
    """True iff b is ≼-minimal against the cell: no stored bound c with
    c ≼ b and c != b. Used for k != 2 (list cells).
    """
    return not any(preceq(c, b) and c != b for (_, c) in cell_items(cell))


def insert_minimal(cell, w, b):
    """Antichain insertion into a list cell: skip if b is already stored
    exactly (equal bound, possibly a different walk) or a stored bound is
    strictly ≼-below b; otherwise add (w, b) and remove every stored
    pair whose bound is strictly ≼-above b. One pair per distinct bound.
    """
    if any(c == b for (_, c) in cell_items(cell)):
        return
    if not minimal_in(b, cell):
        return
    cell[:] = [(w2, c) for (w2, c) in cell
               if not (preceq(b, c) and b != c)]
    cell.append((w, b))


def insert_minimal_k2(cell, w, b):
    """O(1) insertion into a k=2 slot dict: occupy the empty slot, or
    replace the slot's pair only if b is strictly ≼-below its bound.
    """
    key = slot_key(b)                          # L1
    if key not in cell:                        # L2
        cell[key] = (w, b)
    elif preceq(b, cell[key][1]) and b != cell[key][1]:   # L3
        cell[key] = (w, b)


def insert(cell, w, b, k):
    """Insert (w, b) into a cell if w is simple; non-simple walks are
    dropped. Dispatches on k: slot rule for k=2, antichain rule otherwise.
    """
    if not simple(w):                          # L1
        return
    if k == 2:                                 # L2
        insert_minimal_k2(cell, w, b)
    elif minimal_in(b, cell):                  # L3
        insert_minimal(cell, w, b)


def bound_generation(a):
    """Build the bound table T of a Robinson matrix a (trusted input).

    Returns an n×n list of cells: T[i][j] holds the ≼-minimal simple
    (i,j)-upper-bound walks as (walk, bound) pairs; diagonal cells hold
    simple upper-bound cycles (empty until the first cycle is found).
    Cells are slot dicts for k=2, lists otherwise. Every stored bound
    has sum |b[j]| <= n.
    """
    n = len(a)                                 # L1
    k = a[0][0]
    T = [[{} if k == 2 else [] for _ in range(n)] for _ in range(n)]   # L2
    for i in range(n):                         # L3  seed direct walks, i < j only
        for j in range(i + 1, n):              # L4
            if a[i][j] != 0:                   # L5
                insert(T[i][j], (i, j), beta_plus(a, k, i, j), k)
            insert(T[j][i], (j, i), beta_plus(a, k, j, i), k)   # L6
    for s in range(n):                         # L7  sweep ALL ordered pairs
        for i in range(n):                     # L8
            for j in range(n):                 # L9
                for (w1, b1) in list(cell_items(T[i][s])):   # L10  snapshot
                    for (w2, b2) in list(cell_items(T[s][j])):   # L11  snapshot
                        w = w1 + w2[1:]        # L12
                        b = vec_add(b1, b2)
                        insert(T[i][j], w, b, k)   # L13
    return T                                   # L14


# ---- cycle bounds (SPEC.md §5) ----

def cycle_bounds(T):
    """Distinct minimal cycle bounds of the table: a list of (bound, walk)
    pairs, one walk per bound, taken from the diagonal cells. Every walk
    is already a simple cycle (bound generation stores nothing else).
    """
    n = len(T)                                 # L1
    seen = {}
    for v in range(n):                         # L2
        for (w, b) in cell_items(T[v][v]):     # L3
            if b not in seen:                  # L4
                seen[b] = w
    return [(b, seen[b]) for b in seen]        # L5


# ---- given thresholds (SPEC.md §5a) ----

def solve_with_thresholds(a, d):
    """Verify a candidate threshold vector d (length-k Fraction tuple,
    strictly decreasing, positive) and construct the embedding.

    Returns (d, pi) with d unchanged if every stored cycle bound has
    b·d > 0; otherwise (None, (w,)) where w is a simple cycle whose
    bound fails for this d.
    """
    n = len(a)                                 # L1
    T = bound_generation(a)
    for v in range(n):                         # L2
        for (w, b) in cell_items(T[v][v]):     # L3
            if dot(b, d) <= 0:                 # L4
                return None, (w,)
    return d, embedding_from_bounds(a, d, T)   # L5


# ---- k = 2, combinatorial (SPEC.md §5b) ----

def ratio_bounds(T):
    """Intersect the cycle inequalities as constraints on r = d[1]/d[0],
    seeded with the window (0, 1) required by d in D^2. k = 2 only.

    Returns ("FEASIBLE", (lo, hi)): an open rational window (Fractions)
    such that any r in (lo, hi) satisfies every cycle inequality.
    Returns ("INFEASIBLE", (w,)): one cycle unsatisfiable for every r
    in (0, 1); or ("INFEASIBLE", (w_lo, w_hi)): two cycles whose ratio
    intervals are disjoint. All witnesses are simple cycles.
    """
    n = len(T)                                 # L1  window seeded to (0, 1)
    lo = Fraction(0)
    hi = Fraction(1)
    w_lo = None
    w_hi = None
    for v in range(n):                         # L2
        for (w, b) in cell_items(T[v][v]):     # L3
            (a1, a2) = b                       # L4
            if a2 > 0:                         # L5
                r = Fraction(-a1, a2)          # L6
                if r >= 1:                     # L7  contradicts r < 1
                    return "INFEASIBLE", (w,)
                if r > lo:                     # L8
                    lo = r
                    w_lo = w
            elif a2 < 0:                       # L9
                r = Fraction(-a1, a2)          # L10
                if r <= 0:                     # L11  contradicts r > 0
                    return "INFEASIBLE", (w,)
                if r < hi:                     # L12
                    hi = r
                    w_hi = w
            else:                              # L13
                if a1 <= 0:                    # L14
                    return "INFEASIBLE", (w,)
    if lo >= hi:                               # L15
        return "INFEASIBLE", (w_lo, w_hi)
    return "FEASIBLE", (lo, hi)                # L16


def solve_k2(a):
    """k = 2 automatic search, exact rational arithmetic.

    Returns (d, pi) with d = (1, (lo+hi)/2) (Fractions) and pi the
    constructed embedding, or (None, witnesses) as from ratio_bounds.
    """
    T = bound_generation(a)                    # L1
    status, info = ratio_bounds(T)             # L2
    if status == "INFEASIBLE":                 # L3
        return None, info
    (lo, hi) = info                            # L4
    d = (Fraction(1), (lo + hi) / 2)
    return d, embedding_from_bounds(a, d, T)   # L5


# ---- k > 2, LP feasibility (SPEC.md §5c); search is injected, never imported ----

def solve_general(a, search):
    """k > 2 automatic search via the injected callable search(cycles, k)
    (the public layer passes _lp.search_thresholds; the core imports no
    third-party module).

    Returns (d, pi) with d a float tuple from the LP and pi the
    constructed embedding, or (None, witnesses): the deletion-filter
    set of cycles, jointly contradictory with d in D^k.
    """
    k = a[0][0]                                # L1
    T = bound_generation(a)
    d, witnesses = search(cycle_bounds(T), k)  # L2
    if d is None:                              # L3
        return None, witnesses
    return d, embedding_from_bounds(a, d, T)   # L4


# ---- k = 1, trivial (SPEC.md §5d) ----

def solve_k1(a):
    """k = 1: every binary Robinson matrix is feasible (indifference
    graph), so no search is needed. Always returns (d, pi) with
    d = (Fraction(1),); never fails.
    """
    T = bound_generation(a)                    # L1
    d = (Fraction(1),)                         # L2
    return d, embedding_from_bounds(a, d, T)   # L3


# ---- construction (SPEC.md §6) ----

def min_bound(cell, d):
    """Smallest b·d over the cell's bounds (cell nonempty, trusted).
    For an up-cell T[i][v] this is the tightest upper bound on
    pi[v] - pi[i]; for a down-cell T[v][i] its negation is the tightest
    lower bound (lower bounds are read by reversal).
    """
    return min(dot(b, d) for (_, b) in cell_items(cell))   # L1


def embedding_from_bounds(a, d, T):
    """Construct pi from a feasible threshold vector d and the table T.

    Returns a length-n list with pi[0] = 0, same numeric type as d;
    each pi[v] is the midpoint of the tightest lower and upper bounds
    (lb + d[0] if no upper bound exists). For feasible d the result is
    strictly increasing and satisfies every pair inequality.
    """
    n = len(a)                                 # L1
    pi = [0] * n
    for v in range(1, n):                      # L2
        ubs = [pi[i] + min_bound(T[i][v], d)   # L3  skip empty cells
               for i in range(v) if cell_items(T[i][v])]
        ub = min(ubs) if ubs else INF
        lb = max(pi[i] - min_bound(T[v][i], d) for i in range(v))   # L4
        pi[v] = (ub + lb) / 2 if ub < INF else lb + d[0]   # L5
    return pi                                  # L6


def verify_embedding(a, d, pi):
    """Exact boolean check: True iff lower(t) < pi[v] - pi[u] < upper(t)
    for every u < v with t = a[u][v]. Since lower(t) >= 0, this forces
    pi strictly increasing. Raises nothing; the public layer runs it on
    every successful return.
    """
    n = len(a)                                 # L1
    k = a[0][0]
    for u in range(n):                         # L2
        for v in range(u + 1, n):              # L3
            t = a[u][v]                        # L4
            lo = d[t] if t <= k - 1 else 0     # L5
            hi = d[t - 1] if t >= 1 else INF   # L6
            x = pi[v] - pi[u]                  # L7
            if not (lo < x) or not (x < hi):   # L8
                return False
    return True                                # L9


# ---- top-level dispatch (SPEC.md §7) ----

def uniform_embedding_core(a, d=None, search=None):
    """Dispatch on the mode (SPEC.md §7); a is a trusted Robinson matrix.

    d given: verify it and construct (any k). d None: k=1 trivial,
    k=2 combinatorial, k>2 via the injected search callable.
    Returns (d, pi) on success — pi strictly increasing, d scale-free
    (d[0] = 1) except a supplied d, which is returned unchanged — or
    (None, witnesses), a tuple of simple witness cycles. Never raises.
    """
    k = a[0][0]                                # L1
    if d is not None:                          # L2
        return solve_with_thresholds(a, d)
    if k == 1:                                 # L3
        return solve_k1(a)
    if k == 2:                                 # L4
        return solve_k2(a)
    return solve_general(a, search)            # L5