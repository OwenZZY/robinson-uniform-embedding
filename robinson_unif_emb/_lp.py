"""LP backend for the k>2 threshold search (SPEC.md §5c; SPEC-engineering.md §4).

SciPy is imported lazily inside _solve; this module may raise ImportError
(SciPy absent) or RuntimeError (unexpected solver status). Inputs are trusted.
"""


def _solve(cycles, k):
    """The max-slack LP of SPEC.md §5c L1–L7 with d[0] = 1 substituted (L6).

    cycles: list of (bound, walk) pairs from _core.cycle_bounds; k >= 3.
    Variables x = (d[1..k-1], t) (L1); minimize -t (L2); t is x[k-1].
    Returns (feasible, d): feasible iff status 0 and t* > EPS; d is the
    optimal threshold tuple (1.0, d[1], ..., d[k-1]) on status 0, else None.
    Raises ImportError if SciPy is absent, RuntimeError on any solver
    status other than 0 (optimal) or 2 (infeasible).
    """
    try:
        from scipy.optimize import linprog
    except ImportError:
        raise ImportError(
            "install scipy for k>2 automatic search: pip install scipy") from None

    c = [0] * (k - 1) + [-1]                   # L2  minimize -t
    A, rhs = [], []                            # A_ub x <= b_ub (negated >= rows)
    for (b, _) in cycles:                      # L3  b·d - t >= 0, d[0] = 1
        A.append([-b[j] for j in range(1, k)] + [1])
        rhs.append(b[0])
    for i in range(k - 1):                     # L4  d[i] - d[i+1] - t >= 0
        row = [0] * k
        row[k - 1] = 1                         # +t
        if i == 0:                             #     d[0] = 1:  1 - d[1] - t >= 0
            row[0] = 1
            rhs.append(1)
        else:
            row[i - 1] = -1                    # -d[i]
            row[i] = 1                         # +d[i+1]
            rhs.append(0)
        A.append(row)
    row = [0] * k                              # L5  d[k-1] - t >= 0
    row[k - 1] = 1
    row[k - 2] = -1
    rhs.append(0)
    A.append(row)

    res = linprog(c, A_ub=A, b_ub=rhs, bounds=[(None, None)] * k, method="highs")
    if res.status not in (0, 2):
        raise RuntimeError(f"LP solver failure: status={res.status}")
    scale = 1 + max((sum(abs(bj) for bj in b) for (b, _) in cycles), default=0)   # L7
    eps = 1e-9 * scale
    if res.status == 2:
        return False, None
    d = (1.0,) + tuple(float(v) for v in res.x[:k - 1])
    return float(res.x[k - 1]) > eps, d


def search_thresholds(cycles, k):
    """Decide feasibility of {b·d > 0 for each cycle bound; d in D^k},
    by homogeneity with d[0] = 1 (SPEC.md §5c).

    cycles: list of (bound, walk) pairs from _core.cycle_bounds; k >= 3.
    Returns (d, None) with d a length-k float tuple, d[0] = 1, if the
    optimal slack t* > EPS; otherwise (None, witnesses), where witnesses
    is the deletion-filter set of walks, jointly contradictory with
    d in D^k (no single one need fail, no count bound).
    Raises ImportError (SciPy absent) or RuntimeError (solver failure).
    """
    feasible, d = _solve(cycles, k)            # L1–L7
    if not feasible:                           # L8
        return None, minimal_witnesses(cycles, k)
    return d, None                             # L9


def minimal_witnesses(cycles, k):
    """Deletion filter (SPEC.md §5c L1–L5): an irreducible infeasible
    subset of the cycle system.

    cycles: list of (bound, walk) pairs; k >= 3. Tries removing each
    pair in turn and keeps the removal whenever the subset stays
    infeasible (one LP per trial, same EPS rule). Returns the walks of
    the surviving subset as a tuple; they are contradictory together
    with d in D^k, and the count is not bounded by k.
    """
    keep = list(cycles)                        # L1
    for (b, w) in cycles:                      # L2
        trial = [p for p in keep if p != (b, w)]   # L3  keep without (b, w)
        feasible, _ = _solve(trial, k)         # L4  the same LP on the subset
        if not feasible:
            keep = trial
    return tuple(w for (b, w) in keep)         # L5