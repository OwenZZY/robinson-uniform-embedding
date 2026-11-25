#!/usr/bin/env python3
"""Recorded acceptance cases for robinson_unif_emb (round 12).

Run from the repository root:

    python tests/acceptance.py

Requires Python >= 3.10. numpy is used only if installed. The k>2 automatic
cases run only if SciPy is importable; otherwise they are reported SKIP.
No third-party package is needed for the rest.

Groups G1-G4 encode the round-13 hardening target (see REPORT.md T12-T15 and
the review notes). They FAILED on the round-11 tree; as of round 13 all
groups pass, with and without SciPy. The recorded expectations are unchanged.

This file is deliberately self-contained: it re-implements the edge-bound
algebra and the embedding check independently of the package, so a broken
core cannot mask itself.
"""

import importlib.util
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from robinson_unif_emb import uniform_embedding          # noqa: E402
from robinson_unif_emb import _core                      # noqa: E402  (table + spy tests)

HAVE_NUMPY = importlib.util.find_spec("numpy") is not None
HAVE_SCIPY = importlib.util.find_spec("scipy") is not None


# ---------------------------------------------------------------- fixtures ---

A = ((2, 2, 1, 0, 0),
     (2, 2, 2, 1, 1),
     (1, 2, 2, 2, 1),
     (0, 1, 2, 2, 2),
     (0, 1, 1, 2, 2))

B = ((2, 2, 1, 0, 0, 0),
     (2, 2, 2, 1, 1, 1),
     (1, 2, 2, 2, 1, 1),
     (0, 1, 2, 2, 2, 1),
     (0, 1, 1, 2, 2, 2),
     (0, 1, 1, 1, 2, 2))

A3 = tuple(tuple(3 if i == j else A[i][j] for j in range(5)) for i in range(5))
B3 = tuple(tuple(3 if i == j else B[i][j] for j in range(6)) for i in range(6))

K1 = ((1, 1, 0, 0),
      (1, 1, 1, 0),
      (0, 1, 1, 1),
      (0, 0, 1, 1))

# k=3, infeasible only jointly: cycles (0,1,2,0) and (2,4,3,2) are each
# satisfiable, the pair is not (round-11 recorded).
JOINT = ((3, 3, 0, 0, 0),
         (3, 3, 2, 0, 0),
         (0, 2, 3, 1, 1),
         (0, 0, 1, 3, 2),
         (0, 0, 1, 2, 3))


# ------------------------------------------- independent algebra (test-only) ---

def _unit(k, j):
    b = [0] * k
    b[j] = 1
    return tuple(b)


def _beta_plus(a, k, u, v):
    if u < v:
        t = a[u][v]
        return _unit(k, t - 1) if t >= 1 else None
    b = _beta_minus(a, k, v, u)
    return tuple(-x for x in b)


def _beta_minus(a, k, u, v):
    if u < v:
        t = a[u][v]
        return _unit(k, t) if t <= k - 1 else (0,) * k
    b = _beta_plus(a, k, v, u)
    return None if b is None else tuple(-x for x in b)


def _walk_bound(a, k, w):
    b = (0,) * k
    for i in range(1, len(w)):
        step = _beta_plus(a, k, w[i - 1], w[i])
        assert step is not None, f"witness {w} uses a non-edge up-step"
        b = tuple(x + y for x, y in zip(b, step))
    return b


def _simple(w):
    if len(w) < 2:
        return False
    if w[0] == w[-1]:
        return len(w) >= 3 and len(set(w[:-1])) == len(w) - 1
    return len(set(w)) == len(w)


def _is_ubw(a, w):
    return all(not (w[i - 1] < w[i] and a[w[i - 1]][w[i]] == 0)
               for i in range(1, len(w)))


def _check_embedding(a, d, pi, k):
    for u in range(len(a)):
        for v in range(u + 1, len(a)):
            t = a[u][v]
            lo = d[t] if t <= k - 1 else 0
            hi = d[t - 1] if t >= 1 else None
            x = pi[v] - pi[u]
            assert lo < x, f"lower bound violated at ({u},{v}) t={t}: {lo} < {x}"
            assert hi is None or x < hi, f"upper bound violated at ({u},{v}) t={t}: {x} < {hi}"


def _ratio(b):
    a1, a2 = b
    if a2 > 0:
        return "lo", Fraction(-a1, a2)
    if a2 < 0:
        return "hi", Fraction(-a1, a2)
    return "zero", Fraction(a1)


def _k2_witness_contradiction(a, cycles):
    """Each returned cycle must be a simple UBW; together they must contradict
    r in (0,1): one individually dead, or disjoint ratio intervals."""
    assert cycles, "infeasible result without witness cycles"
    lo, hi = [], []
    dead = False
    for w in cycles:
        assert _simple(w) and w[0] == w[-1], f"witness {w} is not a simple cycle"
        assert _is_ubw(a, w), f"witness {w} is not an upper-bound walk"
        kind, r = _ratio(_walk_bound(a, 2, w))
        if kind == "lo":
            dead |= r >= 1
            lo.append(r)
        elif kind == "hi":
            dead |= r <= 0
            hi.append(r)
        else:
            dead |= r <= 0
    assert dead or (lo and hi and max(lo) >= min(hi)), \
        f"witnesses {cycles} do not contradict r in (0,1)"


def _frac_result(res, a, k):
    d = tuple(Fraction(x) for x in res.thresholds)
    pi = [Fraction(v) for _, v in sorted(res.positions.items())]
    _check_embedding(a, d, pi, k)


def _expect_value_error(fn, what):
    try:
        fn()
    except ValueError:
        return
    except Exception as exc:                       # noqa: BLE001 - report the leak
        raise AssertionError(f"{what}: expected ValueError, got {type(exc).__name__}: {exc}")
    raise AssertionError(f"{what}: expected ValueError, got a result")


# ------------------------------------------------------------- recorded data ---

# k=2 feasible: (matrix, lo, hi) with the exact ratio window of the full
# simple-cycle system (exhaustive enumeration, round 11). The library must
# return d = (1, (lo+hi)/2) and a valid embedding.
K2_FEASIBLE = [
    (((2, 0, 0), (0, 2, 2), (0, 2, 2)), "0", "1"),
    (((2, 1, 1, 0, 0), (1, 2, 1, 0, 0), (1, 1, 2, 0, 0),
      (0, 0, 0, 2, 0), (0, 0, 0, 0, 2)), "0", "1/2"),
    (((2, 0, 0, 0, 0), (0, 2, 1, 1, 1), (0, 1, 2, 1, 1),
      (0, 1, 1, 2, 1), (0, 1, 1, 1, 2)), "0", "1/3"),
    (((2, 2, 0), (2, 2, 2), (0, 2, 2)), "1/2", "1"),
    (((2, 2, 2, 1, 0, 0), (2, 2, 2, 1, 0, 0), (2, 2, 2, 2, 1, 0),
      (1, 1, 2, 2, 2, 1), (0, 0, 1, 2, 2, 1), (0, 0, 0, 1, 1, 2)), "1/3", "1"),
]

# k=2 infeasible: the witness pair/single cycle must contradict r in (0,1).
K2_INFEASIBLE = [
    ((2, 0, 0, 0, 0, 0), (0, 2, 0, 0, 0, 0), (0, 0, 2, 2, 0, 0),
     (0, 0, 2, 2, 1, 1), (0, 0, 0, 1, 2, 1), (0, 0, 0, 1, 1, 2)),
    ((2, 2, 1, 0, 0), (2, 2, 2, 0, 0), (1, 2, 2, 1, 1),
     (0, 0, 1, 2, 1), (0, 0, 1, 1, 2)),
    ((2, 1, 1, 0), (1, 2, 1, 0), (1, 1, 2, 2), (0, 0, 2, 2)),
]

K1_CASES = [
    ((1, 1, 1), (1, 1, 1), (1, 1, 1)),
    ((1, 0), (0, 1)),
    ((1, 0, 0), (0, 1, 1), (0, 1, 1)),
]

# given-d k=3: (matrix, d, feasible)
GIVEN_D_K3 = [
    (A3, (8, 6, Fraction(1, 2)), True),
    (A3, (8, 6, 5), False),
    (B3, (8, 6, Fraction(1, 2)), False),
]

# validation cases that must raise ValueError on any tree
INVALID = [
    ("non-square", lambda: uniform_embedding(((1, 2), (3,)))),
    ("ragged", lambda: uniform_embedding(((1, 2), (2,)))),
    ("float entries", lambda: uniform_embedding(((2.0, 1.0), (1.0, 2.0)))),
    ("asymmetric", lambda: uniform_embedding(((2, 1), (0, 2)))),
    ("diag not max", lambda: uniform_embedding(((1, 2), (2, 1)))),
    ("not Robinson", lambda: uniform_embedding(((2, 1, 2), (1, 2, 1), (2, 1, 2)))),
    ("labels length", lambda: uniform_embedding(A, labels=["a", "b"])),
    ("labels duplicate", lambda: uniform_embedding(A, labels=["a", "a", "b", "c", "d"])),
    ("thresholds length", lambda: uniform_embedding(A, thresholds=(8,))),
    ("thresholds order", lambda: uniform_embedding(A, thresholds=(6, 8))),
]


# ------------------------------------------------------------------- runner ---

RESULTS = []


class Skip(Exception):
    """Raised by a test that cannot run in this environment."""


def run(group, name, fn):
    try:
        fn()
    except Skip as exc:
        RESULTS.append((group, name, None, str(exc)))
    except AssertionError as exc:
        RESULTS.append((group, name, False, str(exc)))
    except Exception as exc:                        # noqa: BLE001
        RESULTS.append((group, name, False, f"{type(exc).__name__}: {exc}"))
    else:
        RESULTS.append((group, name, True, ""))


# -- fixtures -----------------------------------------------------------------

def test_A():
    res = uniform_embedding(A)
    assert res.feasible and res.cycles is None
    assert res.thresholds == (1.0, 0.6666666666666666)
    _frac_result(res, A, 2)


def test_B():
    res = uniform_embedding(B)
    assert not res.feasible and res.positions is None and res.thresholds is None
    _k2_witness_contradiction(B, res.cycles)


def test_K1():
    res = uniform_embedding(K1)
    assert res.feasible and res.thresholds == (1.0,)
    _frac_result(res, K1, 1)


def test_A3_given_d():
    res = uniform_embedding(A3, thresholds=(8, 6, Fraction(1, 2)))
    assert res.feasible and res.cycles is None
    _frac_result(res, A3, 3)


def test_B3_given_d():
    res = uniform_embedding(B3, thresholds=(8, 6, Fraction(1, 2)))
    assert not res.feasible
    assert len(res.cycles) == 1
    w = res.cycles[0]
    assert _simple(w) and _is_ubw(B3, w)
    assert all(x == 0 for x in _walk_bound(B3, 3, w))


def test_A3_auto():
    if not HAVE_SCIPY:
        raise Skip("SciPy absent (install the lp extra)")
    res = uniform_embedding(A3)
    assert res.feasible and res.cycles is None
    assert res.thresholds == (1.0, 0.5, 0.25)
    _frac_result(res, A3, 3)


def test_B3_auto():
    if not HAVE_SCIPY:
        raise Skip("SciPy absent (install the lp extra)")
    res = uniform_embedding(B3)
    assert not res.feasible and res.cycles
    for w in res.cycles:
        assert _simple(w) and _is_ubw(B3, w)
    from robinson_unif_emb import _lp
    cyc = [(_walk_bound(B3, 3, w), w) for w in res.cycles]
    feasible, _ = _lp._solve(cyc, 3)
    assert not feasible, "B3 witness set re-solves feasible"


def test_JOINT_auto():
    if not HAVE_SCIPY:
        raise Skip("SciPy absent")
    res = uniform_embedding(JOINT)
    assert not res.feasible and len(res.cycles) >= 2, \
        "a joint contradiction should need at least two cycles"
    from robinson_unif_emb import _lp
    cyc = [(_walk_bound(JOINT, 3, w), w) for w in res.cycles]
    feasible, _ = _lp._solve(cyc, 3)
    assert not feasible, "returned witness set is not contradictory"
    for c in cyc:                                   # each one alone is satisfiable
        single, _ = _lp._solve([c], 3)
        assert single, "the witness set is not minimal (a single cycle suffices)"


# -- recorded k=2 / k=1 / given-d --------------------------------------------

def test_k2_feasible():
    for a, lo_s, hi_s in K2_FEASIBLE:
        lo, hi = Fraction(lo_s), Fraction(hi_s)
        res = uniform_embedding(a)
        assert res.feasible, f"expected feasible for {a}"
        d = tuple(Fraction(x) for x in res.thresholds)
        expected_r = float((lo + hi) / 2)          # the core computes the exact
        assert d[0] == 1 and d[1] == Fraction(expected_r), \
            f"expected r=float((lo+hi)/2)={expected_r}, got {d[1]}"
        assert lo < d[1] < hi, f"{d[1]} outside the oracle window ({lo},{hi})"
        pi = [Fraction(v) for _, v in sorted(res.positions.items())]
        _check_embedding(a, d, pi, 2)


def test_k2_infeasible():
    for a in K2_INFEASIBLE:
        res = uniform_embedding(a)
        assert not res.feasible, f"expected infeasible for {a}"
        _k2_witness_contradiction(a, res.cycles)


def test_k1_cases():
    for a in K1_CASES:
        res = uniform_embedding(a)
        assert res.feasible and res.thresholds == (1.0,)
        _frac_result(res, a, 1)


def test_given_d_k3():
    for a, d_in, feasible in GIVEN_D_K3:
        res = uniform_embedding(a, thresholds=d_in)
        assert res.feasible == feasible, f"{a} d={d_in}: expected feasible={feasible}"
        if feasible:
            _frac_result(res, a, 3)
        else:
            d = tuple(Fraction(x) for x in d_in)
            assert res.cycles, "infeasible without a witness"
            for w in res.cycles:
                assert _simple(w) and _is_ubw(a, w)
                b = _walk_bound(a, 3, w)
                assert sum(x * y for x, y in zip(b, d)) <= 0, \
                    f"witness {w} does not violate the given d"


# -- edge cases ---------------------------------------------------------------

def test_tiny():
    for a, k, d_expected in ((((1,),), 1, (1.0,)),
                             (((2,),), 2, (1.0, 0.5)),
                             (((2, 0), (0, 2)), 2, (1.0, 0.5))):
        res = uniform_embedding(a)
        assert res.feasible and res.thresholds == d_expected, (a, res.thresholds)
        _frac_result(res, a, k)


def test_labels():
    names = list("abcdef")
    res = uniform_embedding(B, labels=names)
    assert not res.feasible and res.cycles
    assert all(isinstance(x, str) for w in res.cycles for x in w)
    idx = {c: i for i, c in enumerate(names)}
    _k2_witness_contradiction(B, [tuple(idx[x] for x in w) for w in res.cycles])
    res = uniform_embedding(A, labels=list("pqrst"))
    assert res.feasible and sorted(res.positions) == list("pqrst")


def test_numpy():
    if not HAVE_NUMPY:
        raise Skip("numpy absent")
    import numpy as np
    res = uniform_embedding(np.array(A, dtype=int))
    assert res.feasible
    res = uniform_embedding(np.array(A, dtype=np.int64),
                            labels=np.array([10, 11, 12, 13, 14]))
    assert res.feasible and sorted(int(k) for k in res.positions) == [10, 11, 12, 13, 14]
    _expect_value_error(lambda: uniform_embedding(np.array(A, dtype=float)),
                        "numpy float array")


# -- validation / G1 / G2 -----------------------------------------------------

def test_invalid_inputs():
    for name, fn in INVALID:
        _expect_value_error(fn, name)


def test_G1_zero_k():
    _expect_value_error(lambda: uniform_embedding(((0, 0), (0, 0))), "k=0 all-zero matrix")


def test_G2_sequences_and_finite():
    _expect_value_error(lambda: uniform_embedding(A, labels=(i for i in range(5))),
                        "labels generator")
    _expect_value_error(lambda: uniform_embedding(A, thresholds=(x for x in (8, 6))),
                        "thresholds generator")
    _expect_value_error(lambda: uniform_embedding(A, thresholds=(float("inf"), 1)),
                        "infinite threshold")


# -- G3 / G4 (internal target checks) ----------------------------------------

def test_G3_no_duplicate_bounds():
    for a in (A3, B3):
        T = _core.bound_generation(a)
        n = len(a)
        for i in range(n):
            for j in range(n):
                pairs = list(_core.cell_items(T[i][j]))
                bounds = [b for _, b in pairs]
                assert len(bounds) == len(set(bounds)), \
                    f"T[{i}][{j}] has duplicate bounds ({len(bounds)} pairs, {len(set(bounds))} distinct)"


def test_G4_exact_verification():
    if not HAVE_SCIPY:
        raise Skip("SciPy absent")
    seen = []
    orig = _core.verify_embedding

    def spy(a, d, pi):
        seen.append((type(d[0]), type(pi[0])))
        return orig(a, d, pi)

    _core.verify_embedding = spy
    try:
        uniform_embedding(A3)
    finally:
        _core.verify_embedding = orig
    assert seen, "verify_embedding was not called"
    assert all(td is Fraction and tp is Fraction for td, tp in seen), \
        f"k>2 verification ran on {seen} (expected Fractions)"


def main():
    run("fixtures", "A feasible + exact embedding", test_A)
    run("fixtures", "B infeasible + (0,0) witness", test_B)
    run("fixtures", "K1 feasible", test_K1)
    run("fixtures", "A3 given d=(8,6,1/2)", test_A3_given_d)
    run("fixtures", "B3 given d=(8,6,1/2)", test_B3_given_d)
    run("fixtures", "A3 automatic (k=3)", test_A3_auto)
    run("fixtures", "B3 automatic (k=3)", test_B3_auto)
    run("fixtures", "joint-witness k=3 (deletion filter)", test_JOINT_auto)
    run("recorded k=2", "feasible windows / exact d / pi", test_k2_feasible)
    run("recorded k=2", "infeasible witnesses", test_k2_infeasible)
    run("recorded k=1", "feasible embeddings", test_k1_cases)
    run("recorded k=3", "given-d feasible/infeasible", test_given_d_k3)
    run("edge", "n=1 / n=2 null", test_tiny)
    run("edge", "labels mapping", test_labels)
    run("edge", "numpy input", test_numpy)
    run("validation", "nine invalid inputs -> ValueError", test_invalid_inputs)
    run("G1", "k=0 all-zero matrix -> ValueError", test_G1_zero_k)
    run("G2", "generators / inf -> ValueError", test_G2_sequences_and_finite)
    run("G3", "no duplicate bounds in k!=2 cells", test_G3_no_duplicate_bounds)
    run("G4", "k>2 verification on Fractions", test_G4_exact_verification)

    width = max(len(name) for _, name, _, _ in RESULTS)
    failures, skipped = [], []
    group = None
    for grp, name, ok, msg in RESULTS:
        if grp != group:
            print(f"\n[{grp}]")
            group = grp
        if ok is None:
            skipped.append((grp, name))
            print(f"  SKIP  {name:<{width}}  {msg}")
        elif ok:
            print(f"  PASS  {name}")
        else:
            failures.append((grp, name, msg))
            print(f"  FAIL  {name:<{width}}  {msg}")

    print("\n" + "-" * 72)
    print(f"{sum(1 for *_, ok, _ in RESULTS if ok is True)} passed, "
          f"{len(failures)} failed, {len(skipped)} skipped")
    if failures:
        print("\nFailures:")
        for grp, name, msg in failures:
            print(f"  [{grp}] {name}: {msg}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
