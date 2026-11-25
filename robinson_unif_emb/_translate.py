"""Translation layer (SPEC-engineering.md §3): the only module touching
external formats. Converts and validates inputs, owns the bijection
labels[i] <-> i, and maps core outputs back to external identifiers.

Validation predicates are boolean checks evaluated in order (1..7);
the first failure raises ValueError immediately. The core never sees
external formats. numpy is duck-typed via tolist(); no third-party imports.
"""

from collections.abc import Hashable
from fractions import Fraction
from numbers import Real


# ---- validation predicates (SPEC-engineering.md §3, in order) ----
# Each returns True iff the condition holds; to_internal raises
# ValueError on the first False.

def _square(A):
    """Predicate 1: A is a nonempty square nested sequence (n >= 1).
    Non-sequences give False, not TypeError.
    """
    if not hasattr(A, "__len__") or len(A) < 1:
        return False
    n = len(A)
    return all(hasattr(row, "__len__") and len(row) == n for row in A)


def _integral(A):
    """Predicate 2a: every entry is integer-typed (int or a numpy
    integer, via __index__); floats such as 2.0 are rejected.
    """
    return all(hasattr(x, "__index__") for row in A for x in row)


def _in_range(A, k):
    """Predicate 2b: every entry lies in {0, ..., k}, k = A[0][0]."""
    return all(0 <= x <= k for row in A for x in row)


def _symmetric(A):
    """Predicate 3: A[u][v] == A[v][u] for all u, v."""
    n = len(A)
    return all(A[u][v] == A[v][u] for u in range(n) for v in range(u + 1, n))


def _diagonal_max(A, k):
    """Predicate 4: the diagonal is constant and equal to the maximum
    entry k.
    """
    n = len(A)
    return (all(A[i][i] == k for i in range(n))
            and k == max(x for row in A for x in row))


def _robinson(A):
    """Predicate 5: A[u][v] >= A[u][w] and A[v][w] >= A[u][w] for all
    u < v < w (entries do not increase away from the diagonal).
    """
    n = len(A)
    return all(A[u][v] >= A[u][w] and A[v][w] >= A[u][w]
               for u in range(n) for v in range(u + 1, n)
               for w in range(v + 1, n))


def _labels_ok(labels, n):
    """Predicate 6 (labels given): exactly n hashable, pairwise distinct
    identifiers, so labels[i] <-> i is a bijection.
    """
    if len(labels) != n or not all(isinstance(x, Hashable) for x in labels):
        return False
    try:
        return len(set(labels)) == n
    except TypeError:      # pathological hashables
        return False


def _thresholds_ok(thresholds, k):
    """Predicate 7 (thresholds given): k real numbers, strictly
    decreasing and positive — a candidate d in D^k.
    """
    if len(thresholds) != k or not all(isinstance(x, Real) for x in thresholds):
        return False
    return (all(x > y for x, y in zip(thresholds, thresholds[1:]))
            and all(x > 0 for x in thresholds))


# ---- conversion ----

def to_internal(A, labels, thresholds):
    """Validate and convert external input to the core's data types.

    A: square nested sequence or numpy array of ints (duck-typed via
    tolist). labels: optional sequence of n unique hashables. thresholds:
    optional candidate d, k strictly decreasing positive reals.
    Returns (a, d0, labels_tuple): a is a tuple of tuples of int;
    d0 a tuple of Fraction (floats by exact binary conversion) or None;
    labels_tuple the external identifiers in matrix order (default
    tuple(range(n))). Raises ValueError on the first failed predicate.
    """
    if hasattr(A, "tolist"):                   # numpy duck-typing
        A = A.tolist()
    if not _square(A):
        raise ValueError("A must be a square nested sequence of size n >= 1")
    if not _integral(A):
        raise ValueError("entries of A must be integer-typed (floats such as 2.0 are rejected)")
    k = A[0][0]
    if not _in_range(A, k):
        raise ValueError("entries of A must lie in {0, ..., k} with k = A[0][0]")
    if not _symmetric(A):
        raise ValueError("A must be symmetric")
    if not _diagonal_max(A, k):
        raise ValueError("the diagonal of A must be constant and equal to the maximum entry k")
    if not _robinson(A):
        raise ValueError("A must be Robinson: a[u][v] >= a[u][w] and a[v][w] >= a[u][w] for all u < v < w")
    n = len(A)
    if labels is not None and not _labels_ok(labels, n):
        raise ValueError("labels must be n unique hashables")
    if thresholds is not None and not _thresholds_ok(thresholds, k):
        raise ValueError("thresholds must be k real numbers, strictly decreasing and positive")
    a = tuple(tuple(int(x) for x in row) for row in A)
    # Fraction(float) is exact binary conversion; no rounding, no stringify
    d0 = tuple(Fraction(x) for x in thresholds) if thresholds is not None else None
    labels_tuple = tuple(labels) if labels is not None else tuple(range(n))
    return a, d0, labels_tuple


def map_positions(pi, labels):
    """Embedding back to external identifiers: {labels[i]: float(pi[i])}.
    pi is the core's length-n list; labels the bijection in matrix order.
    """
    return {labels[i]: float(pi[i]) for i in range(len(pi))}


def map_cycles(witnesses, labels):
    """Witness walks back to external identifiers: a tuple of label
    tuples, first == last preserved. witnesses is the core's tuple of
    0-based vertex walks.
    """
    return tuple(tuple(labels[i] for i in w) for w in witnesses)