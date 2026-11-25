"""Public interface (SPEC-engineering.md §2-§3): one function, one result type.

Wiring: translation -> core -> verification -> result. Infeasibility is a
returned EmbeddingResult, never an exception. _lp is imported lazily, only
when the k>2 automatic threshold search is needed.
"""

from dataclasses import dataclass

from . import _core, _translate

__all__ = ["uniform_embedding", "EmbeddingResult"]


@dataclass(frozen=True)
class EmbeddingResult:
    """Outcome of uniform_embedding; infeasibility is a result, not an error.

    feasible: True iff a uniform embedding was found and verified.
    positions: {label: float} embedding Pi, or None iff infeasible.
    thresholds: the threshold vector d as floats (d[0] = 1 on the
        automatic paths; a supplied d is returned at its own scale),
        or None iff infeasible.
    cycles: simple witness cycles as tuples of labels (first == last),
        or None iff feasible. k=2: one individually unsatisfiable cycle
        or two with disjoint ratio intervals; k>2: a deletion-filter
        set, jointly contradictory with d in D^k (no count bound).
    """

    feasible: bool
    positions: dict | None      # label -> float (Pi); None iff infeasible
    thresholds: tuple | None    # (d_1, ..., d_k) floats; None iff infeasible
    cycles: tuple | None        # witness cycles, tuples of labels with
    #                             first == last; None iff feasible


def uniform_embedding(A, labels=None, thresholds=None) -> EmbeddingResult:
    """Uniform embedding of a Robinson similarity matrix (SPEC.md §1).

    A: square nested sequence or numpy array of ints, entries {0,...,k},
        diagonal all k (the maximum), Robinson ordering.
    labels: optional n unique hashables in matrix order (default
        0..n-1); all output identifiers are mapped back through them.
    thresholds: optional candidate d (length k, strictly decreasing,
        positive). Given: verify it and construct pi (any k, exact).
        Omitted: k=1 trivial, k=2 combinatorial (both exact, no
        dependencies), k>2 LP search (requires the optional SciPy extra).

    Returns an EmbeddingResult: on success positions/thresholds are
    floats and cycles is None; on infeasibility positions/thresholds
    are None and cycles holds the witness cycles. Every successful
    embedding is verified before it is reported feasible.
    Raises ValueError on invalid input, ImportError if the k>2 search
    is requested without SciPy, RuntimeError if verification fails.
    """
    a, d0, labs = _translate.to_internal(A, labels, thresholds)
    k = a[0][0]
    search = None
    if d0 is None and k > 2:
        from . import _lp                      # lazy: only the k>2 search
        search = _lp.search_thresholds
    d, pi = _core.uniform_embedding_core(a, d0, search)
    if d is None:                              # pi holds witnesses here
        return EmbeddingResult(False, None, None, _translate.map_cycles(pi, labs))
    if not _core.verify_embedding(a, d, pi):   # every success is verified
        raise RuntimeError("numerical failure: the constructed embedding did not pass verification "
                           "(k>2: try a different d or report the pruned-table limitation; "
                           "k=1/k=2: report the table bug)")
    return EmbeddingResult(True, _translate.map_positions(pi, labs),
                           tuple(float(x) for x in d), None)