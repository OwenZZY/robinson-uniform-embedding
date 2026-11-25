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
    feasible: bool
    positions: dict | None      # label -> float (Pi); None iff infeasible
    thresholds: tuple | None    # (d_1, ..., d_k) floats; None iff infeasible
    cycles: tuple | None        # witness cycles, tuples of labels with
    #                             first == last; None iff feasible


def uniform_embedding(A, labels=None, thresholds=None) -> EmbeddingResult:
    """Uniform embedding of a Robinson similarity matrix (SPEC.md §1).

    A: square nested sequence or numpy array of ints, entries {0,...,k},
    diagonal all k (the maximum). labels: n unique hashables in matrix
    order (default 0..n-1). thresholds: optional candidate d (length k,
    strictly decreasing, positive); omitted -> k=1 trivial, k=2
    combinatorial, k>2 LP search (requires the optional SciPy extra).
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