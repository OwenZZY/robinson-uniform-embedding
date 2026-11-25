# SPEC — engineering (software design)

Companion to [`SPEC.md`](SPEC.md), the standalone core spec (definitions, data types, numbered pseudocode). This document covers everything that is **not** the algorithm: layers, the public interface, input handling, the LP-backend mechanics, numerics policy, naming, packaging, errors, and the implementation plan. It is self-contained; no paper is needed. A term glossary with source provenance is in [`docs/definitions.md`](docs/definitions.md) (optional reference).

## 1. Layers and modules

```
outer (user / Jupyter)
   │  nested lists, numpy arrays, optional labels
   ▼
translation  ──  robinson_unif_emb/_translate.py
   │  validation predicates (boolean + immediate halt), format conversion,
   │  explicit bijection labels ↔ 0..n-1, Fraction → float on output
   ▼
core  ──  robinson_unif_emb/_core.py
   │  the math pipeline (SPEC.md); pure Python stdlib; trusts its input;
   │  zero try/except/raise; no validation; no third-party imports; any k
   ▼
LP backend  ──  robinson_unif_emb/_lp.py   (optional; only the k>2 search)
      lazy import of scipy.optimize.linprog; consumes core's cycle bounds;
      may raise ImportError if SciPy is absent
```

- **`_core.py`** — the algorithm only. Assumes a suitable internal input (a Robinson matrix, and a feasible `d` where required). No defensive checks, no exceptions, no external formats, no labels, no third-party imports. Comments follow the numbered pseudocode in `SPEC.md`.
- **`_lp.py`** — the only module importing SciPy, lazily inside the search function so the package imports cleanly without it. Used only by the $k>2$ automatic search.
- **`_translate.py`** — the only layer touching external formats: conversion, validation predicates, the label bijection, output conversion.
- **`__init__.py`** — exports the single public function and the result type; wires the layers and dispatches on `k` / `thresholds` (the branch in `SPEC.md` §5).
- **Intermediate layer** — not needed; matrices/vectors translate directly. A NetworkX front end, if ever wanted, is a separate adapter.

## 2. Public interface (single function)

```python
__all__ = ["uniform_embedding", "EmbeddingResult"]

def uniform_embedding(A, labels=None, thresholds=None) -> EmbeddingResult: ...
```

| Parameter | Type | Meaning |
|---|---|---|
| `A` | square nested sequence or numpy array of ints | Robinson similarity matrix, entries $\{0,\dots,k\}$, diagonal all $k$ (the maximum). |
| `labels` | sequence of $n$ unique hashables, optional | External identifiers in matrix order; default `tuple(range(n))`. Bijection `labels[i] ↔ i`. |
| `thresholds` | sequence of $k$ numbers, optional | A candidate $\mathbf d\in\mathbb D^k$. Given → verify + construct (any $k$). Omitted → $k=2$ combinatorial search, $k\ne2$ LP search. |

**Dispatch** (`SPEC.md` §5): `thresholds` given → §5a; else `k=2` → §5b; else → §5c.

**Result** — one frozen dataclass, always returned (infeasibility is an outcome, not an error):

```python
@dataclass(frozen=True)
class EmbeddingResult:
    feasible: bool
    positions: dict | None        # label → float (Π); None iff infeasible
    thresholds: tuple | None      # (d_1,...,d_k) floats; None iff infeasible
    cycles: tuple | None          # infeasible: 1–2 (k=2) or ≤k (k>2) witness cycles,
                                  # each a tuple of labels with first == last; None iff feasible
```

Usage sketch:

```python
from robinson_unif_emb import uniform_embedding
res = uniform_embedding([[2, 2, 1], [2, 2, 2], [1, 2, 2]], labels=["a", "b", "c"])
res.feasible, res.positions, res.thresholds
```

## 3. Translation layer

**Input conversion.** numpy → `tolist()` (duck-typed via `hasattr`); nested sequences → tuple of tuples of `int`. Core indexes `a[u][v]`.

**Validation predicates** (boolean, evaluated in order; first `False` → `ValueError` with a plain message). All here, none in core:

1. square, $n\ge1$;
2. entries integral (numpy integer types accepted) in $\{0,\dots,k\}$, $k=A[0][0]$;
3. symmetric;
4. diagonal constant and equal to the maximum $k$;
5. Robinson: $a_{u,v}\ge a_{u,w}$ and $a_{v,w}\ge a_{u,w}$ for all $u<v<w$;
6. `labels` (if given): length $n$, hashable, unique;
7. `thresholds` (if given): length $k$, strictly decreasing, positive, real.

Any $k\ge1$ is accepted; mode selection needs no input predicate. SciPy availability is checked in `_lp.py`, not here.

**Bijection.** Build `index_of = {label: i}`; keep `labels` for the inverse. Core sees indices only; `positions` and witness `cycles` are mapped back through `labels` before return.

**Output conversion.** `Fraction` → `float` by default (see §5).

## 4. LP backend (`_lp.py`)

```python
def search_thresholds(cycles, k):     # cycles: list[(bound, walk)] from core.cycle_bounds
    # → (d, None), d a feasible tuple with d_1 = 1; or (None, witnesses), ≤ k walks
```

- **Import.** `from scipy.optimize import linprog` inside the function; absent → `ImportError("install scipy for k>2 automatic search: pip install scipy")`.
- **Formulation.** The LP itself (variables, constraints, `d[0]=1` normalization, max-slack `t`) is specified in `SPEC.md` §5c; this module only supplies the SciPy call.
- **Decision.** Strict-feasible iff optimum `t* > 0`; return `d`. If infeasible or `t* ≤ 0`: recover ≤ `k` witness cycles from the Farkas certificate / dual multipliers; if the backend exposes no certificate, fall back to the cycle constraints active at the optimum.
- **Dependency.** SciPy is optional (extra `lp`); the `k=2` and given-`d` modes never import it.

## 5. Numerics

- `k=2` pipeline: `fractions.Fraction` end to end (ratio bounds `-a1/a2` are rational → every strict inequality decided exactly). Convert to `float` only at the public boundary.
- Given-thresholds mode: inherits the user's numeric type (`Fraction` → exact; `float` → float).
- LP search: returns floats (SciPy); the construction then runs in floats.
- An `exact=True` output option (keep `Fraction`s) is a possible later addition; not in v1.

## 6. Internal data representations

The core data types (`a`, `n`, `k`, walk, bound, cell, `T`, `d`, `pi`) are specified in `SPEC.md` §2 — that is the algorithm's contract. Translation's job is to build `a` (tuple of tuples of `int`) from the external input and to map `pi`/witnesses back to labels (§3). The only public-facing types are `EmbeddingResult` and the `labels`/`thresholds` parameters (§2). Lower bounds are never stored; read by reversal from `T[v][i]` (`SPEC.md` §6).

## 7. Naming

- Package `robinson_unif_emb`; distribution/GitHub `robinson-uniform-embedding`.
- Public surface exactly `uniform_embedding` + `EmbeddingResult` (`__all__`).
- Private modules `_core`, `_translate`, `_lp`; internal names are short and mathematical (`pi`, `d`, `ub`, `lb`, `beta_plus`, `preceq`, `T`), matching the notation in `SPEC.md`.

| Notation | Code |
|---|---|
| $\Pi$ | `pi` (core), `positions` (public) |
| $\mathbf d$ | `d` (core), `thresholds` (public) |
| $\beta^\pm(u,v)$ | `beta_plus`, `beta_minus` |
| $\beta^+(W)$ | `walk_bound` |
| $\preceq$ | `preceq` |
| UBW table | `T`, cells `T[i][j]` |
| $\mathrm{ub}_v, \mathrm{lb}_v$ | `ub`, `lb` (locals) |
| bound generation | `bound_generation` |
| $d_2/d_1$ window | `ratio_bounds` → `(lo, hi)` |
| cycle inequalities | `cycle_bounds` (core) → `search_thresholds` (`_lp`) |
| witness cycles | `cycles` (public), `(w1, w2)` (core) |

## 8. Error model

- `ValueError` — translation, on the first failed input predicate (§3).
- `ImportError` — `_lp.py`, if the $k>2$ search is requested without SciPy.
- **Infeasibility is a returned `EmbeddingResult(feasible=False, cycles=...)`, never an exception.**
- The core never raises (zero try/except/raise).

## 9. Dependencies and packaging

- Core: Python stdlib only (`math`, `fractions`).
- Optional: SciPy (extra `lp`) for the $k>2$ search; numpy accepted on input but not required (duck-typed).
- `pyproject.toml` with `[project.optional-dependencies] lp = ["scipy"]`; package imports and modes 2a/2b run with no third-party packages installed.

## 10. Implementation roadmap

API is fixed by `SPEC.md` + this document; modules have disjoint write scopes, so the code can be built by **parallel agents**:

1. `_core.py` primitives + `bound_generation` + `cycle_bounds` (`SPEC.md` §3–§5). Regression: on fixture B below, the diagonal of `T` must contain the bound `(0,0)` — this catches the down-direction under-generation bug.
2. `_core.py` `ratio_bounds`, `embedding_from_bounds`, `solve_k2`, `solve_with_thresholds` (`SPEC.md` §5b, §6). Fixtures: A feasible; B infeasible with witness `(0,1,5,3,0)` or equivalent.
3. `_lp.py` `search_thresholds` + `solve_general` (`SPEC.md` §5c). Verify a `k=3` feasible and an infeasible fixture; confirm the `k=2` and given-`d` modes run with SciPy absent.
4. `_translate.py` predicates + bijection + conversions; `__init__.py` public function + result type + dispatch.
5. Offer (do not build unasked): test suite with the fixtures, README, plotting helper.

**Fixtures** (0-based):
- **A** (5×5, feasible): rows `(2,2,1,0,0) (2,2,2,1,1) (1,2,2,2,1) (0,1,2,2,2) (0,1,1,2,2)`. A valid embedding: `pi=(0,5,6.5,11.75,12.75)`, `d=(8,6)`.
- **B** (6×6, infeasible): rows `(2,2,1,0,0,0) (2,2,2,1,1,1) (1,2,2,2,1,1) (0,1,2,2,2,1) (0,1,1,2,2,2) (0,1,1,1,2,2)`. Witness cycle `(0,1,5,3,0)` has bound `(0,0)`.

## 11. Out of scope for v1

Seriation (recovering order from an unordered matrix); graphon inputs; NetworkX adapters; a pure-Python LP fallback; performance engineering beyond the chain-slot trick; proofs and complexity analysis in comments. A combinatorial feasibility test for `k≥3` is an open problem (none is known); not an implementation target.