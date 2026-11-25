# SPEC — engineering (software design)

Companion to [`SPEC.md`](SPEC.md), the lean data/interface contract (data types, indexing, per-function input/output). The step-by-step pseudocode (P1–P11) and conventions R1–R12 are in [`docs/pseudocode.md`](docs/pseudocode.md). This document covers everything that is **not** the algorithm: layers, the public interface, input handling, the LP-backend mechanics, numerics policy, naming, packaging, errors, and the implementation plan. Glossary: [`docs/definitions.md`](docs/definitions.md).

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

- **`_core.py`** — the algorithm only. Assumes a suitable internal input (a Robinson matrix, and a feasible $\mathbf d$ where required). No defensive checks, no exceptions, no external formats, no labels, no third-party imports. Comments cite the paper line by line.
- **`_lp.py`** — the only module importing SciPy, lazily inside the search function so the package imports cleanly without it. Used only by the $k>2$ automatic search.
- **`_translate.py`** — the only layer touching external formats: conversion, validation predicates, the label bijection, output conversion.
- **`__init__.py`** — exports the single public function and the result type; wires the layers and dispatches on $k$ / `thresholds` (the Stage-2 branch of `SPEC.md`).
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

**Dispatch** (the math branch, `SPEC.md` Stage 2): `thresholds` given → 2a; else $k=2$ → 2b; else → 2c.

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
- **Formulation.** The LP itself (variables, constraints, $d_1=1$ normalization, max-slack $t$) is specified in `docs/pseudocode.md` P10 (R11); this module only supplies the SciPy call.
- **Decision.** Strict-feasible iff optimum $t^*>0$; return $\mathbf d=(d_1,\dots,d_k)$. If infeasible or $t^*\le0$: recover $\le k$ witness cycles from the Farkas certificate / dual multipliers (L558); if the backend exposes no certificate, fall back to the cycle constraints active at the optimum.
- **Dependency.** SciPy is optional (extra `lp`); modes 2a/2b never import it.

## 5. Numerics

- $k=2$ pipeline: `fractions.Fraction` end to end (ratio bounds $-a_1/a_2$ are rational → every strict inequality decided exactly, R10). Convert to `float` only at the public boundary.
- Given-thresholds mode: inherits the user's numeric type (`Fraction` → exact; `float` → float).
- LP search: returns floats (SciPy); the construction then runs in floats.
- An `exact=True` output option (keep `Fraction`s) is a possible later addition; not in v1 (report T2).

## 6. Internal data representations

The core data types (`a`, `n`, `k`, walk, bound, cell, `T`, `d`, `pi`) are specified in `SPEC.md` §2 — that is the algorithm's contract. Translation's job is to build `a` (tuple of tuples of `int`) from the external input and to map `pi`/witnesses back to labels (§3). The only public-facing types are `EmbeddingResult` and the `labels`/`thresholds` parameters (§2). Lower bounds are never stored; read by reversal from `T[v][i]` (`SPEC.md` §5, pseudocode P6).

## 7. Naming

- Package `robinson_unif_emb`; distribution/GitHub `robinson-uniform-embedding`.
- Public surface exactly `uniform_embedding` + `EmbeddingResult` (`__all__`).
- Private modules `_core`, `_translate`, `_lp`; internal names follow the paper where short (`pi`, `d`, `ub`, `lb`, `beta_plus`, `preceq`, `T`) so a mathematician can read against the source. No paper macros leak into identifiers.

| Paper | Code |
|---|---|
| $\Pi$ | `pi` (core), `positions` (public) |
| $\mathbf d$, $\mathbb D^k$ | `d` (core), `thresholds` (public) |
| $\beta^\pm(u,v)$ | `beta_plus`, `beta_minus` |
| $\beta^+(W)$ | `walk_bound` |
| $\preceq$ | `preceq` |
| UBW table | `T`, cells `T[i][j]` |
| $\mathrm{ub}_v, \mathrm{lb}_v$ | `ub`, `lb` (locals) |
| FW-var | `bound_generation` |
| $d_2/d_1$ window | `ratio_bounds` → `(lo, hi)` |
| cycle inequalities | `cycle_bounds` (core) → `search_thresholds` (`_lp`) |
| $C_1, C_2$ | `cycles` (public), `(w1, w2)` (core) |

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

API is fixed by `SPEC.md` + this document; modules have disjoint write scopes, so the code can be built by **parallel agents** (per AGENTS.md):

1. `_core.py` primitives + `bound_generation` + `cycle_bounds` (P1–P4, P9). Regression: Example 1B diagonal must contain bound `(0,0)` (F1).
2. `_core.py` `ratio_bounds`, `embedding_from_bounds`, `solve_k2`, `solve_with_thresholds` (P5–P8). Fixtures: 1A feasible, 1B infeasible with witness `⟨0,1,5,3,0⟩` (0-based) or equivalent.
3. `_lp.py` `search_thresholds` (P10) + general orchestration (P11). Verify a $k=3$ feasible and an infeasible fixture; confirm modes 2a/2b run with SciPy absent.
4. `_translate.py` predicates + bijection + conversions; `__init__.py` public function + result type + dispatch.
5. Offer (do not build unasked): test suite with the fixtures, README, plotting helper.

## 11. Out of scope for v1

Seriation (recovering order from an unordered matrix); graphon inputs; NetworkX adapters; a pure-Python LP fallback; performance engineering beyond the chain-slot trick; proofs and complexity analysis in comments. A combinatorial feasibility test for $k\ge3$ is an open problem of the paper (L602), not an implementation target.