# SPEC — engineering (software design)

Companion to [`SPEC.md`](SPEC.md), the standalone core spec (definitions, data types, numbered pseudocode). This document covers everything that is **not** the algorithm: layers, the public interface, input handling, the LP-backend mechanics, numerics policy, naming, packaging, errors, and the implementation plan. It is self-contained; no paper is needed. A term glossary with source provenance is in [`docs/definitions.md`](docs/definitions.md) (optional reference). The round-7 review and the list of changes from the previous draft are in [`REVIEW-spec.md`](REVIEW-spec.md) and [`HANDOFF.md`](HANDOFF.md).

## 1. Layers and modules

```
outer (user / Jupyter)
   │  nested lists, numpy arrays, optional labels
   ▼
translation  ──  robinson_unif_emb/_translate.py
   │  validation predicates (boolean + immediate halt), format conversion,
   │  explicit bijection labels ↔ 0..n-1, thresholds → Fraction
   ▼
core  ──  robinson_unif_emb/_core.py
   │  the math pipeline (SPEC.md); pure Python stdlib; trusts its input;
   │  zero try/except/raise; no validation; no third-party imports; any k;
   │  the k>2 threshold search is passed in as a callable (no import of _lp)
   ▼
LP backend  ──  robinson_unif_emb/_lp.py   (optional; only the k>2 search)
      lazy import of scipy.optimize.linprog; consumes core's cycle bounds;
      may raise ImportError (SciPy absent) or RuntimeError (solver failure)
```

- **`_core.py`** — the algorithm only. Assumes a suitable internal input. No defensive checks, no exceptions, no external formats, no labels, no third-party imports. `solve_general(a, search)` takes the search function as an argument; the core never imports `_lp`. Comments follow the numbered pseudocode in `SPEC.md`.
- **`_lp.py`** — the only module importing SciPy, lazily inside `search_thresholds`. Also owns the deletion-filter witness procedure (`minimal_witnesses`, `SPEC.md` §5c). Used only by the $k>2$ automatic search.
- **`_translate.py`** — the only layer touching external formats: conversion, validation predicates, the label bijection, output conversion.
- **`__init__.py`** — exports the single public function and the result type; wires translation → core → verification → result. It supplies `_lp.search_thresholds` to the core only when the $k>2$ search is needed; it performs the `pi` verification (§3) and raises `RuntimeError` on failure.
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
| `thresholds` | sequence of $k$ numbers, optional | A candidate $\mathbf d\in\mathbb D^k$. Given → verify + construct (any $k$); converted to `Fraction`. Omitted → $k=1$ trivial, $k=2$ combinatorial search, $k>2$ LP search. |

**Dispatch** (`SPEC.md` §7): `thresholds` given → §5a; else $k=1$ → §5d; else $k=2$ → §5b; else → §5c.

**Result** — one frozen dataclass, always returned (infeasibility is an outcome, not an error):

```python
@dataclass(frozen=True)
class EmbeddingResult:
    feasible: bool
    positions: dict | None        # label → float (Π); None iff infeasible
    thresholds: tuple | None      # (d_1,...,d_k) floats; None iff infeasible
    cycles: tuple | None          # infeasible: simple witness cycles, each a tuple
                                  # of labels with first == last; None iff feasible
```

`cycles` is 1–2 cycles for $k=2$; for $k>2$ it is the deletion-filter set (no count bound); for $k=1$ infeasibility cannot occur. Usage sketch:

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
7. `thresholds` (if given): length $k$, strictly decreasing, positive, real; **converted to `Fraction`** (floats by exact binary conversion).

Any $k\ge1$ is accepted; mode selection needs no input predicate. SciPy availability is checked in `_lp.py`, not here.

**Bijection.** Build `index_of = {label: i}`; keep `labels` for the inverse. Core sees indices only; `positions` and witness `cycles` are mapped back through `labels` before return.

**Output conversion and verification.** `Fraction` → `float`. Run `core.verify_embedding(a, d, pi)` on every successful result before reporting feasibility; on `False` raise

```python
RuntimeError("numerical failure: the constructed embedding did not pass verification "
             "(k>2: try a different d or report the pruned-table limitation; "
             "k=1/k=2: report the table bug)")
```

`k=1` and `k=2` run in exact `Fraction` arithmetic, so a failure there is a table bug, not a numerical one — but it must not be reported as `feasible`.

```python
def _run(a, d0, labels):
    k = a[0][0]
    search = None
    if d0 is None and k > 2:
        from . import _lp
        search = _lp.search_thresholds
    out = _core.uniform_embedding_core(a, d0, search)
    d, pi = out
    if d is None:
        return EmbeddingResult(False, None, None, _map_cycles(out[1], labels))
    if not _core.verify_embedding(a, d, pi):
        raise RuntimeError("numerical failure ...")
    return EmbeddingResult(True, _map_positions(pi, labels), tuple(float(x) for x in d), None)
```

## 4. LP backend (`_lp.py`)

```python
def search_thresholds(cycles, k):
    # cycles: list[(bound, walk)] from core.cycle_bounds; walks are simple cycles.
    # -> (d, None)  with d a length-k tuple of floats, d[0] = 1
    # -> (None, witnesses) with witnesses a tuple of walks (simple cycles)
```

- **Import.** `from scipy.optimize import linprog` inside the function; absent → `ImportError("install scipy for k>2 automatic search: pip install scipy")`.
- **Formulation.** The LP (variables, constraints, `d[0]=1`, max-slack `t`) is specified in `SPEC.md` §5c. Implement it with `d[0]=1` substituted: variables `x = (d[1..k-1], t)`, minimize `-t`, `method="highs"`. Constraint rows are `Σ_{j≥1} b[j]·d[j] - t ≥ -b[0]` for each cycle and the ordering rows `1 - d[1] - t ≥ 0`, `d[i] - d[i+1] - t ≥ 0`, `d[k-1] - t ≥ 0`.
- **Decision.** `status == 0` → optimal; feasible iff `t* > EPS`; `status == 2` → infeasible. Any other status → `RuntimeError("LP solver failure: status=...")`. `EPS = 1e-9 · (1 + max over cycles of Σ_j |b[j]|)`, and `EPS = 1e-9` when `cycles` is empty.
- **Witnesses.** On infeasibility (or `t* ≤ EPS`) run the deletion filter (`SPEC.md` §5c): repeatedly solve the same LP on a trial subset; keep the subset when it stays infeasible; return the surviving walks. One LP per cycle. The result is contradictory together with `d ∈ D^k`; no count bound.
- **Dependency.** SciPy is optional (extra `lp`); the $k=1$, $k=2$ and given-`d` modes never import it.

## 5. Numerics

- `k=1`, `k=2`, given-`thresholds`: `fractions.Fraction` end to end. Ratio bounds `-a1/a2` and every midpoint are rational, so every strict inequality is decided exactly. Supplied thresholds are converted to `Fraction` (floats by exact binary conversion). Convert to `float` only at the public boundary.
- LP search: returns floats; the construction then runs in floats; the public layer verifies `pi` against the chosen `d` (converting floats to exact rationals, `Fraction(d)`) and raises `RuntimeError` on failure. This catches both a float boundary verdict (`t*` at round-off) and a missed cycle bound from the pruned table. The same check runs on every path; for `k=1` and `k=2` it is exact.
- An `exact=True` output option (keep `Fraction`s) remains a possible later addition; not in v1.

## 6. Internal data representations

The core data types (`a`, `n`, `k`, walk, bound, simple walk, cell, `T`, `d`, `pi`) are in `SPEC.md` §2 — that is the algorithm's contract. Translation builds `a` (tuple of tuples of `int`) and maps `pi`/witnesses back to labels (§3). The only public-facing types are `EmbeddingResult` and the `labels`/`thresholds` parameters (§2).

- **One cell interface.** A cell is a list of `(w, b)` pairs for `k≠2`, a slot dict `slot_key → (w, b)` for `k=2`. Every loop uses `cell_items(cell)`; never iterate a cell directly.
- **No lower-bound table.** Lower bounds are read by reversal from `T[v][i]` (`SPEC.md` §6).
- **Invariants.** All stored walks are simple ($|w|\le n$, bounds in $\mathbb Z^k_n$); `T[v][i]` is always nonempty for $i<v$ (the direct down-walk is seeded); `T[i][v]` may be empty; the diagonal is empty until the first cycle is generated (no seed).

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
| witness cycles | `cycles` (public), tuple of walks (core) |
| simple walk test | `simple` |
| cell iteration | `cell_items` |
| embedding check | `verify_embedding` |
| $k=1$ path | `solve_k1` |
| deletion filter | `minimal_witnesses` (`_lp`) |

## 8. Error model

- `ValueError` — translation, on the first failed input predicate (§3).
- `ImportError` — `_lp.py`, if the $k>2$ search is requested without SciPy.
- `RuntimeError` — `_lp.py` on an unexpected solver status; `__init__.py` when the verification of `pi` fails. Neither occurs inside the core.
- **Infeasibility is a returned `EmbeddingResult(feasible=False, cycles=...)`, never an exception.**
- The core never raises (zero try/except/raise).

## 9. Dependencies and packaging

- Core: Python stdlib only (`math`, `fractions`).
- Optional: SciPy (extra `lp`) for the $k>2$ search; numpy accepted on input but not required (duck-typed).
- `pyproject.toml` with `[project.optional-dependencies] lp = ["scipy"]`; package imports and modes `k=1`, `k=2`, given-`d` run with no third-party packages installed.

## 10. Implementation roadmap

The API is fixed by `SPEC.md` + this document. Steps 1–2 touch the same file (`_core.py`) and are sequential; steps 3–4 can be done in parallel with them once the core API is fixed:

1. `_core.py` primitives + `simple` + `cell_items` + `bound_generation` + `cycle_bounds` (`SPEC.md` §3–§5). Acceptance: A and B both complete; on B, a `(0,0)` cycle is present in some `T[v][v]`; table cells stay small (max `|b|₁ ≤ n`).
2. `_core.py` `embedding_from_bounds` + `verify_embedding` + `ratio_bounds` + `solve_k2` + `solve_k1` + `solve_with_thresholds` + `solve_general` (`SPEC.md` §5a–§5d, §6). Fixtures: A feasible; B infeasible with `(0,1,5,3,0)`; A₃ feasible; B₃ infeasible; k=1 example.
3. `_lp.py` `search_thresholds` + `minimal_witnesses` (`SPEC.md` §5c). Acceptance: A₃ feasible via the LP; B₃ infeasible with witnesses; confirm the `k=1`, `k=2` and given-`d` modes run with SciPy absent.
4. `_translate.py` predicates + bijection + conversions + the output verification; `__init__.py` public function + result type + wiring.
5. Offer (do not build unasked): test suite with the fixtures, README, plotting helper.

**Fixtures** (0-based):

- **A** (5×5, $k=2$, feasible): rows `(2,2,1,0,0) (2,2,2,1,1) (1,2,2,2,1) (0,1,2,2,2) (0,1,1,2,2)`. A valid embedding: `pi=(0,5,6.5,11.75,12.75)`, `d=(8,6)`.
- **B** (6×6, $k=2$, infeasible): rows `(2,2,1,0,0,0) (2,2,2,1,1,1) (1,2,2,2,1,1) (0,1,2,2,2,1) (0,1,1,2,2,2) (0,1,1,1,2,2)`. Witness cycle `(0,1,5,3,0)` has bound `(0,0)`.
- **A₃** (= A with diagonal 3, $k=3$, feasible): same rows with `a[0][0]=…=a[4][4]=3`. A valid embedding: the same `pi`, `d=(8,6,0.5)`.
- **B₃** (= B with diagonal 3, $k=3$, infeasible): same rows with `a[0][0]=…=a[5][5]=3`. Witness cycle `(0,1,5,3,0)` has bound `(0,0,0)`.
- **K1** ($k=1$, feasible): rows `(1,1,0,0) (1,1,1,0) (0,1,1,1) (0,0,1,1)`; any valid embedding is accepted (check with `verify_embedding`).

**Acceptance checks.** A: `feasible`, `verify_embedding` true, `pi` strictly increasing. B: `feasible=False`, one or two simple cycles, the `(0,0)` cycle present or an equivalent pair. A₃: feasible, verified. B₃: `feasible=False` through the LP, witnesses jointly contradictory with $d\in\mathbb D^3$. K1: feasible, verified, with no SciPy installed.

## 11. Out of scope for v1

Seriation (recovering order from an unordered matrix); graphon inputs; NetworkX adapters; a pure-Python LP fallback; performance engineering beyond the chain-slot trick; proofs and complexity analysis in comments. A combinatorial feasibility test for $k\ge3$ is an open problem (none is known); not an implementation target. A completeness proof for the pruned $k>2$ table (`SPEC.md` §4) is also out of scope; the verification step is the v1 safeguard.
