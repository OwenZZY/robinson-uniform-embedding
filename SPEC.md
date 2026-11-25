# SPEC — `robinson-unif-emb` (v1)

Implementation spec for the algorithm of Janssen–Zhang, *Uniform Embeddings for Robinson Similarity Matrices* (arXiv:2105.09197v1). Supersedes `spec0-by-hand.md`. This document is the hand-off to the implementation chat: it fixes the layers, the interface, the names, and the data representations. The mathematics lives in `docs/definitions.md` (terms) and `docs/pseudocode.md` (transcription P1–P8 with resolved conventions R1–R10); read those first. Known source issues and the adopted fixes are in `REPORT.md` (findings F1–F10).

**Scope of v1.** The complete $k=2$ combinatorial algorithm (paper §5.3.1): return a uniform embedding or witness cycles proving none exists. Plus the general-$k$ construction from *supplied* thresholds (paper §4.2 + cycle-system check). No LP search for $k>2$, no proofs, no complexity analysis, no tests/README/examples unless offered and accepted.

## 1. Layers

```
outer (user / Jupyter)
   │  nested lists, numpy arrays, optional labels
   ▼
translation  ──  robinson_unif_emb/_translate.py
   │  validation predicates (boolean + immediate halt), format conversion,
   │  explicit bijection labels ↔ 0..n-1, float/Fraction conversion
   ▼
core  ──  robinson_unif_emb/_core.py
      paper algorithms only; trusts its input; pure Python stdlib;
      zero try/except/raise; no validation; paper notation in comments
```

- **Core** (`_core.py`): assumes a perfectly suitable input (a Robinson matrix in the internal representation, a feasible threshold vector when one is required). Correctness is the paper's; the code is a straight transcription of `docs/pseudocode.md`. No defensive checks, no exceptions, no external formats, no labels.
- **Translation** (`_translate.py`): the *only* layer touching external formats. Converts input to the internal representation, runs boolean validation predicates and halts on the first failure (`ValueError` with a plain message), builds the label bijection, and maps results back to external identifiers.
- **Outer**: the user's notebook. Nothing to build beyond a friendly signature and a readable result object.
- **Intermediate layer**: not needed. There is no graph-object input; matrices and vectors translate directly. (If a NetworkX front end is ever wanted, it belongs in a separate adapter, not in core.)

## 2. Public interface (single function)

```python
# robinson_unif_emb/__init__.py
__all__ = ["uniform_embedding", "EmbeddingResult"]

def uniform_embedding(A, labels=None, thresholds=None) -> EmbeddingResult: ...
```

**Parameters**

| Name | Type | Meaning |
|---|---|---|
| `A` | square nested sequence or numpy array of ints | Robinson similarity matrix, entries $\{0,\dots,k\}$, diagonal all $k$ (the maximum). $n=$ size, $k=A[0][0]$. |
| `labels` | sequence of $n$ unique hashables, optional | External identifiers for the items, in matrix order. Default: `tuple(range(n))`. The bijection is `labels[i] ↔ i`. |
| `thresholds` | sequence of $k$ numbers, optional | A candidate $\mathbf d\in\mathbb D^k$ ($d_1>\dots>d_k>0$). If given (any $k$): test it against the cycle system and, if feasible, construct the embedding for it (pseudocode P8). If omitted: $k$ must be $2$ and the combinatorial pipeline finds $\mathbf d$ itself (P7). |

**Modes**

1. `k == 2`, `thresholds is None` → full decision algorithm; finds thresholds or proves infeasibility.
2. `thresholds` given (any `k >= 1`) → feasibility check of that vector + construction. Exact if the supplied numbers are `Fraction`s.
3. `k > 2`, `thresholds is None` → `ValueError` (LP search for $k>2$ is out of scope for v1).

**Return value** — one frozen dataclass, always returned (infeasibility is a mathematical outcome, not an error):

```python
@dataclass(frozen=True)
class EmbeddingResult:
    feasible: bool
    positions: dict | None        # label → float, the embedding Π; None iff infeasible
    thresholds: tuple | None      # (d_1, ..., d_k) as floats; None iff infeasible
    cycles: tuple | None          # infeasible: 1 or 2 witness cycles, each a tuple of labels
                                  # (closed walks: first label == last label); None iff feasible
```

- Success: `feasible=True`, `positions` and `thresholds` populated, `cycles=None`. `positions` is strictly increasing in matrix order and satisfies $d_{t+1}<|\Pi(v)-\Pi(u)|<d_t$ whenever $A_{u,v}=t$ (with $d_0=\infty$, $d_{k+1}=0$).
- Failure: `feasible=False`, `positions=thresholds=None`, `cycles` holds the witness cycle(s) $C_1,C_2$ (paper step 3). A single cycle is returned when one cycle alone is contradictory (e.g. bound $(0,0)$, or $a_2=0,a_1\le0$).
- Thresholds are scale-free; the $k=2$ pipeline normalizes $d_1=1$ (convention R8).

**Errors** — `ValueError` from the translation layer only, on the first failed predicate (§5). The core never raises.

**Usage sketch** (interface illustration, not implementation):

```python
from robinson_unif_emb import uniform_embedding
res = uniform_embedding([[2, 2, 1], [2, 2, 2], [1, 2, 2]], labels=["a", "b", "c"])
res.feasible, res.positions, res.thresholds
```

## 3. Core internal API (`_core.py`)

All functions below are project-internal (private module; not exported). They operate on the internal representations of §4, 0-based, and trust their inputs. Comments cite the paper line by line (see `docs/pseudocode.md`).

| Function | Paper source | Input → Output |
|---|---|---|
| `beta_plus(a, k, u, v)` | def:bounds L186–205; def:beta vu L228–230 | ordered pair → bound tuple or `None` (undefined) |
| `beta_minus(a, k, u, v)` | same | ordered pair → bound tuple or `None` |
| `walk_bound(a, k, w)` | eqn:beta(W) L233 | UBW tuple → bound tuple (sum of step bounds) |
| `preceq(x, y)` | def:boundleq L465–470 | two bound tuples → bool (prefix-sum order) |
| `bound_generation(a)` | alg:FW-mod L508–553, variant R5 | matrix → table `T` of minimal upper-bound walks, all ordered pairs |
| `ratio_bounds(T)` | §5.3 L565–569, R7/R8 | table ($k=2$) → `(FEASIBLE, (lo, hi))` or `(INFEASIBLE, (w1, w2))`, rational window + witness walks |
| `embedding_from_bounds(a, d, T)` | alg:embedding algorithm L343–353, R6 | matrix, feasible `d`, table → `pi` (list of positions) |
| `solve_k2(a)` | §5.3.1 L574–583 | matrix → `(d, pi)` or `(NO_SOLUTION, (w1, w2))` |
| `solve_with_thresholds(a, d)` | cycle system L450–453 + construction | matrix, `d` → same shape |

Implementation notes carried over from the pseudocode (do not re-litigate):

- **R5**: FW sweeps run over *all* ordered pairs `(i, j)`; the reversal update of L546 is dropped (it under-generates down-direction cells; see `REPORT.md` F1). Init covers both directions.
- **R3/R4**: minimality is on $\beta^+$ vectors under `preceq`; inserting a new minimum prunes strictly dominated stored bounds.
- **R2**: `beta_plus(i, i) = χ_k`, `beta_minus(i, i) = 0`; the loop `⟨i, i⟩` seeds each diagonal cell.
- **R6**: `min ∅ = +∞`; if `ub_v` is infinite, `pi[v] = lb_v + d[0]`.
- **R7**: cycle bound with `a2 == 0`: infeasible iff `a1 <= 0`.
- **R8**: ratio window initialized to `(0, 1)`; final `r = (lo + hi) / 2`, `d = (1, r)`.
- **R9**: $k=2$ cells are chain-slot dicts keyed by `(abs(a1)+abs(a2), 1 if a2 > 0 else 2)`; one (tightest) element per slot; `t = 0` allowed.
- **R10**: $k=2$ pipeline uses `fractions.Fraction` end to end; `float` conversion happens only in translation.
- Snapshot cells before iterating when `i == s` or `j == s` (in-place updates during the sweep).
- Witness walks may be non-simple closed walks; reduce to a simple cycle by splitting at repeated vertices and keeping the extremal summand (`simple_cycle` in P5).

## 4. Internal data representations

| Object | Representation | Notes |
|---|---|---|
| matrix `a` | `tuple[tuple[int, ...], ...]`, $n\times n$, 0-based | built by translation; core indexes `a[u][v]` like $a_{u,v}$ |
| `n`, `k` | `len(a)`, `a[0][0]` | diagonal is all $k$ (trusted) |
| walk `w` | `tuple[int, ...]` | vertex indices; concatenation `w1 + w2[1:]`; reversal `w[::-1]` |
| bound `b` | `tuple[int, ...]` of length $k$ | $\beta^+$ values; `None` = undefined |
| cell | $k=2$: `dict[(t, chain), (w, b)]`; general $k$: `list[(w, b)]` antichain | minimal upper-bound walks per ordered pair |
| table `T` | `list[list[cell]]`, $n\times n$ | `T[i][j]` = $(i,j)$-UBW minima; diagonal holds cycles |
| thresholds `d` | `tuple[Fraction, ...]` (auto) or tuple of the user's numeric type (given) | core arithmetic is generic over the type |
| positions `pi` | `list` of the same numeric type as `d` | `pi[0] = 0` |

Lower bounds are never stored: `LowerBounds(i, v) = {−b : b ∈ β+(T[v][i])}` by reversal, so `lb_v = max_i (pi[i] − min_bound(T[v][i], d))` (P6).

## 5. Translation layer: validation and bijection

Boolean predicates, evaluated in order; first `False` → `ValueError("<predicate message>")`. All live in `_translate.py`; none in core.

1. square, $n\ge1$;
2. entries are integers (numpy integer types accepted) in $\{0,\dots,k\}$ with $k=A[0][0]$;
3. symmetric;
4. diagonal constant and equal to the maximum entry $k$;
5. Robinson: for all $u<v<w$, $a_{u,v}\ge a_{u,w}$ and $a_{v,w}\ge a_{u,w}$;
6. `labels`: length $n$, hashable, unique (if given);
7. `thresholds` (if given): length $k$, strictly decreasing, positive, real numbers;
8. mode check: `k > 2` requires `thresholds`.

Bijection: translation builds `index_of = {label: i}` and keeps `labels` for the inverse. Core receives indices only; `positions` and witness `cycles` are mapped back through `labels` before return. Input conversion: numpy → `tolist()` (duck-typed), nested sequences → tuple of tuples of `int`. Output conversion: `Fraction` → `float` (default; see to-decide T2 in `REPORT.md`).

## 6. Naming

- Package `robinson_unif_emb` (import name); distribution/GitHub name `robinson-uniform-embedding`.
- Public: `uniform_embedding` (the paper's term, generic enough for reuse), `EmbeddingResult`.
- Project-internal: modules `_core`, `_translate`; anything not exported is implementation detail. Inside `_core`, names follow the paper where short (`pi`, `d`, `ub`, `lb`, `beta_plus`, `preceq`) so a mathematician can read against the source; docstrings/comments carry the paper notation and line refs.
- Convention: a name without a leading underscore in a private module is still private by module; only `__all__` in `__init__.py` defines the public surface. No paper macros (`\bmD`, `\upperbd`) leak into identifiers.

| Paper | Code |
|---|---|
| $\Pi$ | `pi` (core), `positions` (public) |
| $\mathbf d$, $\mathbb D^k$ | `d` (core), `thresholds` (public) |
| $\beta^+(u,v)$, $\beta^-(u,v)$ | `beta_plus`, `beta_minus` |
| $\beta^+(W)$ | `walk_bound` |
| $\preceq$ | `preceq` |
| UBW table | `T` (core), cells `T[i][j]` |
| $ub_v$, $lb_v$ | `ub`, `lb` (locals in `embedding_from_bounds`) |
| FW-var / Bound-Generation | `bound_generation` |
| $d_2/d_1$ window | `ratio_bounds` → `(lo, hi)` |
| $C_1, C_2$ | `cycles` (public), `(w1, w2)` (core) |
| $\mathcal S^n[k]$ | no runtime object; the validation predicates define it |

## 7. Implementation roadmap (for the next chat)

1. `_core.py` primitives: `beta_plus`, `beta_minus`, `walk_bound`, `preceq` (P1–P2). Check by hand against `docs/definitions.md`.
2. `bound_generation` with general-$k$ antichain cells and $k=2$ chain slots (P3–P4, R2–R5, R9). Verify on Example 1B that the diagonal contains the bound `(0, 0)` (this is the regression test for F1).
3. `ratio_bounds` + `solve_k2` (P5, P7; R7, R8, R10). Fixtures: Example 1A → feasible; Example 1B → infeasible with witness `⟨0,1,5,3,0⟩` (0-based) or an equivalent cycle.
4. `embedding_from_bounds` (P6, R6). Verify the returned `pi` against the strict inequalities for every pair on the fixtures.
5. `_translate.py` predicates + bijection + conversions; `__init__.py` exporting the single function and result type (§2, §5).
6. `solve_with_thresholds` (P8) wired to the `thresholds` parameter.
7. Offer (do not build unasked): a small test suite with the two fixtures, a README, a plotting helper.

## 8. Out of scope for v1

LP search of a feasible $\mathbf d$ for $k>2$ (paper §5, Vaidya discussion); seriation (recovering the order from an unordered matrix); graphon inputs; NetworkX adapters; performance engineering beyond the paper's chain-slot trick; proofs and complexity analysis in code comments.