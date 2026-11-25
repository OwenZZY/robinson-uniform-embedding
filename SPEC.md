# SPEC — uniform embedding for Robinson matrices (standalone core)

This specification is **self-contained**: a software engineer can implement the core algorithm from this document alone, without the paper. All indices are **0-based**. Pseudocode lines are numbered **L1, L2, …** per function; prose refers to those internal line numbers. The public API, input validation, LP-solver call details, naming and packaging are in [`SPEC-engineering.md`](SPEC-engineering.md) (also self-contained). A term glossary with source provenance is in [`docs/definitions.md`](docs/definitions.md) (not needed to implement). The round-7 review and the list of changes from the previous draft are in [`REVIEW-spec.md`](REVIEW-spec.md) and [`HANDOFF.md`](HANDOFF.md); **this document is the contract**.

## 1. Problem and definitions

**Robinson matrix.** `a` is an `n×n` symmetric integer matrix, entries in `{0,…,k}`, every diagonal entry equal to `k` (the maximum), and for all `u<v<w`: `a[u][v] ≥ a[u][w]` and `a[v][w] ≥ a[u][w]`. (`k = a[0][0]`, `n = len(a)`.)

**Threshold vector.** `d` is a length-`k` tuple with `d[0] > d[1] > … > d[k-1] > 0`.

**Uniform embedding.** A strictly increasing list `pi[0..n-1]` such that for every `u<v` with `a[u][v] = t`:

```
lower(t) < pi[v] - pi[u] < upper(t)
upper(t) = d[t-1] if t ≥ 1, else +∞
lower(t) = d[t]   if t ≤ k-1, else 0
```

**Goal.** Given `a`, return `(d, pi)` (a threshold vector and a uniform embedding), or `(None, witnesses)` proving none exists. Optionally a candidate `d` is supplied; then only verify it and construct `pi`.

**Edge bounds.** For `u≠v`, let `t = a[u][v]`:
- `u<v`: `beta_plus(u,v)` = unit vector at index `t-1` if `t≥1`, **undefined** if `t=0`; `beta_minus(u,v)` = unit vector at index `t` if `t≤k-1`, else `zero(k)`.
- `u>v`: `beta_plus(u,v) = -beta_minus(v,u)`; `beta_minus(u,v) = -beta_plus(v,u)` (undefined iff `a[v][u]=0`).
- There is no bound for `u=v`: the one-vertex walk has bound `zero(k)` and is never stored (see §4).

`unit(k,j)` = length-`k` tuple with 1 at index `j`; `zero(k)` = all zeros.

**Walks.** A walk `w` is a tuple of vertex indices. An **upper-bound walk** (UBW) requires every up-step (`w[i-1] < w[i]`) to satisfy `a[w[i-1]][w[i]] > 0`; down-steps are unrestricted. A **lower-bound walk** (LBW) requires every down-step to satisfy `a > 0`; up-steps unrestricted. The bound of a UBW is `beta_plus(W) = Σ_i beta_plus(w[i-1], w[i])`. Reversing a UBW gives an LBW and negates the bound.

**Simple walks.** A walk of length ≥ 1 is **simple** if all its vertices are distinct when it is open, or `len(w) ≥ 3`, `w[0] = w[-1]`, and no other vertex repeats when it is closed. A **simple cycle** is a closed simple walk (at least two distinct vertices). The bound table stores **simple walks only** (upper-bound paths off the diagonal, simple upper-bound cycles on the diagonal). Simple walks have at most `n` vertices, so every stored bound has `Σ_j |b[j]| ≤ n`.

**Prefix order.** `x ≼ y` iff for every `j = 1..k`, `sum(x[0..j-1]) ≤ sum(y[0..j-1])`. This is equivalent to `x·d ≤ y·d` for every threshold vector `d`, so only `≼`-minimal bounds matter.

**Feasibility condition.** A uniform embedding exists iff some threshold vector `d` satisfies `beta_plus(C)·d > 0` for every upper-bound cycle `C`. It suffices to check the `≼`-minimal cycle bounds.

## 2. Data types

| Object | Type | Meaning |
|---|---|---|
| `a` | `tuple[tuple[int,…],…]`, `n×n` | Robinson matrix |
| `n`, `k` | `int` | `len(a)`, `a[0][0]` |
| walk `w` | `tuple[int,…]` | vertex indices; concat `w1 + w2[1:]`, reverse `w[::-1]` |
| bound `b` | `tuple[int,…]` length `k`; `None` = undefined | integer coefficient vector |
| `cell` | `k=2`: `dict[slot_key,(w,b)]`; else `list[(w,b)]` | minimal **simple** UBWs of one ordered pair; diagonal cells hold simple cycles |
| `T` | `list[list[cell]]`, `n×n` | `T[i][j]` = minimal simple `(i,j)`-UBWs |
| `d` | length-`k` tuple: `Fraction` for `k=1`, `k=2` and given-`d`; `float` on the `k>2` LP path | threshold vector |
| `pi` | length-`n` list, same numeric type as `d` | embedding; `pi[0]=0` |

Helpers: `unit(k,j)`, `zero(k)`, `vec_add`, `neg`, `dot`, `simple(w)`, `preceq`, `cell_items(cell)`, `INF = math.inf`.

`cell_items(cell)` returns the `(w,b)` pairs of a cell: the list itself for `k≠2`, `cell.values()` for the `k=2` slot dict. **Every loop over a cell uses `cell_items`**; never iterate a cell directly.

## 3. Primitives

```text
beta_plus(a, k, u, v):                 # bound tuple or None; u ≠ v
L1  if u < v:
L2      t ← a[u][v]
L3      return unit(k, t-1) if t ≥ 1 else None
L4  return neg(beta_minus(a, k, v, u))

beta_minus(a, k, u, v):                # bound tuple or None; u ≠ v
L1  if u < v:
L2      t ← a[u][v]
L3      return unit(k, t) if t ≤ k-1 else zero(k)
L4  b ← beta_plus(a, k, v, u)
L5  return neg(b) if b ≠ None else None

walk_bound(a, k, w):                   # w is a UBW
L1  b ← zero(k)
L2  for i in 1 .. len(w)-1:
L3      b ← vec_add(b, beta_plus(a, k, w[i-1], w[i]))
L4  return b

simple(w):
L1  if len(w) < 2: return False
L2  if w[0] = w[-1]: return len(w) ≥ 3 and len(set(w[:-1])) = len(w)-1
L3  return len(set(w)) = len(w)

preceq(x, y):                          # x ≼ y
L1  for j in 1 .. k:
L2      if sum(x[0..j-1]) > sum(y[0..j-1]): return False
L3  return True

cell_items(cell):
L1  return cell.values() if cell is a dict else cell
```

## 4. Bound generation

`bound_generation(a)` builds `T`: for every ordered pair `(i,j)`, the `≼`-minimal simple `(i,j)`-UBWs, each stored as `(walk, bound)`. Diagonal cells hold simple upper-bound cycles through `v`.

```text
bound_generation(a):
L1  n ← len(a); k ← a[0][0]
L2  T ← n×n table of empty cells
L3  for i in 0 .. n-1:                          # seed direct walks (i < j only)
L4      for j in i+1 .. n-1:
L5          if a[i][j] ≠ 0: insert(T[i][j], (i,j), beta_plus(a,k,i,j), k)
L6          insert(T[j][i], (j,i), beta_plus(a,k,j,i), k)
L7  for s in 0 .. n-1:                          # sweep ALL ordered pairs
L8      for i in 0 .. n-1:
L9          for j in 0 .. n-1:
L10             for (w1,b1) in snapshot(cell_items(T[i][s])):
L11                 for (w2,b2) in snapshot(cell_items(T[s][j])):
L12                     w ← w1 + w2[1:]; b ← vec_add(b1,b2)
L13                     insert(T[i][j], w, b, k)
L14 return T

insert(cell, w, b, k):
L1  if not simple(w): return
L2  if k = 2: insert_minimal_k2(cell, w, b)
L3  elif minimal_in(b, cell): insert_minimal(cell, w, b)
```

Cell operations:
- `minimal_in(b, cell)`: no stored bound `c` with `preceq(c,b)` and `c≠b`.
- `insert_minimal(cell, w, b)`: skip if not `minimal_in(b, cell)`; otherwise add `(w,b)` and remove every stored `(w',c)` with `preceq(b,c)` and `b≠c`.
- `k=2` cells are chain slots (insertion O(1)):

```text
slot_key(b) = (|b[0]|+|b[1]|, 1 if b[1]>0 else 2)

insert_minimal_k2(cell, w, b):
L1  key ← slot_key(b)
L2  if key not in cell: cell[key] ← (w,b)
L3  elif preceq(b, cell[key][1]) and b ≠ cell[key][1]: cell[key] ← (w,b)
```

**Notes (why the code is shaped this way).**
- **Simple walks only.** The paper's `UBW(i,j)` is a set of upper-bound *paths*; the diagonal is a set of simple *cycles*. Storing walks with repeated vertices lets bounds leave `Z^k_n` (e.g. a bound that winds around a cycle twice is incomparable to the one that winds once), so cells grow without bound — with such walks Example 1B does not finish (measured). With L1 of `insert`, every stored bound has `Σ_j |b[j]| ≤ n`; for `k=2`, `t = |b[0]|+|b[1]| ≤ n`, so a cell has at most `2n+2` slots, and for general `k` it is a `≼`-antichain inside `Z^k_n` (the paper's `M`).
- **No diagonal seed.** The one-vertex walk `(v,v)` is not a simple cycle (L2 requires `len(w) ≥ 3`) and is never inserted. Cycles enter through L13 with `i=j`: a simple cycle is split at a vertex `s` into two shorter simple paths, both already present at stage `s`, and their concatenation is inserted. An empty diagonal before the first cycle is expected.
- L7–L9 sweep **all** ordered pairs, including `i>j`. A down-direction cell `T[j][i]` cannot be filled by reversing `T[i][j]`: reversal turns a UBW into an LBW. Both directions must be generated directly. (Restricting to `i≤j` misses witnesses — Example 1B is reported feasible.)
- L10–L11 iterate over a **snapshot**: when `i=s` or `j=s` the cell is mutated in place during the sweep.
- Pruning is safe because `≼` is additive: `x≼y ⇒ x+z≼y+z`, so any concatenation using a pruned bound is dominated by the one using the kept bound.
- **Known limitation (`k>2`).** Pruning is by bounds, not by the walks' vertex sets: a pruned simple path can be needed to assemble a simple cycle. The `k=2` ratio test was checked against exhaustive simple-cycle enumeration (600 random matrices, exact); for `k>2` the pruned table can in principle miss a cycle bound (a 6×6, `k=4` instance was found). The public layer therefore verifies the constructed `pi` before returning (`SPEC-engineering.md` §3) and reports a numerical/table failure rather than a false feasible.

## 5. Thresholds (the branch)

`cycle_bounds(T)` extracts the distinct minimal cycle bounds:

```text
cycle_bounds(T):
L1  n ← len(T); seen ← empty map bound → walk
L2  for v in 0 .. n-1:
L3      for (w,b) in cell_items(T[v][v]):
L4          if b not in seen: seen[b] ← w
L5  return [(b, seen[b]) for b in seen]
```

Every `w` is already a simple cycle (by §4), so no post-processing is needed.

### 5a. Candidate `d` supplied (any `k`) — verify

```text
solve_with_thresholds(a, d):
L1  n ← len(a); T ← bound_generation(a)
L2  for v in 0 .. n-1:
L3      for (w,b) in cell_items(T[v][v]):
L4          if dot(b,d) ≤ 0: return None, (w,)
L5  return d, embedding_from_bounds(a, d, T)
```

### 5b. `k = 2` — combinatorial

Each cycle bound `(a1,a2)` constrains the ratio `r = d[1]/d[0]`: `a2>0 ⇒ r > -a1/a2`; `a2<0 ⇒ r < -a1/a2`; `a2=0 ⇒` require `a1>0`. Use exact rational arithmetic (`Fraction`).

```text
ratio_bounds(T):                       # k = 2
L1  n ← len(T); lo ← 0; hi ← 1; w_lo ← None; w_hi ← None     # r ∈ (0,1)
L2  for v in 0 .. n-1:
L3      for (w,b) in cell_items(T[v][v]):
L4          (a1,a2) ← b
L5          if a2 > 0:
L6              r ← Fraction(-a1,a2)
L7              if r ≥ 1: return INFEASIBLE, (w,)       # contradicts r < 1
L8              if r > lo: lo ← r; w_lo ← w
L9          elif a2 < 0:
L10             r ← Fraction(-a1,a2)
L11             if r ≤ 0: return INFEASIBLE, (w,)       # contradicts r > 0
L12             if r < hi: hi ← r; w_hi ← w
L13         else:
L14             if a1 ≤ 0: return INFEASIBLE, (w,)
L15 if lo ≥ hi: return INFEASIBLE, (w_lo, w_hi)
L16 return FEASIBLE, (lo, hi)

solve_k2(a):
L1  T ← bound_generation(a)
L2  status, info ← ratio_bounds(T)
L3  if status = INFEASIBLE: return None, info
L4  (lo,hi) ← info; d ← (Fraction(1), (lo+hi)/2)
L5  return d, embedding_from_bounds(a, d, T)
```

Notes: L1 seeds the trivial window `(0,1)`. L7 and L11 handle bounds that are individually unsatisfiable inside the window, so L15 always returns two genuine cycles. L13–L14 handle `a2=0` (the inequality is `a1·d[0]>0`); a bound `(0,0)` is immediately infeasible. The witnesses are simple cycles by §4.

### 5c. `k > 2` — LP feasibility

Feasibility of `{b·d > 0 for each minimal cycle bound b; d[0]>…>d[k-1]>0}` is a linear program. By homogeneity fix `d[0]=1` and maximize the minimum slack `t`; the strict system is feasible iff the optimum `t* > EPS`.

```text
search_thresholds(cycles, k):          # lives in _lp.py; see SPEC-engineering.md §4
L1  variables x = (d[0..k-1], t)
L2  maximize t subject to:
L3      dot(b,d) - t ≥ 0            for each (b,_) in cycles
L4      d[i] - d[i+1] - t ≥ 0       for i in 0 .. k-2
L5      d[k-1] - t ≥ 0
L6      d[0] = 1
L7  scale ← 1 + max over cycles of Σ_j |b[j]|;  EPS ← 1e-9 · scale
L8  if infeasible or optimal t ≤ EPS: return None, minimal_witnesses(cycles, k)
L9  return d, None

minimal_witnesses(cycles, k):          # lives in _lp.py
L1  keep ← cycles
L2  for (b,w) in cycles:
L3      trial ← keep without (b,w)
L4      if infeasible(trial, k): keep ← trial
L5  return tuple(w for (b,w) in keep)
```

`infeasible(trial, k)` re-runs the same LP on the subset `trial`; it is true iff the LP is infeasible or its optimum is `≤ EPS`. `minimal_witnesses` is a deletion filter: it returns a subset of cycles that is contradictory **together with `d ∈ D^k`**. The count is not bounded by `k` (in practice it is small).

```text
solve_general(a, search):              # search = _lp.search_thresholds
L1  k ← a[0][0]; T ← bound_generation(a)
L2  d, witnesses ← search(cycle_bounds(T), k)
L3  if d = None: return None, witnesses
L4  return d, embedding_from_bounds(a, d, T)
```

Note: `search` is passed in by the public layer; `_core.py` imports no third-party module. The LP path runs in floats; the public layer verifies `pi` (§6) before returning and raises `RuntimeError` on failure.

### 5d. `k = 1` — trivial

For `k=1` every binary Robinson matrix is feasible (an indifference graph), so no search is needed; any positive `d[0]` is a valid threshold vector.

```text
solve_k1(a):
L1  T ← bound_generation(a)
L2  d ← (Fraction(1),)
L3  return d, embedding_from_bounds(a, d, T)
```

## 6. Construction

Given a feasible `d` and `T`, build `pi`. Lower bounds are read by reversal: the best lower bound on `(i,v)` is `-min{beta_plus(W)·d : W in T[v][i]}`.

```text
min_bound(cell, d):
L1  return min(dot(b,d) for (_,b) in cell_items(cell))

embedding_from_bounds(a, d, T):
L1  n ← len(a); pi ← [0]*n
L2  for v in 1 .. n-1:
L3      ub ← min(pi[i] + min_bound(T[i][v], d) for i in 0..v-1 if cell_items(T[i][v]) nonempty); INF if none
L4      lb ← max(pi[i] - min_bound(T[v][i], d) for i in 0..v-1)
L5      pi[v] ← (ub+lb)/2 if ub < INF else lb + d[0]
L6  return pi

verify_embedding(a, d, pi):
L1  n ← len(a); k ← a[0][0]
L2  for u in 0 .. n-1:
L3      for v in u+1 .. n-1:
L4          t ← a[u][v]
L5          lo ← d[t] if t ≤ k-1 else 0
L6          hi ← d[t-1] if t ≥ 1 else INF
L7          x ← pi[v] - pi[u]
L8          if not (lo < x) or not (x < hi): return False
L9  return True
```

Notes: L3 skips empty cells (no upper bound exists for that pair); if every cell is empty, `ub=INF` and L5 places `pi[v]` just above `lb`. `T[v][i]` is always nonempty (the direct down-walk is seeded), so `lb` is finite. For feasible `d`, `lb < ub` and `pi` satisfies the inequalities (paper's lemmas). `verify_embedding` is an exact boolean check; the public layer runs it on every successful return (and raises `RuntimeError` if it fails) and it raises nothing itself.

## 7. Top-level dispatch and output

```text
uniform_embedding_core(a, d=None, search=None):
L1  k ← a[0][0]
L2  if d ≠ None: return solve_with_thresholds(a, d)
L3  if k = 1: return solve_k1(a)
L4  if k = 2: return solve_k2(a)
L5  return solve_general(a, search)
```

**Output.** Success: `(d, pi)`, `pi` strictly increasing, `d` scale-free (`d[0]=1` for the `k=1`, `k=2` and LP paths; the supplied `d` is returned unchanged). Failure: `(None, witnesses)`.
- `k=2`: one witness — an individually unsatisfiable cycle (bound `(a1,0)` with `a1≤0`, or a ratio bound outside `(0,1)`) — or two cycles whose ratio intervals are disjoint.
- `k>2`: the cycles returned by the deletion filter. They are **jointly** contradictory together with `d ∈ D^k`; no single one need be. The count is not promised to be `≤ k`.
- `k=1`: no failure.

All witnesses are simple upper-bound cycles. The public wrapper (`SPEC-engineering.md`) maps vertex indices to external labels, converts `Fraction` to `float`, and verifies `pi` before returning.
