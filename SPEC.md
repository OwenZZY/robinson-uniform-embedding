# SPEC — uniform embedding for Robinson matrices (standalone core)

This specification is **self-contained**: a software engineer can implement the core algorithm from this document alone, without the paper. All indices are **0-based** and already correct in the pseudocode — there is nothing to convert. Pseudocode lines are numbered **L1, L2, …** per function; prose refers to those internal line numbers. The public API, input validation, LP-solver call details, naming and packaging are in [`SPEC-engineering.md`](SPEC-engineering.md) (also self-contained). A term glossary with source provenance is in [`docs/definitions.md`](docs/definitions.md) (not needed to implement).

## 1. Problem and definitions

**Robinson matrix.** `a` is an `n×n` symmetric integer matrix, entries in `{0,…,k}`, every diagonal entry equal to `k` (the maximum), and for all `u<v<w`: `a[u][v] ≥ a[u][w]` and `a[v][w] ≥ a[u][w]`. (`k = a[0][0]`, `n = len(a)`.)

**Threshold vector.** `d` is a length-`k` tuple with `d[0] > d[1] > … > d[k-1] > 0`.

**Uniform embedding.** A strictly increasing list `pi[0..n-1]` such that for every `u<v` with `a[u][v] = t`:

```
lower(t) < pi[v] - pi[u] < upper(t)
upper(t) = d[t-1] if t ≥ 1, else +∞
lower(t) = d[t]   if t ≤ k-1, else 0
```

**Goal.** Given `a`, return `(d, pi)` (a threshold vector and a uniform embedding), or witness cycles proving none exists. Optionally a candidate `d` is supplied; then only verify it and construct `pi`.

**Edge bounds.** For `u<v`, let `t = a[u][v]`. Define two length-`k` integer vectors:
- `beta_plus(u,v)` = unit vector with 1 at index `t-1` if `t≥1`; **undefined** if `t=0`. Its dot product with `d` is `upper(t)`.
- `beta_minus(u,v)` = unit vector with 1 at index `t` if `t≤k-1`; the zero vector if `t=k`. Its dot product with `d` is `lower(t)`.

For `u>v`: `beta_plus(u,v) = -beta_minus(v,u)` and `beta_minus(u,v) = -beta_plus(v,u)`. For `u=v`: `beta_plus(u,u)` = unit vector at index `k-1`, `beta_minus(u,u)` = zero. (`unit(k,j)` = length-`k` tuple, 1 at index `j`.)

**Walks.** A walk `w` is a tuple of vertex indices. An **upper-bound walk** (UBW) requires every up-step (`w[i-1] < w[i]`) to satisfy `a[w[i-1]][w[i]] > 0`; down-steps are unrestricted. A **lower-bound walk** (LBW) requires every down-step to satisfy `a > 0`; up-steps unrestricted. The bound of a UBW is `beta_plus(W) = Σ_i beta_plus(w[i-1], w[i])`. Reversing a UBW gives an LBW and negates the bound. An **upper-bound cycle** is a closed UBW (`w[0] = w[-1]`); for witnesses it is simple (no other repeated vertex).

**Prefix order.** `x ≼ y` iff for every `j = 1..k`, `sum(x[0..j-1]) ≤ sum(y[0..j-1])`. This is equivalent to `x·d ≤ y·d` for every threshold vector `d`, so only `≼`-minimal bounds matter.

**Feasibility condition.** A uniform embedding exists iff some threshold vector `d` satisfies `beta_plus(C)·d > 0` for every upper-bound cycle `C`. It suffices to check the `≼`-minimal cycle bounds.

## 2. Data types

| Object | Type | Meaning |
|---|---|---|
| `a` | `tuple[tuple[int,…],…]`, `n×n` | Robinson matrix |
| `n`, `k` | `int` | `len(a)`, `a[0][0]` |
| walk `w` | `tuple[int,…]` | vertex indices; concat `w1 + w2[1:]`, reverse `w[::-1]` |
| bound `b` | `tuple[int,…]` length `k`; `None` = undefined | integer coefficient vector |
| cell | `k=2`: `dict[(t,chain),(w,b)]`; else `list[(w,b)]` antichain | minimal UBWs of one ordered pair |
| table `T` | `list[list[cell]]`, `n×n` | `T[i][j]` = minimal `(i,j)`-UBWs; diagonal = cycles |
| `d` | length-`k` tuple (`Fraction` for `k=2`; `float` from LP) | threshold vector |
| `pi` | length-`n` list, same numeric type as `d` | embedding; `pi[0]=0` |

Helpers: `unit(k,j)`, `zero(k)`, `vec_add`, `neg`, `dot`, `INF = math.inf`.

## 3. Primitives

```text
beta_plus(a, k, u, v):                 # bound tuple or None
L1  if u < v:
L2      t ← a[u][v]
L3      return unit(k, t-1) if t ≥ 1 else None
L4  if u = v: return unit(k, k-1)
L5  return neg(beta_minus(a, k, v, u))

beta_minus(a, k, u, v):                # bound tuple or None
L1  if u < v:
L2      t ← a[u][v]
L3      return unit(k, t) if t ≤ k-1 else zero(k)
L4  if u = v: return zero(k)
L5  b ← beta_plus(a, k, v, u)
L6  return neg(b) if b ≠ None else None

walk_bound(a, k, w):                   # w is a UBW
L1  b ← zero(k)
L2  for i in 1 .. len(w)-1:
L3      b ← vec_add(b, beta_plus(a, k, w[i-1], w[i]))
L4  return b

preceq(x, y):                          # x ≼ y
L1  for j in 1 .. k:
L2      if sum(x[0..j-1]) > sum(y[0..j-1]): return False
L3  return True
```

## 4. Bound generation

`bound_generation(a)` builds `T`: for every ordered pair `(i,j)`, the `≼`-minimal UBWs from `i` to `j`, each stored as `(walk, bound)`. Diagonal cells hold cycles.

```text
bound_generation(a):
L1  n ← len(a); k ← a[0][0]
L2  T ← n×n table of empty cells
L3  for i in 0 .. n-1:                          # seed direct walks
L4      for j in i .. n-1:
L5          if a[i][j] ≠ 0: insert(T[i][j], (i,j), beta_plus(a,k,i,j))
L6          insert(T[j][i], (j,i), beta_plus(a,k,j,i))
L7  for s in 0 .. n-1:                          # sweep ALL ordered pairs
L8      for i in 0 .. n-1:
L9          for j in 0 .. n-1:
L10             for (w1,b1) in snapshot(T[i][s]):
L11                 for (w2,b2) in snapshot(T[s][j]):
L12                     w ← w1 + w2[1:]; b ← vec_add(b1,b2)
L13                     if minimal_in(b, T[i][j]):
L14                         insert_minimal(T[i][j], w, b)
L15 return T
```

Cell operations:
- `minimal_in(b, cell)`: no stored bound `c` with `preceq(c,b)` and `c≠b`.
- `insert_minimal(cell, w, b)`: add `(w,b)` (skip if `b` already present); remove every stored `(w',c)` with `preceq(b,c)` and `b≠c`.
- `k=2` cells are chain slots: `slot_key(b) = (|b[0]|+|b[1]|, 1 if b[1]>0 else 2)`; each slot keeps its `≼`-tightest element (within a slot all bounds are comparable), so insertion is O(1):

```text
insert_minimal_k2(cell, w, b):
L1  key ← slot_key(b)
L2  if key not in cell: cell[key] ← (w,b)
L3  elif preceq(b, cell[key][1]) and b ≠ cell[key][1]: cell[key] ← (w,b)
```

**Notes (why the code is shaped this way).**
- L7–L9 sweep **all** ordered pairs, including `i>j`. A down-direction cell `T[j][i]` cannot be filled by reversing `T[i][j]`: reversal turns a UBW into an LBW. Both directions must be generated directly. (Restricting to `i≤j` and filling the rest by reversal misses witnesses — e.g. a matrix whose only contradiction is a cycle bound `(0,0)`.)
- L10–L11 iterate over a **snapshot**: when `i=s` or `j=s` the cell is mutated in place during the sweep.
- L13–L14 keep an antichain. Pruning dominated bounds is safe because `≼` is additive: `x≼y ⇒ x+z≼y+z`, so any concatenation using a pruned bound is matched by one using the tighter bound.
- `T` stores walks (possibly with repeated vertices), not only paths. For a feasible `d`, the minimum bound value over walks equals that over paths, so the construction (§6) is unaffected.

## 5. Thresholds (the branch)

`cycle_bounds(T)` extracts the distinct minimal cycle bounds:

```text
cycle_bounds(T):
L1  seen ← empty map bound → walk
L2  for v in 0 .. n-1:
L3      for (w,b) in T[v][v]:
L4          if b not in seen: seen[b] ← simple_cycle(w)
L5  return [(b, seen[b]) for b in seen]
```

`simple_cycle(w)`: if the closed walk repeats a vertex, split it into simple cycles and keep the summand with the most extreme ratio bound (a sum's ratio bound lies between its summands').

### 5a. Candidate `d` supplied (any `k`) — verify

```text
solve_with_thresholds(a, d):
L1  T ← bound_generation(a)
L2  for v in 0 .. n-1:
L3      for (w,b) in T[v][v]:
L4          if dot(b,d) ≤ 0: return NO_SOLUTION, (simple_cycle(w),)
L5  return d, embedding_from_bounds(a, d, T)
```

### 5b. `k = 2` — combinatorial

Each cycle bound `(a1,a2)` constrains the ratio `r = d[1]/d[0]`: `a2>0 ⇒ r > -a1/a2`; `a2<0 ⇒ r < -a1/a2`; `a2=0 ⇒` require `a1>0`. Use exact rational arithmetic (`Fraction`).

```text
ratio_bounds(T):                       # k = 2
L1  lo ← 0; hi ← 1; w_lo ← None; w_hi ← None     # r ∈ (0,1)
L2  for v in 0 .. n-1:
L3      for (w,b) in T[v][v]:
L4          (a1,a2) ← b
L5          if a2 > 0: r ← Fraction(-a1,a2); if r > lo: lo ← r; w_lo ← w
L6          elif a2 < 0: r ← Fraction(-a1,a2); if r < hi: hi ← r; w_hi ← w
L7          else: if a1 ≤ 0: return INFEASIBLE, (simple_cycle(w), simple_cycle(w))
L8  if lo ≥ hi: return INFEASIBLE, (simple_cycle(w_lo), simple_cycle(w_hi))
L9  return FEASIBLE, (lo, hi)

solve_k2(a):
L1  T ← bound_generation(a)
L2  status, info ← ratio_bounds(T)
L3  if status = INFEASIBLE: return NO_SOLUTION, info
L4  (lo,hi) ← info; d ← (Fraction(1), (lo+hi)/2)
L5  return d, embedding_from_bounds(a, d, T)
```

Notes: L1 seeds the trivial window `(0,1)` (since `d[0]>d[1]>0`). L7 handles `a2=0` (the inequality is `a1·d[0]>0`); a bound `(0,0)` is immediately infeasible. (`k=1` is trivially feasible and may use the general path.)

### 5c. `k > 2` — LP feasibility

Feasibility of `{b·d > 0 for each minimal cycle bound b; d[0]>…>d[k-1]>0}` is a linear program. By homogeneity fix `d[0]=1` and maximize the minimum slack `t`; the strict system is feasible iff the optimum `t* > 0`.

```text
search_thresholds(cycles, k):
L1  variables x = (d[0..k-1], t)
L2  maximize t subject to:
L3      dot(b,d) - t ≥ 0            for each (b,_) in cycles
L4      d[i] - d[i+1] - t ≥ 0       for i in 0 .. k-2
L5      d[k-1] - t ≥ 0
L6      d[0] = 1
L7  if infeasible or optimal t ≤ 0: return None, farkas_witnesses(cycles)
L8  return d, None
```

`farkas_witnesses`: on infeasibility return ≤ `k` cycles whose inequalities are jointly contradictory (from the LP's Farkas certificate / dual multipliers; if unavailable, the cycle constraints active at the optimum). Solver call details: `SPEC-engineering.md` §4.

```text
solve_general(a):                    # k ≠ 2, no d supplied
L1  T ← bound_generation(a)
L2  d, witnesses ← search_thresholds(cycle_bounds(T), k)
L3  if d = None: return NO_SOLUTION, witnesses
L4  return d, embedding_from_bounds(a, d, T)
```

## 6. Construction

Given a feasible `d` and `T`, build `pi`. Lower bounds are read by reversal: the best lower bound on `(i,v)` is `-min{beta_plus(W)·d : W in T[v][i]}`.

```text
min_bound(cell, d):
L1  return min(dot(b,d) for (_,b) in cell)

embedding_from_bounds(a, d, T):
L1  n ← len(a); pi ← [0]*n
L2  for v in 1 .. n-1:
L3      ub ← min(pi[i] + min_bound(T[i][v], d) for i in 0..v-1 if T[i][v] nonempty); INF if none
L4      lb ← max(pi[i] - min_bound(T[v][i], d) for i in 0..v-1)
L5      pi[v] ← (ub+lb)/2 if ub < INF else lb + d[0]
L6  return pi
```

Notes: L3 skips empty cells (no upper bound exists for that pair); if every cell is empty, `ub=INF` and L5 places `pi[v]` just above `lb`. `T[v][i]` is always nonempty (the direct down-walk), so `lb` is finite. The cases `t=0` (no upper bound) and `t=k` (lower bound 0) are already encoded by `beta_plus` undefined / `beta_minus` zero — no sentinel array is stored. For feasible `d`, `lb < ub` and `pi` satisfies the uniform-embedding inequalities.

## 7. Top-level dispatch and output

```text
uniform_embedding_core(a, d=None):
L1  k ← a[0][0]
L2  if d supplied: return solve_with_thresholds(a, d)
L3  if k = 2: return solve_k2(a)
L4  return solve_general(a)
```

**Output.** Feasible: `(d, pi)`, `pi` strictly increasing, `d` scale-free (`d[0]=1` for the `k=2` and LP paths). Infeasible: `NO_SOLUTION` plus witness cycles — 1–2 for `k=2`, up to `k` for `k>2` — each a closed walk whose bound inequality fails for every threshold vector. The public wrapper (`SPEC-engineering.md`) maps vertex indices to external labels and converts numbers for display.