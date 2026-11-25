# SPEC — the algorithm (data and interface contract)

Janssen–Zhang, *Uniform Embeddings for Robinson Similarity Matrices*, arXiv:2105.09197v1.

This spec states the **data, data types, and input/output** of the algorithm — the contract an implementer builds to. It contains no code. The step-by-step **pseudocode is in [`docs/pseudocode.md`](docs/pseudocode.md)** (blocks P1–P11, conventions R1–R12); each function below points to its block. Software concerns (public API, translation, LP-solver mechanics, naming, packaging) are in [`SPEC-engineering.md`](SPEC-engineering.md). Glossary: [`docs/definitions.md`](docs/definitions.md). `L###` = line in `sources/arXiv-2105.09197v1/ArXiv_ver.tex`.

## 1. Top-level input and output

| | Data | Type |
|---|---|---|
| **Input** | `A`, a Robinson similarity matrix | $n\times n$ integers in $\{0,\dots,k\}$, diagonal all $k$; internally `tuple[tuple[int,...],...]` |
| Input (optional) | `thresholds` $\mathbf d$ | length-$k$ tuple, $d_1>\dots>d_k>0$ |
| Input (optional) | `labels` | $n$ unique hashables (external identifiers); default `range(n)` |
| **Output (feasible)** | $(\mathbf d,\Pi)$ | thresholds: length-$k$ tuple; embedding $\Pi$: length-$n$ list of positions, `pi[0]=0` |
| **Output (infeasible)** | witness cycles | 1–2 closed walks ($k=2$) or $\le k$ ($k>2$), each a tuple of vertex indices (first = last) |

## 2. Data types (the contract)

Internal representations are 0-based (§3). The core trusts these shapes; validation is in engineering.

| Object | Concrete type | Mathematical meaning |
|---|---|---|
| `a` | `tuple[tuple[int, ...], ...]`, $n\times n$ | $A\in\mathcal S^n[k]$ |
| `n`, `k` | `int = len(a)`, `int = a[0][0]` | size; largest level (= diagonal = max) |
| walk `w` | `tuple[int, ...]` | $W=\langle w_0,\dots,w_p\rangle$; concat `w1+w2[1:]`, reverse `w[::-1]` |
| bound `b` | `tuple[int, ...]` length $k$; `None` = undefined | $\beta\in\mathbb Z^k$ |
| cell | $k=2$: `dict[(t,chain),(w,b)]`; general: `list[(w,b)]` antichain | $\preceq$-minimal upper-bound walks of one ordered pair |
| table `T` | `list[list[cell]]`, $n\times n$ | `T[i][j]` = minimal $(i,j)$-UBWs; diagonal = cycles |
| thresholds `d` | length-$k$ tuple (`Fraction` for $k=2$ auto; user type if given; `float` from LP) | $\mathbf d\in\mathbb D^k$ |
| positions `pi` | length-$n$ list, same numeric type as `d` | $\Pi:[n]\to\mathbb R$ |

Helpers: `unit(k,j)` = length-$k$ tuple, $1$ at 0-based index $j$; `zero(k)`; `vec_add`, `neg`, `dot`; `INF = math.inf`.

## 3. Indexing convention (the 1-based → 0-based trap)

The paper is 1-based ($[n]=\{1,\dots,n\}$); Python is 0-based. Every index shifts by one; the unit-vector shift is the easiest off-by-one.

| Paper | Code | Note |
|---|---|---|
| vertex $u$ | index `u-1` | `a[u][v]` reads $a_{u+1,v+1}$ |
| $\chi_t$ | `unit(k, t-1)` | $\beta^-$ uses $\chi_{t+1}$ → `unit(k, t)` |
| $\Pi(1)=0$ | `pi[0]=0` | |
| $\min_{i\in[v-1]}$ | `range(v)` | |
| $d_t$ | `d[t-1]` | $d_0=\infty$, $d_{k+1}=0$ are sentinels, not array entries |
| $\langle1,2,6,4,1\rangle$ | `(0,1,5,3,0)` | subtract 1 per vertex |

## 4. Pipeline

```
A ─► bound_generation ─► T
                          │
        ┌─────────────────┼─────────────────┐
   d supplied?          k ≤ 2              k > 2
        │                 │                   │
   verify d        ratio_bounds(T)      cycle_bounds(T) → LP
        │                 │                   │
        └────────► feasible d ◄───────────────┘
                          │
            embedding_from_bounds(A, d, T) ─► Π
                          │
                   (d, Π)  or  witness cycles
```

Upper-bound-walks traverse null-edges only downward; lower-bound-walks only upward. Lower bounds are read by reversal from `T[v][i]`, never stored.

## 5. Function input/output contract

Types are as in §2. "Pseudocode" links to the block in `docs/pseudocode.md`.

| Function | Input | Output | Paper | Pseudocode |
|---|---|---|---|---|
| `beta_plus` | `a, k, u, v` | bound tuple or `None` | def:bounds L186 | P1 |
| `beta_minus` | `a, k, u, v` | bound tuple or `None` | def:bounds L198 | P1 |
| `walk_bound` | `a, k, w` (walk) | bound tuple | eqn:beta(W) L233 | P2 |
| `preceq` | `x, y` (bound tuples) | `bool` ($\mathbf x\preceq\mathbf y$) | def:boundleq L465 | P2 |
| `bound_generation` | `a` | `T` (table) | alg:FW-mod L508 | P3–P4 |
| `cycle_bounds` | `T` | `list[(bound, walk)]`, distinct minimal cycle bounds | step 1 L578 | P9 |
| `ratio_bounds` ($k=2$) | `T` | `(FEASIBLE, (lo, hi))` or `(INFEASIBLE, (w1, w2))` | §5.3 L565 | P5 |
| `search_thresholds` ($k>2$) | `cycles` (from `cycle_bounds`), `k` | `(d, None)` or `(None, witnesses)` | §5 L435 | P10 |
| `embedding_from_bounds` | `a, d, T` | `pi` (positions) | §4.2 L343 | P6 |
| `solve_k2` | `a` | `(d, pi)` or `(NO_SOLUTION, (w1, w2))` | §5.3.1 L574 | P7 |
| `solve_with_thresholds` | `a, d` | `(d, pi)` or `(NO_SOLUTION, (w,))` | cycle system L450 | P8 |
| `solve_general` ($k\ne2$) | `a` | `(d, pi)` or `(NO_SOLUTION, witnesses)` | §5 | P11 |

**Construction rule** (the only formula needed to read the contract; full form in P6): $\Pi(1)=0$ and, for $v\ge2$, $\Pi(v)=(\mathrm{ub}_v+\mathrm{lb}_v)/2$, where $\mathrm{ub}_v$ is the smallest $\Pi(i)+{}$best upper bound on $(i,v)$ and $\mathrm{lb}_v$ the largest $\Pi(i)+{}$best lower bound, over $i<v$. If no upper bound exists, $\mathrm{ub}_v=+\infty$ and $\Pi(v)=\mathrm{lb}_v+d_1$ (R6).

**Threshold rules.** $k=2$: each cycle bound $(a_1,a_2)$ gives a lower ($a_2>0$) or upper ($a_2<0$) bound on $r=d_2/d_1$; $a_2=0$ requires $a_1>0$ (R7); intersect with $(0,1)$ (R8); feasible iff largest lower < smallest upper. $k>2$: LP feasibility of $\beta^+(C)^\top\mathbf d>0$ over minimal cycles, $d_1>\dots>d_k>0$ (R11).

## 6. Tricks and traps (summary)

Full statements and pseudocode in `docs/pseudocode.md` §"Resolved conventions". ⚠ = easy to get wrong.

- **R1 ⚠** Indexing 0-based; $\chi_t\to$ index $t-1$ (§3).
- **R2** Loop init $\beta^+(i,i)=\chi_k$, $\beta^-(i,i)=\mathbf0$.
- **R3** Minimality is on $\beta^+$ vectors under $\preceq$, not on walks.
- **R4** Insert-with-prune: delete stored bounds the new minimum dominates ($\preceq$ is additive).
- **R5 ⚠** FW sweeps **all ordered pairs**; the published reversal update (L546) under-generates and misses Example 1B's witness — drop it (see `REPORT.md` F1).
- **R6 ⚠** Empty upper set: $\min\emptyset=+\infty$; then $\Pi(v)=\mathrm{lb}_v+d_1$.
- **R7 ⚠** $k=2$ cycle bound with $a_2=0$: infeasible iff $a_1\le0$ (Example 1B is $(0,0)$).
- **R8** Ratio window intersected with $(0,1)$; $d_1=1$.
- **R9** $k=2$ chain slots: $t=|a_1|+|a_2|$, chain 1 if $a_2>0$ else 2; $t=0$ allowed.
- **R10** $k=2$ uses `Fraction` (exact); `float` only at the public boundary.
- **R11** LP strictness via homogeneous max-slack, $d_1=1$; feasible iff optimum $t^*>0$.
- **R12** SciPy optional, lazy import in `_lp.py`; core pure stdlib.
- **⚠ Snapshot** cells before iterating when $i=s$ or $j=s$ (in-place updates).
- **Walks vs paths:** `T` stores walks; for feasible $\mathbf d$, $\min$ over walks $=\min$ over paths; reduce witnesses to simple cycles.

## 7. Output and witnesses

Feasible: $(\mathbf d,\Pi)$, $\Pi$ strictly increasing in matrix order, $\mathbf d$ scale-free ($d_1=1$ for $k=2$/LP). Infeasible: witness cycles whose bound inequality fails for every $\mathbf d\in\mathbb D^k$; mapped to external labels by the engineering bijection.