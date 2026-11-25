# SPEC — the algorithm (self-contained)

Janssen–Zhang, *Uniform Embeddings for Robinson Similarity Matrices*, arXiv:2105.09197v1.

**This is the single document the core implementer reads.** It specifies the data types (§1), the indexing convention (§2), the pipeline (§3), and every core function with its **explicit pseudocode inline** (§4–§7). Implementation tricks and traps are collected in §8. You should not need to open another file to implement the core.

Companions (reference only, not required to implement the core): [`SPEC-engineering.md`](SPEC-engineering.md) — public API, translation/validation, LP-solver mechanics, numerics policy, naming, packaging, errors, roadmap. [`docs/definitions.md`](docs/definitions.md) — full glossary with source lines. [`docs/overview.md`](docs/overview.md) — paper map.

Pseudocode blocks are labelled **P1–P11**; resolved conventions **R1–R12** are stated in §8 and referenced inline. `L###` = line in `sources/arXiv-2105.09197v1/ArXiv_ver.tex`.

## 1. Data types (the contract)

Internal representations are 0-based (§2). The core trusts these shapes; validation lives in engineering.

| Object | Concrete type | Mathematical meaning |
|---|---|---|
| `a` | `tuple[tuple[int, ...], ...]`, $n\times n$ | $A\in\mathcal S^n[k]$, Robinson, diagonal all $k$ |
| `n` | `int = len(a)` | number of items |
| `k` | `int = a[0][0]` | largest similarity level (= diagonal = max entry) |
| walk `w` | `tuple[int, ...]` of vertex indices | $W=\langle w_0,\dots,w_p\rangle$; concat `w1 + w2[1:]`, reverse `w[::-1]` |
| bound `b` | `tuple[int, ...]` of length $k$; `None` = undefined | $\beta\in\mathbb Z^k$ |
| cell | $k=2$: `dict[(t, chain), (w, b)]`; general $k$: `list[(w, b)]` antichain | $\preceq$-minimal upper-bound walks for one ordered pair |
| table `T` | `list[list[cell]]`, $n\times n$ | `T[i][j]` = minimal $(i,j)$-UBWs; diagonal = cycles |
| thresholds `d` | `tuple` of length $k$ (`Fraction` for $k=2$ auto; user type if given; `float` from LP) | $\mathbf d\in\mathbb D^k$ |
| positions `pi` | `list` of length $n$, same numeric type as `d` | $\Pi:[n]\to\mathbb R$, `pi[0]=0` |

Helpers: `unit(k, j)` = length-$k$ tuple with $1$ at **0-based** index $j$; `zero(k)` = all-zero tuple; `vec_add`, `neg`, `dot` = componentwise add / negate / dot product; `INF = math.inf`.

## 2. Indexing — the 1-based → 0-based trap

The paper uses $[n]=\{1,\dots,n\}$; Python is 0-based. **Every index shifts by one.** The unit-vector shift is the easiest off-by-one: $\beta^-$ uses $\chi_{t+1}$, which lands at index $t$, not $t+1$.

| Paper (1-based) | Code (0-based) | Note |
|---|---|---|
| vertex $u\in[n]$ | index `u-1` | `a[u][v]` reads $a_{u+1,v+1}$ |
| $\chi_t$ (1 at position $t$) | `unit(k, t-1)` | 1 at index $t-1$ |
| $\beta^-=\chi_{t+1}$ | `unit(k, t)` | index $(t+1)-1=t$ |
| $\Pi(1)=0$ | `pi[0]=0` | |
| $\min_{i\in[v-1]}$ | `for i in range(v)` | $i=0,\dots,v-1$ |
| $d_t$ | `d[t-1]` | threshold array index $t-1$ |
| $d_0=\infty$, $d_{k+1}=0$ | sentinels, **not** array entries | upper bound for $t=0$ is $\infty$; lower bound for $t=k$ is $0$ |
| witness $\langle1,2,6,4,1\rangle$ | `(0, 1, 5, 3, 0)` | subtract 1 from each vertex |

Check any implementation against the fixtures (Example 1A/1B, `docs/overview.md`); an off-by-one in $\chi$ shows up immediately as a wrong bound.

## 3. The pipeline

**Input.** $A\in\mathcal S^n[k]$; optionally $\mathbf d\in\mathbb D^k$. **Output.** $(\mathbf d,\Pi)$, or witness cycles proving none exists.

```
A ∈ S^n[k] ─► bound_generation ─► T          (minimal upper-bound walks; cycles on the diagonal)
                                    │
                 ┌──────────────────┼──────────────────┐
            d supplied?           k ≤ 2               k > 2
                 │                  │                    │
           verify d          ratio_bounds(T)       cycle_bounds(T) → LP
                 │                  │                    │
                 └────────► feasible d ◄─────────────────┘
                                    │
                  embedding_from_bounds(A, d, T) ─► Π
                                    │
                            (d, Π)  or  witness cycles
```

An **upper-bound-walk** (UBW) traverses null-edges ($a=0$) only downward; a **lower-bound-walk** (LBW) only upward (def:bound-walk L215). Reversal exchanges the two and negates the bound. Lower bounds are never stored — read them by reversal from `T[v][i]` (§7).

## 4. Primitives (P1–P2)

**$\beta^\pm$ edge bounds** (def:bounds L186–205; def:beta vu L228). For $u<v$: $\beta^+(u,v)=\chi_t$ if $a_{u,v}=t\ge1$, undefined if $0$; $\beta^-(u,v)=\chi_{t+1}$ if $t\le k-1$, $\mathbf0$ if $t=k$. Reversed pairs negate and swap. Diagonal convention R2.

```python
# P1
def beta_plus(a, k, u, v):            # β+(u,v); None = undefined
    if u < v:
        t = a[u][v]
        return unit(k, t - 1) if t >= 1 else None      # χ_t  → index t-1
    if u == v:
        return unit(k, k - 1)                          # R2: χ_k
    return neg(beta_minus(a, k, v, u))                 # β+(v,u) = −β−(u,v)

def beta_minus(a, k, u, v):           # β−(u,v); None = undefined
    if u < v:
        t = a[u][v]
        return unit(k, t) if t <= k - 1 else zero(k)   # χ_{t+1} → index t ; 0 if t=k
    if u == v:
        return zero(k)                                 # R2
    b = beta_plus(a, k, v, u)                          # β−(v,u) = −β+(u,v)
    return neg(b) if b is not None else None
```

**Walk bound and prefix order** (eqn:beta(W) L233; def:boundleq L465; thm:bound_equivalence L472). $\beta^+(W)=\sum_i\beta^+(w_{i-1},w_i)$; $\mathbf a\preceq\mathbf b$ iff every prefix sum of $\mathbf a$ is $\le$ that of $\mathbf b$ — equivalent to $\mathbf a^\top\mathbf d\le\mathbf b^\top\mathbf d$ for all $\mathbf d\in\mathbb D^k$.

```python
# P2
def walk_bound(a, k, w):              # β+(W), W a UBW
    return vec_sum(beta_plus(a, k, w[i - 1], w[i]) for i in range(1, len(w)))

def preceq(x, y):                     # x ≼ y
    return all(sum(x[:t]) <= sum(y[:t]) for t in range(1, k + 1))
```

## 5. Stage 1 — `bound_generation`: walk finding (P3–P4)

- **Input:** `a`. **Output:** `T` (§1). **Paper:** `alg:FW-mod` L508–553, all-ordered-pairs variant (R5; trap §8/F1).

Seed each pair with its direct walk (up-step needs an edge; down-step always allowed). Sweep intermediate $s$ over **all ordered pairs** $(i,j)$, concatenate, keep only $\preceq$-minimal bounds and prune dominated ones (R3, R4). Minimality is on $\beta^+$ vectors, not walks.

```python
# P3
def bound_generation(a):
    n, k = len(a), a[0][0]
    T = [[new_cell(k) for _ in range(n)] for _ in range(n)]
    # Init L528–535 (both directions):
    for i in range(n):
        for j in range(i, n):
            if a[i][j] != 0:
                insert(T[i][j], (i, j), beta_plus(a, k, i, j))
            insert(T[j][i], (j, i), beta_plus(a, k, j, i))   # j=i gives the loop (R2)
    # Sweeps L536–549; R5: j over ALL n; reversal update L546 dropped:
    for s in range(n):
        for i in range(n):
            for j in range(n):
                for (w1, b1) in snapshot(T[i][s]):     # snapshot: cell may update in place
                    for (w2, b2) in snapshot(T[s][j]): # when i==s or j==s (trap §8)
                        w = w1 + w2[1:]                # concatenation L237
                        b = vec_add(b1, b2)
                        if minimal_in(b, T[i][j]):     # no stored bound ≺ b (R3)
                            insert_minimal(T[i][j], w, b)   # add + prune ≻ b (R4)
    return T
```

`new_cell`, `insert`, `minimal_in`, `insert_minimal`: general $k$ = antichain list scanned with `preceq`; $k=2$ = chain slots (P4), $O(1)$ per insertion.

```python
# P4 — k=2 chain slots (lemma:k2minimal L487–500; rule L585; R9)
def slot_key(b):                      # b = (a1, a2)
    t = abs(b[0]) + abs(b[1])
    return (t, 1 if b[1] > 0 else 2)  # a2 = 0 → chain 2 (R9); t = 0 allowed

def insert_minimal_k2(cell, w, b):
    key = slot_key(b)
    if key not in cell:
        cell[key] = (w, b)
    elif preceq(b, cell[key][1]) and b != cell[key][1]:
        cell[key] = (w, b)            # same chain ⇒ comparable; keep tighter
    # else stored chain-minimum dominates b: discard
```

## 6. Stage 2 — choosing thresholds (the branch)

Each door consumes `T` and yields a feasible `d` (→ §7) or witness walks (NO SOLUTION).

### 6a. `verify_thresholds` — $\mathbf d$ supplied, any $k$ (P8)

- **Input:** `a`, `d`. **Output:** `(d, pi)` or `(NO_SOLUTION, (w,))`. **Paper:** cycle system L450–453.

```python
# P8
def solve_with_thresholds(a, d):
    T = bound_generation(a)
    for v in range(n):
        for (w, b) in T[v][v]:
            if dot(b, d) <= 0:                       # some cycle fails for this d
                return NO_SOLUTION, (simple_cycle(w),)
    return d, embedding_from_bounds(a, d, T)
```

### 6b. `ratio_bounds` — $k\le2$, combinatorial (P5, P7)

- **Input:** `T` ($k=2$). **Output:** feasible `d=(1,r)` or `(NO_SOLUTION, (w1,w2))`. **Paper:** §5.3 L565–569; §5.3.1 L574–583.

Each cycle bound $(a_1,a_2)$ constrains $r=d_2/d_1$: $a_2>0\Rightarrow r>-a_1/a_2$; $a_2<0\Rightarrow r<-a_1/a_2$; $a_2=0\Rightarrow$ need $a_1>0$ else infeasible (R7). Intersect with $(0,1)$ (R8).

```python
# P5
def ratio_bounds(T):                  # k = 2
    lo, hi = Fraction(0), Fraction(1)        # R8: trivial window r ∈ (0,1)
    w_lo, w_hi = None, None
    for v in range(n):
        for (w, b) in T[v][v]:
            a1, a2 = b
            if a2 > 0:
                r = Fraction(-a1, a2)
                if r > lo: lo, w_lo = r, w
            elif a2 < 0:
                r = Fraction(-a1, a2)
                if r < hi: hi, w_hi = r, w
            else:                            # a2 = 0 (R7)
                if a1 <= 0:
                    return INFEASIBLE, (simple_cycle(w), simple_cycle(w))
    if lo >= hi:
        return INFEASIBLE, (simple_cycle(w_lo), simple_cycle(w_hi))
    return FEASIBLE, (lo, hi)

# P7
def solve_k2(a):
    T = bound_generation(a)
    status, info = ratio_bounds(T)
    if status == INFEASIBLE:
        return NO_SOLUTION, info
    lo, hi = info
    d = (Fraction(1), (lo + hi) / 2)         # R8, R10
    return d, embedding_from_bounds(a, d, T)
```

($k=1$ is trivially feasible for Robinson input; handle via 6c or a direct trivial path.)

### 6c. `cycle_bounds` + LP search — $k>2$ (P9–P11)

- **Input:** minimal cycle bounds from `T`. **Output:** feasible `d` ($d_1=1$) or `(NO_SOLUTION, witnesses)` ($\le k$ cycles). **Paper:** §5 L435; certificate L558; combinatorial $k\ge3$ open (L602). **Backend mechanics:** `SPEC-engineering.md` §4.

Feasibility of $\{\beta^+(C)^\top\mathbf d>0\ \forall C,\ d_1>\dots>d_k>0\}$ is an LP in $k$ variables. By homogeneity fix $d_1=1$, maximize the minimum slack $t$; strict-feasible iff $t^*>0$ (R11).

```python
# P9
def cycle_bounds(T):                  # any k; dedup diagonal cycle bounds
    seen = {}
    for v in range(n):
        for (w, b) in T[v][v]:
            if b not in seen:
                seen[b] = simple_cycle(w)
    return [(b, seen[b]) for b in seen]

# P10 — in _lp.py; SciPy imported lazily (R12)
def search_thresholds(cycles, k):
    from scipy.optimize import linprog        # absent → ImportError (engineering §4)
    # x = (d_1..d_k, t); minimize -t; d_1 = 1; constraints ≥ 0:
    #   b·d - t ≥ 0  per cycle;  d_i - d_{i+1} - t ≥ 0;  d_k - t ≥ 0
    res = linprog(c=[0]*k + [-1], A_ub=..., b_ub=...,
                  bounds=[(0, 1)]*k + [(None, 1)],
                  constraints={"type": "eq", ... d_1 = 1 ...}, method="highs")
    if res.status == INFEASIBLE or res.x[-1] <= 0:
        return None, farkas_cycles(res, cycles)   # ≤ k witnesses (L558)
    return tuple(res.x[:k]), None

# P11 — general orchestration (wired in the public function, not in core)
def solve_general(a):                 # k != 2, no thresholds
    T = bound_generation(a)
    d, witnesses = search_thresholds(cycle_bounds(T), k)
    if d is None:
        return NO_SOLUTION, witnesses
    return d, embedding_from_bounds(a, d, T)
```

## 7. Stage 3 — `embedding_from_bounds`: construction (P6)

- **Input:** `a`, feasible `d`, `T`. **Output:** `pi`. **Paper:** §4.2 L343–353.

$\Pi(1)=0$; for $v\ge2$, with lower bounds read by reversal ($\operatorname{LowerBounds}(i,v)=-\beta^+(T[v][i])$):

$$\mathrm{ub}_v=\min_{i<v}\Big(\Pi(i)+\min_{W\in T[i,v]}\beta^+(W)^\top\mathbf d\Big),\qquad
\mathrm{lb}_v=\max_{i<v}\Big(\Pi(i)-\min_{W\in T[v,i]}\beta^+(W)^\top\mathbf d\Big),$$

$$\Pi(v)=\frac{\mathrm{ub}_v+\mathrm{lb}_v}{2}\quad\big(\text{or }\mathrm{lb}_v+d_1\text{ if }\mathrm{ub}_v=+\infty,\ \text{R6}\big).$$

```python
# P6
def min_bound(cell, d):
    return min(dot(b, d) for (_, b) in cell)

def embedding_from_bounds(a, d, T):
    n = len(a)
    pi = [0] * n                             # Π(1) = 0
    for v in range(1, n):
        ub = min((pi[i] + min_bound(T[i][v], d)
                  for i in range(v) if T[i][v]), default=INF)      # R6
        lb = max(pi[i] - min_bound(T[v][i], d) for i in range(v))  # T[v][i] ≠ ∅
        pi[v] = (ub + lb) / 2 if ub < INF else lb + d[0]           # R6
    return pi
```

For feasible $\mathbf d$, lemmas L357–426 give $\mathrm{lb}_v<\mathrm{ub}_v$ and $d_{t+1}<\Pi(v)-\Pi(u)<d_t$ whenever $a_{u,v}=t$ ($u<v$). $k$-generic.

## 8. Implementation tricks and traps (R1–R12)

The dangerous ones are marked **⚠**.

- **R1 ⚠ — Indexing.** 0-based internally; $\chi_t\to$ index $t-1$, $\chi_{t+1}\to$ index $t$; sentinels $d_0=\infty,d_{k+1}=0$ are not array entries. Full table in §2.
- **R2 — Loop init.** $\beta^+(i,i)=\chi_k$, $\beta^-(i,i)=\mathbf0$; the loop $\langle i,i\rangle$ seeds each diagonal cell. Concatenations that stutter at $i$ gain a $\chi_k\succeq\mathbf0$ term and are dominated/filtered.
- **R3 — Minimality on bounds.** "Minimal" (L544) is under $\preceq$ on $\beta^+$ vectors, not on walks (author note L554).
- **R4 — Antichain pruning.** Inserting a new minimum also deletes stored bounds it strictly dominates. Safe because $\preceq$ is additive ($\mathbf a\preceq\mathbf b\Rightarrow\mathbf a+\mathbf c\preceq\mathbf b+\mathbf c$).
- **R5 ⚠ — FW over all ordered pairs.** The published sweep targets only $i\le j$ and fills `T[j][i]` by reversal (L546), but reversal maps UBW→LBW (L241), so down-direction cells are under-generated and Example 1B's witness $(0,0)$ is missed. **Fix:** sweep all $(i,j)$, drop L546; init already covers both directions. See `REPORT.md` F1 for the trace.
- **R6 ⚠ — Empty upper set.** $\operatorname{UpperBounds}(i,v)$ can be empty (e.g. $n=2$, $a_{1,2}=0$). $\min\emptyset=+\infty$; if $\mathrm{ub}_v=+\infty$, set $\Pi(v)=\mathrm{lb}_v+d_1$. $\mathrm{lb}_v$ is always finite (down-walk $\langle v,i\rangle$ always present).
- **R7 ⚠ — $a_2=0$ cycle branch ($k=2$).** Inequality $a_1d_1>0$: vacuous if $a_1>0$, else NO SOLUTION. Required by Example 1B (bound $(0,0)$).
- **R8 — Ratio window.** Intersect cycle bounds with $(0,1)$ (since $\mathbf d\in\mathbb D^2$); $d_1=1$, $d_2=r$ midpoint.
- **R9 — Chain membership ($k=2$).** $t=|a_1|+|a_2|$, chain $1$ if $a_2>0$ else $2$; $t=0$ allowed. Cross-chain domination is not checked (as in the paper); a few extra non-minimal elements may remain — harmless.
- **R10 — Exact arithmetic.** $k=2$ uses `Fraction` end to end (ratio bounds rational → strict inequalities exact); `float` only at the public boundary.
- **R11 — LP strictness.** Homogeneous max-slack with $d_1=1$; strict-feasible iff optimum $t^*>0$; infeasible → Farkas certificate of $\le k$ cycles.
- **R12 — Optional backend.** SciPy imported lazily in `_lp.py`; core stays pure stdlib; $k=2$ and given-$\mathbf d$ modes run without SciPy.

Additional traps (not numbered conventions):

- **⚠ In-place snapshot.** When $i=s$ or $j=s$, `T[i][s]`/`T[s][j]` mutate during the sweep — iterate over a snapshot, not the live cell.
- **Walks vs paths.** `T` stores walks (may repeat vertices). For feasible $\mathbf d$, $\min$ over walks $=\min$ over paths (lemma:loopless L318), so the construction is unaffected. Diagonal cells may hold non-simple closed walks; `simple_cycle(w)` splits at a repeated vertex and keeps the extremal summand for a witness.
- **Lower bounds by reversal.** Never store LBWs; $\operatorname{LowerBounds}(i,v)=-\beta^+(T[v][i])$, hence the minus sign in $\mathrm{lb}_v$ (§7).
- **Strict convention.** From `thm:monotone` onward $d_{k+1}=0$ (not $-\infty$), $\Pi$ strictly increasing, all inequalities strict (L159). `\dist` is the signed $\Pi(v)-\Pi(u)$, not an absolute value.

## 9. Output and witnesses

- **Feasible:** $(\mathbf d,\Pi)$; $\Pi$ strictly increasing in matrix order; $\mathbf d$ scale-free ($d_1=1$ for $k=2$ and LP paths).
- **Infeasible:** witness cycles. $k=2$: two cycles $C_1,C_2$ (or one if a single cycle is contradictory, e.g. $(0,0)$); $k>2$: up to $k$ Farkas cycles. Each is a closed walk whose bound inequality fails for every $\mathbf d\in\mathbb D^k$. Returned to the user as label tuples via the engineering bijection.