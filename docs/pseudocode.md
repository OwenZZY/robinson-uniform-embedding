# Pseudocode transcription (Python-like)

Source: `sources/arXiv-2105.09197v1/ArXiv_ver.tex` (`L###` = line number there; labels are the paper's).
This file transcribes the paper's algorithms into Python-like pseudocode, line by line. It is **design documentation, not implementation**. Identifiers are provisional; final names are fixed in `SPEC.md`.

Indexing: the paper uses $[n]=\{1,\dots,n\}$; all pseudocode below is **0-based** ($\chi_t$ becomes a unit vector with index $t-1$). Walks are tuples of vertex indices; bounds are tuples of $k$ integers; `None` encodes "undefined" for $\beta^\pm$.

## Resolved conventions

The paper's pseudocode and formulas leave a few points implicit. Each is resolved here and referenced from the pseudocode as `R#`. (Typos and one substantive pseudocode gap are recorded in `REPORT.md`, findings bucket.)

- **R1 — Indexing.** 0-based internally; the public bijection maps external labels to `0..n-1`.
- **R2 — Diagonal / loop initialization.** $\beta^\pm(u,u)$ is not defined in the paper, but the init (L528–535) puts $\langle i,i\rangle$ into `UBW(i,i)` since $a_{i,i}=k\ne0$. Adopt $\beta^+(i,i)=\chi_k$, $\beta^-(i,i)=\mathbf0$ (the $u<v$ formulas applied at equality). The loop is a trivial upper-bound-cycle with bound $\chi_k\succ\mathbf0$; concatenations that stutter at $i$ gain a $\chi_k\succeq\mathbf0$ term and are dominated, hence filtered.
- **R3 — Minimality is on bound vectors.** "Minimal" in L544 means minimal under $\preceq$ applied to $\beta^+$, not to walks (confirmed by the author's note at L554).
- **R4 — Antichain maintenance.** When inserting a new minimal bound, also delete existing bounds it strictly dominates. $\preceq$ is additive ($\mathbf a\preceq\mathbf b\Rightarrow\mathbf a+\mathbf c\preceq\mathbf b+\mathbf c$), so pruning never loses a needed concatenation. The paper only tests the newcomer; pruning keeps the sets at the paper's claimed sizes.
- **R5 — Floyd–Warshall over all ordered pairs (deviation from L536–549).** As written, the sweep targets only $i\le j$ and fills `UBW(j,i)` with reversals (L546). But reversal maps upper-bound-walks to *lower*-bound-walks (L241), so down-direction cells receive walks that are not UBWs and are never concatenation targets. Concrete consequence: on Example 1B the witness cycle bound $(0,0)$ is never generated (trace in `REPORT.md`, F1). **Adopted fix:** loop $j$ over all $n$ and drop L546; the init (L528–535) already covers both directions. Lower bounds are still obtained by reversal, at read time (P6).
- **R6 — Empty upper-bound sets.** $\operatorname{UpperBounds}(i,v)$ can be empty (e.g. $a_{i,v}=0$ and no up-edge path exists; simplest case $n=2$, $a_{1,2}=0$). Take $\min\emptyset=+\infty$; if $ub_v=+\infty$, set $\Pi(v)=lb_v+d_1$ (any value $>lb_v$ is valid: no upper constraint involves $v$). $lb_v$ is always finite because the down-walk $\langle v,i\rangle$ is always in the table.
- **R7 — Cycle bounds with $a_2=0$ ($k=2$).** §5.3 splits only $a_2>0$ / $a_2<0$ (L565–567). If $a_2=0$ the inequality is $a_1d_1>0$, i.e. $a_1>0$: vacuous if $a_1>0$, otherwise NO SOLUTION with that cycle as witness. This branch is required by the paper's own Example 1B, whose witness bound is $(0,0)$.
- **R8 — Ratio window.** $r=d_2/d_1$ must lie in $(0,1)$ by the definition of $\mathbb D^2$, independently of any cycle. Intersect the cycle-derived window with $(0,1)$; if empty, NO SOLUTION. Scale-free choice $d_1=1$, $d_2=r$.
- **R9 — Chain membership boundary.** The rule of L585 assigns $a_2=0$ to chain 2. Although $(-t,0)$ is listed in $S_t^1$ (L495), it is $\preceq$-below every element of $S_t^2$, so the assignment is consistent for per-chain-minimum bookkeeping. Allow $t=0$ (the zero bound) as a slot.
- **R10 — Exact arithmetic.** All $k=2$ ratio bounds are rational ($-a_1/a_2$); use `fractions.Fraction` internally so every strict inequality is decided exactly. Convert to `float` only at the public boundary.

**Walks vs paths.** The table stores walks, possibly with repeated vertices (concatenation can create cycles). For a feasible $\mathbf d$, the loop-removal lemma (`lemma:loopless`, L318–323) gives $\min$ over walks $=\min$ over paths, so the embedding formula is unaffected. Diagonal cells may contain non-simple closed walks; their bounds are sums of cycle bounds, so the extra inequalities are harmless, and a witness can be reduced to a simple cycle (P5).

## P1 — Edge bounds (def:bounds L186–205; def:beta vu L228–230)

```python
def beta_plus(a, k, u, v):            # β+(u,v); None = undefined
    if u < v:
        t = a[u][v]
        return unit(k, t) if t >= 1 else None        # χ_t ; undefined on null-edge
    if u == v:
        return unit(k, k)                            # R2
    return neg(beta_minus(a, k, v, u))               # β+(v,u) = −β−(u,v)

def beta_minus(a, k, u, v):           # β−(u,v); None = undefined
    if u < v:
        t = a[u][v]
        return unit(k, t + 1) if t <= k - 1 else zero(k)   # χ_{t+1} ; 0 if t = k
    if u == v:
        return zero(k)                               # R2
    b = beta_plus(a, k, v, u)                        # β−(v,u) = −β+(u,v)
    return neg(b) if b is not None else None
```

Note: the sentence at L207 calls $\beta^+(u,v)$ a "lower bound"; the definition and all later use (e.g. L213) make it the **upper** bound. Transcribed accordingly.

## P2 — Walk bound and prefix order (eqn:beta(W) L233; def:boundleq L465–470)

```python
def walk_bound(a, k, w):              # β+(W) = Σ_i β+(w_{i-1}, w_i), W a UBW
    return vec_sum(beta_plus(a, k, w[i - 1], w[i]) for i in range(1, len(w)))

def preceq(x, y):                     # x ≼ y  iff every prefix sum of x ≤ that of y
    return all(sum(x[:t]) <= sum(y[:t]) for t in range(1, k + 1))
```

## P3 — Bound-Generation, FW-var (Algorithm alg:FW-mod, L508–553; adopted variant per R5)

```python
def bound_generation(a):
    n, k = len(a), a[0][0]                               # diagonal = k (assumed)
    T = [[new_cell() for _ in range(n)] for _ in range(n)]
    # T[i][j] holds minimal (walk, β+-bound) pairs for (i,j)-upper-bound-walks.
    # k = 2: cell = dict slot→(walk,bound) (P4); general k: cell = antichain list.

    # Init, L528–535: up-steps need an edge; down-steps always allowed.
    for i in range(n):
        for j in range(i, n):
            if a[i][j] != 0:
                insert(T[i][j], (i, j), beta_plus(a, k, i, j))
            insert(T[j][i], (j, i), beta_plus(a, k, j, i))   # j = i gives the loop (R2)

    # Sweeps, L536–549. R5: j runs over ALL n; the reversal update L546 is dropped.
    for s in range(n):
        for i in range(n):
            for j in range(n):
                for (w1, b1) in snapshot(T[i][s]):     # snapshot: cells may update
                    for (w2, b2) in snapshot(T[s][j]): # in place when i == s or j == s
                        w = w1 + w2[1:]                # concatenation, L237
                        b = vec_add(b1, b2)            # β+(W1+W2) = β+(W1)+β+(W2)
                        if minimal_in(b, T[i][j]):     # L544, on β+ vectors (R3)
                            insert_minimal(T[i][j], w, b)   # add + prune dominated (R4)
    return T
```

`minimal_in(b, cell)` = no stored bound is strictly $\prec b$. `insert_minimal` adds $(w,b)$ (deduplicating equal bounds) and removes stored bounds strictly $\succ b$ (R4). For $k=2$ both operations use chain slots (P4) instead of scanning.

## P4 — $k=2$ chain slots (lemma:k2minimal L487–500; membership rule L585; R9)

```python
def slot_key(b):                      # b = (a1, a2)
    t = abs(b[0]) + abs(b[1])
    chain = 1 if b[1] > 0 else 2      # L585; a2 = 0 → chain 2 (R9)
    return (t, chain)

def insert_minimal_k2(cell, w, b):    # cell: dict slot → (walk, bound)
    key = slot_key(b)
    if key not in cell:
        cell[key] = (w, b)
    elif preceq(b, cell[key][1]) and b != cell[key][1]:
        cell[key] = (w, b)            # same chain ⇒ comparable; keep the tighter
    # else the stored chain-minimum dominates b: discard.
```

Each chain holds at most one element, its $\preceq$-minimum; the $O(1)$ test of L585. Cross-chain domination is not checked (as in the paper); a few extra non-minimal elements may remain — harmless for both the cycle test and the embedding minima.

## P5 — Ratio bounds from cycles (§5.3 L565–569; steps 2–3 L579–580; R7, R8)

```python
def ratio_bounds(T):                  # k = 2; reads diagonal cells UBW(v,v)
    lo, hi = Fraction(0), Fraction(1)        # R8: the trivial window r ∈ (0,1)
    w_lo, w_hi = None, None
    for v in range(n):
        for (w, b) in T[v][v]:               # minimal closed walks at v
            a1, a2 = b
            if a2 > 0:                       # lower bound r > −a1/a2
                r = Fraction(-a1, a2)
                if r > lo: lo, w_lo = r, w
            elif a2 < 0:                     # upper bound r < −a1/a2
                r = Fraction(-a1, a2)
                if r < hi: hi, w_hi = r, w
            else:                            # a2 = 0 (R7): need a1 > 0
                if a1 <= 0:
                    return INFEASIBLE, (simple_cycle(w), simple_cycle(w))
    if lo >= hi:                             # step 3, L580
        return INFEASIBLE, (simple_cycle(w_lo), simple_cycle(w_hi))
    return FEASIBLE, (lo, hi)
```

`simple_cycle(w)`: if the stored closed walk repeats a vertex, split it into simple cycles and keep the summand with the extremal ratio bound (a sum's ratio bound is a mediant of its summands', so an extremum is attained at a simple cycle). Duplicates of the same cycle across different $v$ are harmless (L588).

## P6 — Embedding construction (alg:embedding algorithm L343–348; equi:ubv and lbv L350–353; R6)

```python
def min_bound(cell, d):
    return min(dot(b, d) for (_, b) in cell)

def embedding_from_bounds(a, d, T):
    n = len(a)
    pi = [0] * n                             # Π(1) = 0
    for v in range(1, n):
        # ub_v = min_{i<v} ( Π(i) + min{ b·d : b ∈ UpperBounds(i,v) } )
        ub = min((pi[i] + min_bound(T[i][v], d)
                  for i in range(v) if T[i][v]), default=INF)      # R6
        # LowerBounds(i,v) = { −β+(W) : W ∈ UBW(v,i) }  (reversal, L241):
        # lb_v = max_{i<v} ( Π(i) − min{ b·d : b ∈ T[v][i] } )
        lb = max(pi[i] - min_bound(T[v][i], d) for i in range(v))  # T[v][i] ≠ ∅
        pi[v] = (ub + lb) / 2 if ub < INF else lb + d[0]           # R6
    return pi
```

For feasible $\mathbf d$ the lemmas at L357–426 guarantee $lb_v<ub_v$ and that `pi` satisfies the strict system; the core trusts this and does not re-check.

## P7 — $k=2$ pipeline (§5.3.1, L574–583)

```python
def solve_k2(a):
    T = bound_generation(a)                          # step 1, L578
    status, info = ratio_bounds(T)                   # steps 2–3, L579–580
    if status == INFEASIBLE:
        return NO_SOLUTION, info                     # witness cycles C1, C2
    lo, hi = info
    r = (lo + hi) / 2                                # step 4, L581 (Fraction, R10)
    d = (Fraction(1), r)
    pi = embedding_from_bounds(a, d, T)              # step 5, L582
    return d, pi
```

## P8 — General $k$ with supplied thresholds (design hook; see SPEC, to-decide T1)

For any $k$, a *given* $\mathbf d\in\mathbb D^k$ can be tested combinatorially against the cycle system (no LP needed), then used for construction:

```python
def solve_with_thresholds(a, d):
    T = bound_generation(a)
    for v in range(n):
        for (w, b) in T[v][v]:
            if dot(b, d) <= 0:                       # cycle system, L450–453
                return NO_SOLUTION, (simple_cycle(w),)
    return d, embedding_from_bounds(a, d, T)
```

Finding a feasible $\mathbf d$ for $k>2$ (LP over the cycle inequalities, §5 L435) is deferred; see `SPEC.md`.