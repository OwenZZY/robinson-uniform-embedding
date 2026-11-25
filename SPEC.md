# SPEC — the algorithm (mathematics)

Janssen–Zhang, *Uniform Embeddings for Robinson Similarity Matrices*, arXiv:2105.09197v1.

This document is the **mathematical pipeline**: what goes in, what comes out, and which function to read for each case. It is written to be followed by a mathematician. Every software concern — input formats, validation, the LP-solver mechanics, numerics, naming, packaging — lives in [`SPEC-engineering.md`](SPEC-engineering.md). Full notation: [`docs/definitions.md`](docs/definitions.md). Executable transcription: [`docs/pseudocode.md`](docs/pseudocode.md) (P1–P11, resolved conventions R1–R12).

## The pipeline

**Input.** A Robinson matrix $A\in\mathcal S^n[k]$; optionally a threshold vector $\mathbf d\in\mathbb D^k$.
**Output.** Either $(\mathbf d,\Pi)$ — thresholds and a uniform embedding — or witness cycles proving no uniform embedding exists.

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

Read the trunk top to bottom. At Stage 2, take the door matching your case; each door leads either to a feasible $\mathbf d$ (continue to Stage 3) or to witness cycles (stop: NO SOLUTION). Every function below is a page with a fixed **Input / Output / Paper / Pseudocode** header.

## Primitives used throughout

| Primitive | Meaning | Paper |
|---|---|---|
| $\beta^+(u,v)$ | For $u<v$: $\bm\chi_t$ if $a_{u,v}=t\ge1$; undefined if $a_{u,v}=0$. Reversed: $\beta^+(v,u)=-\beta^-(u,v)$. | def:bounds L186–205; def:beta vu L228 |
| $\beta^-(u,v)$ | For $u<v$: $\bm\chi_{t+1}$ if $a_{u,v}=t\le k-1$; $\mathbf0$ if $t=k$. Reversed: $\beta^-(v,u)=-\beta^+(u,v)$. | def:bounds L198–204 |
| $\beta^+(W)$ | $\sum_{i}\beta^+(w_{i-1},w_i)$ along an upper-bound-walk $W$. | eqn:beta(W) L233 |
| $\mathbf a\preceq\mathbf b$ | Prefix-sum order: $\sum_{i\le t}a_i\le\sum_{i\le t}b_i$ for all $t\in[k]$; equivalent to $\mathbf a^\top\mathbf d\le\mathbf b^\top\mathbf d$ for all $\mathbf d\in\mathbb D^k$. | def:boundleq L465; thm:bound_equivalence L472 |

Pseudocode P1–P2. An **upper-bound-walk** may traverse null-edges ($a=0$) only downward; a **lower-bound-walk** only upward (def:bound-walk L215). Reversal exchanges the two and negates the bound.

## Stage 1 — `bound_generation`: walk finding

- **Input.** $A$ (an $n\times n$ Robinson matrix, 0-based).
- **Output.** A table $T$: $T[i][j]$ is the set of $\preceq$-minimal $(i,j)$-upper-bound-walks, each stored with its bound $\beta^+(W)$. The diagonal cells $T[v][v]$ hold the upper-bound-cycles.
- **Paper.** Algorithm `alg:FW-mod` (L508–553), with the all-ordered-pairs variant (convention R5; see `REPORT.md` F1).
- **Pseudocode.** P3 (general $k$ antichains; $k=2$ chain slots P4).

Floyd–Warshall-style enumeration. Seed each pair with its direct walk (an up-step requires an edge; a down-step is always allowed). Sweep an intermediate vertex $s$ over all ordered pairs $(i,j)$, concatenate $W_1\in T[i,s]$ with $W_2\in T[s,j]$, and keep the concatenation only if its bound is $\preceq$-minimal, pruning any bounds it dominates (R3, R4). Minimality is on the bound vectors $\beta^+$, not on the walks. The table is $k$-generic.

## Stage 2 — choosing thresholds (the branch)

All three doors consume $T$ and produce either a feasible $\mathbf d$ (go to Stage 3) or witness cycles (NO SOLUTION).

### 2a — `verify_thresholds`: $\mathbf d$ supplied (any $k$)

- **Input.** $T$ and a candidate $\mathbf d\in\mathbb D^k$.
- **Output.** "feasible" (use this $\mathbf d$), or a witness cycle.
- **Paper.** Cycle system `eqn:cycle_system` (L450–453).
- **Pseudocode.** P8.

Check $\beta^+(C)^\top\mathbf d>0$ for every minimal cycle bound on the diagonal. If any is $\le0$, that cycle witnesses infeasibility of this $\mathbf d$. No search, no LP — the task is only verify-then-construct.

### 2b — `ratio_bounds`: $k\le2$, combinatorial

- **Input.** $T$ (with $k=2$).
- **Output.** A feasible $\mathbf d=(1,r)$, or NO SOLUTION with two witness cycles.
- **Paper.** §5.3 (L565–569); pipeline §5.3.1 (L574–583).
- **Pseudocode.** P5, P7.

For $k=2$ each cycle bound $(a_1,a_2)$ gives a constraint on the single ratio $r=d_2/d_1$:
$a_2>0\Rightarrow r>-a_1/a_2$ (lower); $a_2<0\Rightarrow r<-a_1/a_2$ (upper); $a_2=0\Rightarrow$ require $a_1>0$, else infeasible (R7). Intersect all bounds with the trivial window $(0,1)$ (R8). If the largest lower bound $\ge$ the smallest upper bound, NO SOLUTION, returning the two extremal cycles; otherwise pick $r$ strictly between (e.g. the midpoint) and set $\mathbf d=(1,r)$. ($k=1$ is trivially feasible for Robinson input.)

### 2c — LP search: $k>2$

- **Input.** The minimal cycle bounds extracted from $T$ (`cycle_bounds`, P9).
- **Output.** A feasible $\mathbf d$ (normalized $d_1=1$), or NO SOLUTION with $\le k$ Farkas witness cycles.
- **Paper.** §5 opening (L435); certificate (L558). A combinatorial test for $k\ge3$ is open (Conclusion L602).
- **Pseudocode.** P9–P11. Backend mechanics: `SPEC-engineering.md` §4.

Feasibility of the cycle system $\{\beta^+(C)^\top\mathbf d>0\ \forall C,\ d_1>\dots>d_k>0\}$ is a linear program in $k$ variables. By homogeneity, fix $d_1=1$ and maximize the minimum slack $t$:
$\beta^+(C)^\top\mathbf d\ge t$ per minimal cycle, $d_i-d_{i+1}\ge t$, $d_k\ge t$.
The strict system is feasible iff the optimum $t^*>0$ (R11); then $\mathbf d$ is feasible. If $t^*\le0$, the Farkas certificate gives $\le k$ cycles whose inequalities are jointly contradictory — the "bottleneck" of L558.

## Stage 3 — `embedding_from_bounds`: construction

- **Input.** $A$, a feasible $\mathbf d$, and $T$.
- **Output.** $\Pi:[n]\to\mathbb R$, a uniform embedding of $A$ with respect to $\mathbf d$.
- **Paper.** §4.2, `alg:embedding algorithm` (L343–348) and `equi:ubv and lbv` (L350–353).
- **Pseudocode.** P6.

Set $\Pi(1)=0$. For $v\ge2$, with lower bounds read from $T$ by reversal ($\operatorname{LowerBounds}(i,v)=-\beta^+(T[v][i])$):

$$ub_v=\min_{i<v}\Big(\Pi(i)+\min_{W\in T[i,v]}\beta^+(W)^\top\mathbf d\Big),\qquad
lb_v=\max_{i<v}\Big(\Pi(i)-\min_{W\in T[v,i]}\beta^+(W)^\top\mathbf d\Big),$$

$$\Pi(v)=\frac{ub_v+lb_v}{2}\quad\big(\text{or }lb_v+d_1\text{ if }ub_v=+\infty,\ \text{R6}\big).$$

For feasible $\mathbf d$, the lemmas at L357–426 guarantee $lb_v<ub_v$ and that $\Pi$ satisfies the strict system $d_{t+1}<\Pi(v)-\Pi(u)<d_t$ whenever $a_{u,v}=t$ ($u<v$). The construction is $k$-generic.

## Output

- **Feasible:** $(\mathbf d,\Pi)$. $\Pi$ is strictly increasing in matrix order; $\mathbf d$ is scale-free (the $k=2$ and LP paths normalize $d_1=1$).
- **Infeasible:** witness cycle(s). For $k=2$, two cycles $C_1,C_2$ (or one if a single cycle is already contradictory, e.g. bound $(0,0)$); for $k>2$, up to $k$ cycles from the Farkas certificate. Each is a closed walk whose bound inequality fails for every $\mathbf d\in\mathbb D^k$.