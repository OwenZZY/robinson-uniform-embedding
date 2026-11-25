# Paper overview — Uniform Embeddings for Robinson Similarity Matrices

Janssen–Zhang, arXiv:2105.09197v1 (May 2021). LaTeX source: `sources/arXiv-2105.09197v1/ArXiv_ver.tex`.
Companion files: `SPEC.md` (the data/interface contract), `docs/pseudocode.md` (Python-like transcription P1–P11 with conventions R1–R12), and `docs/definitions.md` (every term, one per line, with source lines). This overview is the reading map.

## The problem

Items $1,\dots,n$ carry pairwise similarities in $\{0,1,\dots,k\}$ ("$k$ = very similar", "$0$ = not at all"), forming a **Robinson matrix** $A$: symmetric, entries increasing toward the diagonal, diagonal all $k$. Seriation only recovers the order; this paper asks for positions on the line that respect the *values*: find thresholds $d_1>\dots>d_k>0$ and $\Pi:[n]\to\mathbb R$ such that similarity $t$ corresponds exactly to distance in $(d_{t+1},d_t]$ — a **uniform embedding** (Definition `def:uniform embedding`).

Equivalently: $A=\sum_{t=1}^k A^{(t)}$ decomposes into nested indifference (unit interval) graphs $G^{(t)}$; a uniform embedding is a *simultaneous* indifference embedding of the whole nested family. Not every Robinson matrix admits one (Example 1B).

## The three results

1. **Normalization** (`thm:monotone`, §2). If a uniform embedding exists, one exists that is strictly increasing with all inequalities strict: for $u<v$, $a_{u,v}=t \iff d_{t+1}<\Pi(v)-\Pi(u)<d_t$, where $d_{k+1}=0$, $d_0=\infty$. Everything after §2 uses this strict form.
2. **Existence criterion** (`thm:main theorem`, §4). Walks in the complete graph on $[n]$ generate linear bounds on $\Pi(v)-\Pi(u)$ with integer coefficient vectors in $\mathbb Z^k$ (§3). $A$ has a uniform embedding iff some threshold vector $\mathbf d$ makes every lower-bound-path value smaller than every upper-bound-path value, for every pair. Equivalently (§5.1): $\beta^+(C)^\top\mathbf d>0$ for every **upper-bound-cycle** $C$ — a finite system in the $k$ variables $d_i$ only.
3. **Algorithms** (§5). A Floyd–Warshall-style procedure (**FW-var**, Algorithm `alg:FW-mod`) enumerates the minimal bounds; minimality is under the prefix-sum order $\preceq$, which exactly captures domination for all $\mathbf d\in\mathbb D^k$ (`thm:bound_equivalence`). For **$k=2$** the cycle system collapses to lower/upper bounds on the single ratio $r=d_2/d_1$, giving a fully combinatorial decision-plus-construction algorithm (§5.3.1) that returns either $(\mathbf d,\Pi)$ or two witness cycles $C_1,C_2$ proving non-existence. For general $k$, feasibility of the cycle system can be handed to an LP solver; once a feasible $\mathbf d$ is known, the embedding is built by an explicit iterative formula (§4.2).

## The construction formula (§4.2)

Given feasible $\mathbf d$ and the bound tables: $\Pi(1)=0$, and for $v\ge2$,
$\Pi(v)=(\mathrm{ub}_v+\mathrm{lb}_v)/2$ where $\mathrm{ub}_v$ ($\mathrm{lb}_v$) is the smallest (largest) value of $\Pi(i)+{}$best upper (lower) bound on $(i,v)$, over $i<v$. The midpoint lands strictly inside every required interval. See `docs/pseudocode.md` P6 for the executable form, including the empty-upper-set convention.

## The $k=2$ pipeline (§5.3.1)

1. FW-var → minimal upper-bound walks for all ordered pairs; cycles live in the diagonal cells.
2. Each cycle bound $(a_1,a_2)$ gives $r>-a_1/a_2$ (if $a_2>0$) or $r<-a_1/a_2$ (if $a_2<0$); track the extreme bounds and their witness cycles.
3. Largest lower $\ge$ smallest upper → **NO SOLUTION**, return the two witnesses.
4. Otherwise pick $r$ strictly between (and in $(0,1)$), e.g. $d_1=1,d_2=r$.
5. Build $\Pi$ by the formula; return $(\mathbf d,\Pi)$.

The $k=2$ speed-up: $\mathbb Z^2_n$ partitions into $2n$ chains $S_t^1,S_t^2$ (`lemma:k2minimal`), so each pair's minimal-bound set is maintained in $O(1)$ per insertion via chain slots.

## Reading order for implementers

1. `docs/definitions.md` — notation and terms.
2. `SPEC.md` — data types (§2), indexing trap (§3), and the per-function input/output contract (§5).
3. `docs/pseudocode.md` — pseudocode P1–P11 and conventions R1–R12.
4. `SPEC-engineering.md` — public interface, translation, LP backend, naming, packaging.
5. `REPORT.md` — findings, especially F1 (a gap in the published FW-var pseudocode and the adopted fix).

## Section map of the source

| § | Content | tex lines |
|---|---|---|
| 1 | Introduction, Robinson matrices, seriation, nested-graph view, Example 1 (matrices $A$, $B$) | L43–128 |
| 2 | Uniform embedding definition; strict-monotone normalization | L131–159 |
| 3 | Bounds, edges/null-edges, walks, $\beta^\pm$, concatenation, reversal | L163–251 |
| 4 | Path condition; main theorem; cycles→paths; embedding formula | L257–428 |
| 5 | Cycle system; prefix order $\preceq$; FW-var; $k=2$ combinatorial algorithm | L433–594 |
| 6 | Conclusions; open question for $k\ge3$ | L598–602 |
| A | Proofs (monotonicity, concatenation, loop removal, bound equivalence) | L608–1260 |

## Reference fixtures (Example 1, L72–106)

- $A$ (5×5, embeddable): rows `(2,2,1,0,0) (2,2,2,1,1) (1,2,2,2,1) (0,1,2,2,2) (0,1,1,2,2)`; witness embedding $\Pi=(0,5,6.5,11.75,12.75)$, $d_1=8,d_2=6$.
- $B$ (6×6, not embeddable): rows `(2,2,1,0,0,0) (2,2,2,1,1,1) (1,2,2,2,1,1) (0,1,2,2,2,1) (0,1,1,2,2,2) (0,1,1,1,2,2)`; upper-bound-cycle $\langle1,2,6,4,1\rangle$ (1-based) has $\beta^+=(0,0)$, i.e. the impossible inequality $0>0$.

Both are verified against the definitions and serve as the first test cases for the implementation round.