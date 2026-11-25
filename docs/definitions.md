# Definitions and mathematical terms

Source: Janssen–Zhang, *Uniform Embeddings for Robinson Similarity Matrices*, arXiv:2105.09197v1.
Reference convention: `L###` = line number in `sources/arXiv-2105.09197v1/ArXiv_ver.tex`; `label` = the paper's LaTeX label; `§` = section. Macro meanings come from `sources/arXiv-2105.09197v1/preambles/alias.tex`.
One term per line. Terms marked **[proof-only]** appear only in the Appendix and are not needed by the algorithm.

## Core objects

| Term | Definition | Source |
|---|---|---|
| Robinson matrix | Symmetric matrix $A=(a_{u,v})$ with $a_{u,v}\ge a_{u,w}$ and $a_{v,w}\ge a_{u,w}$ for all $u<v<w$; entries increase toward the diagonal in every row and column. | eq. `robinson:intro`, §1, L47–50 |
| $k$ | Largest similarity level; entries lie in $\{0,1,\dots,k\}$; all diagonal entries are assumed equal to $k$ throughout. | §1, L57–58 |
| $[t]$ | Index set $\{1,2,\dots,t\}$ for $t\in\mathbb Z_+$; $[n]$ is the set of items/vertices. | §2, L134 |
| $\mathcal S^n[k]$ | Set of all $n\times n$ Robinson matrices with entries in $\{0,1,\dots,k\}$. | §2, L134 |
| $\mathbb D^k$ (macro `\bmD`) | Set of threshold vectors $\mathbf d=(d_1,\dots,d_k)\in\mathbb R^k$ with $d_1>d_2>\dots>d_k>0$. | §2, L134 |
| $\Pi$ | Embedding map $\Pi:[n]\to\mathbb R$; $\Pi(v)$ is the position of item $v$ on the line. | def:uniform embedding, L137 |
| $\operatorname{dist}(u,v)$ (macro `\dist`) | Signed difference $\Pi(v)-\Pi(u)$ (not an absolute value); positive for $u<v$ under the standing monotonicity convention. | alias.tex L206 |
| Uniform embedding | For $A\in\mathcal S^n[k]$, $\mathbf d\in\mathbb D^k$: a map $\Pi$ with $a_{u,v}=t \iff d_{t+1}<|\Pi(v)-\Pi(u)|\le d_t$ for all pairs, $t\in\{0,\dots,k\}$, with $d_0=\infty$, $d_{k+1}=-\infty$. | def:uniform embedding, L136–142 |
| Standing convention (strict form) | From Theorem `thm:monotone` onward: $\Pi$ is strictly increasing and all inequalities are strict: for $u<v$, $a_{u,v}=t \iff d_{t+1}<\Pi(v)-\Pi(u)<d_t$, with $d_{k+1}=0$ (changed from $-\infty$) and $d_0=\infty$. | thm:monotone L147–154; convention stated L159 |
| $d_0$, $d_{k+1}$ | Sentinel thresholds: $d_0=\infty$ always; $d_{k+1}=-\infty$ in the original definition, $d_{k+1}=0$ under the strict convention (forces distinct positions). | L141, L153, L159 |
| Seriation problem | Find an ordering of items making the similarity matrix Robinson; this paper strengthens it to a metric embedding respecting the similarity values. | §1, L54–55 |
| Level matrix $A^{(t)}$ | Binary matrix with $a^{(t)}_{u,v}=1$ iff $a_{u,v}\ge t$; $A=\sum_{t=1}^k A^{(t)}$. | §1, L66–68 |
| $G^{(t)}$ | Indifference graph with adjacency matrix $A^{(t)}-I$; the family is nested: $G^{(t+1)}\subseteq G^{(t)}$. | §1, L68 |
| Indifference / unit interval / proper interval graph | Graph admitting a vertex embedding into $\mathbb R$ with adjacency iff distance $\le d$ for one threshold $d>0$; the three classes coincide. | §1, L60–62 |
| Graphon | Symmetric function $w:[0,1]^2\to[0,1]$; related-work context only (continuous-boundary results do not apply to matrices). | §1.1, L118 |

## Bounds, walks, paths, cycles

| Term | Definition | Source |
|---|---|---|
| $\bm\chi_i$ | Unit vector of $\mathbb Z^k$ with $1$ in position $i$, $0$ elsewhere. | §3, L184 |
| $\mathbf 0$ | All-zero vector of $\mathbb Z^k$; a lower bound on every $(u,v)$, $u<v$. | §3, L177 |
| Upper bound on $(u,v)$ | Vector $\mathbf b\in\mathbb Z^k$ such that $\Pi(v)-\Pi(u)<\mathbf b^\top\mathbf d$ is implied by the strict inequality system for every uniform embedding $\Pi$ and threshold vector $\mathbf d$. | def:bound, L169–173 |
| Lower bound on $(u,v)$ | Vector $\mathbf b\in\mathbb Z^k$ such that $\mathbf b^\top\mathbf d<\Pi(v)-\Pi(u)$ is implied by the strict system. | def:bound, L173 |
| Bound sign flip | If $\mathbf b$ is an upper bound on $(u,v)$ then $-\mathbf b$ is a lower bound on $(v,u)$, and vice versa. | §3, L178 |
| Edge | Pair $\{u,v\}$ with $a_{u,v}>0$. | §3, L213 |
| Null-edge | Pair $\{u,v\}$ with $a_{u,v}=0$. | §3, L213 |
| $\beta^+(u,v)$ (edge upper bound) | For $u<v$: $\bm\chi_t$ if $a_{u,v}=t\ge1$; undefined if $a_{u,v}=0$. | def:bounds, L190–196 |
| $\beta^-(u,v)$ (edge lower bound) | For $u<v$: $\bm\chi_{t+1}$ if $a_{u,v}=t\le k-1$; $\mathbf 0$ if $a_{u,v}=k$. | def:bounds, L198–204 |
| $\beta^\pm$ on reversed pairs | For $u<v$: $\beta^+(v,u)=-\beta^-(u,v)$ and $\beta^-(v,u)=-\beta^+(u,v)$. | eq. `def:beta vu`, L228–230 |
| $(u,v)$-walk | Sequence $W=\langle w_0,\dots,w_p\rangle$ in the complete graph on $[n]$ with $w_0=u$, $w_p=v$; repetitions allowed. | def:bound-walk, L215–216 |
| Upper-bound-walk (UBW) | Walk whose every up-step ($w_{i-1}<w_i$) is an edge; down-steps may be edges or null-edges (null-edges only traversed downward). | def:bound-walk, L216–224 |
| Lower-bound-walk (LBW) | Walk whose every down-step is an edge; up-steps may be edges or null-edges (null-edges only traversed upward). | def:bound-walk, L225 |
| $\beta^+(W)$, $\beta^-(W)$ (walk bounds) | $\beta^+(W)=\sum_{i=1}^p\beta^+(w_{i-1},w_i)$ for a UBW; $\beta^-(W)=\sum_{i=1}^p\beta^-(w_{i-1},w_i)$ for a LBW. | eqs. `eqn:beta(W)`, `eqn:betamin(W)`, L231–235 |
| Concatenation $W_1+W_2$ | For $W_1=\langle u_0,\dots,u_s\rangle$, $W_2=\langle u_s,\dots,u_p\rangle$: the walk $\langle u_0,\dots,u_p\rangle$ (junction vertex appears once). | def:bound-walk, L237 |
| Reverse $W^\leftarrow$ | $\langle w_p,\dots,w_0\rangle$; reverses UBW↔LBW; $\beta^+(W)=-\beta^-(W^\leftarrow)$. | §3, L240–241 |
| Concatenation lemma | Walk bounds are valid bounds: $\beta^+(W)$ upper on $(u,v)$ for a $(u,v)$-UBW; $\beta^-(W)$ lower for a LBW. | lemma:concatenation, L244–248 |
| $(u,v)$-upper/lower-bound-path | A UBW/LBW with no repeated vertices; sets $\mathcal U_{u,v}$, $\mathcal L_{u,v}$ (finite, unlike walks). | §4, L261–262 |
| $\operatorname{UpperBounds}(u,v)$, $\operatorname{LowerBounds}(u,v)$ | Sets of bound vectors $\{\beta^+(W):W\in\mathcal U_{u,v}\}$, $\{\beta^-(W):W\in\mathcal L_{u,v}\}$. | §4.2, L342 |
| Upper-bound-cycle | $(u,u)$-UBW with no repeated vertex other than $u=v$; traversal direction matters (the reverse is a lower-bound-cycle). | §4.1, L294 |
| $\mathcal C$ | Set of all upper-bound-cycles. | §5.1, L450 |
| $\mathbb Z^k_n$ | $\{\mathbf a\in\mathbb Z^k:\sum_{i=1}^k|a_i|\le n\}$; contains every cycle bound (cycle length $\le n$, each step contributes $\pm\bm\chi_i$). | §5.1, L454–457 |
| $M$ | Number of minimal elements of $\mathbb Z^k_n$ under $\preceq$ (used in the paper's complexity statements). | §5.2, L556 |
| $N$ | Input size $N=n(n-1)/2$, the number of upper-diagonal entries. | §5.3, L594 |

## Conditions and theorems (statements only)

| Term | Definition | Source |
|---|---|---|
| Path condition (system `eqn:condition_d`) | For all $u<v$, all $\mathbf b\in\operatorname{UpperBounds}(u,v)$, $\mathbf a\in\operatorname{LowerBounds}(u,v)$: $\mathbf a^\top\mathbf d<\mathbf b^\top\mathbf d$. | §4, L264–268 |
| Main theorem | $A\in\mathcal S^n[k]$ has a uniform embedding iff some $\mathbf d\in\mathbb D^k$ satisfies the path condition. | thm:main theorem, L271–273 |
| Cycle system (system `eqn:cycle_system`) | For all $C\in\mathcal C$: $\beta^+(C)^\top\mathbf d>0$; equivalent to the path condition. | §5.1, L449–453 |
| Cycle-positivity lemma | If $\mathbf d$ satisfies the path condition then $\beta^+(C)^\top\mathbf d>0$ for every upper-bound-cycle $C$. | lemma:ub cycle geq 0, L297–301 |
| Loop-removal lemma | Under a feasible $\mathbf d$, every walk bound is matched or tightened by a path bound (remove cycles iteratively). | lemma:loopless, L318–323 |
| $\preceq$ (prefix-sum order) | $\mathbf a\preceq\mathbf b$ iff $\sum_{i=1}^t a_i\le\sum_{i=1}^t b_i$ for all $t\in[k]$. | def:boundleq, L465–470 |
| Tighter / wider | $\mathbf a\preceq\mathbf b$: $\mathbf a$ is tighter, $\mathbf b$ wider; otherwise incomparable (terminology from the commented source). | commented L469 |
| Bound-equivalence theorem | $\mathbf a\preceq\mathbf b \iff \mathbf a^\top\mathbf d\le\mathbf b^\top\mathbf d$ for all $\mathbf d\in\mathbb D^k$; hence $\preceq$ is a partial order. | thm:bound_equivalence, L472–475 |
| Minimal-element reduction | It suffices that the cycle system hold for the $\preceq$-minimal cycle bounds (an antichain). | corollary, L481–483 |
| $k=2$ size lemma | For $k=2$ the minimal-bound cycle system has at most $2n$ inequalities. | lemma:k2minimal, L487–489 |
| Chains $S_t$, $S_t^1$, $S_t^2$ | $S_t=\{\mathbf a\in\mathbb Z^2:|a_1|+|a_2|=t\}$; $S_t^1=\{(-t+i,i):0\le i\le t\}\cup\{(i,t-i):1\le i<t\}$; $S_t^2=\{(-t+i,-i):1\le i\le t\}\cup\{(i,-t+i):1\le i\le t\}$; each is a $\preceq$-chain, and they partition $S_t$. | lemma:k2minimal proof, L493–499 |
| Chain membership rule ($k=2$) | Bound $(a_1,a_2)$ lies in chain $S_t^i$ with $t=|a_1|+|a_2|$, $i=1$ if $a_2>0$, $i=2$ otherwise. | §5.3, L585 |
| Ratio test ($k=2$) | Cycle inequality $a_1d_1+a_2d_2>0$ gives, for $r=d_2/d_1$: lower bound $r>-a_1/a_2$ if $a_2>0$, upper bound $r<-a_1/a_2$ if $a_2<0$; feasible iff largest lower < smallest upper. | §5.3, L565–567 |

## Algorithmic quantities

| Term | Definition | Source |
|---|---|---|
| Bound-Generation / FW-var | Floyd–Warshall-style enumeration of minimal upper-bound walks between all pairs; lower bounds obtained by reversal. | Algorithm `alg:FW-mod`, L508–553 |
| UBW table | Lookup $\operatorname{UBW}(i,j)$ = minimal elements (under $\preceq$ on $\beta^+$) of upper-bound walks/paths from $i$ to $j$. | alg:FW-mod output, L520–524 |
| $\mathrm{ub}_v$ | $\mathrm{ub}_v=\min_{i\in[v-1]}\big(\Pi(i)+\min\{\mathbf b^\top\mathbf d:\mathbf b\in\operatorname{UpperBounds}(i,v)\}\big)$. | eq. `equi:ubv and lbv`, L350–353 |
| $\mathrm{lb}_v$ | $\mathrm{lb}_v=\max_{i\in[v-1]}\big(\Pi(i)+\max\{\mathbf a^\top\mathbf d:\mathbf a\in\operatorname{LowerBounds}(i,v)\}\big)$. | eq. `equi:ubv and lbv`, L350–353 |
| Embedding formula | $\Pi(1)=0$; $\Pi(v)=(\mathrm{ub}_v+\mathrm{lb}_v)/2$ for $2\le v\le n$; yields a uniform embedding for any feasible $\mathbf d$. | eq. `alg:embedding algorithm`, L343–348; lemmas L357–426 |
| Uniform Embedding Algorithm ($k=2$) | Steps: (1) generate minimal bounds, extract cycles from $\operatorname{UBW}(v,v)$; (2) convert each cycle bound to a ratio bound on $d_2/d_1$, keep witnesses $C_1,C_2$; (3) if largest lower $\ge$ smallest upper: NO SOLUTION, return $C_1,C_2$; (4) else pick $d_1,d_2$ with ratio strictly between; (5) build $\Pi$ by the embedding formula, return it. | §5.3.1, L574–583 |
| Stated complexities | FW-var $O(n^3M^3)$ general; $k=2$ bound generation $O(n^5)$ via chain slots, cycle inequalities $O(n^2)$, embedding $O(n^3)$, overall $O(N^{2.5})$; LP route $O(n^4)$ (Vaidya). Recorded as stated, not re-derived. | L556, L585–594 |

## Appendix / proof-only terms

| Term | Definition | Source |
|---|---|---|
| Repeating rows | Rows $i,j$ with $a_{i,w}=a_{j,w}$ for all $w$; handled by the monotonicity theorem's proof. **[proof-only]** | lemma:no repeating, L617–618 |
| $A[I]$, index set $I$ | Induced submatrix keeping the first appearance of each distinct row. **[proof-only]** | L769–771 |
| $\epsilon_i$ | Half the minimum slack of $\Pi_0$ at row $i$; used to space repeating rows. **[proof-only]** | eq. `eqn:thm:monotone`, L775–781 |
| $[x]_+$ | Positive part: $x$ if $x>0$, else $0$. **[proof-only]** | commented L992–999, used L1008 |
| Buffer vector $\mathbf c$ | Vector constructed from $\mathbf a\preceq\mathbf b$ with $\mathbf a^\top\mathbf d\le\mathbf c^\top\mathbf d\le\mathbf b^\top\mathbf d$; proof device for bound-equivalence. **[proof-only]** | def:buffer vector c, L1002–1017 |
| $e_t$, $e_{t,i}$, $f_{t,i}$, $c_{t,i}$ | Auxiliary redistribution sequences in the buffer-vector construction. **[proof-only]** | L1007–1016 |

## Worked examples (fixtures)

| Term | Definition | Source |
|---|---|---|
| Example 1 matrix $A$ (embeddable) | $5\times5$, rows: $(2,2,1,0,0),(2,2,2,1,1),(1,2,2,2,1),(0,1,2,2,2),(0,1,1,2,2)$; uniform embedding $\Pi=(0,5,6.5,11.75,12.75)$ with $d_1=8,d_2=6$. | example:counter-example, L72–98 |
| Example 1 matrix $B$ (not embeddable) | $6\times6$, rows: $(2,2,1,0,0,0),(2,2,2,1,1,1),(1,2,2,2,1,1),(0,1,2,2,2,1),(0,1,1,2,2,2),(0,1,1,1,2,2)$; contradiction from cycle $\langle1,4,6,2,1\rangle$ (lower-bound orientation); the upper-bound orientation is $\langle1,2,6,4,1\rangle$ with $\beta^+=(0,0)$. | example:counter-example, L83–104 |