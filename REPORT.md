# REPORT — round 1 (documentation & design)

Consolidated report, kept current each round. Four buckets: **I know / Decided / To-decide / Findings**.

## Reading guide

1. This report's **Findings** first — F1 is a substantive gap in the published FW-var pseudocode, with a concrete trace on the paper's own Example 1B and the adopted fix.
2. `SPEC.md` — the v1 design and interface (the hand-off document for implementation).
3. `docs/pseudocode.md` — Python-like transcription P1–P8 and resolved conventions R1–R10.
4. `docs/definitions.md` — every mathematical term, one per line, with source lines.
5. `docs/overview.md` — the paper map for readers and other models.

## I know

- The task (from `spec0-by-hand.md`): analyze the paper, document all mathematical terms and Python-like pseudocode in `docs/`, then write the design and interface into `SPEC.md`. **No implementation code this round.** Done.
- The paper was read in full: `sources/arXiv-2105.09197v1/ArXiv_ver.tex` (1264 lines, including the Appendix) plus `preambles/alias.tex` for macro meanings. About 60 terms are recorded with line-precise references.
- The algorithmic content transcribed: edge bounds $\beta^\pm$, walk bounds, the prefix order $\preceq$, FW-var bound generation, the $k=2$ chain-slot trick and ratio test, the embedding formula, the $k=2$ pipeline, and a general-$k$ given-thresholds hook.
- The two Example 1 matrices are recorded as fixtures and hand-checked: $A$'s published embedding satisfies every pair inequality; $B$'s witness cycle $\langle1,2,6,4,1\rangle$ (upper-bound orientation) has bound $(0,0)$.
- Constraints honored: nothing in `sources/` modified; no proofs or complexity analysis added (the paper's stated complexities are recorded as statements only); no tests/README/examples built.
- Round 2 scope clarification (with the author): the problem is general $\mathcal S^n[k]$; the core (bound generation + construction) is $k$-generic and was never dropped. Only the combinatorial *threshold search* is $k=2$-specific; $k>2$ feasibility is an LP per the paper (L435), with combinatorial $k\ge3$ left open (L602). v1 now implements the full general-$k$ pipeline (D3, D11).

## Decided

- **D1 — Deliverables.** `docs/definitions.md`, `docs/pseudocode.md`, `docs/overview.md`; `SPEC.md` (v1, supersedes spec0); this report.
- **D2 — Architecture.** Three layers: core `_core.py` (paper algorithms, trusts input, zero try/except/raise, pure stdlib), translation `_translate.py` (validation predicates + immediate halt, format conversion, label bijection), outer (user notebook). No intermediate layer needed. Single public function `uniform_embedding(A, labels=None, thresholds=None)` returning a frozen `EmbeddingResult` (`feasible`, `positions`, `thresholds`, `cycles`).
- **D3 — v1 scope (revised round 2).** Complete for all of $\mathcal S^n[k]$. $k=2$: combinatorial algorithm, no dependencies. Any $k$ with supplied thresholds: combinatorial cycle-system check + construction. $k\ne2$ without thresholds: automatic search via LP feasibility over the cycle system, SciPy as an optional dependency; infeasibility returns Farkas witness cycles. The problem is general $k$; only the *combinatorial threshold search* is $k=2$-specific (ratio collapse, §5.3). For $k>2$ the paper itself routes feasibility through LP (L435) and leaves a combinatorial test for $k\ge3$ open (L602).
- **D4 — Identifiers.** Internal indices 0-based; external labels via an explicit bijection built in translation; default labels `range(n)`.
- **D5 — Numerics.** $k=2$ pipeline computes with `fractions.Fraction` (all ratio bounds are rational, so every strict inequality is exact); conversion to `float` only at the public boundary. Given-thresholds mode inherits the user's numeric type.
- **D6 — FW fix (from F1).** Bound generation sweeps all ordered pairs; the reversal update of tex L546 is dropped. Init already covers both directions.
- **D7 — Minimality bookkeeping.** Minimality on $\beta^+$ vectors under $\preceq$; insert-with-prune (antichain); $k=2$ cells are chain-slot dicts with the membership rule of L585.
- **D8 — Boundary conventions.** R2 (loop init $\beta^+(i,i)=\chi_k$), R6 ($ub_v=+\infty \Rightarrow \Pi(v)=lb_v+d_1$), R7 ($a_2=0$ cycle branch), R8 (ratio window intersected with $(0,1)$, $d_1=1$). All recorded in `docs/pseudocode.md`.
- **D9 — Naming.** Package `robinson_unif_emb`; public surface exactly `uniform_embedding` + `EmbeddingResult`; private modules and paper-flavored internal names (`pi`, `d`, `beta_plus`, `preceq`, `T`), mapping table in `SPEC.md` §6.
- **D10 — Version control.** Repo initialized and this round committed (per AGENTS.md "commit per round"); no push. `.gitignore` added for `.DS_Store`/`__pycache__`. `sources/` committed unmodified (see T3).
- **D11 — LP backend layering (round 2).** New optional module `_lp.py` holds the SciPy import (lazy) and `search_thresholds`; `_core.py` stays pure stdlib for every $k$ and gains `cycle_bounds`. The general-$k$ automatic path is wired in the public function (core → `_lp` → core), not inside core. LP formulation: homogeneous max-slack with $d_1=1$ (R11); strict-feasible iff optimum $t^*>0$.

## To-decide (defaults stated; skipping is safe)

- **T1 — RESOLVED (round 2, author's choice).** v1 includes the automatic general-$k$ threshold search via the cycle-system LP with a SciPy backend; $k=2$ stays the dependency-free fast path. See D3, D11.
- **T2 — Exact-output option?** Default: floats only in v1; an `exact=True` flag returning `Fraction`s can be added later if floating-point collapse of strict gaps ever matters in practice.
- **T3 — Track `sources/` in git?** Default: yes (your own arXiv source; nothing is pushed without your say-so). Alternative: keep it untracked/local.
- **T4 — Richer failure payload?** Default: witness cycles only (as the paper returns $C_1,C_2$). Alternative: also expose the two critical ratio values / bound vectors.
- **T5 — Tests, README, plotting helper?** Default: offer after the core works, build only on your go-ahead (per AGENTS.md).

## Findings

- **F1 — The published FW-var pseudocode under-generates down-direction cells; Example 1B's witness is missed.** In `alg:FW-mod`, the sweep targets only $i\le j$ (L538) and fills `UBW(j,i)` by inserting $(W_1+W_2)^\leftarrow$ (L546). But reversal maps upper-bound-walks to *lower*-bound-walks (L241), and down-direction cells are never concatenation targets. Concrete trace on Example 1B (1-based): the witness cycle $\langle1,2,6,4,1\rangle$ has $\beta^+=(0,0)$ and decomposes at $s=6$ into $\langle1,2,6\rangle$ (bound $(1,1)$) and $\langle6,4,1\rangle$ (bound $(-1,-1)$). The latter is a $(6,1)$-UBW, but it is not the init walk $\langle6,1\rangle$, and it is the reversal of $\langle1,4,6\rangle$, which is *not* a UBW (up-step $1\to4$ is a null-edge, $b_{1,4}=0$). Since $(6,1)$ is never a concatenation target, $\langle6,4,1\rangle$ never enters the table; the diagonal then yields only trivial ratio bounds (e.g. $(0,1)$, $(1,-1)$) and the algorithm would report feasibility for $B$. The same holds for every rotation of the witness. **Fix adopted (R5/D6):** sweep all ordered pairs $(i,j)$ and drop L546. Then $\langle6,4\rangle+\langle4,1\rangle$ produces $\langle6,4,1\rangle$ at $s=4$, and the $(0,0)$ cycle appears at $s=6$. The init (L528–535) is already correct for both directions, so the fix is local. The mathematical claims are unaffected; this is a pseudocode-level gap.
- **F2 — The $a_2=0$ cycle branch is missing from §5.3 but is needed by Example 1B.** The ratio test splits only $a_2>0$ / $a_2<0$ (L565–567). $B$'s witness bound is $(0,0)$: the inequality $0>0$ is immediately false. Resolved (R7): $a_2=0 \Rightarrow$ require $a_1>0$, else NO SOLUTION with that cycle as witness.
- **F3 — The embedding formula needs an empty-set convention.** $\operatorname{UpperBounds}(i,v)$ can be empty (simplest: $n=2$, $a_{1,2}=0$; the up-step is a null-edge and no other vertex exists). Then $ub_v=\min\emptyset$ and $(ub_v+lb_v)/2$ is undefined as written. Resolved (R6): $\min\emptyset=+\infty$; if $ub_v=+\infty$, set $\Pi(v)=lb_v+d_1$. Any value $>lb_v$ is valid since no upper constraint involves $v$; $lb_v$ is always finite (the down-walk $\langle v,i\rangle$ is always present).
- **F4 — Minimality in L544 is on bound vectors, not walks.** Confirmed by the author's own note at L554 ("line 10 should it be $\beta^+()$…? YEs"). Transcribed accordingly (R3); insert-with-prune is safe because $\preceq$ is additive (R4).
- **F5 — The table stores walks, not only paths.** Concatenation can create repeated vertices; off-diagonal cells may contain walks with cycles, and diagonal cells may contain non-simple closed walks. For feasible $\mathbf d$ this is harmless: the loop-removal lemma gives $\min$ over walks $=\min$ over paths, and extra closed-walk inequalities are sums of cycle inequalities. Witness walks are reduced to simple cycles when returned (P5).
- **F6 — The ratio window must include $(0,1)$.** $\mathbf d\in\mathbb D^2$ requires $0<d_2/d_1<1$ independently of any cycle inequality; the cycle system alone need not encode it. Resolved (R8): initialize the window to $(0,1)$ and intersect.
- **F7 — Typos resolved in transcription** (none affect the mathematics): L51 "asymmetric" → symmetric; L207 calls $\beta^+(u,v)$ a "lower bound" → upper (definition and L213 agree); L519 input "$\mathcal S^k$" → $\mathcal S^n[k]$; L579 "largest upper bound and the smallest upper bound" → largest *lower* bound and smallest upper bound (step 3 confirms); L435 "Theoremr".
- **F8 — The $d_{k+1}$ switch is deliberate, not a typo.** Definition uses $d_{k+1}=-\infty$; from `thm:monotone` onward $d_{k+1}=0$ with strict inequalities and $\Pi$ strictly increasing (explained at L159). The implementation uses the strict form throughout.
- **F9 — `\dist` is signed.** The macro expands to $\Pi(v)-\Pi(u)$ (alias.tex L206), not $|\Pi(v)-\Pi(u)|$; absolute values appear only in the original definition. All transcriptions use the signed form under the monotonicity convention.
- **F10 — Chain-membership boundary ($k=2$).** L585 assigns $a_2=0$ to chain 2, while $(-t,0)$ is listed in $S_t^1$ (L495). Consistent: $(-t,0)$ is $\preceq$-below every element of $S_t^2$, so per-chain-minimum bookkeeping is unaffected. The zero bound ($t=0$) is allowed as a slot (R9).