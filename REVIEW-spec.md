# Spec readiness review — round 7 (review only; no spec edits)

> **Status (round 8): applied.** The author approved the defaults; `SPEC.md` and `SPEC-engineering.md` now implement B1–B8/N1–N8 (T6–T9 resolved). See `HANDOFF.md` for what changed and what the next chat must notice. This file remains the evidence for the review and the known `k>2` pruning caveat.

**Verdict.** Not ready to feed as-is. The interfaces, data types, definitions, fixtures and the exact `k=2` arithmetic are implementation-ready and were verified. The core of `SPEC.md` is not: §4's bound table is defined over walks with repeated vertices, and as specified it does not finish on the spec's own Example 1B (measured: stage 4 already has a cell of 419 bounds with `|b|₁ = 1760`; stage 5 did not finish in 150 s). Two headline claims tied to that table — "≤ 2n" cells and O(1) chain-slot insertion — are false for the table actually defined. The witness/error contract (`simple_cycle`, `farkas_witnesses`, `NO_SOLUTION`) and the core↔LP wiring are also not implementable as written. A coding agent fed these two files will produce code that passes Example 1A, hangs on Example 1B, and invents the `k>2` witnesses.

**Fix list:** 8 blocking items (§2, B1–B8), 8 small ones (§3, N1–N8), and 4 decisions that are yours (§4, T6–T9). Everything else in the two files is good; §1 says what not to touch.

## Reading guide

1. §1 — what is ready; do not rewrite these parts.
2. §2 — **B1** is the real blocker; B2–B4 are contract holes; B5–B8 are cheap but must be fixed before feeding.
3. §3 — wording/consistency fixes; safe to batch.
4. §4 — algorithmic decisions only you can make (defaults stated).
5. §5 — how each claim above was checked (reproducible; nothing in `sources/` touched).

## 1. What is ready (leave alone)

- §1 definitions (Robinson matrix, threshold vector, strict uniform-embedding inequalities, β±, walks, ≼) are correct and 0-based throughout; I verified the index conventions on both fixtures.
- §2 data types are clear, except the k=2 cell/iteration mismatch in B6.
- §5b ratio logic is exactly right: `a2>0 ⇒ r > -a1/a2`, `a2<0 ⇒ r < -a1/a2`, `a2=0 ⇒ a1>0`, window seeded to (0,1), exact `Fraction`s. Verified against exhaustive simple-cycle enumeration on 600 random Robinson matrices (n = 3…7): status and ratio window always matched, and every constructed `π` satisfied every pair inequality exactly.
- §6 construction formula (including the empty-upper-set convention `ub=∞ ⇒ π[v]=lb+d[0]`) is correct; it is the formula I used in all verification runs.
- §7 dispatch on `d` / `k=2` / `k≠2` and the `EmbeddingResult` shape (engineering §2) are a good public contract.
- Fixtures A and B are correct as given: A's `π=(0,5,6.5,11.75,12.75)`, `d=(8,6)` satisfies every pair inequality; B's witness `(0,1,5,3,0)` is a simple upper-bound cycle with bound `(0,0)` and its bound really is `(0,0)` (checked step by step).
- The layer split (core / `_lp` / `_translate`), the one-public-function rule, validation-in-translation, exact-`Fraction` policy for k=2, and the packaging/naming sections are sound.

## 2. Blocking (fix before feeding)

### B1 — §4's table is unbounded on Example 1B; the 2n / O(1) claims are false

Locations: `SPEC.md` §4 lines 78–117 (especially the note at 117), cell rules 101–111; acceptance test `SPEC-engineering.md` §10 step 1 (line 146).

What the spec says: each cell holds the `≼`-minimal `(walk, bound)` pairs, walks may repeat vertices ("T stores walks … not only paths"), every ordered pair can be a concatenation target (the F1 fix), and for `k=2` each slot keeps its tightest element so insertion is O(1) and a cell has at most `2n` elements.

What happens: concatenating a walk with a closed walk adds a bound that is usually *incomparable* to the existing one (e.g. `(-1,3)` vs `(-2,7)`), so it is not pruned; bounds then leave `Z^k_n` and the table grows superlinearly. Measured on the spec's own fixtures (pure-stdlib prototype, exact arithmetic):

| run | result |
|---|---|
| A (5×5), all walks | completes, but 265 bounds total, largest cell 29, max `|b|₁ = 115` (n = 5; claim is `≤ 2n = 10`) |
| B (6×6), all walks | stage 4: 1292 bounds, largest cell 419, max `|b|₁ = 1760`; stage 5 did not finish in 150 s |
| A, k=2 chain slots without the global `minimal_in` gate | did not finish in 120 s (this is the literal "insertion is O(1)" reading; the list version finishes in < 0.01 s) |
| B, simple walks/cycles only | infeasible in < 0.01 s, witness `(0,1,5,3,0)`, max `|b|₁ = 6` |

Why: the note at line 117 ("the minimum bound value over walks equals that over paths") is about the value of the minimum at a *feasible* `d`; it does not bound the number of *minimal elements of the walk closure*, which is what the table stores. The paper's `M`/`2n` bounds count minimal elements of `Z^k_n`, i.e. bounds of simple (cycle-free) paths and cycles. The k>2 path has the same issue — the bound machinery is k-generic.

Consequence: an implementer following §4 literally ships a module that passes A and hangs on B; the roadmap's own B regression (diagonal contains `(0,0)`) never runs.

Required edit: state the walk class and a size invariant explicitly in §2/§4; delete or rewrite the line-117 note; add to §10 an acceptance test "B completes and `(0,0)` is in some `T[v][v]`". Which walk class to use is T6 below.

### B2 — `simple_cycle` is undefined for k>2 and can return a non-witness

Location: `SPEC.md` line 132; call sites 128, 141, 157–158.

Three defects: (i) "the summand with the most extreme ratio bound" has no meaning for k>2 (`ratio` is defined only in §5b); (ii) no direction (largest vs smallest ratio) is fixed; (iii) the chosen summand need not violate the inequality at the supplied `d`, so `solve_with_thresholds` (line 141) can return a "witness" with `dot(b,d) > 0`. Concrete arithmetic: a closed walk decomposing into bounds `(-4,9)` and `(3,-10)` sums to `(-1,-1)`; at `d=(1,1/2)` the dots are `+0.5` and `-2`; the rule picks the max-`|r|` summand `(-4,9)` (4/9 > 3/10), which is not violating. Also, `simple_cycle(w_lo)` at line 158 assumes the extreme ratio is preserved by the chosen decomposition — that is asserted, not proved.

Required edit: either (a) make the table hold simple cycles only, so `simple_cycle` becomes the identity and can be deleted (this is the cleanest; see T6), or (b) define `simple_cycle(w, d)` = a summand with `β·d ≤ 0` (exists because the sum is `≤ 0`), and for §5b pass a direction (`lo`/`hi`) and verify the returned cycle attains the recorded extreme.

### B3 — `farkas_witnesses` is not an algorithm; the `≤ k` claim is unjustified

Location: `SPEC.md` line 187 (the whole procedure), 183; `SPEC-engineering.md` §4 lines 89–95, and the result contract line 54.

The text promises "≤ k cycles whose inequalities are jointly contradictory (from the LP's Farkas certificate / dual multipliers; if unavailable, the cycle constraints active at the optimum)". Problems: (i) turning SciPy `linprog` marginals into a certificate is not specified; (ii) the fallback is impossible in the infeasibility branch — an infeasible LP has no optimum; (iii) `≤ k` is not derived anywhere (the dual certificate lives with the ordering constraints too, and Carathéodory gives a bound on the *constraint* count, not on cycles); (iv) it is not said whether the cycles are contradictory *alone* or together with `d ∈ D^k` (the ordering constraints are part of the system, so a certificate may need them); (v) `_lp` is said to return walks while the public contract says cycles — no one owns the simplification.

Required edit: replace with a concrete, verifiable procedure. Simplest defensible one: a deletion filter — start from the LP's cycle constraints (or the dual support, if extractable), drop one at a time while the remaining system (cycles + `d ∈ D^k`, `d₀=1`) stays infeasible (one LP call per trial); return the remaining cycles; state that the certificate is "these cycles together with `d ∈ D^k`", and drop the `≤ k` promise or prove it. Alternatively, expose the dual support and verify it before returning.

### B4 — the k>2 dispatch is placed in two contradictory layers

Locations: `SPEC.md` §5c `solve_general` lines 190–195 and §7 `uniform_embedding_core` lines 219–224 (the LP call is inside core); `SPEC-engineering.md` §1 lines 16–24 and §2 line 27 ("`__init__` … wires the layers and dispatches"; D11: "wired in the public function, not inside core"); §10 step 3 puts `solve_general` in `_lp.py`.

If the agent follows `SPEC.md` §7, `_core.py` contains `solve_general` and therefore calls the LP, contradicting "core: pure stdlib, no third-party imports" (it would import `_lp`, which is at least a reversed dependency). If it follows the engineering doc, `SPEC.md`'s `solve_general` is dead code and the LP branch must live in `__init__.py`. Also note `solve_general` is not a pure LP function: it calls `bound_generation` and `embedding_from_bounds` (core).

Required edit: pick one call graph and write it down once, in both files. The clean one: `_core.solve_general(a, search)` takes the LP search as an argument (dependency injection), `__init__.py` passes `_lp.search_thresholds`; or `__init__.py` sequences `bound_generation` → `search_thresholds` → `embedding_from_bounds` itself. Either is fine; silence is not.

### B5 — the k=2 cell path never says whether `minimal_in` still runs

Locations: `SPEC.md` lines 96–97 (L13–L14 call `minimal_in` then `insert_minimal`), 101–111 (the k=2 chain-slot rule "insertion is O(1)"), and `SPEC-engineering.md` §10 step 1.

If the k=2 branch replaces the whole L13–L14 block with `insert_minimal_k2` (the reading that makes "O(1)" true), the table explodes even on Example 1A (measured: did not finish in 120 s). If it keeps the global `minimal_in` gate, A finishes but insertion is not O(1) and the "≤2n" claim still fails once bounds leave `Z^k_n` (B1). The pseudocode has no `if k = 2` dispatch at all; the data-type table only implies one.

Required edit: add the dispatch line to `bound_generation` and state per `k` exactly which test runs. After the T6 decision, the honest statement is probably "k=2: per-slot comparison plus the simple-walk filter; O(1) because all bounds lie in `Z^2_n`" (true only with T6(a)).

### B6 — the cell type and the iteration helpers disagree

Locations: `SPEC.md` §2 line 41 (`k=2`: `dict[(t,chain),(w,b)]`), vs §4 L10–L11 lines 93–94, §5 `cycle_bounds` lines 125–129, §6 `min_bound` lines 202–203: all written as `for (w,b) in cell` / `for (_,b) in cell`.

Iterating a dict yields keys, so for a k=2 cell these loops either unpack a 2-tuple key into `(w,b)` (binding `w=t`, `b=chain`, then `dot(b,d)` fails on an int) or yield ints. `solve_with_thresholds` and `cycle_bounds` are the given-`d` and feasibility paths, so k=2 with supplied thresholds breaks too. This is exactly the kind of defect that survives a spec read and dies in the first run.

Required edit: define one cell interface (e.g. a `Cell` with `.items()`, or "a cell is a list of `(w,b)` plus an optional slot index; iteration always yields `(w,b)`") and use it in every function.

### B7 — the strict LP has no tolerance policy and the final `π` is not re-checked

Locations: `SPEC.md` §5c lines 173–187; `SPEC-engineering.md` §5 lines 98–103.

"Strict-feasible iff the optimum `t* > 0`" is a decision rule for exact arithmetic, not for HiGHS floats: a boundary case (`t* = 0` up to round-off) is exactly the infeasible case, and a tiny positive `t*` can produce a `π` whose strict inequalities round to equality. The spec also says (correctly) that the construction is in floats on this path but never says what to do if the constructed `π` fails verification.

Required edit: state the decision rule (`t* > ε·scale`, with `ε` and `scale` named), and add "verify `π` against all pair inequalities with the chosen `d`; if it fails, treat as a numerical failure (retry at higher precision / report)". Also say whether `d` and `π` are returned as the LP floats or re-solved exactly.

### B8 — the output contract misdescribes joint witnesses and leaves `NO_SOLUTION` undefined

Locations: `SPEC.md` line 226 ("each a closed walk whose bound inequality fails for every threshold vector"), line 157 (the `a2=0` branch returns the same cycle twice), lines 141/164/193 (`NO_SOLUTION` sentinel); `SPEC-engineering.md` line 54 ("1–2 (k=2) or ≤k (k>2)").

For the k=2 pair `(w_lo, w_hi)` and for the k>2 Farkas set, each cycle *individually* is satisfiable — it is the **pair/set together with `d ∈ D^k`** that is contradictory. Only the `a2=0, a1≤0` branch is individually dead. The line-226 sentence would push an implementer to look for impossible individual witnesses. Also `NO_SOLUTION` is used as a value but never defined (sentinel? exception? pair shape?).

Required edit: "infeasible is returned as `(None, witnesses)` / `EmbeddingResult(feasible=False, cycles=...)`; the witnesses are jointly contradictory together with `d ∈ D^k`; a single cycle suffices only when its bound is `(a₁,0)` with `a₁ ≤ 0`" — plus one definition of the sentinel and the exact tuple shape for each branch.

## 3. Non-blocking (batch these)

- **N1 — the diagonal seed's convention is a landmine.** §1 line 25 and §3 L4 line 55 define `β⁺(u,u) = unit(k,k−1)`, i.e. the diagonal is treated as an edge of value k, while `walk_bound(a,k,(u,u)) = 0` (empty sum). The seed is harmless *as specified* (its bound makes the trivially true `d_k > 0` and can only prune trivially true bounds), but an implementer who "fixes" it to the trivial walk's bound 0 makes every matrix infeasible (`0·d > 0`). Say why the edge convention is used and that the trivial walk is excluded, or drop the seed and exclude length-1 walks from `cycle_bounds`.
- **N2 — k=1 needs SciPy.** §5b note line 169 says k=1 "may use the general path"; the general path is the LP, so a trivially feasible input needs the optional dependency, against `SPEC-engineering.md` §9 line 140. Route k=1 to a trivial construction (`π[v] = v·d₀`-style, or `d=(1,)` + the §6 formula) or document the dependency.
- **N3 — integer thresholds lose exactness silently.** `SPEC-engineering.md` §5 line 101 says given-`d` mode "inherits the user's numeric type", but `(ub+lb)/2` on ints produces floats. Convert supplied thresholds to `Fraction` (consistent with D5) or state that ints become floats.
- **N4 — §2's `d` type conflicts with given-`d` mode.** Line 43 says "Fraction for k=2; float from LP"; given-`d` keeps whatever the user supplied. One sentence.
- **N5 — the roadmap's "disjoint write scopes" is false.** `SPEC-engineering.md` §10 steps 1, 2 and 3 all write `_core.py`, and `solve_general` is core in `SPEC.md` but `_lp` in step 3. Fix the scopes or drop the parallelism claim.
- **N6 — `REPORT.md`'s reading guide was stale** (D8 pointed to "`SPEC.md` §8", D9 to a mapping table "in `SPEC.md` §6", D13 to `docs/pseudocode.md`). Fixed in this round's `REPORT.md`; listed here so the next agent does not inherit the old pointers.
- **N7 — the `_lp` ⇄ core type boundary is undefined.** `SPEC-engineering.md` line 89 says `_lp` receives `(bound, walk)` pairs and returns `≤ k` walks; the public contract (`SPEC.md` line 226) says cycles; who simplifies and who maps to labels is unclear. One line.
- **N8 — `min_bound` on an empty cell and the `INF` comparison.** §6 line 208 guards emptiness, but `min_bound` itself is total only on nonempty cells; `min_bound(T[v][i])` at line 209 relies on a note rather than an invariant. Also `pi[v] = (ub+lb)/2 if ub < INF else ...` mixes `Fraction` and `float('inf')` — harmless in CPython but worth stating the numeric-type contract once.

## 4. To decide (defaults stated; all are your call)

- **T6 — which walks may enter the table?** This decides B1/B2/B5 and the complexity story.
  - (a) *Simple walks only* (open paths, simple cycles): matches the paper's own definition of `U_{u,v}` ("upper-bound-paths"), restores `Z^k_n`, the `≤2n`/`M` bounds, O(1) slots, makes `simple_cycle` unnecessary, and passes A, B and 600 random k=2 instances against exhaustive enumeration. **Caveat found:** pruning can drop a path bound whose representative is needed to assemble a cycle; on one random k=4 matrix a cycle bound `(-2,1,2,1)` was missed (verdict unchanged there, but completeness is not proved). Recommend (a) plus a brute-force cross-check for k=3/k≥4 before promising the k>2 witness.
  - (b) *All walks* (status quo): complete but exponential — fails the B acceptance test.
  - (c) *Walks with a length/norm cap*: needs a completeness argument; none is in the text.
  - My default: (a), with the caveat recorded as a known limitation of the k>2 path and a brute-force test that flags missed bounds.
- **T7 — k>2 witness contract.** Default: verified irreducible-infeasible-subset (deletion filter) returned as cycles + "together with `d ∈ D^k`", count not promised as `≤ k`. Alternative: keep `≤ k` only if you can prove/derive it from the dual (then say how).
- **T8 — strict-LP tolerance.** Default: `t* > 10⁻⁹ · Σ|b|₁` and an exact re-verification of the final `π`; on failure report a numerical failure rather than a false feasible.
- **T9 — k=1.** Default: trivially feasible, combinatorial, no SciPy.

## 5. Evidence (reproducible; `sources/` untouched)

Harness: a throwaway pure-stdlib prototype of `SPEC.md` in `/tmp` (not added to the repo), `fractions.Fraction` end to end, plus an exhaustive simple-cycle enumerator used as the oracle (all directed simple upper-bound cycles, exact ratio window). Nothing in `sources/` was modified; the only repo changes this round are this review and the `REPORT.md` pointer.

- Both fixtures: A feasible with a valid `π`; B infeasible with witness `(0,1,5,3,0)` / `(0,0)` — both re-derived independently, not copied from the spec.
- The paper's `i≤j` loop (state before the F1 fix) returns **feasible** for B, and for A yields a `π` that violates the lower bound at `(0,2)` (`dist = 2/3 = d₂`, the lower threshold for `t=1`; equality is not allowed) — so the F1 fix is genuinely needed, and the spec is right about that.
- `SPEC.md` §4 as written (all walks): A 265 bounds / max cell 29 / max `|b|₁ 115`; B explodes at stage 4 (1292 / 419 / 1760) and stage 5 did not finish in 150 s.
- Simple-walk table (T6(a)): A and B both correct in < 0.01 s; random k=2, n = 3…7, 600 matrices: status and exact window always equal to the oracle; every constructed `π` verified against every pair inequality.
- Given-`d` mode on the simple-walk table: A with `d=(8,6)` yields a valid `π`; B with `d=(8,6)` returns the `(0,0)` cycle.
- k=3: 120 random matrices, every oracle cycle bound is `≽` some diagonal bound (no misses). k=4: 1 miss in 120 (matrix `[[4,4,0,0,0,0],[4,4,3,1,1,0],[0,3,4,4,2,0],[0,1,4,4,4,0],[0,1,2,4,4,3],[0,0,0,0,3,4]]`, missed `(-2,1,2,1)`), verdict unchanged (that table's system is already infeasible) — the reason for the T6 caveat.

## 6. Minimal edit sequence if you accept the defaults

1. T6(a) + B1 + B5 + B6 in §2/§4 (one coherent "cell = simple walk bound + slot index" contract, explicit `k=2` dispatch, one iteration interface).
2. B2 becomes "delete `simple_cycle`; the table has only simple cycles".
3. B3 + T7 in §5c/§4 + engineering §4.
4. B4 (call graph) + B7 (tolerance) + B8 (`NO_SOLUTION`, witness wording) across both files.
5. N1–N8 (N6 was done in this round's `REPORT.md`), then a final read-through against the fixtures.
