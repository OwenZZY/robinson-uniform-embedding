# HANDOFF — what to hand to the implementing chat

This file is for the **next chat** (the coding agent). Read [`SPEC.md`](SPEC.md) first as the contract, then [`SPEC-engineering.md`](SPEC-engineering.md). This handoff explains what was wrong with the previous draft, what changed, and what not to "fix". The detailed evidence for the changes is in [`REVIEW-spec.md`](REVIEW-spec.md); that file is history, **the two spec files are the contract**.

Status: round 8, applied to `SPEC.md` and `SPEC-engineering.md` after the round-7 review. The fixtures below were re-verified against the patched pseudocode before it was committed.

## 1. What was wrong with the previous draft

1. **The bound table was defined over walks with repeated vertices.** The previous draft said that was harmless (true only for the value of the minimum at a feasible `d`), but it makes the table grow without bound: on the spec's own Example 1B the cell at stage 4 had 419 bounds with `|b|₁ = 1760`, and the algorithm did not finish. The paper's `≤ 2n` / `M` size bounds are for bounds of *simple* paths and cycles (`Z^k_n`).
2. **The `k=2` chain-slot rule conflicted with the "minimality" test.** The literal "insertion is O(1)" reading (compare only within the slot) already hung on Example 1A; keeping the global test made the O(1)/`2n` claims false.
3. **`simple_cycle` was not implementable.** It was defined via "the most extreme ratio", which means nothing for `k>2`, has no direction, and could return a cycle that does not violate the inequality — an invalid witness for the given-`d` mode.
4. **`farkas_witnesses` was a promise, not an algorithm.** No procedure to extract or verify the `k>2` witnesses, an impossible fallback ("active constraints at the optimum" when the LP is infeasible), and an unproved `≤ k` bound.
5. **The `k>2` dispatch lived in three places at once** (core in `SPEC.md`, public function in the engineering doc, `_lp` in the roadmap), contradicting "core: no third-party imports".
6. **The `k=2` cell type and the loops disagreed** (dict cells iterated as if they were lists), so `cycle_bounds`, `min_bound` and the given-`d` mode broke for `k=2`.
7. **The strict LP had no tolerance policy or verification**, and the output contract misdescribed joint witnesses ("each fails for every threshold vector" is false for a `k=2` pair and for a `k>2` Farkas set); `NO_SOLUTION` was undefined.

## 2. What changed (the patched contract)

**`SPEC.md`**

- §1: `beta_plus`/`beta_minus` are defined only for `u≠v` (no diagonal edge convention); added the **simple walk** definition (open: all distinct; closed: `len(w) ≥ 3`, only `w[0]=w[-1]` repeats; a simple cycle has ≥ 2 distinct vertices); every stored bound has `Σ|b[j]| ≤ n`.
- §2: cell interface fixed — `cell_items(cell)` is the only way to iterate a cell (`k=2`: slot dict; otherwise: list); `d`/`pi` types stated per path.
- §3: primitives; `simple(w)`; `cell_items`; no `u=v` cases.
- §4: **init only for `i<j`** (no diagonal seed); `insert` filters non-simple walks and dispatches on `k` (`insert_minimal_k2` for `k=2`, `minimal_in` + `insert_minimal` otherwise); all ordered pairs are swept; empty diagonal before the first cycle is expected; size invariant (`|b|₁ ≤ n`, `k=2` slots `t ≤ n`, at most `2n+2`); a **known limitation note** for `k>2` (pruning can in principle miss a cycle bound; the public layer verifies `pi`).
- §5: `cycle_bounds` no longer simplifies walks (they are simple cycles by §4); §5a given-`d`; §5b `k=2` with two new single-witness branches (`a2>0` and `r ≥ 1`; `a2<0` and `r ≤ 0`) so the ratio window never returns `None` as a witness; §5c LP + `EPS = 1e-9·(1+max Σ|b|)` + the **deletion-filter** `minimal_witnesses` (one LP per trial, witnesses jointly contradictory with `d ∈ D^k`, no `≤ k` claim); §5d `k=1` (always feasible, `d=(1,)`, no SciPy).
- §6: construction + `verify_embedding` (exact boolean check, used by the wrapper on every successful return).
- §7: dispatch `thresholds` → `k=1` → `k=2` → general; output is `(d, pi)` or `(None, witnesses)`; witness semantics: single cycle only for individually-unsatisfiable bounds (`k=2`), otherwise a jointly contradictory set.

**`SPEC-engineering.md`**

- §1: `_core.solve_general(a, search)` takes the search callable; the core never imports `_lp`; `__init__.py` supplies `_lp.search_thresholds` only when the `k>2` search is needed.
- §2–§3: thresholds are converted to `Fraction`; **every** successful result is verified (`core.verify_embedding`) and a failure raises `RuntimeError` (not a false feasible); the `_run` sketch shows the wiring.
- §4: `_lp` signature, `linprog` formulation with `d[0]=1` substituted, status handling (`0`/`2`/other), `EPS`, deletion filter, lazy import.
- §5: numerics (exact for `k=1`/`k=2`/given-`d`; floats + verification for the LP).
- §8: error model (`ValueError`, `ImportError`, `RuntimeError`; infeasibility is a result; core never raises).
- §10: accurate roadmap (steps 1–2 sequential in `_core.py`; 3–4 parallel after the core API is fixed) and fixtures **A, B, A₃, B₃, K1** with acceptance checks.
- §11: out of scope, including a completeness proof for the pruned `k>2` table.

## 3. What to notice while implementing

- **0-based everywhere**; pseudocode lines `L#` are the contract. Do not transpose indices.
- **The table stores simple walks only.** Do not add a diagonal seed, do not store repeated-vertex walks, do not reintroduce `simple_cycle`. The `k=2`/`k=1` completeness and the `2n` bound depend on this.
- **Cells:** always `cell_items`; `k=2` cells are dicts (`slot_key → (w,b)`), everything else is a list. `bound_generation`'s `insert` dispatches on `k`; `k=2` uses the slot rule only (no global `minimal_in`).
- **`k=2`:** the ratio window is seeded to `(0,1)`; a bound outside the window is an immediate single-cycle witness. `d = (Fraction(1), (lo+hi)/2)`. Exact `Fraction`s.
- **`k>2`:** the LP is separate and optional; witnesses are a deletion-filter set, **not** promised `≤ k`, and are contradictory *together with* `d ∈ D^k`. The LP path is floats; the wrapper verifies `pi`.
- **`k=1`:** no search, no SciPy; `d=(1,)` and the §6 construction.
- **Sentinel:** failure is `(None, witnesses)`; do not invent an exception for infeasibility.
- **The wrapper verifies every result** and raises `RuntimeError` on failure; that is the v1 safeguard for the pruned `k>2` table. It must not silently return `feasible=True` after a failed check.
- **Where the spec overrides the paper:** the paper's `UBW(i,j)` are upper-bound *paths* — the spec makes that explicit; the `k=2` individually-unsatisfiable branches, the `k=1` path, the deletion filter and the verification are spec-level additions. Follow the spec, not the paper, when they differ.

## 4. Acceptance tests (run these before reporting done)

Use the fixtures in `SPEC-engineering.md` §10 (0-based):

- **A** (5×5, `k=2`): `feasible=True`; `verify_embedding` true; `pi` strictly increasing.
- **B** (6×6, `k=2`): `feasible=False`; 1–2 simple cycles; `(0,1,5,3,0)` with bound `(0,0)` or an equivalent pair. **B must complete quickly** — if it does not, the table is storing non-simple walks.
- **A₃** (A with diagonal 3, `k=3`): feasible via the LP; verified.
- **B₃** (B with diagonal 3, `k=3`): infeasible via the LP; the `(0,0,0)` cycle or an equivalent jointly contradictory set.
- **K1** (4×4, `k=1`): feasible; verified; **no SciPy installed**.
- Optional but recommended (not a deliverable): compare the `k=2` ratio window against exhaustive simple-cycle enumeration on small random Robinson matrices (this is how the patches were checked; 600 random matrices, exact match).

If you have SciPy, also confirm that `k=1`, `k=2` and given-`thresholds` modes run with SciPy absent (simulate with an import guard or a fresh environment).

## 5. Known limitation (record it, do not try to fix it in v1)

The `k>2` bound table is pruned by bounds, not by the walks' vertex sets; a pruned simple path can in principle be needed to assemble a simple cycle. A concrete 6×6, `k=4` matrix where a cycle bound was missed exists (in `REVIEW-spec.md` §5). The `k=2` pipeline was checked exhaustively against the oracle; the `k>2` failure mode is a **false feasible**, which the `verify_embedding` check turns into a `RuntimeError`. A completeness proof or a different enumeration is out of scope for v1; do not silently weaken the check or claim completeness.

## 6. What not to do

- Do not modify anything in `sources/`.
- Do not re-add the diagonal seed, the walk table, `simple_cycle`, or the global `minimal_in` for `k=2`.
- Do not "simplify" the `k=2` ratio branches back to the old two-branch version (the single-witness branches are needed for bounds outside `(0,1)`).
- Do not turn infeasibility into an exception, and do not let a failed verification pass as `feasible`.
- Do not promise `≤ k` witnesses for `k>2` or "each witness fails for every threshold vector".
- Do not add tests/README/examples beyond what the owner asked for; offer them.

## 7. Suggested order

1. `_core.py`: primitives, `simple`, `cell_items`, `bound_generation`, `cycle_bounds` (`SPEC.md` §3–§5). Gate: A and B complete; B's `(0,0)` cycle present; cells small.
2. `_core.py`: `embedding_from_bounds`, `verify_embedding`, `ratio_bounds`, `solve_k2`, `solve_k1`, `solve_with_thresholds`, `solve_general` (`SPEC.md` §5a–§5d, §6). Gate: A, B, A₃, B₃, K1.
3. `_lp.py`: `search_thresholds`, `minimal_witnesses` (`SPEC.md` §5c, engineering §4). Gate: A₃/B₃ through the LP; no SciPy needed elsewhere.
4. `_translate.py` + `__init__.py`: predicates, bijection, conversions, verification, result. Gate: the public examples; error cases.

If something in the spec is genuinely ambiguous while implementing, ask the owner rather than guessing, and record the resolution in `REPORT.md`.
