# Acceptance cases (recorded round 12)

Runnable suite: [`tests/acceptance.py`](../tests/acceptance.py). It is self-contained: it re-implements the edge-bound algebra, the embedding verifier and the k=2 witness check independently of the package, so a broken core cannot mask itself. It uses the package's private `_core` only for the two internal target checks G3/G4.

Run from the repository root:

```bash
python tests/acceptance.py           # needs Python >= 3.10; numpy optional
# with the optional SciPy extra (the k>2 automatic cases run only then):
uv venv --python 3.13 && source .venv/bin/activate
uv pip install -e ".[lp]"
python tests/acceptance.py
```

Without SciPy the four `k>2` automatic/joint cases and G4 are reported `SKIP`; without numpy the numpy-input case is skipped too. The suite exits non-zero on any failure.

## Current status (round 13, all groups pass)

With SciPy: **20 passed, 0 failed, 0 skipped**. Without SciPy (with numpy): **16 passed, 0 failed, 4 skipped**; without numpy one more skip (15/0/5). The round-13 hardening targets G1–G4 are fixed (see `REPORT.md`, T12–T15 resolved); the recorded expectations are unchanged. The suite exits non-zero on any failure.

## What is covered

### Fixtures (recorded from the paper/review)

| Case | Input | Expected |
|---|---|---|
| A | 5×5, k=2 | feasible, `d=(1, 2/3)`, exact `π` verified |
| B | 6×6, k=2 | infeasible, witness `(0,1,5,3,0)` or an equivalent contradiction on `r=d₂/d₁` |
| K1 | 4×4, k=1 | feasible, `d=(1,)`, exact `π` verified |
| A₃ | A with diagonal 3, k=3 | feasible via the LP, `d=(1, 0.5, 0.25)`; given `d=(8,6,1/2)` also feasible |
| B₃ | B with diagonal 3, k=3 | infeasible via the LP (the `(0,0,0)` cycle); given `d=(8,6,1/2)` also infeasible |
| JOINT | 5×5, k=3, recorded matrix | infeasible **only jointly**: `(0,1,2,0)` and `(2,4,3,2)` are each satisfiable, the pair is not; the deletion filter must return ≥2 cycles, the set re-solves infeasible, and every single cycle alone re-solves feasible |

The JOINT matrix (0-based):

```
(3,3,0,0,0)
(3,3,2,0,0)
(0,2,3,1,1)
(0,0,1,3,2)
(0,0,1,2,3)
```

### Recorded k=2 cases (exact windows from exhaustive simple-cycle enumeration)

The five feasible matrices and their exact windows `(lo, hi)` for `r=d₂/d₁`; the library must return `d=(1, float((lo+hi)/2))` and a valid embedding:

| n | window |
|---|---|
| 3 | (0, 1) |
| 5 | (0, 1/2) |
| 5 | (0, 1/3) |
| 3 | (1/2, 1) |
| 6 | (1/3, 1) |

The three recorded infeasible k=2 matrices must return 1–2 simple cycles whose bounds contradict `r ∈ (0,1)` (one cycle individually dead, or disjoint ratio intervals). The matrices are in `K2_FEASIBLE` / `K2_INFEASIBLE` in the suite.

### Recorded k=1 cases

Three binary Robinson matrices (the all-ones 3×3, the 2×2 identity, and a 3-node path graph) must be feasible with `d=(1,)` and a verified embedding.

### Given-`d` k=3

`A₃` with `d=(8,6,1/2)` (feasible), `A₃` with `d=(8,6,5)` (infeasible; the witness's bound must have negative dot with `d`), `B₃` with `d=(8,6,1/2)` (infeasible; the `(0,0,0)` cycle).

### Edge cases

`n=1` at k=1 and k=2; `n=2` with a null edge at k=2; labels mapping (B with labels `a..f`, the witness must be a tuple of labels with first == last); numpy integer input accepted, numpy float input rejected; nine invalid inputs must raise `ValueError` (non-square, ragged, float entries, asymmetric, diagonal not max, not Robinson, wrong/duplicate labels, wrong/unordered thresholds).

### Hardening targets G1–G4 (round 13, all fixed)

- **G1** `((0,0),(0,0))` → `ValueError` (the all-zero matrix is outside the `k≥1` domain). Was `TypeError`; fixed in round 13.
- **G2** generator `labels`, generator `thresholds`, and `thresholds=(inf, 1)` → `ValueError`. Was `TypeError` / `OverflowError`; fixed in round 13.
- **G3** in `k≠2` cells, the distinct-bound count equals the pair count (no duplicate bounds per cell). Was 8 pairs / 4 distinct in A₃'s `T[0][0]`; fixed in round 13.
- **G4** the `k>2` verification runs on `Fraction` `d`/`π` (the wrapper converts before calling `verify_embedding`). Ran on floats; fixed in round 13.

## How the recorded values were produced

- k=2 windows: exhaustive enumeration of all directed simple upper-bound cycles (test-only oracle, not committed) with exact `Fraction` ratio intervals. The 200-random-matrix run in the round-11 review used the same oracle; the five windows above are a subset chosen to cover the window shapes.
- k=3 witnesses: an exact rational max-slack LP (vertex enumeration; test-only) established that B₃'s witness is individually infeasible and that JOINT's two cycles are each feasible while the pair is not; the suite re-checks the structure and, when SciPy is present, re-solves the LP on the returned witness set.
- The A/B/A₃/B₃/K1 expectations are the published/verified fixtures; the exact `π` values are not pinned except where the suite re-verifies inequalities exactly.

## What is not covered (v1)

Large-`n` performance, adversarial float boundaries (the `k>2` LP path is verified exactly), and the completeness of the pruned `k>2` table (a documented limitation; the verification step turns a false feasible into `RuntimeError`).
