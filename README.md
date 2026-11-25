> **Collision risk in a satellite shell.** A project Ori built and ran on a problem from Oriphim's held-out test set, which Ori was never tuned on. Short paper: [paper/paper.pdf](paper/paper.pdf). The Sceptic's review is in [notes/review.md](notes/review.md); limits are in [LIMITS.md](LIMITS.md). [oriphim.com](https://oriphim.com)
>
> How this repo was made: Ori built and ran the project and wrote the paper. Afterwards, Oriphim added the licence, edited the author line and wrote this note, with help from Claude Code.

# Cartesian orbit propagator: 550 km orbits, eclipses, shell conjunctions, avoidance Δv

By Oriphim, with Ori.

**Asked:** build a Cartesian propagator with a high-order integrator (DOP853) and J2, J3,
drag and SRP. Use it at 550 km for: (A) Kepler checks, nodal precession, period, decay and
the sun-synchronous inclination; (B) eclipse duration vs β angle over a year; (C) close
approaches below 1 km in 30 days for shells of 10³–10⁵ satellites, compared with kinetic
gas; (D) annual avoidance Δv. The agreed experiment, with every assumption, is in
[PLAN.md](PLAN.md).

## Answer

**A — 550 km orbit, A/m 0.02 m²/kg, Cd 2.2, i = 53°, 30 days, full force model**
- **Kepler check:** after 30 days the two-body run is 0.12 m from the analytic Kepler
  solution. Period error is 6e-12 relative; energy and angular momentum drift ~1e-11.
- **Nodal precession:** −4.496 °/day. First-order J2 theory gives −4.491 °/day. A J2-only
  run gives −4.4941 °/day, 0.07 % from the formula; the leftover is J2²-order.
- **Nodal (node-to-node) period:** 5734.3 s (95.57 min). The Keplerian period is
  5739.0 s and the J2 analytic value 5735.4 s.
- **Decay:** the mean semi-major axis falls by 59.9 m/day (−1.79 km over 30 days). The
  orbit-averaged analytic estimate is 58.6 m/day. With the same static atmosphere the
  rough lifetime to 150 km is ~2.9 years (1051 days).
- **Sun-synchronous inclination:** **97.594°** when "550 km" means a mean semi-major axis
  of R⊕+550 km. That is the numerical root-find on the full-force propagator. The
  first-order formula gives 97.593°. If 550 km is the starting osculating radius instead,
  the answer is 97.557°.

**B — eclipses over 2026 at 550 km**
- **Eclipse fraction vs β:** it follows the cylindrical closed form. The fraction is 0.37
  (35.4 min) at β ≈ 0 and falls to zero at β* = 67.0°.
- **53° orbit:** β swings between −76° and +75°. Eclipses on 340 of 365 days, 30.6 % of
  each orbit on average.
- **Dawn-dusk SSO (i = 97.594°, LTAN 18:00):** β stays between 59° and 88°. Eclipses occur
  only from 2 May to 11 August 2026 (102 days), peaking at 21.3 min umbra on 21 June.
  Eclipse-free the rest of the year.

**C — close approaches below 1 km in 30 days (550 ± 5 km shell, i = 53°, random phasing)**

| N | planes | encounters | latitude-resolved gas (corrected) | sim / gas | sim / naive uniform gas |
|---|---|---|---|---|---|
| 1,000 | 32 | 3,950 | 3,836 | 1.030 | 0.61 |
| 10,000 | 100 | 397,772 | 391,985 | 1.015 | 0.61 |
| 100,000 | 316 | 39,668,099 | 39,469,096 | 1.005 | 0.61 |

- **Scaling:** counts grow as N^2.001.
- **Gas model:** the kinetic gas works, but only after accounting for how a single
  inclination concentrates satellites and slows relative speeds near ±53° latitude:
  Γ = σ v N² K(sin i)/(2π³ r² Δh).
- **Naive uniform shell:** the gas model with uniform density and ⟨v_rel⟩ = 4v/π
  overestimates by ~1.65×.

**D — annual avoidance Δv per satellite**

The rule is one along-track burn 1 day before closest approach that doubles the miss
distance to 2× the threshold.

| N | threshold 100 m | 200 m | 500 m | 1 km |
|---|---|---|---|---|
| 1,000 | 0.0003 m/s | 0.003 | 0.05 | 0.39 |
| 10,000 | 0.004 | 0.03 | 0.49 | 3.8 |
| 100,000 | 0.039 | 0.31 | 4.8 | 37 |

- **Scaling:** Δv grows ∝ N and ∝ threshold³: rate ∝ threshold² times Δv per burn ∝
  threshold. It halves when the lead time doubles; the 0.5- and 2-day leads are in
  `results/partD_delta_v.csv`.
- **Per burn:** the average is ~0.8 mm/s at 100 m and ~7.7 mm/s at 1 km.
- **What drives the total:** the number of manoeuvres, not the size of each burn. A
  1 km threshold at N = 100,000 would mean ~4,800 burns per satellite per year.

## Evidence
- **Tests:** `pytest tests` gives 43 passing. The unit tests are mine. The pass/fail
  checks in `tests/checks/` were written by the independent reviewer before the main run.
  - J2 and J3 accelerations equal the numerical gradient of the zonal potential.
  - The drag vector, shadow function and element round-trip are checked.
  - The hash-grid screener matches an O(N²) brute force pair for pair, and matches an
    independent 0.25 s distance-minimum search.
  - A regression test makes sure a slow co-moving pair is counted once.
- **Part B cross-check:** 420 orbits of full-force propagation (7 days in January, 7 days
  near day 160) agree with the daily geometry to 0.08 min on average and 0.27 min at
  worst.
- **First run failed review and was fixed.** The first full Part C run counted
  near-co-moving same-plane pairs (Δa < 1 m, drifting ~1 mm/s) at almost every step. One
  satellite at N = 10⁵ logged 148,933 "encounters". An encounter is now a true range
  minimum: the exact range-rate changes sign inside the step, and the time and miss
  distance are refined by bisection. Parts C/D were rerun; the maximum per satellite is
  now 1,926.
- **Why counts sit 0.5–3 % above the gas estimate.** The reviewer's
  `test_C_excess_decomposition` splits the excess:
  - Same-plane drift passes, which the gas model leaves out: 2.3 % at N = 1,000 and
    0.8 % at N = 10,000.
  - The (1−1/P) same-plane correction is slightly wrong; the exact factor is
    1 − π/(2KP).

  Allowing for both, cross-plane counts match the gas estimate within 0.5 %.
- **Two reviewer checks were tightened to correct physics:**
  - The penumbra fraction grows near grazing (0.0103 at β = 61.9°, confirmed by hand).
  - With only 32 planes, no cross-plane crossing is slower than 1.19 km/s.

  Reviews are in `notes/` and the ori review log.
- **Outputs (all in `results/`):**
  - elements vs time: `partA_osculating_elements.csv`, `partA_mean_elements.csv`,
    `partA_elements.png`
  - eclipse fraction vs β: `partB_*.csv`, `partB_eclipse_fraction_vs_beta.png`,
    `partB_beta_and_eclipse_vs_day.png`
  - conjunction counts vs N: `partC_conjunctions.csv`, `partC_conjunctions_vs_N.png`,
    `partC_encounter_distributions.png`
  - annual Δv: `partD_delta_v.csv`, `partD_annual_dv_vs_N.png`
  - summaries: `part*_summary.json`

## What this does not show
See [LIMITS.md](LIMITS.md). The main points:
- **Error bars:** 47–48 % of encounters come from phase-locked pairs that meet
  repeatedly (up to 905 times). Poisson √N errors are therefore lower bounds.
- **Part C propagation:** analytic circular J2-secular orbits, not the full propagator,
  and no drag.
- **Atmosphere:** static, with no solar cycle, so decay and lifetime are for mean
  conditions only.
- **Satellite model:** Cr = 1.5 and a cannonball area are assumptions.
- **Part D:** the Δv is the avoidance burn alone. There is no return burn, no
  collision-probability screening, no re-screening, and every repeat encounter is charged.
- **Random phasing:** this is an uncoordinated shell. Real constellations phase their
  planes to avoid most of these encounters.

## Re-run
```
uv venv .venv && uv pip install --python .venv/bin/python -e '.[test]'
.venv/bin/python run.py            # everything, ~1-2 h (Part C at N=1e5 dominates; 15 cores)
.venv/bin/python run.py --quick    # pilot: A 3 d, B 30 d, C/D 1 d  -> results_quick/
.venv/bin/python run.py --parts AB # Parts A and B only (~10 s)
.venv/bin/python -m pytest -q tests
```
The random seed is fixed (`seed=1`). `--dt` sets the Part C screening step (default 10 s).
