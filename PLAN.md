# PLAN — Cartesian orbit propagator, eclipses, shell conjunctions, avoidance Δv

Requested by: Oriphim. Planned by Ori (no attachments; all inputs are standard values
cited below or labelled ASSUMPTION).

Units inside the code: km, s, kg. Frame: Earth-centred inertial (mean equator/equinox of
J2000, precession/nutation ignored). Epoch for every part: 2026-01-01 00:00 UTC.

## Physics (the propagator)

State `[x y z vx vy vz]` integrated with **DOP853** (Dormand–Prince 8(5,3), scipy
`solve_ivp`), rtol 1e-12, atol 1e-9 km (tolerances are numerical settings, tunable).
Acceleration = sum of:
1. Two-body, μ = 398600.4418 km³/s² (EGM96 / WGS-84).
2. Zonal J2 and J3, R⊕ = 6378.137 km, J2 = −√5·C̄20, J3 = −√7·C̄30 with EGM96
   C̄20 = −4.841653717e-4, C̄30 = 9.572541738e-7 (→ J2 = 1.0826267e-3, J3 = −2.53266e-6).
   Formulas: Vallado (2013) *Fundamentals of Astrodynamics and Applications* 4th ed. §8.7.
3. Drag, a = −½ ρ (Cd A/m) |v_rel| v_rel, v_rel = v − ω⊕×r, ω⊕ = 7.292115e-5 rad/s.
   Density: Vallado (2013) Table 8-4 piecewise exponential model (static, no solar
   cycle) — the "exponential" option the request allows.
4. Solar radiation pressure, cannonball: a = −P☉ Cr (A/m) (AU/d)² ν ŝ, P☉ = 4.56e-6 N/m²
   at 1 AU, ν = conical umbra/penumbra shadow function (Montenbruck & Gill 2000, §3.4.2).
   Sun position: low-precision almanac formula (Vallado 2013 Alg. 29, ~0.01°).

Satellite (Part A): A/m = 0.02 m²/kg, Cd = 2.2 (given). ASSUMPTIONS: Cr = 1.5 for SRP;
cannonball (attitude-averaged) area — a "flat" plate's real A/m depends on attitude.

## Part A — 550 km circular orbit
- A1 Kepler check: two-body only, i = 53°, 30 days. Compare to the analytic Kepler
  solution: period 2π√(a³/μ), position error after 30 days, energy and |h| drift.
- A2 J2-only: measured nodal regression rate and nodal period vs first-order analytic
  (Ω̇ = −1.5 n J2 (R/a)² cos i).
- A3 Full model (J2+J3+drag+SRP), i = 53° (ASSUMPTION: "standard inclined" = Starlink-like
  53°), 30 days: osculating elements vs time (CSV + figure), measured nodal precession
  (linear fit of Ω at ascending nodes), nodal period (node-to-node), semi-major-axis decay
  rate (fit of a at ascending nodes) vs analytic −√(μa) ρ Cd A/m; rough lifetime from
  orbit-averaged integration of the same density model (labelled rough).
- A4 Sun-synchronous inclination at 550 km: first-order analytic from
  Ω̇ = 360°/365.2422 d (tropical year), then refined by a secant root-find on the
  numerically measured Ω̇ of the full-force propagator (5-day arcs).

## Part B — eclipse vs beta angle over 2026
Two 550 km circular orbits: inclined i = 53° (RAAN₀ = 0), and dawn-dusk SSO
(i from A4, LTAN 18:00, i.e. RAAN₀ = α☉ + 90°). Daily over 365 days: RAAN advanced at
the J2 secular rate (validated in A2/A3), β = asin(ŝ·ĥ), eclipse duration from one
sampled orbit (1 s) with the conical shadow model (umbra and penumbra) and compared to
the cylindrical closed form f = (1/π) acos(√(h²+2Rh)/((R+h)cos β)). Cross-check: 7-day
full-force propagation of each orbit, eclipse durations per orbit vs the daily estimate.
Deliver eclipse fraction vs β and vs day for both orbits.

## Part C — conjunctions in a shell of N = 1,000 / 10,000 / 100,000
Shell (ASSUMPTIONS): circular, i = 53°, altitude uniform in 550 ± 5 km (10 km thick),
P = round(√N) planes evenly spaced in RAAN, satellites assigned round-robin to planes,
**random mean argument of latitude** per satellite (randomized phasing), fixed RNG seed.
Propagation: analytic mean-element J2 secular model (Ω̇, u̇ per satellite) — the full
numerical propagator is far too slow for 10⁵ objects × 30 days; this scale-down is part
of the plan and the secular rates are checked against A2/A3. No drag in Part C
(satellites assumed station-kept).
Screening: time step dt (numerical, default 10 s), spatial hash grid with cell ≥
R_s = d + v_rel,max dt/2 (+margin), linear relative-motion closest approach, each
encounter counted once (in the step whose ±dt/2 window holds the TCA). Count events with
miss distance < d = 1 km over 30 days. Validated against brute-force O(N²) on small cases.
Kinetic-gas comparison: (i) naive uniform-density shell, rate = ½ N² σ ⟨v_rel⟩ / V with
σ = πd², V = 4πr²Δ, ⟨v_rel⟩ = 4v/π (isotropic horizontal); (ii) latitude-resolved gas for
a single-inclination shell, Γ = σ v N² K(sin i) / (2π³ r² Δ) (K = complete elliptic
integral), with derived (not fitted) corrections for finite shell thickness and for
same-plane pairs that cannot cross. Deliver counts vs N with Poisson errors and the ratio
to both gas estimates.

## Part D — avoidance Δv
Rule: if predicted miss < threshold d_th, one of the two satellites (the cheaper)
performs a single along-track burn at lead time T before TCA, sized with the
Clohessy–Wiltshire solution so the miss (projected on the plane ⟂ v_rel) becomes
2·d_th. Thresholds d_th ∈ {0.1, 0.2, 0.5, 1.0} km, lead times T ∈ {0.5, 1, 2} days
(ASSUMPTIONS, sweep). Annual Δv per satellite = (Σ Δv over 30-day encounters / N) ×
365.25/30. Avoidance burn only (a return-to-slot burn would ~double it); no
maneuver-induced secondary conjunctions.

## What counts as the answer
- A: Kepler check errors at round-off/tolerance level; numeric Ω̇ and nodal period, da/dt
  for the given satellite; SSO inclination (analytic and numeric), expected ≈ 97.6°.
- B: eclipse fraction vs β curves for both orbits; dawn-dusk eclipse-free most of year.
- C: table of conjunction counts (<1 km, 30 d) vs N and ratio to kinetic-gas estimates;
  counts should scale ~N².
- D: table of annual Δv per satellite vs N, threshold and lead time.

## Run
`python run.py` (all parts, full size); `python run.py --quick` (pilot: shortened arcs,
1 day for Part C). Outputs in results/.
