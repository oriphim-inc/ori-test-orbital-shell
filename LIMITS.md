# LIMITS: what these results do NOT show

Written by the Sceptic after the final run (Parts A/B from the first full run, Parts C/D from
the rerun after commit 5c309cd). Everything below is either an assumption in PLAN.md or a
limit found during review.

## Assumptions made in place of missing inputs (PLAN.md)
- "Standard inclined orbit" means i = 53 deg (Starlink-like). Other inclinations change the
  nodal rate, the beta-angle history and the conjunction latitude pattern.
- SRP uses Cr = 1.5 and a cannonball area: the A/m of 0.02 m^2/kg is treated as
  attitude-averaged. A real flat plate's drag and SRP area changes with attitude by up to
  a factor of several, so da/dt and lifetime scale with whatever the true mean area is.
- Shell (Part C): every satellite circular at i = 53 deg, altitude uniform in 550 +- 5 km,
  P = round(sqrt N) evenly spaced planes, independent random phase per satellite. Avoidance
  rule (Part D): along-track burn on the cheaper satellite, miss raised to 2 x threshold,
  thresholds 0.1 to 1 km, lead times 0.5 to 2 days.

## Force model and propagation (Parts A, B)
- Static exponential atmosphere (Vallado Table 8-4). It has no solar cycle, no day/night
  bulge and no geomagnetic storms. At 550 km the real density varies by about 10x over a
  solar cycle. So da/dt (-60 m/day) and the rough lifetime (~2.9 yr down to 150 km) hold
  for this mean model only. They are not forecasts for 2026, which is near solar maximum,
  when decay would be faster.
- Gravity is J2 and J3 only. Tesseral harmonics, higher zonals, lunisolar third-body
  forces, tides and relativity are left out. Over 30 days these move the node and phase
  much less than J2 does, but they are not zero.
- Frame: precession and nutation are ignored, and the Sun comes from a low-precision
  almanac (~0.01 deg). The almanac Sun is mean-of-date while the frame is called J2000, so
  there is a ~0.36 deg (26-year precession) inconsistency. It does not matter for beta
  angles of the dawn-dusk orbit (its RAAN is set from the same Sun). The 53 deg orbit's
  RAAN0 = 0 is arbitrary anyway.
- Earth is a sphere of radius R_E for the drag altitude and for the shadow. Earth
  oblateness changes the shadow edge (eclipse entry/exit) by seconds, not minutes.
- SSO inclination: 97.594 deg is the value for a MEAN (orbit-averaged) semi-major axis of
  R+550 km. If the same 550 km is taken as an osculating radius at the node, the answer is
  97.557 deg (both are in partA_summary.json). The "J2 analytic" comparisons use
  first-order secular theory with orbit-averaged osculating elements, not Brouwer mean
  elements, so agreement at the 1e-3 level is the expected best.

## Eclipses (Part B)
- The daily eclipse values advance RAAN at the first-order J2 rate on a circular orbit of
  radius R+550 km. A 7-day full-force propagation agrees to within 0.27 min per orbit (at
  days 0 and 160 only).
- The conical model uses a constant solar radius of 0.267 deg in the check. Atmospheric
  refraction and the atmosphere's own shadow layer are ignored (in reality these add a few
  seconds of partial dimming).
- Penumbra is NOT always ~10 s per crossing. It grows as the orbit grazes the shadow
  (penumbra fraction ~0.003 at beta = 0, ~0.010 at beta = 62 deg, ~0.024 at beta = 66 deg).
  This is cone geometry, not a numerical artefact.

## Conjunctions (Part C)
- No numerical propagation. Satellites follow analytic circular mean-element orbits with
  first-order J2 secular rates only. There is no drag, no station-keeping noise, no
  eccentricity, and no short-period J2 terms (osculating radius varies by about +-10 km in
  reality, which is as large as the whole shell).
- Random independent phasing is a statistical model. Real constellations use controlled
  phasing (e.g. Walker patterns), which can remove most intra-shell close approaches, or
  can give strongly non-Poisson recurrences. The counts are an expectation for an
  uncoordinated shell. They are NOT a worst case for any real operator, and not a
  prediction for a designed constellation.
- Recurrence: about 47% of events come from pairs that meet more than once. Some
  phase-locked cross-plane pairs meet at both plane intersections every half orbit for the
  whole 30 days (up to 905 times for one N=1e4 pair). Perturbations and station-keeping
  noise would break this. So Poisson sigma = sqrt(count) understates the real scatter: the
  effective variance is 3-7x larger. The quoted Poisson errors are lower bounds.
- Same-plane pairs with |da| < 1 km drift past each other at m/s speeds and add a class of
  event the kinetic-gas formula leaves out: ~2.3% of events at N=1e3, ~0.8% at N=1e4,
  ~0.3% at N=1e5. Their number depends entirely on the assumed altitude scatter inside
  each plane. A real plane with altitude control would not produce them this way.
- The gas "same-plane" correction (1 - 1/P) overcorrects slightly. The derived factor is
  1 - pi/(2 K P), which adds +0.66% / +0.21% / +0.07% to the estimate for N = 1e3/1e4/1e5.
  With these two effects accounted for, the cross-plane counts match the gas model within
  about 1 effective sigma (checked for N = 1e3 and 1e4).
- With only 32 planes (N=1e3), no cross-plane crossing is slower than ~1.19 km/s. The
  low-speed share of encounters therefore depends on the plane count, not only on N.
- Miss distance is the distance between centres of point masses. There is no covariance
  and no probability of collision (Pc). Operators screen on Pc or ellipsoids, which give
  different maneuver counts.
- 30 days with one random seed. Nothing is shown about other seeds, other shell
  thicknesses, or multiple shells and inclinations interacting.
- partC_encounters_N100000.csv is a truncated sample (per-chunk cap; header says
  complete=False). It is not the full event list.

## Avoidance delta-v (Part D)
- Avoidance burn only. There is no return-to-slot burn (which roughly doubles the cost), no
  re-screening after a maneuver, and no maneuver-induced secondary conjunctions. Phase-locked
  repeat encounters are each charged a burn, even though in reality one burn would
  de-phase the pair.
- Clohessy-Wiltshire linear response about a circular reference with a single impulsive
  along-track burn. Differential J2 from the small da change, drag, and errors in
  predicting the conjunction (covariance growth over a 0.5-2 day lead) are ignored. Real
  operations decide on Pc with uncertain miss vectors, which typically means more and
  larger burns than this deterministic rule.
- Annual per-satellite delta-v = (30-day sum / N) x 365.25/30. This assumes the 30-day
  window is representative and that the burden is shared evenly. The per-satellite
  distribution is very uneven (N=1e4: mean 80 vs max 986 events per satellite in 30 days).
- At N = 1e5 and a 1 km threshold the rule implies ~4800 maneuvers per satellite per
  year (~13 per day). That is operationally impossible, so the delta-v there only measures
  how strongly the cost scales. It is not a feasible operations budget.
