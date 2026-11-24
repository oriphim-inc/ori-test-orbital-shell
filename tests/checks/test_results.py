"""Sceptic's pass/fail checks on results/ (written BEFORE the main run).

Reference values are recomputed here from textbook formulas and standard constants,
independently of src/, so a result cannot pass by agreeing with its own code.
Tolerances come from physics (size of neglected terms), not from pilot numbers.

Run:  .venv/bin/pytest tests/checks            (reads results/)
      ORI_RESULTS=results_quick .venv/bin/pytest tests/checks   (pilot; year-long checks skip)
Every test skips if its result file is absent.
"""
import csv
import json
import math
import os

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ellipk

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(ROOT, os.environ.get("ORI_RESULTS", "results"))

# Independent standard constants (EGM96 / WGS-84; Vallado 2013)
MU = 398600.4418          # km^3/s^2
RE = 6378.137             # km
J2 = 1.0826267e-3
WE = 7.292115e-5          # rad/s
H = 550.0
A0 = RE + H
DAY = 86400.0
SUN_ANG_RAD_DEG = 0.267   # apparent solar radius at 1 AU (deg)
BETA_STAR = math.degrees(math.asin(RE / A0))   # 67.04 deg


def need(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        pytest.skip(f"{p} not produced yet")
    return p


def jload(name):
    with open(need(name)) as f:
        return json.load(f)


def csv_rows(name):
    with open(need(name)) as f:
        return list(csv.DictReader(f))


def csv_array(name):
    return np.genfromtxt(need(name), delimiter=",", names=True)


def j2_raan_rate(a, i_rad):
    n = math.sqrt(MU / a ** 3)
    return -1.5 * n * J2 * (RE / a) ** 2 * math.cos(i_rad)


def vallado_rho(h):
    """Vallado (2013) Table 8-4, only the bands needed between 150 and 600 km (+ edges)."""
    tab = [(150, 2.070e-9, 22.523), (180, 5.464e-10, 29.740), (200, 2.789e-10, 37.105),
           (250, 7.248e-11, 45.546), (300, 2.418e-11, 53.628), (350, 9.518e-12, 53.298),
           (400, 3.725e-12, 58.515), (450, 1.585e-12, 60.828), (500, 6.967e-13, 63.822),
           (600, 1.454e-13, 71.835)]
    h0, r0, sh = [t for t in tab if t[0] <= h][-1]
    return r0 * math.exp(-(h - h0) / sh)


# ------------------------------------------------------------------ hygiene
ALL_FILES = ["partA_summary.json", "partA_osculating_elements.csv", "partA_mean_elements.csv",
             "partB_summary.json", "partB_inclined53.csv", "partB_dawndusk_sso.csv",
             "partB_propagated_check.csv", "partC_conjunctions.csv", "partD_delta_v.csv",
             "partCD_summary.json"]


@pytest.mark.parametrize("name", ALL_FILES)
def test_no_nan_or_inf(name):
    p = need(name)
    if name.endswith(".json"):
        txt = open(p).read()
        assert "NaN" not in txt and "Infinity" not in txt, f"{name} contains NaN/Infinity"
        return
    with open(p) as f:
        rows = list(csv.reader(f))[1:]
    assert rows, f"{name} is empty"
    for row in rows:
        for x in row:
            try:
                v = float(x)
            except ValueError:
                continue  # text field (e.g. orbit name)
            assert math.isfinite(v), f"{name} has non-finite value {x!r}"


# ------------------------------------------------------------------ Part A
def test_A1_kepler_at_integrator_tolerance():
    k = jload("partA_summary.json")["kepler_check"]
    T_ref = 2 * math.pi * math.sqrt(A0 ** 3 / MU)
    assert k["period_kepler_s"] == pytest.approx(T_ref, rel=1e-12)
    assert abs(T_ref - 5738.99) < 0.05          # 550 km two-body period
    # rtol 1e-12 over <= 30 d (~450 revs): global error must stay well below 1 m
    assert k["pos_err_km"] < 1e-2 * max(1.0, k["days"] / 30)   # < 10 m (generous)
    assert k["vel_err_kms"] < 1e-5
    assert k["energy_rel_drift"] < 1e-9
    assert k["angmom_rel_drift"] < 1e-9
    assert k["period_rel_err"] < 1e-8


def test_A2_j2_nodal_rate_vs_first_order():
    j = jload("partA_summary.json")["j2_check"]
    # numeric vs first-order secular theory: neglected terms are O(J2) relative ~1e-3
    assert j["raan_rate_rel_diff"] < 5e-3
    # independent reference at a = R+550 km, i = 53 deg (mean a differs by a few km)
    ref = math.degrees(j2_raan_rate(A0, math.radians(53))) * DAY
    assert ref == pytest.approx(-4.48, abs=0.05)
    assert j["raan_rate_numeric_deg_day"] == pytest.approx(ref, rel=1e-2)
    # nodal period: first-order theory vs numeric, O(J2) -> 0.1 % (~6 s)
    assert j["nodal_period_numeric_s"] == pytest.approx(j["nodal_period_analytic_s"], rel=1e-3)
    # at 53 deg the J2 nodal period is SHORTER than the Kepler period of the same a
    assert j["nodal_period_numeric_s"] < 2 * math.pi * math.sqrt(j["a_mean_at_node_km"] ** 3 / MU)


def test_A3_full_force_rates():
    s = jload("partA_summary.json")["full_force_53deg"]
    # J3, drag, SRP barely move the node: within 1 % of J2 analytic
    assert s["raan_rate_numeric_deg_day"] == pytest.approx(s["raan_rate_J2_analytic_deg_day"], rel=1e-2)
    assert s["nodal_period_numeric_s"] == pytest.approx(s["nodal_period_J2_analytic_s"], rel=1e-3)
    assert 5650 < s["nodal_period_numeric_s"] < 5800
    # mean altitude must really be ~550 km (flag if the osculating start shifted it > 10 km)
    assert abs(s["a_mean_start_km"] - A0) < 10.0


def test_A3_decay_vs_analytic_and_units():
    s = jload("partA_summary.json")["full_force_53deg"]
    assert s["a_decay_numeric_m_day"] < 0
    # numeric vs code analytic (same density): within 10 % (eccentricity / altitude
    # oscillation of +-10 km in an H~64 km atmosphere, SRP, J3)
    assert s["a_decay_numeric_m_day"] == pytest.approx(s["a_decay_analytic_m_day"], rel=0.10)
    # independent SI recomputation at the run's mean a: da/dt = -sqrt(mu a) rho Cd A/m F
    a = s["a_mean_start_km"] * 1e3
    mu = MU * 1e9
    i = math.radians(53)
    F = (1 - WE * a * math.cos(i) / math.sqrt(mu / a)) ** 2
    ref = -math.sqrt(mu * a) * vallado_rho(a / 1e3 - RE) * 2.2 * 0.02 * F * DAY   # m/day
    assert s["a_decay_analytic_m_day"] == pytest.approx(ref, rel=0.02)
    assert s["a_decay_numeric_m_day"] == pytest.approx(ref, rel=0.12)
    # total drop consistent with rate x duration (within 15 %)
    assert s["a_drop_total_km"] * 1e3 == pytest.approx(s["a_decay_numeric_m_day"] * s["days"], rel=0.15)


def test_A3_elements_physical():
    el = csv_array("partA_osculating_elements.csv")
    assert np.all(el["e"] < 5e-3)                     # stays near-circular
    assert np.all(np.abs(el["i_deg"] - 53) < 0.1)     # J2 short-period i oscillation << 0.1 deg
    alt = el["radius_minus_Re_km"]
    assert alt.min() > 520 and alt.max() < 580


def test_A3_lifetime_independent():
    L = jload("partA_summary.json")["lifetime_rough"]
    h_end = L["h_end_km"]
    i = math.radians(53)
    mu = MU * 1e9

    def inv_rate(h):  # s per m of a
        a = (RE + h) * 1e3
        F = (1 - WE * a * math.cos(i) / math.sqrt(mu / a)) ** 2
        return 1.0 / (math.sqrt(mu * a) * vallado_rho(h) * 2.2 * 0.02 * F)

    bps = [x for x in (180, 200, 250, 300, 350, 400, 450, 500) if h_end < x < H]
    t, _ = quad(lambda h: inv_rate(h) * 1e3, h_end, H, points=bps, limit=200)
    assert L["lifetime_days"] == pytest.approx(t / DAY, rel=0.02)
    assert 1.0 * 365 < L["lifetime_days"] < 30 * 365


def test_A4_sso_inclination():
    A = jload("partA_summary.json")
    om = 2 * math.pi / (365.2422 * DAY)
    n = math.sqrt(MU / A0 ** 3)
    i_ref = math.degrees(math.acos(-om / (1.5 * n * J2 * (RE / A0) ** 2)))
    assert i_ref == pytest.approx(97.59, abs=0.02)          # literature ~97.6 deg at 550 km
    for key in ("sso_mean_550", "sso_osculating_550"):
        s = A[key]
        assert s["i_analytic_deg"] == pytest.approx(i_ref, abs=1e-4)
        assert abs(s["residual_rate_deg_day"]) < 1e-4
        assert s["target_rate_deg_day"] == pytest.approx(0.98565, abs=1e-4)
    # mean-a numeric: first-order theory neglects O(J2^2) -> 0.05 deg
    assert A["sso_mean_550"]["i_numeric_deg"] == pytest.approx(i_ref, abs=0.05)
    assert 97.45 < A["sso_mean_550"]["i_numeric_deg"] < 97.75
    # osculating start lowers the mean a by a few km -> i shifts by < 0.2 deg
    assert A["sso_osculating_550"]["i_numeric_deg"] == pytest.approx(i_ref, abs=0.2)


# ------------------------------------------------------------------ Part B
def cyl(beta_deg):
    a = A0
    if abs(beta_deg) >= BETA_STAR:
        return 0.0
    return math.acos(math.sqrt(H * H + 2 * RE * H) / (a * math.cos(math.radians(beta_deg)))) / math.pi


def penumbra_fraction(beta_deg):
    """Orbit fraction in penumbra only, conical shadow, circular orbit of radius A0."""
    b = math.asin(RE / A0)
    a = math.radians(SUN_ANG_RAD_DEG)
    cb = np.cos(np.radians(np.asarray(beta_deg, float)))
    th = lambda psi: np.arccos(np.clip(math.cos(psi) / cb, -1, 1))
    return (th(b + a) - th(b - a)) / math.pi


@pytest.mark.parametrize("name", ["partB_inclined53.csv", "partB_dawndusk_sso.csv"])
def test_B_eclipse_geometry(name):
    b = csv_array(name)
    beta = b["beta_deg"]
    fu, fa, fc = b["frac_umbra"], b["frac_umbra_plus_penumbra"], b["frac_cylindrical"]
    assert np.all(np.abs(beta) <= 90)
    # cylindrical column equals the closed form recomputed here
    ref = np.array([cyl(x) for x in beta])
    assert np.allclose(fc, ref, atol=2e-5)
    # beta* = asin(R/(R+h)): no umbra beyond it; penumbra at most ~2 solar radii further
    assert np.all(fu[np.abs(beta) >= BETA_STAR] == 0)
    assert np.all(fa[np.abs(beta) >= BETA_STAR + 2 * SUN_ANG_RAD_DEG + 0.1] == 0)
    # ordering: umbra <= umbra+penumbra; conical umbra slightly inside the cylinder
    assert np.all(fu <= fa + 1e-12)
    inner = np.abs(beta) < BETA_STAR - 5
    assert np.all(np.abs(fu[inner] - fc[inner]) < 0.01)
    # REVISED after the run (review 2026-11-14). The pre-run assertion was a flat
    # "penumbra fraction < 0.01 (~10 s per crossing)" for |beta| < beta* - 5 deg. That was
    # wrong physics: the penumbra arc widens as the track grazes the shadow, and the
    # cone geometry gives 0.0103 at beta = 61.9 deg (inside the old band), ~0.003 at
    # beta = 0. Replaced by the analytic penumbra width (sun at infinity, apparent radii
    # b = asin(R/r), a = 0.267 deg; point on orbit at angle theta from midnight has
    # cos(psi) = cos(beta) cos(theta)); the umbra/penumbra boundaries are psi = b -/+ a.
    # Tolerance: sampling 5e-5, solar-radius seasonal +-1.7 %, parallax 3e-3 deg.
    near = np.abs(beta) < BETA_STAR - 1
    assert np.all(np.abs((fa - fu)[near] - penumbra_fraction(beta[near]))
                  < 1e-3 + 0.03 * penumbra_fraction(beta[near]))
    assert np.all(fa[inner] - fu[inner] < 0.0125)     # coarse cap (value at beta* - 5)
    assert np.all(fu[inner] <= fc[inner] + 1e-3)
    assert np.all(fa[inner] >= fc[inner] - 1e-3)
    assert np.all(fu <= 0.3724 + 1e-3)                # beta=0 maximum 0.3723


def test_B_beta0_duration():
    s = jload("partB_summary.json")["inclined53"]
    T = 2 * math.pi * math.sqrt(A0 ** 3 / MU) / 60
    assert s["nodal_period_min"] == pytest.approx(T, rel=5e-3)
    if s["days"] < 365:
        pytest.skip("year-long beta coverage only in the full run")
    # beta crosses zero during 2026 -> max umbra ~ 0.372 * T (~35.6 min), conical slightly less
    assert 34.5 < s["max_umbra_min"] < 36.0
    assert 70 < max(abs(s["beta_min_deg"]), abs(s["beta_max_deg"])) <= 53 + 23.45 + 0.5


def test_B_dawn_dusk_year():
    s = jload("partB_summary.json")["dawndusk_sso"]
    if s["days"] < 365:
        pytest.skip("full year only")
    # beta = asin(sin(i + delta_sun)) for an exact 18:00 LTAN; min near June solstice
    # ~ 180 - 97.6 - 23.44 = 59 deg (+- equation-of-time / drift effects)
    assert 55 < s["beta_min_deg"] < 64
    assert s["beta_max_deg"] > 80
    # eclipse only when beta < 67 deg: roughly May-August, i.e. 70-130 days; free most of year
    assert 70 <= s["days_with_eclipse"] <= 130
    b = csv_array("partB_dawndusk_sso.csv")
    ecl_days = b["day"][b["frac_umbra_plus_penumbra"] > 0]
    assert np.all(np.abs(ecl_days - 171) < 80)        # clustered around the June solstice
    assert s["max_umbra_min"] == pytest.approx(cyl(s["beta_min_deg"]) * s["nodal_period_min"], abs=1.0)


def test_B_propagated_vs_daily_geometry():
    rows = csv_rows("partB_propagated_check.csv")
    assert len(rows) > 0
    daily = {n: csv_array(f"partB_{n}.csv") for n in ("inclined53", "dawndusk_sso")}
    diffs_all, diffs_inner = [], []
    for r in rows:
        b = daily[r["orbit"]]
        d = float(r["day"])
        up, ap = float(r["umbra_min_propagated"]), float(r["umbra_plus_penumbra_min_propagated"])
        assert up <= ap + 1e-9
        if d > b["day"][-1]:
            continue
        g = np.interp(d, b["day"], b["min_umbra"])    # daily geometry interpolated in time
        diffs_all.append(abs(up - g))
        k = int(math.floor(d))
        nb = b["beta_deg"][k:k + 2]
        if np.all(np.abs(nb) < BETA_STAR - 3):        # away from the steep shadow edge
            diffs_inner.append(abs(up - g))
    # +-10 km radius oscillation, 2 s sampling, secular-RAAN approximation: < 1 min/orbit
    if diffs_inner:
        assert max(diffs_inner) < 1.0
    assert np.mean(diffs_all) < 0.5


# ------------------------------------------------------------------ Part C
def gas_reference(N, P, days, d=1.0, dh=10.0, inc=53.0):
    r = A0
    v = math.sqrt(MU / r)
    sig = math.pi * d * d
    T = days * DAY
    naive = 0.5 * N * N * sig * (4 * v / math.pi) / (4 * math.pi * r * r * dh) * T
    K = ellipk(math.sin(math.radians(inc)) ** 2)
    lat = N * N * sig * v * K / (2 * math.pi ** 3 * r * r * dh) * T
    corr = lat * (1 - 4 * d / (3 * math.pi * dh)) * (1 - 1.0 / P)
    return naive, lat, corr


def slow_fraction_model(v0=1.0, inc=53.0):
    from scipy.optimize import brentq
    v = math.sqrt(MU / A0)
    s2 = math.sin(math.radians(inc)) ** 2
    vr = lambda p: 2 * v * math.sqrt(max(s2 - math.sin(p) ** 2, 0.0)) / math.cos(p)
    p1 = brentq(lambda p: vr(p) - v0, 0.0, math.radians(inc) - 1e-12)
    f, _ = quad(lambda p: 1 / math.sqrt(s2 - math.sin(p) ** 2), p1, math.radians(inc), limit=200)
    return f / ellipk(s2)


def test_C_gas_estimates_independent():
    for r in csv_rows("partC_conjunctions.csv"):
        N, P, days = int(r["N"]), int(r["planes"]), float(r["days"])
        assert P == round(math.sqrt(N))
        naive, lat, corr = gas_reference(N, P, days)
        assert float(r["gas_naive_uniform"]) == pytest.approx(naive, rel=1e-3)
        assert float(r["gas_latitude"]) == pytest.approx(lat, rel=1e-3)
        assert float(r["gas_latitude_corrected"]) == pytest.approx(corr, rel=1e-3)


def test_C_counts_match_gas():
    """Falsifiable: with uniformly random independent phases, the EXPECTED count equals
    the latitude-resolved gas rate (exact in expectation); thin-shell, discrete-plane and
    linear-TCA effects are O(1/P, d/dh, 1e-3). Pass: within 10 % and within 5 sigma."""
    for r in csv_rows("partC_conjunctions.csv"):
        c = int(r["count"])
        g = float(r["gas_latitude_corrected"])
        if g < 100:
            continue  # too few events to test (quick pilot)
        assert abs(c / g - 1) < 0.10 + 5 / math.sqrt(g), (r["N"], c, g)
        assert float(r["poisson_sigma"]) == pytest.approx(math.sqrt(c), rel=1e-4)
        assert float(r["ratio_to_naive"]) == pytest.approx(c / float(r["gas_naive_uniform"]), rel=1e-3)
        # slow (<1 km/s) encounters: gas model puts them in the turning-latitude band
        # where v_rel = 2 v sqrt(sin^2 i - sin^2 phi)/cos phi < 1 km/s (~4.1 %).
        # REVISED after the run (review 2026-11-14): that continuum-RAAN band only exists
        # if adjacent planes can cross slower than 1 km/s. With P discrete planes the
        # slowest cross-plane crossing is v_min(P) = 2 v sin(gamma_1/2),
        # cos gamma_1 = 1 - sin^2 i (1 - cos(2 pi/P)): 1.19 km/s at P = 32 (N = 1e3), so the
        # band is empty there and only same-plane drift passes (m/s) are slow. Apply the
        # gas-band test only when v_min(P) < 0.5 km/s; otherwise require that the slow
        # share is small (< 5 %) and test the floor with the records (test below).
        if c >= 1000:
            if vmin_cross(int(r["planes"])) < 0.5:
                assert int(r["slow_lt_1kms"]) / c == pytest.approx(slow_fraction_model(), abs=0.01)
            else:
                assert int(r["slow_lt_1kms"]) / c < 0.05


def vmin_cross(P, inc=53.0, days=0.0, dh=10.0):
    """Slowest cross-plane crossing for P evenly spaced planes. With days > 0 the plane
    spacing may have closed by the differential J2 regression of a shell of thickness
    dh: dOmega = 3.5 |Omega_dot| dh / a * T (Omega_dot ~ a^-3.5)."""
    v = math.sqrt(MU / A0)
    s2 = math.sin(math.radians(inc)) ** 2
    dO = 2 * math.pi / P - 3.5 * abs(j2_raan_rate(A0, math.radians(inc))) * dh / A0 * days * DAY
    g = math.acos(1 - s2 * (1 - math.cos(max(dO, 0.0))))
    return 2 * v * math.sin(g / 2)


def test_C_excess_decomposition():
    """Added after the run (review 2026-11-14). Counts exceed the corrected gas estimate by
    3.0 / 1.5 / 0.5 % (N = 1e3/1e4/1e5). Falsifiable explanation, checked on complete
    encounter records: (1) same-plane pairs (|da| < 1 km) drift past each other at m/s
    speeds - a class the gas model omits; they must all be slow (< 10 m/s: 1.5 n da plus
    the plane opening from differential J2 regression over 30 d); (2) the discrete-plane
    correction is 1 - f(0)/(P fbar) = 1 - pi/(2 K P), not 1 - 1/P (per-pair rate for
    rings at angle gamma is prop. to 1/cos(gamma/2); its mean over dOmega is 2K/pi);
    (3) cross-plane events never go below v_min(P); (4) cross-plane count then matches
    the gas estimate within 4 effective sigma, where sigma is inflated by pair repeats
    (phase-locked pairs meet at both plane intersections every half orbit)."""
    s2 = math.sin(math.radians(53)) ** 2
    K = ellipk(s2)
    tested = 0
    for r in csv_rows("partC_conjunctions.csv"):
        N, P, c = int(r["N"]), int(r["planes"]), int(r["count"])
        p = os.path.join(RES, f"partC_encounters_N{N}.csv")
        if not os.path.exists(p):
            continue
        with open(p) as f:
            if "complete=True" not in f.readline():
                continue
        rec = np.loadtxt(p, delimiter=",", skiprows=1, ndmin=2)
        if rec.shape[0] != c or c < 1000:
            continue
        i, j, vr = rec[:, 1].astype(np.int64), rec[:, 2].astype(np.int64), rec[:, 4]
        same = (i % P) == (j % P)
        assert np.all(vr[same] < 0.010)
        assert np.all(vr[~same] > 0.99 * vmin_cross(P, days=float(r["days"])))
        assert np.all(same[vr < 0.05])           # nothing slow except same-plane drift
        g = float(r["gas_latitude_corrected"]) / (1 - 1.0 / P) * (1 - math.pi / (2 * K * P))
        _, cnt = np.unique((i * N + j)[~same], return_counts=True)
        infl = float(np.sum(cnt.astype(float) ** 2) / np.sum(cnt))
        sig = math.sqrt(infl * (~same).sum())
        assert abs((~same).sum() - g) < 4 * sig, (N, int((~same).sum()), g, sig)
        tested += 1
    if tested == 0:
        pytest.skip("no complete encounter records")


def test_C_count_scaling_quadratic():
    rows = csv_rows("partC_conjunctions.csv")
    if len(rows) < 2:
        pytest.skip("need >= 2 populations")
    N = np.array([float(r["N"]) for r in rows])
    C = np.array([float(r["count"]) for r in rows])
    if C.min() < 30:
        pytest.skip("too few events for a slope")
    slope = np.polyfit(np.log(N), np.log(C), 1)[0]
    assert 1.95 < slope < 2.05
    s = jload("partCD_summary.json")
    assert s["count_vs_N_loglog_slope"] == pytest.approx(slope, abs=1e-3)


def test_C_miss_distance_and_latitude_distributions():
    rows = csv_rows("partC_conjunctions.csv")
    Nmax = max(int(r["N"]) for r in rows)
    h = np.load(need(f"partC_hist_N{Nmax}.npz"))
    hd, hl = h["hist_d"].astype(float), h["hist_lat"].astype(float)
    if hd.sum() < 2000:
        pytest.skip("too few events")
    # shell (10 km) >> d (1 km): offsets uniform in the 2-D plane normal to v_rel
    # -> P(miss < d/2) = 1/4 (thickness correction -> 0.255). NOT 1/2.
    frac_half = hd[: hd.size // 2].sum() / hd.sum()
    assert 0.22 < frac_half < 0.29
    # latitude: pdf ~ 1/sqrt(sin^2 i - sin^2 phi), none beyond |phi| > i
    lat = np.arange(-89.5, 90, 1.0)
    assert hl[np.abs(lat) > 53.5].sum() <= 1e-3 * hl.sum()
    s2 = math.sin(math.radians(53)) ** 2
    K = ellipk(s2)
    Fin, _ = quad(lambda p: 1 / math.sqrt(s2 - math.sin(p) ** 2), 0, math.radians(45))
    model_frac_gt45 = 1 - Fin / K
    sim = hl[np.abs(lat) > 45].sum() / hl.sum()
    assert sim == pytest.approx(model_frac_gt45, abs=0.03)


# ------------------------------------------------------------------ Part D
def test_D_bookkeeping_and_bounds():
    rows = csv_rows("partD_delta_v.csv")
    conj = {int(r["N"]): r for r in csv_rows("partC_conjunctions.csv")}
    for r in rows:
        N, thr, T = int(r["N"]), float(r["threshold_km"]), float(r["lead_days"]) * DAY
        days = float(conj[N]["days"])
        m = int(r["maneuvers_30d"])
        assert float(r["target_miss_km"]) == pytest.approx(2 * thr)
        assert float(r["maneuvers_per_sat_per_year"]) == pytest.approx(m * 365.25 / days / N, rel=1e-4)
        assert float(r["annual_dv_per_sat_m_s"]) == pytest.approx(
            float(r["maneuvers_per_sat_per_year"]) * float(r["mean_dv_per_maneuver_mm_s"]) / 1e3,
            rel=1e-3, abs=1e-12)
        c = int(conj[N]["count"])
        if thr == 1.0:
            assert m == c
        # 2-D disk of offsets: maneuvers(thr)/count ~ thr^2
        if c * thr * thr > 400:
            assert m / c == pytest.approx(thr * thr, rel=0.25)
        if m >= 20:
            # CW along-track burn: displacement ~ 3 T dv; the component normal to v_rel is
            # between sin(90-53)=0.6 and 1 of it for i=53 crossings. Required change is
            # between thr (miss at threshold) and 2 thr. Bounds on mean dv * T / thr:
            x = float(r["mean_dv_per_maneuver_mm_s"]) * 1e-6 * T / thr
            n = math.sqrt(MU / A0 ** 3)
            lo = 1.0 * T / (3 * T + 8 / n)
            hi = 1.3 * 2.0 / (0.6 * 3)
            assert lo < x < hi, (N, thr, r["lead_days"], x)


def test_D_scaling():
    rows = csv_rows("partD_delta_v.csv")
    big = [r for r in rows if int(r["maneuvers_30d"]) >= 200]
    by = {(int(r["N"]), float(r["threshold_km"]), float(r["lead_days"])): float(r["annual_dv_per_sat_m_s"])
          for r in big}
    # dv ~ 1/T at fixed (N, thr)
    for (N, thr, L), v in by.items():
        if (N, thr, 2 * L) in by:
            assert by[(N, thr, 2 * L)] / v == pytest.approx(0.5, rel=0.1)
    # per-satellite annual dv ~ N (count ~ N^2)
    Ns = sorted({k[0] for k in by})
    for a, b in zip(Ns, Ns[1:]):
        for (N, thr, L), v in by.items():
            if N == a and (b, thr, L) in by:
                assert by[(b, thr, L)] / v == pytest.approx(b / a, rel=0.2)
