import math
import os

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "numba"))
import numpy as np
import pytest

from orbitprop.constants import MU, RE, J2, J3, JD_EPOCH
from orbitprop.forces import accel, make_params, density, shadow_fraction, sun_position
from orbitprop.orbits import coe2rv, rv2coe
from orbitprop import shell as sh
from orbitprop.part_b import cylindrical_fraction


def potential(r, j2=True, j3=True):
    x, y, z = r
    rn = math.sqrt(x * x + y * y + z * z)
    s = z / rn
    U = MU / rn
    if j2:
        U -= MU / rn * J2 * (RE / rn) ** 2 * 0.5 * (3 * s * s - 1)
    if j3:
        U -= MU / rn * J3 * (RE / rn) ** 3 * 0.5 * (5 * s ** 3 - 3 * s)
    return U


@pytest.mark.parametrize("j2,j3", [(True, False), (False, True), (True, True)])
def test_zonal_accel_is_gradient(j2, j3):
    p = make_params(j2=j2, j3=j3, drag=False, srp=False)
    rng = np.random.default_rng(0)
    for _ in range(5):
        r = rng.normal(size=3)
        r = r / np.linalg.norm(r) * (RE + 550)
        a = accel(0.0, np.r_[r, 0, 0, 0], p)[3:]
        g = np.zeros(3)
        for k in range(3):
            e = np.zeros(3); e[k] = 1e-3
            g[k] = (potential(r + e, j2, j3) - potential(r - e, j2, j3)) / 2e-3
        assert np.allclose(a, g, rtol=1e-7, atol=1e-13)


def test_density_table_points():
    assert density(500.0) == pytest.approx(6.967e-13)
    assert density(550.0) == pytest.approx(6.967e-13 * math.exp(-50 / 63.822))
    # continuity of the piecewise model within ~15 % at the 600 km boundary
    assert density(599.999) / density(600.0) == pytest.approx(1.0, rel=0.15)


def test_drag_direction_and_magnitude():
    r, v = coe2rv(RE + 550, 0, 0.0, 0, 0, 0)  # equatorial prograde
    p = make_params(j2=False, j3=False, drag=True, srp=False)
    p0 = make_params(j2=False, j3=False, drag=False, srp=False)
    ad = accel(0.0, np.r_[r, v], p)[3:] - accel(0.0, np.r_[r, v], p0)[3:]
    vrel = v - np.cross([0, 0, 7.292115e-5], r)
    expect = -0.5 * density(550.0) * 2.2 * 0.02 * np.linalg.norm(vrel) * vrel * 1e3
    assert np.allclose(ad, expect, rtol=1e-12)


def test_shadow_function():
    rs = sun_position(JD_EPOCH)
    s = rs / np.linalg.norm(rs)
    assert shadow_fraction(s * (RE + 550), rs) == 1.0
    assert shadow_fraction(-s * (RE + 550), rs) == 0.0


def test_coe_roundtrip():
    r, v = coe2rv(7000.0, 0.01, 0.9, 1.0, 2.0, 0.5)
    el = rv2coe(r, v)
    assert el["a"][0] == pytest.approx(7000.0)
    assert el["e"][0] == pytest.approx(0.01)
    assert el["i"][0] == pytest.approx(0.9)
    assert el["raan"][0] == pytest.approx(1.0)
    assert el["argp"][0] == pytest.approx(2.0)
    assert el["nu"][0] == pytest.approx(0.5)


def test_cylindrical_beta0():
    h = 550.0
    f = cylindrical_fraction(0.0, h)
    assert f == pytest.approx(math.acos(math.sqrt(h * h + 2 * RE * h) / (RE + h)) / math.pi)
    assert cylindrical_fraction(math.radians(70), h) == 0.0


def test_grid_equals_bruteforce():
    shell = sh.make_shell(150, seed=7)
    r = sh.run_screen(shell, 0.2, dt=10.0, d=25.0, thr=(25.0,), nchunks=4)
    b = sh.brute_force(shell, 0.2, 10.0, 25.0)
    assert r["count"] == len(b) > 20
    gi = {tuple(x) for x in r["records"][:, 1:3].astype(int)}
    bi = {tuple(x) for x in b[:, 1:3].astype(int)}
    assert gi == bi


def test_grid_against_fine_minimum_search():
    """Independent check: 0.25 s sampling of all pair distances, count local minima
    below d (excluding ends), compare with the 10 s linear-TCA screener."""
    shell = sh.make_shell(60, seed=11, thickness=10.0)
    d, days = 40.0, 0.05
    r = sh.run_screen(shell, days, dt=10.0, d=d, thr=(d,), nchunks=2)
    N = shell["N"]
    tt = np.arange(0, days * 86400 + 1e-9, 0.25)
    pos = np.empty((tt.size, N, 3)); vel = np.empty((N, 3)); P = np.empty((N, 3))
    for k, t in enumerate(tt):
        sh._states(t, shell["a"], shell["inc"], shell["raan0"], shell["u0"], shell["raan_dot"],
                   shell["u_dot"], P, vel)
        pos[k] = P
    iu, ju = np.triu_indices(N, 1)
    D = np.linalg.norm(pos[:, ju] - pos[:, iu], axis=2)
    mins = (D[1:-1] < D[:-2]) & (D[1:-1] <= D[2:]) & (D[1:-1] < d)
    n_fine = int(mins.sum())
    assert n_fine > 5
    assert abs(r["count"] - n_fine) <= 1


def test_gas_scaling_is_quadratic():
    g1 = sh.gas_estimates(sh.make_shell(100, seed=1), 30)
    g2 = sh.gas_estimates(sh.make_shell(400, seed=1), 30)
    assert g2["latitude_resolved"] / g1["latitude_resolved"] == pytest.approx(16.0)


def test_state1_matches_vectorised_states():
    shell = sh.make_shell(20, seed=5)
    args = (shell["a"], shell["inc"], shell["raan0"], shell["u0"], shell["raan_dot"], shell["u_dot"])
    P = np.empty((20, 3)); V = np.empty((20, 3)); r = np.empty(3); v = np.empty(3)
    sh._states(1234.5, *args, P, V)
    for k in range(20):
        sh._state1(k, 1234.5, *args, r, v)
        assert np.allclose(r, P[k], atol=1e-9) and np.allclose(v, V[k], atol=1e-12)


def test_comoving_pair_counted_once():
    """Regression: same-plane pair, da = 1 m, 0.5 km along-track apart, drifting through
    each other at ~mm/s. One physical closest approach (~1 m) -> exactly one encounter."""
    from orbitprop.orbits import j2_secular_rates
    s = sh.make_shell(2, seed=1, n_planes=1)
    s["a"][:] = [RE + 550.0, RE + 550.001]
    s["u0"][:] = [0.0, 0.5 / (RE + 550.0)]
    for k in range(2):
        s["raan_dot"][k], _, _, s["u_dot"][k] = j2_secular_rates(s["a"][k], s["inc"])
    r = sh.run_screen(s, 8.0, dt=10.0, nchunks=1, thr=(1.0,))
    assert r["count"] == 1
    # miss = radial 1 m (+ ~0.7 m cross-track from the slightly different J2 nodal rate)
    assert 0.001 <= r["records"][0, 3] < 0.0013
