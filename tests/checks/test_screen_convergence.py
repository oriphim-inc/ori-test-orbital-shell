"""Sceptic's convergence check for the Part C screener (runs project code, ~seconds).

Time-step convergence: the encounter count must not depend on the screening step.
Linear relative motion over +-dt/2 with v_rel <= 12 km/s and differential gravity
n^2 * 60 km gives a TCA-position error ~1 m at dt = 10 s, so counts at dt = 10 s and
dt = 5 s must agree to within boundary/tangency cases (< 1 % + 3 events).
Skips until results/partC_conjunctions.csv exists (i.e. after the main run).
"""
import os

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "numba"))
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(ROOT, os.environ.get("ORI_RESULTS", "results"))


def test_C_count_converged_in_dt():
    if not os.path.exists(os.path.join(RES, "partC_conjunctions.csv")):
        pytest.skip("Part C not run yet")
    from orbitprop import shell as sh
    shell = sh.make_shell(1000, seed=1)
    d, days = 3.0, 2.0       # d = 3 km to get a few thousand events from 1000 sats
    c10 = sh.run_screen(shell, days, dt=10.0, d=d, thr=(d,), nchunks=8)["count"]
    c5 = sh.run_screen(shell, days, dt=5.0, d=d, thr=(d,), nchunks=8)["count"]
    assert c10 > 300
    assert abs(c10 - c5) <= 0.01 * c5 + 3, (c10, c5)
