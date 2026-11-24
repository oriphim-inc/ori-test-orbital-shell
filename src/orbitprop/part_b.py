"""Part B: eclipse duration vs beta angle over a year."""
import math
import numpy as np
from numba import njit

from .constants import MU, RE, DAY, JD_EPOCH
from .forces import sun_position, shadow_fraction, make_params
from .orbits import j2_secular_rates, propagate, ascending_node_event, coe2rv


def orbit_basis(raan, inc):
    P = np.array([math.cos(raan), math.sin(raan), 0.0])
    Q = np.array([-math.cos(inc) * math.sin(raan), math.cos(inc) * math.cos(raan), math.sin(inc)])
    W = np.cross(P, Q)
    return P, Q, W


@njit(cache=True)
def _shadow_over_orbit(P, Q, a, rs, nsamp):
    nu = np.empty(nsamp)
    r = np.empty(3)
    for k in range(nsamp):
        u = 2 * math.pi * k / nsamp
        cu, su = math.cos(u), math.sin(u)
        for j in range(3):
            r[j] = a * (cu * P[j] + su * Q[j])
        nu[k] = shadow_fraction(r, rs)
    return nu


def cylindrical_fraction(beta, h):
    """Closed-form eclipse fraction, cylindrical shadow, circular orbit."""
    a = RE + h
    beta_star = math.asin(RE / a)
    if abs(beta) >= beta_star:
        return 0.0
    return math.acos(math.sqrt(h * h + 2 * RE * h) / (a * math.cos(beta))) / math.pi


def year_sweep(inc_deg, raan0, h=550.0, days=365, nsamp=20000):
    """Daily beta angle and eclipse duration. RAAN advanced at the J2 secular rate."""
    a = RE + h
    inc = math.radians(inc_deg)
    rd = j2_secular_rates(a, inc)[0]
    n = math.sqrt(MU / a ** 3)
    _, _, _, ud = j2_secular_rates(a, inc)
    T = 2 * math.pi / ud            # nodal period (s)
    rows = []
    for d in range(days):
        t = d * DAY
        raan = raan0 + rd * t
        P, Q, W = orbit_basis(raan, inc)
        rs = sun_position(JD_EPOCH + d)
        s = rs / np.linalg.norm(rs)
        beta = math.asin(float(np.dot(s, W)))
        nu = _shadow_over_orbit(P, Q, a, rs, nsamp)
        f_umbra = float(np.mean(nu == 0.0))
        f_any = float(np.mean(nu < 1.0))
        f_cyl = cylindrical_fraction(beta, h)
        rows.append((d, math.degrees(beta), f_umbra, f_any, f_cyl, f_umbra * T / 60, f_any * T / 60,
                     f_cyl * T / 60))
    return np.array(rows), T


HEADER = "day,beta_deg,frac_umbra,frac_umbra_plus_penumbra,frac_cylindrical,min_umbra,min_umbra_plus_penumbra,min_cylindrical"


def dawn_dusk_raan(ltan_h=18.0):
    rs = sun_position(JD_EPOCH)
    alpha = math.atan2(rs[1], rs[0])
    return alpha + math.radians((ltan_h - 12.0) * 15.0)


def propagated_check(inc_deg, raan0, days=7, h=550.0, dt=2.0, start_day=0):
    """Full-force propagation; per-orbit umbra duration from sampled shadow function,
    to compare with the daily geometric estimate."""
    from .part_a import initial_state
    y0 = initial_state(inc_deg, raan=raan0, h=h, mean_a=True)   # mean a = RE + h
    jd0 = JD_EPOCH + start_day
    p = make_params(jd0=jd0)
    t_eval = np.arange(0, days * DAY, dt)
    sol = propagate(y0, days * DAY, p, t_eval=t_eval, events=ascending_node_event,
                    rtol=1e-10, atol=1e-8)
    nu = np.array([shadow_fraction(sol.y[:3, k], sun_position(jd0 + sol.t[k] / DAY))
                   for k in range(sol.t.size)])
    tn = sol.t_events[0]
    out = []
    for k in range(tn.size - 1):
        m = (sol.t >= tn[k]) & (sol.t < tn[k + 1])
        out.append((0.5 * (tn[k] + tn[k + 1]) / DAY + start_day, np.sum(nu[m] == 0.0) * dt / 60,
                    np.sum(nu[m] < 1.0) * dt / 60))
    return np.array(out)
