"""Element conversions, Kepler solution, DOP853 propagation wrapper."""
import math
import numpy as np
from scipy.integrate import solve_ivp

from .constants import MU, RE, J2
from .forces import accel


def coe2rv(a, e, i, raan, argp, nu, mu=MU):
    p = a * (1 - e * e)
    r_pf = np.array([math.cos(nu), math.sin(nu), 0.0]) * p / (1 + e * math.cos(nu))
    v_pf = np.array([-math.sin(nu), e + math.cos(nu), 0.0]) * math.sqrt(mu / p)
    cO, sO, ci, si, cw, sw = (math.cos(raan), math.sin(raan), math.cos(i), math.sin(i),
                              math.cos(argp), math.sin(argp))
    R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                  [sw * si, cw * si, ci]])
    return R @ r_pf, R @ v_pf


def rv2coe(r, v, mu=MU):
    """Vectorised osculating elements. r, v: (..., 3). Returns dict of arrays
    a, e, i, raan, argp, nu, u (argument of latitude), all angles in rad."""
    r = np.atleast_2d(r)
    v = np.atleast_2d(v)
    rn = np.linalg.norm(r, axis=-1)
    vn2 = np.sum(v * v, axis=-1)
    h = np.cross(r, v)
    hn = np.linalg.norm(h, axis=-1)
    n = np.stack([-h[:, 1], h[:, 0], np.zeros_like(hn)], axis=-1)
    nn = np.linalg.norm(n, axis=-1)
    evec = ((vn2 - mu / rn)[:, None] * r - np.sum(r * v, axis=-1)[:, None] * v) / mu
    e = np.linalg.norm(evec, axis=-1)
    a = 1.0 / (2.0 / rn - vn2 / mu)
    i = np.arccos(np.clip(h[:, 2] / hn, -1, 1))
    raan = np.mod(np.arctan2(n[:, 1], n[:, 0]), 2 * np.pi)
    # argument of latitude: angle from node to r
    cu = np.sum(n * r, axis=-1) / (nn * rn)
    su = np.sum(np.cross(n, r) * h, axis=-1) / (nn * rn * hn)
    u = np.mod(np.arctan2(su, cu), 2 * np.pi)
    cw = np.sum(n * evec, axis=-1) / (nn * np.maximum(e, 1e-300))
    sw = np.sum(np.cross(n, evec) * h, axis=-1) / (nn * np.maximum(e, 1e-300) * hn)
    argp = np.mod(np.arctan2(sw, cw), 2 * np.pi)
    nu = np.mod(u - argp, 2 * np.pi)
    return dict(a=a, e=e, i=i, raan=raan, argp=argp, nu=nu, u=u)


def kepler_solve(M, e, tol=1e-15):
    E = M if e < 0.8 else math.pi
    for _ in range(50):
        dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E


def kepler_state(r0, v0, t, mu=MU):
    """Analytic two-body state at time t from (r0, v0) via elements."""
    el = rv2coe(r0, v0, mu)
    a, e = el['a'][0], el['e'][0]
    nu0 = el['nu'][0]
    E0 = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu0 / 2), math.sqrt(1 + e) * math.cos(nu0 / 2))
    M0 = E0 - e * math.sin(E0)
    n = math.sqrt(mu / a ** 3)
    M = math.fmod(M0 + n * t, 2 * math.pi)
    E = kepler_solve(M, e)
    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    return coe2rv(a, e, el['i'][0], el['raan'][0], el['argp'][0], nu, mu)


def propagate(y0, t_end, params, t_eval=None, rtol=1e-12, atol=1e-9, events=None,
              dense_output=False):
    """Integrate with DOP853 (Dormand-Prince 8(5,3))."""
    fun = lambda t, y: accel(t, y, params)
    sol = solve_ivp(fun, (0.0, t_end), np.asarray(y0, float), method='DOP853',
                    t_eval=t_eval, rtol=rtol, atol=atol, events=events,
                    dense_output=dense_output)
    if sol.status < 0:
        raise RuntimeError(f"DOP853 failed: {sol.message}")
    return sol


def ascending_node_event(t, y):
    return y[2]


ascending_node_event.direction = 1.0


def j2_secular_rates(a, i, e=0.0):
    """First-order J2 secular rates (rad/s): raan_dot, argp_dot, M_dot, u_dot (circular)."""
    n = math.sqrt(MU / a ** 3)
    p = a * (1 - e * e)
    k = J2 * (RE / p) ** 2
    ci = math.cos(i)
    raan_dot = -1.5 * n * k * ci
    argp_dot = 0.75 * n * k * (5 * ci * ci - 1)
    M_dot = n * (1 + 0.75 * k * math.sqrt(1 - e * e) * (3 * ci * ci - 1))
    return raan_dot, argp_dot, M_dot, argp_dot + M_dot
