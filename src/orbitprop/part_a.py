"""Part A: Kepler verification, J2 check, full-force 550 km orbit, SSO inclination."""
import math
import numpy as np
from scipy.integrate import solve_ivp

from .constants import MU, RE, J2, DAY, OMEGA_SSO, OMEGA_E
from .forces import make_params, density, accel
from .orbits import (coe2rv, rv2coe, kepler_state, propagate, ascending_node_event,
                     j2_secular_rates)

H0 = 550.0
A0 = RE + H0


def initial_state(inc_deg, raan=0.0, h=H0, mean_a=False):
    """Circular state at the ascending node. mean_a=False: osculating radius R+h.
    mean_a=True: radius adjusted (fixed-point, J2-only 1-day arcs) so the orbit-averaged
    semi-major axis equals R+h."""
    r0 = RE + h
    if mean_a:
        p = make_params(j2=True, j3=False, drag=False, srp=False)
        for _ in range(6):
            r, v = coe2rv(r0, 0.0, math.radians(inc_deg), raan, 0.0, 0.0)
            sol = propagate(np.r_[r, v], DAY, p, events=ascending_node_event, dense_output=True)
            a_m = orbit_means(sol)['a'].mean()
            r0 += (RE + h) - a_m
    r, v = coe2rv(r0, 0.0, math.radians(inc_deg), raan, 0.0, 0.0)
    return np.r_[r, v]


def orbit_means(sol, npts=400):
    """Average osculating elements over each node-to-node revolution (removes the
    short-period J2 terms to first order). Needs dense_output and node events."""
    tn = sol.t_events[0]
    out = dict(t=[], a=[], e=[], i=[], raan=[])
    for k in range(tn.size - 1):
        tt = np.linspace(tn[k], tn[k + 1], npts, endpoint=False)
        Y = sol.sol(tt)
        el = rv2coe(Y[:3].T, Y[3:].T)
        out['t'].append(0.5 * (tn[k] + tn[k + 1]))
        for key in ('a', 'e', 'i'):
            out[key].append(el[key].mean())
        rr = np.unwrap(el['raan'])
        out['raan'].append(rr.mean())
    o = {k: np.array(v) for k, v in out.items()}
    o['raan'] = np.unwrap(o['raan'])
    return o


def kepler_check(days, inc_deg=53.0):
    y0 = initial_state(inc_deg)
    p = make_params(j2=False, j3=False, drag=False, srp=False)
    T = days * DAY
    t_eval = np.linspace(0, T, int(days * 24 * 6) + 1)
    sol = propagate(y0, T, p, t_eval=t_eval, events=ascending_node_event)
    r, v = sol.y[:3].T, sol.y[3:].T
    eps = 0.5 * np.sum(v * v, 1) - MU / np.linalg.norm(r, axis=1)
    h = np.linalg.norm(np.cross(r, v), axis=1)
    rk, vk = kepler_state(y0[:3], y0[3:], T)
    tn = sol.t_events[0]
    period_num = float(np.mean(np.diff(tn))) if tn.size > 1 else float('nan')
    period_kep = 2 * math.pi * math.sqrt(A0 ** 3 / MU)
    return dict(
        days=days, period_kepler_s=period_kep, period_numeric_s=period_num,
        period_rel_err=abs(period_num - period_kep) / period_kep,
        pos_err_km=float(np.linalg.norm(sol.y[:3, -1] - rk)),
        vel_err_kms=float(np.linalg.norm(sol.y[3:, -1] - vk)),
        energy_rel_drift=float(np.max(np.abs(eps - eps[0])) / abs(eps[0])),
        angmom_rel_drift=float(np.max(np.abs(h - h[0])) / h[0]),
        nfev=int(sol.nfev))


def node_fit(sol):
    """Linear fits of orbit-averaged (mean) elements: raan rate, a-decay; nodal period
    from node-to-node intervals."""
    tn = sol.t_events[0]
    m = orbit_means(sol)
    rate = np.polyfit(m['t'], m['raan'], 1)[0]
    da = np.polyfit(m['t'], m['a'], 1)[0]
    return dict(t_nodes=tn, mean=m, raan_rate=rate,
                nodal_period=float(np.mean(np.diff(tn))), a_rate=da)


def j2_check(days, inc_deg=53.0):
    y0 = initial_state(inc_deg, mean_a=True)
    p = make_params(j2=True, j3=False, drag=False, srp=False)
    sol = propagate(y0, days * DAY, p, events=ascending_node_event, dense_output=True)
    f = node_fit(sol)
    a_mean = float(np.mean(f['mean']['a']))
    i_mean = float(np.mean(f['mean']['i']))
    rd, wd, Md, ud = j2_secular_rates(a_mean, i_mean)
    return dict(days=days, raan_rate_numeric_deg_day=math.degrees(f['raan_rate']) * DAY,
                raan_rate_analytic_deg_day=math.degrees(rd) * DAY,
                raan_rate_rel_diff=abs(f['raan_rate'] - rd) / abs(rd),
                nodal_period_numeric_s=f['nodal_period'],
                nodal_period_analytic_s=2 * math.pi / ud,
                kepler_period_s=2 * math.pi * math.sqrt(A0 ** 3 / MU),
                a_mean_at_node_km=a_mean)


def full_run(days, inc_deg=53.0, sample_s=60.0, cd=2.2, am=0.02, cr=1.5):
    y0 = initial_state(inc_deg, mean_a=True)   # mean altitude 550 km (as Parts B, C)
    p = make_params(cd=cd, am=am, cr=cr)
    T = days * DAY
    t_eval = np.arange(0, T + 1e-9, sample_s)
    sol = propagate(y0, T, p, t_eval=t_eval, events=ascending_node_event, dense_output=True)
    el = rv2coe(sol.y[:3].T, sol.y[3:].T)
    f = node_fit(sol)
    a_mean = float(np.mean(f['mean']['a'][:15]))   # first day
    i_rad = float(np.mean(f['mean']['i'][:15]))
    rd, _, _, ud = j2_secular_rates(a_mean, i_rad)
    # analytic circular decay: da/dt = -sqrt(mu a) rho Cd A/m * F, F wind factor
    v = math.sqrt(MU / a_mean)
    F = (1 - OMEGA_E * a_mean * math.cos(i_rad) / v) ** 2
    da_an = -math.sqrt(MU * a_mean) * density(a_mean - RE) * cd * am * 1e3 * F
    return dict(sol_t=sol.t, el=el, alt=np.linalg.norm(sol.y[:3].T, axis=1) - RE, fit=f,
                summary=dict(
                    days=days, inc_deg=inc_deg, a_mean_start_km=a_mean,
                    raan_rate_numeric_deg_day=math.degrees(f['raan_rate']) * DAY,
                    raan_rate_J2_analytic_deg_day=math.degrees(rd) * DAY,
                    nodal_period_numeric_s=f['nodal_period'],
                    nodal_period_J2_analytic_s=2 * math.pi / ud,
                    kepler_period_s=2 * math.pi * math.sqrt(A0 ** 3 / MU),
                    a_decay_numeric_m_day=f['a_rate'] * 1e3 * DAY,
                    a_decay_analytic_m_day=da_an * 1e3 * DAY,
                    a_drop_total_km=float(f['mean']['a'][-1] - f['mean']['a'][0]),
                    nfev=int(sol.nfev)))


def measured_raan_rate(inc_deg, days, mean_a=False, **kw):
    y0 = initial_state(inc_deg, mean_a=mean_a)
    p = make_params(**kw)
    sol = propagate(y0, days * DAY, p, events=ascending_node_event, dense_output=True)
    return node_fit(sol)['raan_rate']


def sso_inclination(days=5.0, tol_deg=1e-5, mean_a=False):
    a = A0
    n = math.sqrt(MU / a ** 3)
    ci = -OMEGA_SSO / (1.5 * n * J2 * (RE / a) ** 2)
    i_an = math.degrees(math.acos(ci))
    # secant root find on numerically measured nodal rate (full force model)
    f = lambda i: measured_raan_rate(i, days, mean_a=mean_a) - OMEGA_SSO
    x0, x1 = i_an, i_an + 0.05
    f0, f1 = f(x0), f(x1)
    hist = [(x0, f0), (x1, f1)]
    for _ in range(20):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0 = x1, f1
        x1, f1 = x2, f(x2)
        hist.append((x1, f1))
        if abs(x1 - x0) < tol_deg:
            break
    return dict(mean_a_init=mean_a, i_analytic_deg=i_an, i_numeric_deg=x1,
                residual_rate_deg_day=math.degrees(f1) * DAY,
                target_rate_deg_day=math.degrees(OMEGA_SSO) * DAY, arc_days=days,
                iterations=[(float(x), float(math.degrees(y) * DAY)) for x, y in hist])


def lifetime_rough(cd=2.2, am=0.02, inc_deg=53.0, h_end=150.0):
    """Orbit-averaged circular decay with the same static exponential atmosphere."""
    i_rad = math.radians(inc_deg)

    def rhs(t, a):
        v = math.sqrt(MU / a[0])
        F = (1 - OMEGA_E * a[0] * math.cos(i_rad) / v) ** 2
        return [-math.sqrt(MU * a[0]) * density(a[0] - RE) * cd * am * 1e3 * F]

    ev = lambda t, a: a[0] - RE - h_end
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 200 * 365.25 * DAY), [A0], method='DOP853', events=ev,
                    rtol=1e-8, dense_output=False, max_step=5 * DAY)
    return dict(lifetime_days=float(sol.t_events[0][0] / DAY) if sol.t_events[0].size else float('nan'),
                h_end_km=h_end)
