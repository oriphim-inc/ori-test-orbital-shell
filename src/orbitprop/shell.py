"""Part C/D: shell generation, conjunction screening (hash grid + linear TCA),
kinetic-gas estimates, CW-sized avoidance delta-v.

Satellites move on circular mean-element orbits with first-order J2 secular rates
(raan_dot, u_dot); see PLAN.md for why the full propagator is not used here.
"""
import math
import numpy as np
from numba import njit, prange
from scipy.special import ellipk

from .constants import MU, RE
from .orbits import j2_secular_rates

NB_D = 100          # miss-distance histogram bins (over [0, d))
NB_V = 160          # relative-speed histogram, 0.1 km/s bins up to 16 km/s
NB_LAT = 180        # latitude histogram, 1 deg bins


def make_shell(N, h=550.0, thickness=10.0, inc_deg=53.0, seed=1, n_planes=None):
    rng = np.random.default_rng(seed)
    P = n_planes or max(1, int(round(math.sqrt(N))))
    plane = np.arange(N) % P
    raan0 = 2 * math.pi * plane / P
    u0 = rng.uniform(0, 2 * math.pi, N)
    a = RE + h + rng.uniform(-thickness / 2, thickness / 2, N)
    inc = math.radians(inc_deg)
    rd = np.empty(N)
    ud = np.empty(N)
    for k in range(N):
        r_, _, _, u_ = j2_secular_rates(a[k], inc)
        rd[k] = r_
        ud[k] = u_
    return dict(N=N, P=P, a=a, inc=inc, raan0=raan0, u0=u0, raan_dot=rd, u_dot=ud,
                plane=plane, h=h, thickness=thickness, seed=seed)


@njit(cache=True, parallel=True)
def _states(t, a, inc, raan0, u0, rd, ud, pos, vel):
    ci, si = math.cos(inc), math.sin(inc)
    for k in prange(a.shape[0]):
        O = raan0[k] + rd[k] * t
        u = u0[k] + ud[k] * t
        cO, sO, cu, su = math.cos(O), math.sin(O), math.cos(u), math.sin(u)
        Px, Py = cO, sO
        Qx, Qy, Qz = -ci * sO, ci * cO, si
        pos[k, 0] = a[k] * (cu * Px + su * Qx)
        pos[k, 1] = a[k] * (cu * Py + su * Qy)
        pos[k, 2] = a[k] * (su * Qz)
        # d/dt: in-plane motion + nodal regression
        dPx, dPy = -sO, cO
        dQx, dQy = -ci * cO, -ci * sO
        vu = a[k] * ud[k]
        vO = a[k] * rd[k]
        vel[k, 0] = vu * (-su * Px + cu * Qx) + vO * (cu * dPx + su * dQx)
        vel[k, 1] = vu * (-su * Py + cu * Qy) + vO * (cu * dPy + su * dQy)
        vel[k, 2] = vu * (cu * Qz)


@njit(cache=True)
def _state1(k, t, a, inc, raan0, u0, rd, ud, r, v):
    ci, si = math.cos(inc), math.sin(inc)
    O = raan0[k] + rd[k] * t
    u = u0[k] + ud[k] * t
    cO, sO, cu, su = math.cos(O), math.sin(O), math.cos(u), math.sin(u)
    r[0] = a[k] * (cu * cO - su * ci * sO)
    r[1] = a[k] * (cu * sO + su * ci * cO)
    r[2] = a[k] * su * si
    vu = a[k] * ud[k]
    vO = a[k] * rd[k]
    v[0] = vu * (-su * cO - cu * ci * sO) + vO * (-cu * sO - su * ci * cO)
    v[1] = vu * (-su * sO + cu * ci * cO) + vO * (cu * cO - su * ci * sO)
    v[2] = vu * cu * si


@njit(cache=True)
def _rdot(i, j, t, a, inc, raan0, u0, rd, ud, ri, vi, rj, vj):
    _state1(i, t, a, inc, raan0, u0, rd, ud, ri, vi)
    _state1(j, t, a, inc, raan0, u0, rd, ud, rj, vj)
    f = 0.0
    for q in range(3):
        f += (rj[q] - ri[q]) * (vj[q] - vi[q])
    return f


@njit(cache=True)
def _refine(i, j, tlo, thi, a, inc, raan0, u0, rd, ud, ri, vi, rj, vj):
    """Exact range minimum of pair (i, j) inside [tlo, thi): requires r.w < 0 at tlo and
    >= 0 at thi (range-rate sign change - to +), then bisection. Returns TCA (and leaves
    the exact states at TCA in ri, vi, rj, vj) or -1 if there is no minimum in the window."""
    if thi <= tlo:
        return -1.0
    flo = _rdot(i, j, tlo, a, inc, raan0, u0, rd, ud, ri, vi, rj, vj)
    fhi = _rdot(i, j, thi, a, inc, raan0, u0, rd, ud, ri, vi, rj, vj)
    if not (flo < 0.0 and fhi >= 0.0):
        return -1.0
    lo, hi = tlo, thi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        fm = _rdot(i, j, mid, a, inc, raan0, u0, rd, ud, ri, vi, rj, vj)
        if fm < 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-6:
            break
    tc = 0.5 * (lo + hi)
    _rdot(i, j, tc, a, inc, raan0, u0, rd, ud, ri, vi, rj, vj)
    return tc


@njit(cache=True, inline='always')
def _hash(cx, cy, cz, mask):
    return ((cx * 73856093) ^ (cy * 19349663) ^ (cz * 83492791)) & mask


@njit(cache=True)
def _cw_dv(m, w, rk, vk, dsafe, T, sign):
    """Minimum |dv| (km/s) of an along-track burn at lead time T on satellite k so the
    miss vector m (from k's partner frame: m = r_other - r_k, sign handles which one)
    projected perpendicular to w reaches dsafe. Clohessy-Wiltshire response."""
    rn = math.sqrt(rk[0] ** 2 + rk[1] ** 2 + rk[2] ** 2)
    n = math.sqrt(MU / rn ** 3)
    # LVLH: radial xr, along-track ya = h x r / |..|
    hx = rk[1] * vk[2] - rk[2] * vk[1]
    hy = rk[2] * vk[0] - rk[0] * vk[2]
    hz = rk[0] * vk[1] - rk[1] * vk[0]
    ax = hy * rk[2] - hz * rk[1]
    ay = hz * rk[0] - hx * rk[2]
    az = hx * rk[1] - hy * rk[0]
    an = math.sqrt(ax * ax + ay * ay + az * az)
    cr = 2.0 / n * (1 - math.cos(n * T))
    ca = 4.0 / n * math.sin(n * T) - 3.0 * T
    g0 = sign * (cr * rk[0] / rn + ca * ax / an)
    g1 = sign * (cr * rk[1] / rn + ca * ay / an)
    g2 = sign * (cr * rk[2] / rn + ca * az / an)
    ww = w[0] ** 2 + w[1] ** 2 + w[2] ** 2
    if ww > 0:
        gw = (g0 * w[0] + g1 * w[1] + g2 * w[2]) / ww
        g0 -= gw * w[0]
        g1 -= gw * w[1]
        g2 -= gw * w[2]
    gg = g0 * g0 + g1 * g1 + g2 * g2
    mg = m[0] * g0 + m[1] * g1 + m[2] * g2
    mm = m[0] ** 2 + m[1] ** 2 + m[2] ** 2
    if gg <= 0:
        return np.inf
    disc = mg * mg + gg * (dsafe * dsafe - mm)
    if disc < 0:
        return 0.0
    return (math.sqrt(disc) - abs(mg)) / gg


@njit(cache=True, parallel=True)
def screen(a, inc, raan0, u0, rd, ud, t_end, dt, d, vrel_max, thr, safe_mult, leads,
           nchunks, cap):
    N = a.shape[0]
    Rs = d + 0.5 * vrel_max * dt + 0.5
    L = Rs
    tsize = 1
    while tsize < 2 * N:
        tsize *= 2
    mask = tsize - 1
    head = np.empty(tsize, np.int64)
    nxt = np.empty(N, np.int64)
    pos = np.empty((N, 3))
    vel = np.empty((N, 3))
    cell = np.empty((N, 3), np.int64)
    nthr = thr.shape[0]
    nlead = leads.shape[0]
    cnt = np.zeros(nchunks, np.int64)
    hist_d = np.zeros((nchunks, NB_D), np.int64)
    hist_v = np.zeros((nchunks, NB_V), np.int64)
    hist_lat = np.zeros((nchunks, NB_LAT), np.int64)
    per_sat = np.zeros((nchunks, N), np.int32)
    dv_sum = np.zeros((nchunks, nthr, nlead))
    dv_n = np.zeros((nchunks, nthr), np.int64)
    rec = np.zeros((nchunks, cap, 6 + nthr * nlead))
    nrec = np.zeros(nchunks, np.int64)
    nsteps = int(math.ceil(t_end / dt))
    csize = (N + nchunks - 1) // nchunks
    for s in range(nsteps + 1):
        t = s * dt
        _states(t, a, inc, raan0, u0, rd, ud, pos, vel)
        head[:] = -1
        for k in range(N):
            cx = int(math.floor(pos[k, 0] / L))
            cy = int(math.floor(pos[k, 1] / L))
            cz = int(math.floor(pos[k, 2] / L))
            cell[k, 0] = cx
            cell[k, 1] = cy
            cell[k, 2] = cz
            hh = _hash(cx, cy, cz, mask)
            nxt[k] = head[hh]
            head[hh] = k
        for c in prange(nchunks):
            m = np.empty(3)
            w = np.empty(3)
            ri_ = np.empty(3)
            vi_ = np.empty(3)
            rj_ = np.empty(3)
            vj_ = np.empty(3)
            for i in range(c * csize, min(N, (c + 1) * csize)):
                for ox in range(-1, 2):
                    for oy in range(-1, 2):
                        for oz in range(-1, 2):
                            nx = cell[i, 0] + ox
                            ny = cell[i, 1] + oy
                            nz = cell[i, 2] + oz
                            j = head[_hash(nx, ny, nz, mask)]
                            while j != -1:
                                if j > i and cell[j, 0] == nx and cell[j, 1] == ny and cell[j, 2] == nz:
                                    rx = pos[j, 0] - pos[i, 0]
                                    ry = pos[j, 1] - pos[i, 1]
                                    rz = pos[j, 2] - pos[i, 2]
                                    r2 = rx * rx + ry * ry + rz * rz
                                    if r2 < Rs * Rs:
                                        wx = vel[j, 0] - vel[i, 0]
                                        wy = vel[j, 1] - vel[i, 1]
                                        wz = vel[j, 2] - vel[i, 2]
                                        w2 = wx * wx + wy * wy + wz * wz
                                        ts = 0.0
                                        if w2 > 1e-18:
                                            ts = -(rx * wx + ry * wy + rz * wz) / w2
                                        lx = rx + wx * ts
                                        ly = ry + wy * ts
                                        lz = rz + wz * ts
                                        # slow (|w| < 0.05 km/s), rotation-dominated pairs: the
                                        # straight-line TCA is unreliable, test them directly
                                        reach = d + 0.1 + math.sqrt(w2) * dt
                                        if (abs(ts) <= dt and lx * lx + ly * ly + lz * lz < (d + 0.1) ** 2) \
                                                or (w2 < 0.0025 and r2 < reach * reach):
                                            tlo = max(0.0, t - 0.5 * dt)
                                            thi = min(t_end, t + 0.5 * dt)
                                            tst = _refine(i, j, tlo, thi, a, inc, raan0, u0, rd, ud,
                                                          ri_, vi_, rj_, vj_)
                                            if tst >= 0.0:
                                                for q3 in range(3):
                                                    m[q3] = rj_[q3] - ri_[q3]
                                                    w[q3] = vj_[q3] - vi_[q3]
                                                dm = math.sqrt(m[0] ** 2 + m[1] ** 2 + m[2] ** 2)
                                                if dm < d:
                                                    vr = math.sqrt(w[0] ** 2 + w[1] ** 2 + w[2] ** 2)
                                                    cnt[c] += 1
                                                    hist_d[c, min(NB_D - 1, int(dm / d * NB_D))] += 1
                                                    hist_v[c, min(NB_V - 1, int(vr / 0.1))] += 1
                                                    rin = math.sqrt(ri_[0] ** 2 + ri_[1] ** 2 + ri_[2] ** 2)
                                                    lat = math.degrees(math.asin(ri_[2] / rin))
                                                    hist_lat[c, min(NB_LAT - 1, int(lat + 90.0))] += 1
                                                    per_sat[c, i] += 1
                                                    per_sat[c, j] += 1
                                                    qrec = -1
                                                    if nrec[c] < cap:
                                                        q = nrec[c]
                                                        qrec = q
                                                        rec[c, q, 0] = tst
                                                        rec[c, q, 1] = i
                                                        rec[c, q, 2] = j
                                                        rec[c, q, 3] = dm
                                                        rec[c, q, 4] = vr
                                                        rec[c, q, 5] = lat
                                                        for q3 in range(nthr * nlead):
                                                            rec[c, q, 6 + q3] = np.nan
                                                        nrec[c] += 1
                                                    for it in range(nthr):
                                                        if dm < thr[it]:
                                                            dv_n[c, it] += 1
                                                            ds = safe_mult * thr[it]
                                                            for il in range(nlead):
                                                                d1 = _cw_dv(m, w, ri_, vi_, ds, leads[il], 1.0)
                                                                d2 = _cw_dv(m, w, rj_, vj_, ds, leads[il], -1.0)
                                                                best = min(d1, d2)
                                                                dv_sum[c, it, il] += best
                                                                if qrec >= 0:
                                                                    rec[c, qrec, 6 + it * nlead + il] = best
                                j = nxt[j]
    return (cnt.sum(), hist_d.sum(0), hist_v.sum(0), hist_lat.sum(0), per_sat.sum(0),
            dv_sum.sum(0), dv_n.sum(0), rec, nrec)


def vrel_bound(shell):
    v = math.sqrt(MU / shell['a'].min())
    return 2 * v * math.sin(shell['inc']) * 1.02 + 0.1


def run_screen(shell, days, dt=10.0, d=1.0, thr=(0.1, 0.2, 0.5, 1.0), safe_mult=2.0,
               leads_days=(0.5, 1.0, 2.0), nchunks=64, cap=20000):
    thr = np.asarray(thr, float)
    leads = np.asarray(leads_days, float) * 86400.0
    out = screen(shell['a'], shell['inc'], shell['raan0'], shell['u0'], shell['raan_dot'],
                 shell['u_dot'], days * 86400.0, dt, d, vrel_bound(shell), thr, safe_mult, leads,
                 nchunks, cap)
    cnt, hd, hv, hl, per_sat, dvs, dvn, rec, nrec = out
    recs = np.concatenate([rec[c, :nrec[c]] for c in range(rec.shape[0])])
    recs = recs[np.argsort(recs[:, 0])] if recs.size else recs.reshape(0, rec.shape[2])
    return dict(count=int(cnt), hist_d=hd, hist_v=hv, hist_lat=hl, per_sat=per_sat,
                dv_sum=dvs, dv_n=dvn, records=recs, records_complete=bool(np.all(nrec < cap)),
                thr=thr, leads_days=np.asarray(leads_days, float), safe_mult=safe_mult,
                days=days, dt=dt, d=d)


def brute_force(shell, days, dt, d):
    """O(N^2) reference: per step all pairs; same prefilter + exact range-minimum rule."""
    a, inc = shell['a'], shell['inc']
    args = (a, inc, shell['raan0'], shell['u0'], shell['raan_dot'], shell['u_dot'])
    N = a.size
    pos = np.empty((N, 3))
    vel = np.empty((N, 3))
    t_end = days * 86400.0
    nsteps = int(math.ceil(t_end / dt))
    found = []
    iu, ju = np.triu_indices(N, 1)
    ri, vi, rj, vj = (np.empty(3) for _ in range(4))
    for s in range(nsteps + 1):
        t = s * dt
        _states(t, *args, pos, vel)
        r = pos[ju] - pos[iu]
        w = vel[ju] - vel[iu]
        w2 = np.sum(w * w, 1)
        ts = np.where(w2 > 1e-18, -np.sum(r * w, 1) / np.maximum(w2, 1e-300), 0.0)
        dm = np.linalg.norm(r + w * ts[:, None], axis=1)
        reach = d + 0.1 + np.sqrt(w2) * dt
        cand = ((np.abs(ts) <= dt) & (dm < d + 0.1)) | ((w2 < 0.0025) & (np.sum(r * r, 1) < reach ** 2))
        for k in np.nonzero(cand)[0]:
            tc = _refine(iu[k], ju[k], max(0.0, t - dt / 2), min(t_end, t + dt / 2), *args, ri, vi, rj, vj)
            if tc >= 0:
                dd = float(np.linalg.norm(rj - ri))
                if dd < d:
                    found.append((tc, iu[k], ju[k], dd))
    return np.array(found).reshape(-1, 4)


# ---------------------------------------------------------------- kinetic gas
def gas_estimates(shell, days, d=1.0):
    N = shell['N']
    r = RE + shell['h']
    Dh = shell['thickness']
    v = math.sqrt(MU / r)
    i = shell['inc']
    T = days * 86400.0
    sigma = math.pi * d * d
    V = 4 * math.pi * r * r * Dh
    naive = 0.5 * N * N * sigma * (4 * v / math.pi) / V * T
    K = float(ellipk(math.sin(i) ** 2))
    lat = sigma * v * N * N * K / (2 * math.pi ** 3 * r * r * Dh) * T
    # derived corrections: finite thickness (uniform altitudes), same-plane pairs
    f_thick = 1 - 4 * d / (3 * math.pi * Dh)
    P = shell['P']
    f_plane = 1 - 1.0 / P
    return dict(naive_uniform=naive, latitude_resolved=lat,
                latitude_resolved_corrected=lat * f_thick * f_plane,
                f_thickness=f_thick, f_same_plane=f_plane, K=K, v=v)


def gas_latitude_pdf(inc, lat_deg):
    """Expected latitude distribution of crossing encounters (per degree), from the
    gas model: rate density per dphi proportional to 1/sqrt(sin^2 i - sin^2 phi)."""
    phi = np.radians(lat_deg)
    s = np.sin(inc) ** 2 - np.sin(phi) ** 2
    f = np.where(s > 0, 1 / np.sqrt(np.maximum(s, 1e-300)), 0.0)
    K = float(ellipk(math.sin(inc) ** 2))
    return f / (2 * K) * math.pi / 180
