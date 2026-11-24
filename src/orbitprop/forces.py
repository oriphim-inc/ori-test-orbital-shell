"""Force model, compiled with numba. Units km, s.

Accelerations: two-body, J2, J3 (Vallado 2013 sec. 8.7), drag with a co-rotating
exponential atmosphere (Vallado 2013 Table 8-4), cannonball SRP with conical shadow
(Montenbruck & Gill 2000 sec. 3.4). Sun: low-precision almanac (Vallado 2013 Alg. 29).
"""
import math
import numpy as np
from numba import njit

from .constants import MU, RE, J2, J3, OMEGA_E, AU, R_SUN, P_SUN, JD_EPOCH

# Vallado (2013) Table 8-4: base altitude h0 [km], nominal density rho0 [kg/m^3],
# scale height H [km].
_ATM = np.array([
    [0, 1.225, 7.249], [25, 3.899e-2, 6.349], [30, 1.774e-2, 6.682],
    [40, 3.972e-3, 7.554], [50, 1.057e-3, 8.382], [60, 3.206e-4, 7.714],
    [70, 8.770e-5, 6.549], [80, 1.905e-5, 5.799], [90, 3.396e-6, 5.382],
    [100, 5.297e-7, 5.877], [110, 9.661e-8, 7.263], [120, 2.438e-8, 9.473],
    [130, 8.484e-9, 12.636], [140, 3.845e-9, 16.149], [150, 2.070e-9, 22.523],
    [180, 5.464e-10, 29.740], [200, 2.789e-10, 37.105], [250, 7.248e-11, 45.546],
    [300, 2.418e-11, 53.628], [350, 9.518e-12, 53.298], [400, 3.725e-12, 58.515],
    [450, 1.585e-12, 60.828], [500, 6.967e-13, 63.822], [600, 1.454e-13, 71.835],
    [700, 3.614e-14, 88.667], [800, 1.170e-14, 124.64], [900, 5.245e-15, 181.05],
    [1000, 3.019e-15, 268.00]])
ATM_H0 = np.ascontiguousarray(_ATM[:, 0])
ATM_RHO0 = np.ascontiguousarray(_ATM[:, 1])
ATM_H = np.ascontiguousarray(_ATM[:, 2])


@njit(cache=True)
def density(h_km):
    """Exponential atmosphere density in kg/m^3 at geometric altitude h (km)."""
    if h_km < 0.0:
        h_km = 0.0
    k = 0
    n = ATM_H0.shape[0]
    for j in range(n):
        if h_km >= ATM_H0[j]:
            k = j
    return ATM_RHO0[k] * math.exp(-(h_km - ATM_H0[k]) / ATM_H[k])


@njit(cache=True)
def sun_position(jd):
    """Geocentric Sun vector (km), mean equator/equinox, Vallado (2013) Alg. 29."""
    T = (jd - 2451545.0) / 36525.0
    lam_m = 280.460 + 36000.771 * T
    M = math.radians(357.5291092 + 35999.05034 * T)
    lam = math.radians(lam_m + 1.914666471 * math.sin(M) + 0.019994643 * math.sin(2 * M))
    r = (1.000140612 - 0.016708617 * math.cos(M) - 0.000139589 * math.cos(2 * M)) * AU
    eps = math.radians(23.439291 - 0.0130042 * T)
    out = np.empty(3)
    out[0] = r * math.cos(lam)
    out[1] = r * math.cos(eps) * math.sin(lam)
    out[2] = r * math.sin(eps) * math.sin(lam)
    return out


@njit(cache=True)
def shadow_fraction(r, rs):
    """Fraction of the solar disc visible from r (0 umbra, 1 sunlit), conical model
    (Montenbruck & Gill 2000, eq. 3.79-3.87). r, rs: geocentric km."""
    d0 = rs[0] - r[0]
    d1 = rs[1] - r[1]
    d2 = rs[2] - r[2]
    dsun = math.sqrt(d0 * d0 + d1 * d1 + d2 * d2)
    rn = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
    a = math.asin(min(1.0, R_SUN / dsun))
    b = math.asin(min(1.0, RE / rn))
    cosc = -(r[0] * d0 + r[1] * d1 + r[2] * d2) / (rn * dsun)
    cosc = max(-1.0, min(1.0, cosc))
    c = math.acos(cosc)
    if c >= a + b:
        return 1.0
    if c <= b - a:
        return 0.0
    if c <= a - b:  # annular (cannot happen for Earth in LEO), keep for completeness
        return 1.0 - (b * b) / (a * a)
    x = (c * c + a * a - b * b) / (2 * c)
    y = math.sqrt(max(0.0, a * a - x * x))
    A = a * a * math.acos(max(-1.0, min(1.0, x / a))) + \
        b * b * math.acos(max(-1.0, min(1.0, (c - x) / b))) - c * y
    return 1.0 - A / (math.pi * a * a)


@njit(cache=True)
def accel(t, y, p):
    """Total acceleration. p = [use_J2, use_J3, use_drag, use_srp, CdAm, CrAm, jd0]
    with CdAm = Cd*A/m and CrAm = Cr*A/m in m^2/kg, jd0 = Julian date at t = 0."""
    x, yy, z = y[0], y[1], y[2]
    r2 = x * x + yy * yy + z * z
    r = math.sqrt(r2)
    out = np.empty(6)
    out[0] = y[3]
    out[1] = y[4]
    out[2] = y[5]
    k = -MU / (r2 * r)
    ax = k * x
    ay = k * yy
    az = k * z
    zr2 = z * z / r2
    if p[0] != 0.0:
        f = -1.5 * J2 * MU * RE * RE / (r2 * r2 * r)
        ax += f * x * (1.0 - 5.0 * zr2)
        ay += f * yy * (1.0 - 5.0 * zr2)
        az += f * z * (3.0 - 5.0 * zr2)
    if p[1] != 0.0:
        # Vallado (2013) eq. 8-30 (J3 terms)
        f = -2.5 * J3 * MU * RE ** 3 / (r2 * r2 * r2 * r)
        ax += f * x * (3.0 * z - 7.0 * z * zr2)
        ay += f * yy * (3.0 * z - 7.0 * z * zr2)
        az += f * (6.0 * z * z - 7.0 * z * z * zr2 - 0.6 * r2)
    if p[2] != 0.0:
        vrx = y[3] + OMEGA_E * yy
        vry = y[4] - OMEGA_E * x
        vrz = y[5]
        vr = math.sqrt(vrx * vrx + vry * vry + vrz * vrz)
        rho = density(r - 6378.137)  # kg/m^3 ; altitude above a sphere of radius RE
        # a [km/s^2] = -0.5 rho[kg/m^3] CdAm[m^2/kg] |v|v [km^2/s^2] * 1e3 (m->km: v^2 1e6 m^2/s^2 ... )
        # rho*CdAm has units 1/m; v^2 in km^2/s^2 -> times 1e3 m/km gives km/s^2
        fd = -0.5 * rho * p[4] * vr * 1e3
        ax += fd * vrx
        ay += fd * vry
        az += fd * vrz
    if p[3] != 0.0:
        rs = sun_position(p[6] + t / 86400.0)
        nu = shadow_fraction(y[:3], rs)
        if nu > 0.0:
            dx = x - rs[0]
            dy = yy - rs[1]
            dz = z - rs[2]
            dn = math.sqrt(dx * dx + dy * dy + dz * dz)
            # N/m^2 * m^2/kg = m/s^2 -> /1e3 km/s^2 ; pushes away from the Sun
            fs = nu * P_SUN * p[5] * (AU / dn) ** 2 / 1e3 / dn
            ax += fs * dx
            ay += fs * dy
            az += fs * dz
    out[3] = ax
    out[4] = ay
    out[5] = az
    return out


def make_params(j2=True, j3=True, drag=True, srp=True, cd=2.2, cr=1.5, am=0.02, jd0=JD_EPOCH):
    return np.array([float(j2), float(j3), float(drag), float(srp), cd * am, cr * am, jd0])
