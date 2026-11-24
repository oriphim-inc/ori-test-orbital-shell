"""Physical constants (km, s, kg). Sources in PLAN.md."""
import math

MU = 398600.4418            # km^3/s^2, EGM96 / WGS-84
RE = 6378.137               # km, WGS-84 equatorial radius
C20_NORM = -4.841653717e-4  # EGM96 normalised C20
C30_NORM = 9.572541738e-7   # EGM96 normalised C30
J2 = -math.sqrt(5.0) * C20_NORM   # 1.0826267e-3
J3 = -math.sqrt(7.0) * C30_NORM   # -2.53266e-6
OMEGA_E = 7.292115e-5       # rad/s, Earth rotation rate
AU = 149597870.700          # km
R_SUN = 696000.0            # km (Montenbruck & Gill 2000)
P_SUN = 4.56e-6             # N/m^2 at 1 AU
DAY = 86400.0
TROPICAL_YEAR_DAYS = 365.2422
OMEGA_SSO = 2.0 * math.pi / (TROPICAL_YEAR_DAYS * DAY)  # rad/s, mean Sun motion

# Epoch for all parts: 2026-01-01 00:00 UTC
JD_EPOCH = 2461041.5
