"""Reproduce every result: python run.py [--quick] [--parts ABCD] [--dt 10]

--quick: pilot sizes (A: 3 d, SSO arcs 2 d, B: 30 d, C/D: 1 d). Writes to results/
(or results_quick/ with --quick).
"""
import argparse
import json
import math
import os
import sys
import time

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "mpl"))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "numba"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from orbitprop.constants import DAY, RE
from orbitprop import part_a, part_b, shell as sh
from orbitprop.orbits import j2_secular_rates


def jdump(obj, path):
    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=conv)


def run_a(cfg, out):
    t0 = time.time()
    res = {}
    res["kepler_check"] = part_a.kepler_check(cfg["a_days"])
    res["j2_check"] = part_a.j2_check(cfg["a_days"])
    full = part_a.full_run(cfg["a_days"])
    res["full_force_53deg"] = full["summary"]
    res["sso_osculating_550"] = part_a.sso_inclination(cfg["sso_days"], mean_a=False)
    res["sso_mean_550"] = part_a.sso_inclination(cfg["sso_days"], mean_a=True)
    res["lifetime_rough"] = part_a.lifetime_rough()
    jdump(res, f"{out}/partA_summary.json")
    el = full["el"]
    t = full["sol_t"]
    np.savetxt(f"{out}/partA_osculating_elements.csv",
               np.c_[t, el["a"], el["e"], np.degrees(el["i"]), np.degrees(el["raan"]),
                     np.degrees(el["argp"]), np.degrees(el["u"]), full["alt"]],
               delimiter=",", fmt="%.9g",
               header="t_s,a_km,e,i_deg,raan_deg,argp_deg,u_deg,radius_minus_Re_km", comments="")
    m = full["fit"]["mean"]
    np.savetxt(f"{out}/partA_mean_elements.csv",
               np.c_[m["t"], m["a"], m["e"], np.degrees(m["i"]), np.degrees(m["raan"])],
               delimiter=",", fmt="%.10g", header="t_s,a_mean_km,e_mean,i_mean_deg,raan_mean_deg_unwrapped",
               comments="")
    fig, ax = plt.subplots(5, 1, figsize=(9, 12), sharex=True)
    td, tm = t / DAY, m["t"] / DAY
    ax[0].plot(td, el["a"], lw=0.3, color="0.7", label="osculating")
    ax[0].plot(tm, m["a"], "k", lw=1.2, label="orbit-averaged")
    ax[0].set_ylabel("a [km]"); ax[0].legend(fontsize=8)
    ax[1].plot(td, el["e"], lw=0.3, color="0.7"); ax[1].plot(tm, m["e"], "k"); ax[1].set_ylabel("e")
    ax[2].plot(td, np.degrees(el["i"]), lw=0.3, color="0.7"); ax[2].plot(tm, np.degrees(m["i"]), "k")
    ax[2].set_ylabel("i [deg]")
    ax[3].plot(tm, np.degrees(m["raan"]), "k"); ax[3].set_ylabel("RAAN (unwrapped) [deg]")
    ax[4].plot(td, full["alt"], lw=0.3, color="0.6"); ax[4].set_ylabel("|r| - Re [km]")
    ax[4].set_xlabel("days since 2026-01-01")
    s = full["summary"]
    ax[0].set_title(f"Part A: 550 km, i=53 deg, J2+J3+drag+SRP, A/m=0.02, Cd=2.2\n"
                    f"RAAN rate {s['raan_rate_numeric_deg_day']:.4f} deg/d, nodal period "
                    f"{s['nodal_period_numeric_s']:.2f} s, da/dt {s['a_decay_numeric_m_day']:.1f} m/d")
    fig.tight_layout(); fig.savefig(f"{out}/partA_elements.png", dpi=130); plt.close(fig)
    print(f"[A] done in {time.time() - t0:.1f} s")
    return res


def run_b(cfg, out, i_sso):
    t0 = time.time()
    orbits = {"inclined53": (53.0, 0.0), "dawndusk_sso": (i_sso, part_b.dawn_dusk_raan(18.0))}
    res = {}
    fig1, ax1 = plt.subplots(figsize=(8, 5))
    fig2, ax2 = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    checks = []
    for name, (inc, raan0) in orbits.items():
        rows, T = part_b.year_sweep(inc, raan0, days=cfg["b_days"])
        np.savetxt(f"{out}/partB_{name}.csv", rows, delimiter=",", fmt="%.6g", header=part_b.HEADER,
                   comments="")
        res[name] = dict(inc_deg=inc, raan0_deg=math.degrees(raan0), nodal_period_min=T / 60,
                         beta_min_deg=rows[:, 1].min(), beta_max_deg=rows[:, 1].max(),
                         max_umbra_min=rows[:, 5].max(), mean_umbra_frac=rows[:, 2].mean(),
                         days_with_eclipse=int((rows[:, 3] > 0).sum()), days=int(rows.shape[0]),
                         max_abs_conical_minus_cylindrical_frac=float(np.abs(rows[:, 2] - rows[:, 4]).max()))
        ax1.plot(rows[:, 1], rows[:, 2], ".", ms=3, label=f"{name} (umbra, conical)")
        ax2[0].plot(rows[:, 0], rows[:, 1], label=name)
        ax2[1].plot(rows[:, 0], rows[:, 5], label=f"{name} umbra")
        ax2[1].plot(rows[:, 0], rows[:, 6], ":", lw=0.8, label=f"{name} umbra+penumbra")
        for sd in cfg["b_check_days"]:
            rd = j2_secular_rates(RE + 550.0, math.radians(inc))[0]
            pc = part_b.propagated_check(inc, raan0 + rd * sd * DAY,
                                         days=cfg["b_check_len"], start_day=sd)
            for row in pc:
                d = int(math.floor(row[0]))
                if d < rows.shape[0]:
                    checks.append((name, row[0], row[1], row[2], rows[d, 5], rows[d, 6]))
    b = np.linspace(-80, 80, 801)
    ax1.plot(b, [part_b.cylindrical_fraction(math.radians(x), 550.0) for x in b], "k-", lw=1,
             label="cylindrical closed form")
    ax1.set_xlabel("beta angle [deg]"); ax1.set_ylabel("eclipse fraction of orbit")
    ax1.set_title("Part B: eclipse fraction vs beta, 550 km circular"); ax1.legend(fontsize=8)
    fig1.savefig(f"{out}/partB_eclipse_fraction_vs_beta.png", dpi=130); plt.close(fig1)
    ax2[0].set_ylabel("beta [deg]"); ax2[0].legend(fontsize=8)
    ax2[1].set_ylabel("eclipse duration [min]"); ax2[1].set_xlabel("day of 2026"); ax2[1].legend(fontsize=7)
    fig2.tight_layout(); fig2.savefig(f"{out}/partB_beta_and_eclipse_vs_day.png", dpi=130); plt.close(fig2)
    with open(f"{out}/partB_propagated_check.csv", "w") as f:
        f.write("orbit,day,umbra_min_propagated,umbra_plus_penumbra_min_propagated,umbra_min_daily_geometry,"
                "umbra_plus_penumbra_min_daily_geometry\n")
        for c in checks:
            f.write(",".join(str(x) for x in c) + "\n")
    diffs = np.array([abs(c[2] - c[4]) for c in checks])
    res["propagated_check"] = dict(n_orbits=len(checks), max_abs_umbra_diff_min=float(diffs.max()),
                                   mean_abs_umbra_diff_min=float(diffs.mean()))
    jdump(res, f"{out}/partB_summary.json")
    print(f"[B] done in {time.time() - t0:.1f} s")
    return res


def run_cd(cfg, out, dt):
    t0 = time.time()
    rows, drows, hists = [], [], {}
    for N in cfg["populations"]:
        t1 = time.time()
        shell = sh.make_shell(N, seed=cfg["seed"])
        r = sh.run_screen(shell, cfg["c_days"], dt=dt, nchunks=64 if N >= 10000 else 8)
        g = sh.gas_estimates(shell, cfg["c_days"])
        slow = int(r["hist_v"][:10].sum())  # v_rel < 1 km/s
        rows.append(dict(N=N, planes=shell["P"], days=cfg["c_days"], count=r["count"],
                         poisson_sigma=math.sqrt(r["count"]), slow_lt_1kms=slow,
                         gas_naive_uniform=g["naive_uniform"], gas_latitude=g["latitude_resolved"],
                         gas_latitude_corrected=g["latitude_resolved_corrected"],
                         ratio_to_naive=r["count"] / g["naive_uniform"],
                         ratio_to_latitude_corrected=r["count"] / g["latitude_resolved_corrected"],
                         per_sat_mean=float(r["per_sat"].mean()), per_sat_max=int(r["per_sat"].max()),
                         runtime_s=time.time() - t1))
        hists[N] = r
        scale = 365.25 / cfg["c_days"] / N
        for it, thr in enumerate(r["thr"]):
            for il, lead in enumerate(r["leads_days"]):
                n_m = r["dv_n"][it]
                drows.append(dict(N=N, threshold_km=thr, target_miss_km=r["safe_mult"] * thr, lead_days=lead,
                                  maneuvers_30d=int(n_m), maneuvers_per_sat_per_year=n_m * scale,
                                  mean_dv_per_maneuver_mm_s=1e6 * r["dv_sum"][it, il] / max(n_m, 1),
                                  annual_dv_per_sat_m_s=1e3 * r["dv_sum"][it, il] * scale))
        np.savez(f"{out}/partC_hist_N{N}.npz", hist_d=r["hist_d"], hist_v=r["hist_v"], hist_lat=r["hist_lat"],
                 per_sat=r["per_sat"])
        if r["records"].shape[0]:
            np.savetxt(f"{out}/partC_encounters_N{N}.csv", r["records"][:, :6], delimiter=",", fmt="%.6f",
                       header=f"t_s,i,j,miss_km,vrel_kms,lat_deg  (complete={r['records_complete']})",
                       comments="")
        print(f"[C] N={N}: {r['count']} encounters < {r['d']} km in {cfg['c_days']} d "
              f"(gas corrected {g['latitude_resolved_corrected']:.0f}) [{time.time() - t1:.0f} s]")
    keys = list(rows[0].keys())
    with open(f"{out}/partC_conjunctions.csv", "w") as f:
        f.write(",".join(keys) + "\n")
        for r_ in rows:
            f.write(",".join(f"{r_[k]:.6g}" if isinstance(r_[k], float) else str(r_[k]) for k in keys) + "\n")
    keys = list(drows[0].keys())
    with open(f"{out}/partD_delta_v.csv", "w") as f:
        f.write(",".join(keys) + "\n")
        for r_ in drows:
            f.write(",".join(f"{r_[k]:.6g}" if isinstance(r_[k], float) else str(r_[k]) for k in keys) + "\n")
    Ns = np.array([r_["N"] for r_ in rows], float)
    C = np.array([r_["count"] for r_ in rows], float)
    slope = float(np.polyfit(np.log(Ns), np.log(np.maximum(C, 1)), 1)[0]) if len(Ns) > 1 else float("nan")
    # figures
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(Ns, C, yerr=np.sqrt(C), fmt="ko", label="simulated (<1 km, 30 d)" if cfg["c_days"] == 30 else "simulated")
    NN = np.logspace(2.8, 5.2, 50)
    s0 = sh.make_shell(1000, seed=cfg["seed"])
    g1 = sh.gas_estimates(s0, cfg["c_days"])
    ax.plot(NN, g1["naive_uniform"] * (NN / 1000) ** 2, "b--", label="kinetic gas, uniform shell")
    ax.plot(NN, g1["latitude_resolved"] * (NN / 1000) ** 2, "r-", label="kinetic gas, latitude-resolved")
    ax.plot(Ns, [r_["gas_latitude_corrected"] for r_ in rows], "rs", mfc="none",
            label="latitude-resolved + thickness/same-plane corr.")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("satellites in shell N")
    ax.set_ylabel(f"close approaches < 1 km in {cfg['c_days']} d")
    ax.set_title(f"Part C: 550±5 km, i=53°, random phasing; slope {slope:.3f}"); ax.legend(fontsize=8)
    fig.savefig(f"{out}/partC_conjunctions_vs_N.png", dpi=130); plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    r = hists[max(hists)]
    lat = np.arange(-89.5, 90, 1.0)
    ax[0].bar(lat, r["hist_lat"] / r["hist_lat"].sum(), width=1, color="0.6", label=f"sim N={max(hists)}")
    ax[0].plot(lat, sh.gas_latitude_pdf(math.radians(53), lat), "r", label="gas model")
    ax[0].set_xlabel("latitude of TCA [deg]"); ax[0].set_ylim(0, 0.08); ax[0].legend()
    vv = np.arange(0.05, 16, 0.1)
    ax[1].bar(vv, r["hist_v"], width=0.1, color="0.6"); ax[1].set_xlabel("relative speed [km/s]"); ax[1].set_yscale("log")
    dd = (np.arange(sh.NB_D) + 0.5) / sh.NB_D
    ax[2].bar(dd, r["hist_d"], width=1 / sh.NB_D, color="0.6"); ax[2].set_xlabel("miss distance [km]")
    ax[2].set_title("miss-distance distribution (<1 km)")
    fig.tight_layout(); fig.savefig(f"{out}/partC_encounter_distributions.png", dpi=130); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 5))
    for thr in hists[max(hists)]["thr"]:
        sel = [d for d in drows if d["threshold_km"] == thr and d["lead_days"] == 1.0]
        ax.plot([d["N"] for d in sel], [d["annual_dv_per_sat_m_s"] for d in sel], "o-", label=f"threshold {thr*1e3:.0f} m")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("satellites in shell N")
    ax.set_ylabel("avoidance Δv per satellite per year [m/s]")
    ax.set_title("Part D: along-track burn 1 d before TCA, miss → 2×threshold"); ax.legend()
    fig.savefig(f"{out}/partD_annual_dv_vs_N.png", dpi=130); plt.close(fig)
    summ = dict(conjunctions=rows, count_vs_N_loglog_slope=slope, delta_v=drows, dt_s=dt)
    jdump(summ, f"{out}/partCD_summary.json")
    print(f"[C/D] done in {time.time() - t0:.1f} s")
    return summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--parts", default="ABCD")
    ap.add_argument("--dt", type=float, default=10.0, help="Part C screening step [s]")
    ap.add_argument("--c-days", type=float, default=None)
    ap.add_argument("--populations", default="1000,10000,100000")
    args = ap.parse_args()
    q = args.quick
    cfg = dict(a_days=3 if q else 30, sso_days=2 if q else 5, b_days=30 if q else 365,
               b_check_days=[0] if q else [0, 160], b_check_len=1 if q else 7,
               c_days=args.c_days or (1 if q else 30), seed=1,
               populations=[int(x) for x in args.populations.split(",")])
    out = "results_quick" if q else "results"
    os.makedirs(out, exist_ok=True)
    jdump(dict(cfg=cfg, dt=args.dt, parts=args.parts), f"{out}/run_config.json")
    i_sso = None
    if "A" in args.parts:
        ra = run_a(cfg, out)
        i_sso = ra["sso_mean_550"]["i_numeric_deg"]
    if "B" in args.parts:
        if i_sso is None:
            i_sso = part_a.sso_inclination(cfg["sso_days"], mean_a=True)["i_numeric_deg"]
        run_b(cfg, out, i_sso)
    if "C" in args.parts or "D" in args.parts:
        run_cd(cfg, out, args.dt)


if __name__ == "__main__":
    main()
