#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fit_rock_h2o_hhe.py
-------------------
smint-style retrieval of the interior composition of a planet with
  * a rocky interior (iron core + silicate mantle, 32% Fe / 68% silicates by mass)
  * a miscible envelope of H2O + H/He

Free composition parameters (mass fractions of the whole planet):
  f_rock + f_h2o + f_hhe = 1

Three cases are produced and overlaid on one corner plot:
  1) radius only           : Rp = 2.4 +/- 0.1 Re, broad uniform mass prior
  2) mass + radius         : Rp = 2.4 +/- 0.1 Re, Mp = 5.5 +/- 0.5 Me
  3) mass + radius + H2O/H : case 2 posterior importance-reweighted by an
                             envelope number ratio H2O/H = 0.25

NOTE ON THE GRID
  The public smint release (v1.0.0) has no grid for a rock + miscible H2O/H-He
  envelope. This script follows smint's structure (emcee + corner, Gaussian
  radius likelihood, mass/Teq priors, NaN = invalid model) but reads the
  structure grid from a CSV you supply (see load_grid_csv). The grid must have
  been computed for a 32% Fe / 68% silicate interior.

  --demo builds a TOY analytic grid (NOT PHYSICAL) only so you can check that
  the pipeline runs. Do not interpret demo results.

Usage
  python fit_rock_h2o_hhe.py --grid my_grid.csv --teq 600 --teq-err 50
  python fit_rock_h2o_hhe.py --demo                       # pipeline test
  python fit_rock_h2o_hhe.py --grid my_grid.csv --run-fit False   # replot only
"""

import argparse
import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import emcee
import corner
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MU_H2O = 18.015   # g/mol
MU_H = 1.008      # g/mol (atomic H)
# Protosolar H and He mass fractions (Asplund et al. 2009), renormalised to an
# H+He-only gas, since the metals are carried by the water component here.
X_PROTO, Y_PROTO = 0.7154, 0.2703
X_H_IN_GAS = X_PROTO / (X_PROTO + Y_PROTO)   # ~0.726


# ---------------------------------------------------------------------------
# Grid handling
# ---------------------------------------------------------------------------
# The grid is a regular 4D grid in
#   Mp   [Earth masses]   (interpolated linearly in log10 Mp, as in smint)
#   fenv [mass fraction of envelope = f_h2o + f_hhe]
#   Zenv [water mass fraction WITHIN the envelope = f_h2o / fenv]
#   Teq  [K]
# fenv x Zenv is a rectangle (unlike f_h2o x f_hhe, which is a triangle),
# which is what RegularGridInterpolator needs.
# Missing / non-converged models should be NaN (or absent from the CSV);
# they are treated as zero likelihood, like smint's validity flag.

DEFAULT_COLS = dict(Mp="Mp_earth", fenv="f_env", Zenv="Z_env",
                    Teq="Teq_K", R="Rp_earth")


def load_grid_csv(path, cols=DEFAULT_COLS):
    """Read a CSV with one row per model and build the regular grid arrays.
    Adapt `cols` (or this function) to your grid's format."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    keys = ["Mp", "fenv", "Zenv", "Teq"]
    raw = {k: np.asarray(data[cols[k]], dtype=float) for k in keys}
    axes = {k: np.unique(raw[k]) for k in keys}

    R = np.full([len(axes[k]) for k in keys], np.nan)
    idx = tuple(np.searchsorted(axes[k], raw[k]) for k in keys)
    R[idx] = np.asarray(data[cols["R"]], dtype=float)

    n_missing = np.isnan(R).sum()
    if n_missing:
        print(f"[grid] {n_missing}/{R.size} grid nodes missing -> treated as invalid")
    return axes, R


def make_demo_grid():
    """TOY grid for pipeline testing only. Not a physical structure model.
    Loosely mimics: rocky M-R scaling, water inflating the interior, and a
    Lopez & Fortney (2014)-like power law for the H/He layer thickness."""
    print("\n*** DEMO MODE: toy analytic grid, results are NOT physical ***\n")
    axes = dict(Mp=np.geomspace(0.5, 30, 40),
                fenv=np.linspace(0.0, 0.8, 41),
                Zenv=np.linspace(0.0, 1.0, 21),
                Teq=np.linspace(300, 1500, 7))
    M, F, Z, T = np.meshgrid(*[axes[k] for k in ["Mp", "fenv", "Zenv", "Teq"]],
                             indexing="ij")
    f_h2o, f_hhe = F * Z, F * (1 - Z)
    r_int = (1.0 + 0.55 * f_h2o) * M ** 0.27
    r_gas = 2.06 * M ** -0.21 * (np.clip(f_hhe, 1e-6, None) / 0.05) ** 0.59 \
        * (T / 278.) ** (4 * 0.044) * (1 + 0.3 * Z)
    R = r_int + np.where(f_hhe > 1e-4, r_gas, 0.0)
    return axes, R


def make_interpolator(axes, R):
    pts = (np.log10(axes["Mp"]), axes["fenv"], axes["Zenv"], axes["Teq"])
    return RegularGridInterpolator(pts, R, method="linear",
                                   bounds_error=False, fill_value=np.nan)


def grid_limits(axes):
    return {k: (axes[k].min(), axes[k].max()) for k in axes}


def radius_from_comp(interp, Mp, f_h2o, f_hhe, Teq):
    fenv = f_h2o + f_hhe
    Zenv = f_h2o / fenv if fenv > 0 else 0.0
    return float(interp((np.log10(Mp), fenv, Zenv, Teq)))


# ---------------------------------------------------------------------------
# Posterior
# theta = (f_h2o, x_hhe, Mp, Teq), x_hhe = f_hhe or log10(f_hhe)
# ---------------------------------------------------------------------------
def unpack(theta, p):
    f_h2o, x_hhe, Mp, Teq = theta
    f_hhe = 10 ** x_hhe if p["log_fhhe"] else x_hhe
    return f_h2o, f_hhe, Mp, Teq


def lnprior(theta, p):
    f_h2o, f_hhe, Mp, Teq = unpack(theta, p)
    lim = p["grid_lim"]
    fenv = f_h2o + f_hhe
    # composition: inside the simplex and inside the grid
    if f_h2o < 0 or f_hhe < 0 or fenv > lim["fenv"][1]:
        return -np.inf
    if p["log_fhhe"] and not (p["log_fhhe_min"] <= theta[1] <= np.log10(lim["fenv"][1])):
        return -np.inf
    # Teq: Gaussian, truncated to the grid
    if not (lim["Teq"][0] <= Teq <= lim["Teq"][1]):
        return -np.inf
    lp = -0.5 * ((Teq - p["Teq"]) / p["err_Teq"]) ** 2
    # mass
    if not (max(lim["Mp"][0], p["Mp_min"]) <= Mp <= min(lim["Mp"][1], p["Mp_max"])):
        return -np.inf
    if p["mass_constrained"]:
        lp += -0.5 * ((Mp - p["Mp"]) / p["err_Mp"]) ** 2
    # (otherwise uniform in Mp between Mp_min and Mp_max)
    return lp


def lnlike(theta, p, interp):
    f_h2o, f_hhe, Mp, Teq = unpack(theta, p)
    R = radius_from_comp(interp, Mp, f_h2o, f_hhe, Teq)
    if not np.isfinite(R):
        return -np.inf      # outside grid or invalid model
    return -0.5 * ((R - p["Rp"]) / p["err_Rp"]) ** 2


def lnprob(theta, p, interp):
    lp = lnprior(theta, p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, p, interp)


def initial_positions(p, interp, rng):
    """Draw walkers from the prior until they all have finite posterior."""
    pos = []
    lim = p["grid_lim"]
    while len(pos) < p["nwalkers"]:
        f = rng.dirichlet([1, 1, 1])                  # flat on the simplex
        f_h2o, f_hhe = f[1] * lim["fenv"][1], f[2] * lim["fenv"][1]
        x_hhe = np.log10(max(f_hhe, 10 ** p["log_fhhe_min"])) if p["log_fhhe"] else f_hhe
        if p["mass_constrained"]:
            Mp = rng.normal(p["Mp"], p["err_Mp"])
        else:
            Mp = rng.uniform(p["Mp_min"], p["Mp_max"])
        Teq = rng.normal(p["Teq"], p["err_Teq"])
        th = np.array([f_h2o, x_hhe, Mp, Teq])
        if np.isfinite(lnprob(th, p, interp)):
            pos.append(th)
    return np.array(pos)


def run_fit(p, interp, label, rng):
    print(f"\n=== Running MCMC: {label} ===")
    p0 = initial_positions(p, interp, rng)
    sampler = emcee.EnsembleSampler(p["nwalkers"], 4, lnprob, args=(p, interp))
    sampler.run_mcmc(p0, p["nsteps"], progress=True)
    try:
        print("  autocorrelation times:", np.round(sampler.get_autocorr_time(tol=0), 1))
    except Exception:
        pass
    print(f"  mean acceptance fraction: {np.mean(sampler.acceptance_fraction):.2f}")
    chain = sampler.get_chain(discard=p["nburnin"], thin=p["thin"], flat=True)
    f_h2o, f_hhe, Mp, Teq = (np.array(a) for a in zip(*[unpack(t, p) for t in chain]))
    f_rock = 1.0 - f_h2o - f_hhe
    return np.column_stack([f_rock, f_h2o, f_hhe, Mp, Teq])


# ---------------------------------------------------------------------------
# Case 3: post-processing with the envelope H2O/H number ratio
# ---------------------------------------------------------------------------
def h2o_to_h_ratio(f_h2o, f_hhe, include_water_H=False):
    """Number ratio n(H2O)/n(H) in the envelope.
    By default n(H) counts H atoms in the H/He gas only.
    With include_water_H, the 2 H atoms per water molecule are added too."""
    n_h2o = f_h2o / MU_H2O
    n_h = f_hhe * X_H_IN_GAS / MU_H
    if include_water_H:
        n_h = n_h + 2 * n_h2o
    with np.errstate(divide="ignore", invalid="ignore"):
        return n_h2o / n_h


def reweight_by_ratio(samples, p, rng):
    ratio = h2o_to_h_ratio(samples[:, 1], samples[:, 2], p["include_water_H"])
    good = np.isfinite(ratio) & (ratio > 0)
    lnw = np.full(len(samples), -np.inf)
    lnw[good] = -0.5 * ((np.log10(ratio[good]) - np.log10(p["ratio"]))
                        / p["ratio_sig_dex"]) ** 2
    w = np.exp(lnw - lnw.max())
    w /= w.sum()
    ess = 1.0 / np.sum(w ** 2)
    print(f"\n=== Case 3: importance reweighting by H2O/H = {p['ratio']} "
          f"(+/- {p['ratio_sig_dex']} dex) ===")
    print(f"  effective sample size: {ess:.0f} of {len(samples)}")
    if ess < 500:
        print("  WARNING: low ESS. Run case 2 longer, widen ratio_sig_dex, or fit "
              "case 3 directly (fix Zenv from the ratio and sample f_env).")
    idx = rng.choice(len(samples), size=len(samples), replace=True, p=w)
    return samples[idx], ratio


def zenv_from_ratio(r, include_water_H=False):
    """Water mass fraction of the envelope implied by a given H2O/H ratio
    (useful for checking / for a direct case-3 fit)."""
    if include_water_H:
        if r >= 0.5:
            return np.nan        # H2O/H can never exceed 0.5 if water H is counted
        r = r / (1 - 2 * r)      # convert to ratio relative to gas H only
    a = r * MU_H2O * X_H_IN_GAS / MU_H     # = f_h2o / f_hhe
    return a / (1 + a)


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
LABELS = [r"$f_{\rm rock}$", r"$f_{\rm H_2O}$", r"$f_{\rm H/He}$",
          r"$M_p\,[M_\oplus]$", r"$T_{\rm eq}$ [K]"]
NAMES = ["f_rock", "f_h2o", "f_hhe", "Mp", "Teq"]


def summarize(results, path):
    lines = ["case,param,p16,p50,p84"]
    for case, s in results.items():
        print(f"\n{case}:")
        for i, n in enumerate(NAMES):
            lo, med, hi = np.percentile(s[:, i], [15.9, 50, 84.1])
            print(f"  {n:7s} = {med:.4f} +{hi - med:.4f} -{med - lo:.4f}")
            lines.append(f"{case},{n},{lo:.6g},{med:.6g},{hi:.6g}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def plot_corner(results, colors, path, include_mass=False):
    cols = [0, 1, 2, 3] if include_mass else [0, 1, 2]
    stacked = np.vstack([s[:, cols] for s in results.values()])
    ranges = [(max(0, np.percentile(stacked[:, j], 0.1) - 0.02),
               min(1 if j < 3 else np.inf, np.percentile(stacked[:, j], 99.9) + 0.02))
              for j in range(len(cols))]
    fig = None
    for (case, s), c in zip(results.items(), colors):
        fig = corner.corner(s[:, cols], fig=fig, labels=[LABELS[j] for j in cols],
                            range=ranges, color=c, bins=40, smooth=1.0,
                            levels=(0.393, 0.865), plot_datapoints=False,
                            plot_density=False, fill_contours=False,
                            hist_kwargs=dict(density=True, color=c, lw=1.5),
                            contour_kwargs=dict(colors=[c], linewidths=1.5))
    handles = [Line2D([], [], color=c, lw=2, label=k)
               for k, c in zip(results.keys(), colors)]
    fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=11)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    print(f"\nCorner plot saved to {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def str2bool(v):
    return str(v).lower() in ("1", "true", "yes", "y")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", help="CSV structure grid (rock 32%%Fe + H2O/H-He envelope)")
    ap.add_argument("--demo", action="store_true", help="use TOY grid (pipeline test)")
    ap.add_argument("--outdir", default="results_rock_h2o_hhe")
    ap.add_argument("--run-fit", type=str2bool, default=True,
                    help="False = reload saved chains and only post-process/plot")
    # planet
    ap.add_argument("--rp", type=float, default=2.4)
    ap.add_argument("--rp-err", type=float, default=0.1)
    ap.add_argument("--mp", type=float, default=5.5)
    ap.add_argument("--mp-err", type=float, default=0.5)
    ap.add_argument("--mp-min", type=float, default=0.5, help="uniform mass prior, case 1")
    ap.add_argument("--mp-max", type=float, default=20.0, help="uniform mass prior, case 1")
    ap.add_argument("--teq", type=float, default=500.0, help="SET THIS for your planet")
    ap.add_argument("--teq-err", type=float, default=50.0)
    # case 3
    ap.add_argument("--ratio", type=float, default=0.25, help="envelope H2O/H number ratio")
    ap.add_argument("--ratio-sig-dex", type=float, default=0.05,
                    help="width of the ratio constraint in dex (a delta function "
                         "would give zero weight to every sample)")
    ap.add_argument("--include-water-H", action="store_true",
                    help="count H atoms in H2O in the denominator of H2O/H")
    # priors / sampler
    ap.add_argument("--log-fhhe", action="store_true",
                    help="log-uniform prior on f_HHe instead of flat-on-simplex")
    ap.add_argument("--log-fhhe-min", type=float, default=-4.0)
    ap.add_argument("--nwalkers", type=int, default=64)
    ap.add_argument("--nsteps", type=int, default=6000)
    ap.add_argument("--nburnin", type=int, default=2000)
    ap.add_argument("--thin", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--corner-with-mass", action="store_true")
    a = ap.parse_args()

    if not a.demo and not a.grid:
        ap.error("provide --grid path/to/grid.csv, or --demo to test the pipeline")

    rng = np.random.default_rng(a.seed)
    os.makedirs(a.outdir, exist_ok=True)

    axes, R = make_demo_grid() if a.demo else load_grid_csv(a.grid)
    interp = make_interpolator(axes, R)
    lim = grid_limits(axes)
    print("Grid limits:", {k: tuple(np.round(v, 3)) for k, v in lim.items()})

    base = dict(Rp=a.rp, err_Rp=a.rp_err, Mp=a.mp, err_Mp=a.mp_err,
                Mp_min=a.mp_min, Mp_max=a.mp_max, Teq=a.teq, err_Teq=a.teq_err,
                grid_lim=lim, log_fhhe=a.log_fhhe, log_fhhe_min=a.log_fhhe_min,
                nwalkers=a.nwalkers, nsteps=a.nsteps, nburnin=a.nburnin, thin=a.thin,
                ratio=a.ratio, ratio_sig_dex=a.ratio_sig_dex,
                include_water_H=a.include_water_H)

    c1 = f"1) R = {a.rp}±{a.rp_err} R_E"
    c2 = f"2) + M = {a.mp}±{a.mp_err} M_E"
    c3 = f"3) + H2O/H = {a.ratio}"
    f1, f2 = os.path.join(a.outdir, "case1_samples.npy"), os.path.join(a.outdir, "case2_samples.npy")

    if a.run_fit:
        s1 = run_fit({**base, "mass_constrained": False}, interp, c1, rng)
        s2 = run_fit({**base, "mass_constrained": True}, interp, c2, rng)
        np.save(f1, s1)
        np.save(f2, s2)
    else:
        s1, s2 = np.load(f1), np.load(f2)

    s3, ratio2 = reweight_by_ratio(s2, base, rng)
    np.save(os.path.join(a.outdir, "case3_samples.npy"), s3)
    print(f"  H2O/H over case-2 posterior: median {np.nanmedian(ratio2):.3g}, "
          f"68% [{np.nanpercentile(ratio2, 15.9):.3g}, {np.nanpercentile(ratio2, 84.1):.3g}]")
    print(f"  -> implied envelope water mass fraction Z_env = "
          f"{zenv_from_ratio(a.ratio, a.include_water_H):.3f}")

    results = {c1: s1, c2: s2, c3: s3}
    summarize(results, os.path.join(a.outdir, "constraints.csv"))
    plot_corner(results, ["#999999", "#1f77b4", "#d62728"],
                os.path.join(a.outdir, "corner_fractions.png"), a.corner_with_mass)


if __name__ == "__main__":
    main()
