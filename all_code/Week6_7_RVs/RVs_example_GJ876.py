# %%
#######################
# GJ 876 RADIAL-VELOCITY ANALYSIS (EXOPLANET ANALYSIS CLASS)
#######################

"""
End-to-end RV analysis of GJ 876 using public data from DACE.

Workflow
    Data access with dace_query (multi-instrument)
    Data checks and cleaning (units, QC, error bars, outliers, jitter)
    Periodograms: Lomb-Scargle, BGLS (l1 periodogram: see l1_periodogram_gj876.py)
    Activity indicators: rotation and RV correlations
    Keplerian fits with radvel and exoplanet (Kipping 2013 eccentricity prior)
    Periodogram of the residuals
    Bayesian model comparison (lnL, lnZ, BIC, WAIC)
    Full trans-dimensional analysis with kima

Install (suggested)
    uv add dace-query astropy radvel "exoplanet>=0.6" pymc arviz dynesty corner kima
    uv add git+https://github.com/jvines/BayesianGLS

Physics note for students
    GJ 876 b and c are in a 2:1 mean-motion resonance and interact strongly.
    Over a ~20 yr baseline their orbits precess, so a sum of Keplerians is only an
    approximation: expect residuals larger than the error bars and spurious
    'extra' signals. A full treatment needs an N-body (dynamical) fit. Setting
    TIME_WINDOW below to a shorter span is a good exercise to see this.
"""

# %%

#######################
# IMPORTS AND GLOBAL SETTINGS
#######################

import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats, optimize
from astropy.timeseries import LombScargle

from RVs_functions import *

warnings.filterwarnings("ignore", category=RuntimeWarning)
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

# %%

STAR = "GJ876"
MSTAR, MSTAR_ERR = 0.37, 0.01            # Msun (von Braun et al. 2014)
KNOWN_PERIODS = {"d": 1.9378, "c": 30.08, "b": 61.08, "e": 123.8}   # days
PROT_RANGE = (80.0, 100.0)               # approximate literature rotation period range
PMIN, PMAX, OFAC = 1.2, 3000.0, 10       # periodogram period range and oversampling
MIN_POINTS = 10                          # drop instruments with fewer points
KIPPING_A, KIPPING_B = 0.867, 3.03       # Kipping (2013) Beta prior on eccentricity
TIME_WINDOW = None                       # e.g. (55000, 58000) in RJD to limit the baseline

# Switches for the slow parts
RUN_RADVEL_MCMC = True
RUN_JAXOPLANET = True
RUN_DYNESTY = True
RUN_XO_COMPARISON = True
RUN_KIMA = True

OUTDIR = Path("gj876_outputs")
OUTDIR.mkdir(exist_ok=True)

# %%
#######################
# DATA ACCESS
#######################

# %%
#######################
# DOWNLOAD THE GJ 876 TIME SERIES FROM DACE
#######################
from dace_query.spectroscopy import Spectroscopy

raw = Spectroscopy.get_timeseries(STAR, sorted_by_instrument=False, output_format="pandas")
raw = pd.DataFrame(raw)
if 'instrument_name' in raw.columns:
    raw=raw.rename(columns={'instrument_name':'ins_name','instrument_id':'ins_id','instrument_group':'ins_group'})
raw.to_csv(OUTDIR / "gj876_dace_raw.csv", index=False)
print("Raw table shape:", raw.shape)
print("Columns:", list(raw.columns))

for c in ["ins_name", "ins_mode", "drs_version"]:
    if c not in raw.columns:
        raw[c] = "unknown"
raw[["ins_name", "ins_mode", "drs_version"]] = raw[["ins_name", "ins_mode", "drs_version"]].fillna("unknown").astype(str)
print(raw.groupby(["ins_name", "drs_version", "ins_mode"]).size())
# %%
#######################
# KEEP ONE DRS REDUCTION PER INSTRUMENT AND BUILD A TIDY TABLE
#######################
# DACE often stores the same spectra reduced with several pipeline versions.
# Keeping all of them would double-count observations.


def pick_drs(group):
    counts = group.groupby("drs_version").size().sort_index()
    best = counts.index[counts.values == counts.max()][-1]
    return group[group["drs_version"] == best]


drs_counts = (raw.groupby(["ins_name", "drs_version"]).size().reset_index(name="n")
              .sort_values(["ins_name", "n", "drs_version"]))
best_drs = drs_counts.groupby("ins_name").tail(1)[["ins_name", "drs_version"]]
sel = raw.merge(best_drs, on=["ins_name", "drs_version"], how="inner").reset_index(drop=True)

col_t = find_col(sel, ["rjd", "bjd", "time"])
col_rv = find_col(sel, ["rv", "vrad"])
col_err = find_col(sel, ["rv_err", "svrad", "rv_error"])
print(f"Using columns: t={col_t}, rv={col_rv}, err={col_err}")

df = pd.DataFrame({
    "t": pd.to_numeric(sel[col_t], errors="coerce"),
    "rv": pd.to_numeric(sel[col_rv], errors="coerce"),
    "err": pd.to_numeric(sel[col_err], errors="coerce"),
})
# Separate instrument modes only when an instrument has more than one
modes_per_inst = sel.groupby("ins_name")["ins_mode"].nunique()
df["inst"] = [n if modes_per_inst[n] == 1 else f"{n}_{m}"
              for n, m in zip(sel["ins_name"], sel["ins_mode"])]
df["inst"] = df["inst"].str.replace(r"[^A-Za-z0-9]+", "", regex=True)  # safe for radvel/kima names

# Quality-control flag (missing flag means 'not assessed', so keep the point)
col_qc = find_col(sel, ["drs_qc"])
if col_qc:
    q = sel[col_qc]
    df["qc"] = q.isna() | q.astype(str).str.lower().isin(["true", "1", "1.0"])
else:
    df["qc"] = True

# Activity indicators (names vary between DACE versions, so search for them)
INDICATORS = {"fwhm": ["fwhm","ccf_fwhm"], "bis": ["bispan", "bis_span", "bis", "ccf_bispan"], "contrast": ["contrast","ccf_contrast"],
              "rhk": ["rhk","spectro_rhk"], "smw": ["s_mw", "smw","spectra_smw"], "halpha": ["ha", "halpha","spectra_halpha"],
              "naD": ["na","spectro_na"], "caII": ["ca","spectro_ca"]}
available_ind = []
for key, cands in INDICATORS.items():
    c = find_col(sel, cands)
    if c is None:
        continue
    v = pd.to_numeric(sel[c], errors="coerce")
    v[v <= -9999] = np.nan                      # sentinel values
    df[key] = v
    ce = find_col(sel, [c + "_err"])
    if ce:
        ve = pd.to_numeric(sel[ce], errors="coerce")
        ve[(ve <= 0) | (ve <= -9999)] = np.nan
        df[key + "_err"] = ve
    else:
        df[key + "_err"] = np.nan
    if df[key].notna().sum() > 0:
        available_ind.append(key)
print("Activity indicators found:", available_ind)

if TIME_WINDOW is not None:
    df = df[(df["t"] >= TIME_WINDOW[0]) & (df["t"] <= TIME_WINDOW[1])]
df = df.sort_values("t").reset_index(drop=True)
print(df.groupby("inst").size())

# %%
#######################
# DATA CHECKS AND CLEANING
#######################

# %%
#######################
# BASIC VALIDITY: FINITE VALUES, POSITIVE ERRORS, QC FLAG
#######################
clean_log = []
n0 = len(df)
bad = ~np.isfinite(df[["t", "rv", "err"]]).all(axis=1) | (df["err"] <= 0) | (~df["qc"])
for inst, n in df[bad].groupby("inst").size().items():
    clean_log.append((inst, "invalid / failed QC", n))
df = df[~bad].reset_index(drop=True)
print(f"Removed {n0 - len(df)} invalid points")

# %%
#######################
# UNIT SANITY CHECK PER INSTRUMENT
#######################
# DACE RVs should be in m/s. A km/s-labelled dataset shows up as tiny errors (~0.001)
# and a tiny scatter compared with GJ 876's ~200 m/s planetary signal.
summary = df.groupby("inst").agg(N=("rv", "size"), med_rv=("rv", "median"),
                                 rms_rv=("rv", "std"), med_err=("err", "median"),
                                 min_err=("err", "min"))
print(summary)

for inst, row in summary.iterrows():
    if row["med_err"] < 0.05 and row["rms_rv"] < 5:
        print(f"  {inst}: looks like km/s -> converting to m/s")
        m = df["inst"] == inst
        df.loc[m, ["rv", "err"]] *= 1000.0
        clean_log.append((inst, "converted km/s -> m/s", int(m.sum())))

# %%
#######################
# ANOMALOUS ERROR BARS (TOO LARGE OR TOO SMALL)
#######################
for inst in df["inst"].unique():
    m = df["inst"] == inst
    e = df.loc[m, "err"]
    med, sd = np.median(e), mad_std(e)
    too_big = m & (df["err"] > med + 5 * max(sd, 0.05 * med))
    if too_big.sum():
        clean_log.append((inst, "error bar > median + 5 MAD", int(too_big.sum())))
        df = df[~too_big]
    m = df["inst"] == inst
    too_small = m & (df["err"] < 0.25 * med)
    if too_small.sum():
        # Unrealistically small errors would dominate the fit: raise them to the median
        clean_log.append((inst, "error bar < 0.25 x median (inflated)", int(too_small.sum())))
        df.loc[too_small, "err"] = med
df = df.reset_index(drop=True)

# %%
#######################
# DROP INSTRUMENTS WITH TOO FEW POINTS
#######################
counts = df["inst"].value_counts()
for inst in counts[counts < MIN_POINTS].index:
    clean_log.append((inst, f"fewer than {MIN_POINTS} points (instrument dropped)", int(counts[inst])))
df = df[df["inst"].isin(counts[counts >= MIN_POINTS].index)].reset_index(drop=True)
print("Instruments kept:", sorted(df["inst"].unique()))

# %%
#######################
# RV OUTLIERS RELATIVE TO A QUICK MULTI-SINUSOID MODEL
#######################
# Sigma-clipping the raw RVs would remove real planetary maxima/minima (K_b ~ 200 m/s).
# Instead, clip on residuals after removing the dominant periodic signals.
for it in range(2):
    pw = prewhiten(df, n_signals=2, verbose=False)
    res = pw["fit"]["resid"]
    outlier = np.zeros(len(df), bool)
    for inst in df["inst"].unique():
        m = (df["inst"] == inst).to_numpy()
        r = res[m]
        outlier[m] = np.abs(r - np.median(r)) > 5 * mad_std(r)
    for inst, n in df[outlier].groupby("inst").size().items():
        clean_log.append((inst, f"5-sigma residual outlier (pass {it + 1})", n))
    df = df[~outlier].reset_index(drop=True)
    if outlier.sum() == 0:
        break

# %%
#######################
# SCATTER VERSUS ERROR BARS (IS EXTRA JITTER NEEDED?)
#######################
pw = prewhiten(df, n_signals=2, verbose=False)
df["resid_quick"] = pw["fit"]["resid"]
rows = []
for inst, g in df.groupby("inst"):
    rms = np.std(g["resid_quick"])
    rows.append(dict(inst=inst, N=len(g), med_err=np.median(g["err"]), rms_resid=rms,
                     ratio=rms / np.median(g["err"]),
                     jitter_est=np.sqrt(max(rms**2 - np.mean(g["err"]**2), 0))))
scatter_tab = pd.DataFrame(rows).set_index("inst")
print(scatter_tab.round(2))
print("\nratio >> 1: extra noise (activity, instrument systematics, or unmodelled planets).")
print("ratio << 1: error bars probably overestimated or mislabelled.")

print("\nCleaning log:")
print(pd.DataFrame(clean_log, columns=["inst", "action", "N"]).to_string(index=False))
df.to_csv(OUTDIR / "gj876_rv_clean.csv", index=False)


# %%
#######################
# OVERVIEW PLOT OF THE CLEANED DATA
#######################
offsets0 = linear_fit(df, pw["periods"])["offsets"]
fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
for inst, g in df.groupby("inst"):
    ax[0].errorbar(g["t"], g["rv"] - offsets0[inst], g["err"], fmt="o", ms=3, label=inst)
    ax[1].errorbar(g["t"], g["resid_quick"], g["err"], fmt="o", ms=3)
ax[0].set_ylabel("RV - offset [m/s]")
ax[1].set_ylabel("Resid. (2 sines) [m/s]")
ax[1].set_xlabel("RJD [d]")
ax[0].legend(ncol=4, fontsize=8)
fig.tight_layout()
plt.show()

# %%
#######################
# PERIODOGRAMS
#######################

# %%
#######################
# COMMON INPUTS AND THE SPECTRAL WINDOW
#######################
t = df["t"].to_numpy()
err = df["err"].to_numpy()
freq = freq_grid(t)
rv0 = offset_subtracted(df, linear_fit(df)["offsets"])   # weighted-mean offsets only

window = LombScargle(t, np.ones_like(t), fit_mean=False, center_data=False).power(freq)
fig, ax = plt.subplots(figsize=(11, 3))
plot_periodogram(ax, freq, window, marks=None, prot=False, ylabel="Window power")
ax.set_title("Spectral window (look for 1 d, ~29.5 d lunar, 365 d aliases)")
plt.show()

# %%
#######################
# GENERALISED LOMB-SCARGLE (ASTROPY) ON THE COMBINED DATA
#######################
ls, power, fal = ls_periodogram(t, rv0, err, freq)
fig, ax = plt.subplots(figsize=(11, 4))
plot_periodogram(ax, freq, power, fal=fal)
ax.set_title("GLS periodogram, all instruments (offsets removed)")
ax.legend(fontsize=7, loc="upper left")
plt.show()

top = np.argsort(power)[::-1][:5]
print("Top peaks (not de-duplicated):", np.round(1 / freq[top], 3))

# %%
#######################
# GLS PER INSTRUMENT (CONSISTENCY CHECK)
#######################
insts = sorted(df["inst"].unique())
fig, axes = plt.subplots(len(insts), 1, figsize=(11, 2.2 * len(insts)), sharex=True, squeeze=False)
for ax, inst in zip(axes[:, 0], insts):
    g = df[df["inst"] == inst]
    if len(g) < 15:
        ax.text(0.5, 0.5, f"{inst}: too few points", transform=ax.transAxes, ha="center")
        continue
    f_i = freq_grid(g["t"].to_numpy())
    _, p_i, fal_i = ls_periodogram(g["t"], g["rv"], g["err"], f_i)
    plot_periodogram(ax, f_i, p_i, fal=fal_i)
    ax.set_title(inst, fontsize=9)
fig.tight_layout()
plt.show()

# %%
#######################
# ITERATIVE PREWHITENING WITH GLS
#######################
print("Prewhitening the combined data set:")
pw = prewhiten(df, n_signals=4)
fig, axes = plt.subplots(len(pw["powers"]), 1, figsize=(11, 2.4 * len(pw["powers"])), sharex=True)
for i, (ax, p) in enumerate(zip(axes, pw["powers"])):
    plot_periodogram(ax, pw["freq"], p)
    ax.set_title(f"After removing {i} signal(s); peak at {pw['periods'][i]:.3f} d "
                 f"(FAP {pw['faps'][i]:.1e})", fontsize=9)
fig.tight_layout()
plt.show()

# %%
#######################
# BAYESIAN GLS (MORTIER ET AL. 2015)
#######################
# Uses the jvines/BayesianGLS package if importable, otherwise an equivalent local implementation.


def bgls_local(t, y, err, freq, chunk=1500):
    """log of the (unnormalised) BGLS probability, Mortier et al. (2015) eq. in sec. 2."""
    t, y, err = map(np.asarray, (t, y, err))
    w = 1.0 / err**2
    W, bigY = w.sum(), (w * y).sum()
    logp = np.empty(freq.size)
    for i0 in range(0, freq.size, chunk):
        om = 2 * np.pi * freq[i0:i0 + chunk][:, None]
        theta = 0.5 * np.arctan2((w * np.sin(2 * om * t)).sum(1), (w * np.cos(2 * om * t)).sum(1))
        x = om * t - theta[:, None]
        cx, sx = np.cos(x), np.sin(x)
        C, S = (w * cx).sum(1), (w * sx).sum(1)
        YC, YS = (w * y * cx).sum(1), (w * y * sx).sum(1)
        CC, SS = (w * cx**2).sum(1), (w * sx**2).sum(1)
        K = (C**2 * SS + S**2 * CC - W * CC * SS) / (2 * CC * SS)
        L = (bigY * CC * SS - C * YC * SS - S * YS * CC) / (CC * SS)
        M = (YC**2 * SS + YS**2 * CC) / (2 * CC * SS)
        logp[i0:i0 + chunk] = -0.5 * np.log(CC * SS * np.abs(K)) + (M - L**2 / (4 * K))
    return logp - np.nanmax(logp)


def run_bgls(t, y, err, freq):
    try:
        import bgls as bgls_pkg
        periods, prob = bgls_pkg.bgls(np.asarray(t), np.asarray(y), np.asarray(err),
                                      plow=1 / freq.max(), phigh=1 / freq.min(), ofac=OFAC)
        f_out = 1.0 / np.asarray(periods)
        prob = np.asarray(prob, float)
        logp = np.log(np.clip(prob, 1e-300, None)) if prob.min() >= 0 else prob
        print("BGLS: using the BayesianGLS package")
        order = np.argsort(f_out)
        return f_out[order], (logp - logp.max())[order]
    except Exception as exc:  # package missing or different API
        print(f"BGLS: package not used ({type(exc).__name__}); local implementation")
        return freq, bgls_local(t, y, err, freq)


f_b, logp_b = run_bgls(t, rv0, err, freq)
fig, ax = plt.subplots(figsize=(11, 4))
plot_periodogram(ax, f_b, logp_b, ylabel="log P(f|D) - max")
ax.set_ylim(max(logp_b.min(), -200), 5)
ax.set_title("BGLS, combined data (offsets removed)")
plt.show()

# %%
#######################
# L1 PERIODOGRAM NOTE
#######################
# The l1 periodogram (Hara et al. 2017) is not pip-installable, so it lives in a
# separate script: l1_periodogram_gj876.py. It reads the cleaned data written
# above (gj876_outputs/gj876_rv_clean.csv).
print("l1 periodogram: run l1_periodogram_gj876.py separately (see its header for setup).")

# %%
#######################
# ACTIVITY INDICATORS
#######################

# %%
#######################
# TIME SERIES OF THE INDICATORS
#######################
ind_ok = [k for k in available_ind if df[k].notna().sum() >= 20]
print("Indicators with enough data:", ind_ok)
if ind_ok:
    fig, axes = plt.subplots(len(ind_ok), 1, figsize=(11, 2.0 * len(ind_ok)), sharex=True, squeeze=False)
    for ax, key in zip(axes[:, 0], ind_ok):
        for inst, g in df.groupby("inst"):
            g = g[g[key].notna()]
            if len(g):
                ax.plot(g["t"], g[key] - np.median(g[key]), "o", ms=2.5, label=inst)
        ax.set_ylabel(key)
    axes[0, 0].legend(ncol=4, fontsize=7)
    axes[-1, 0].set_xlabel("RJD [d]")
    fig.tight_layout()
    plt.show()

# %%
#######################
# PERIODOGRAMS OF THE INDICATORS (LOOK FOR THE ~80-100 D ROTATION)
#######################
# Shaded band: literature rotation period range. A planet candidate whose period
# (or its harmonics Prot/2, Prot/3) appears in an indicator is suspicious.
if ind_ok:
    fig, axes = plt.subplots(len(ind_ok), 1, figsize=(11, 2.3 * len(ind_ok)), sharex=True, squeeze=False)
    for ax, key in zip(axes[:, 0], ind_ok):
        g = df[df[key].notna()].copy()
        g[key] = g[key] - g.groupby("inst")[key].transform("median")
        e_k = g[key + "_err"].to_numpy()
        e_k = e_k if np.all(np.isfinite(e_k)) else None
        f_k = freq_grid(g["t"].to_numpy(), pmin=PMIN, pmax=min(PMAX, np.ptp(g["t"])))
        _, p_k, fal_k = ls_periodogram(g["t"], g[key], e_k, f_k)
        plot_periodogram(ax, f_k, p_k, fal=fal_k)
        P_best = 1 / f_k[np.argmax(p_k)]
        ax.set_title(f"{key}: highest peak {P_best:.1f} d", fontsize=9)
    fig.tight_layout()
    plt.show()

# %%
#######################
# CORRELATIONS BETWEEN RVS AND INDICATORS
#######################
# Planetary signals dwarf activity in GJ 876, so correlate the indicators with the
# RV residuals after prewhitening, per instrument (offsets differ between instruments).
df["resid_pw"] = pw["fit"]["resid"]
rows = []
for key in ind_ok:
    for inst, g in df.groupby("inst"):
        g = g[g[key].notna()]
        if len(g) < 15:
            continue
        for which in ["rv", "resid_pw"]:
            r, p = stats.pearsonr(g[key], g[which])
            rho, prho = stats.spearmanr(g[key], g[which])
            rows.append(dict(indicator=key, inst=inst, y=which, N=len(g),
                             pearson_r=r, p_value=p, spearman_rho=rho, p_spearman=prho))
corr_tab = pd.DataFrame(rows)
if len(corr_tab):
    print(corr_tab.round(3).to_string(index=False))

    main_inst = df["inst"].value_counts().index[0]
    g = df[df["inst"] == main_inst]
    keys = [k for k in ind_ok if g[k].notna().sum() >= 15]
    if keys:
        fig, axes = plt.subplots(1, len(keys), figsize=(3.2 * len(keys), 3), squeeze=False)
        for ax, key in zip(axes[0], keys):
            gg = g[g[key].notna()]
            ax.plot(gg[key], gg["resid_pw"], "o", ms=3)
            ax.set_xlabel(key)
        axes[0, 0].set_ylabel("RV resid. [m/s]")
        fig.suptitle(f"{main_inst}: residual RVs vs indicators")
        fig.tight_layout()
        plt.show()

# %%
#######################
# KEPLERIAN FITS
#######################

# %%
#######################
# INITIAL GUESSES FROM THE PREWHITENING
#######################
N_PL_FIT = 2   # GJ 876 b and c dominate; extend after the residual analysis
guess_fit = linear_fit(df, pw["periods"][:N_PL_FIT])
guesses = guess_fit["sines"]
offsets_guess = guess_fit["offsets"]
for i, gs in enumerate(guesses, 1):
    print(f"planet {i}: P = {gs['P']:.4f} d, K = {gs['K']:.1f} m/s, tc = {gs['tc']:.2f}")

# %%
#######################
# RADVEL MODEL BUILDER (KIPPING ECCENTRICITY PRIOR)
#######################
import radvel


def kipping_logprior(x):
    e = np.clip(x[0] ** 2 + x[1] ** 2, 1e-8, 1 - 1e-8)
    return stats.beta.logpdf(e, KIPPING_A, KIPPING_B)


def build_radvel_posterior(df, guesses, offsets, kipping=True):
    """
    Fitting basis 'per tc secosw sesinw k'. A uniform prior on (sqrt(e)cos w, sqrt(e)sin w)
    inside the unit disk is uniform in e and w, so adding log Beta(e) gives exactly the
    Kipping (2013) prior on e.
    """
    n = len(guesses)
    params = radvel.Parameters(n, basis="per tc secosw sesinw k")
    for i, gs in enumerate(guesses, 1):
        params[f"per{i}"] = radvel.Parameter(value=gs["P"])
        params[f"tc{i}"] = radvel.Parameter(value=gs["tc"])
        params[f"secosw{i}"] = radvel.Parameter(value=0.1)
        params[f"sesinw{i}"] = radvel.Parameter(value=0.1)
        params[f"k{i}"] = radvel.Parameter(value=max(gs["K"], 1.0))
    params["dvdt"] = radvel.Parameter(value=0.0, vary=False)
    params["curv"] = radvel.Parameter(value=0.0, vary=False)
    model = radvel.RVModel(params, time_base=np.median(df["t"]))

    likes = []
    for inst, g in df.groupby("inst"):
        like = radvel.likelihood.RVLikelihood(model, g["t"].to_numpy(), g["rv"].to_numpy(),
                                              g["err"].to_numpy(), suffix="_" + inst)
        like.params["gamma_" + inst] = radvel.Parameter(value=offsets[inst])
        like.params["jit_" + inst] = radvel.Parameter(value=2.0)
        likes.append(like)
    post = radvel.posterior.Posterior(radvel.likelihood.CompositeLikelihood(likes))

    post.priors += [radvel.prior.EccentricityPrior(n), radvel.prior.PositiveKPrior(n)]
    for i, gs in enumerate(guesses, 1):
        post.priors += [radvel.prior.HardBounds(f"per{i}", 0.9 * gs["P"], 1.1 * gs["P"])]
        if kipping:
            post.priors += [radvel.prior.UserDefinedPrior([f"secosw{i}", f"sesinw{i}"],
                                                          kipping_logprior, f"Kipping e{i}")]
    for inst in df["inst"].unique():
        post.priors += [radvel.prior.HardBounds("jit_" + inst, 0.0, 50.0)]
    return post


def radvel_map(post, n_restarts=2):
    for _ in range(n_restarts):
        for method in ["Powell", "Nelder-Mead"]:
            res = optimize.minimize(post.neglogprob_array, post.get_vary_params(), method=method,
                                    options=dict(maxiter=200000, maxfev=200000, xatol=1e-8, fatol=1e-8)
                                    if method == "Nelder-Mead" else dict(maxiter=200000))
            post.set_vary_params(res.x)
    return post


# %%
#######################
# RADVEL MAXIMUM A POSTERIORI FIT
#######################
post = build_radvel_posterior(df, guesses, offsets_guess)
post = radvel_map(post)
print(post)
print("ln L (max) =", post.likelihood.logprob())

try:
    from radvel.plot import orbit_plots
    mp = orbit_plots.MultipanelPlot(post, saveplot=str(OUTDIR / "radvel_map_fit.pdf"))
    mp.plot_multipanel()
    plt.show()
except Exception as exc:
    print("radvel multipanel plot failed:", exc)

# %%
#######################
# RADVEL MCMC SAMPLING OF THE POSTERIOR
#######################
if RUN_RADVEL_MCMC:
    chains = radvel.mcmc(post, nwalkers=40, nrun=4000, ensembles=4, savename=None)
    for i in range(1, len(guesses) + 1):
        chains[f"e{i}"] = chains[f"secosw{i}"] ** 2 + chains[f"sesinw{i}"] ** 2
        chains[f"w{i}"] = np.degrees(np.arctan2(chains[f"sesinw{i}"], chains[f"secosw{i}"]))
        mstar = np.random.normal(MSTAR, MSTAR_ERR, len(chains))
        chains[f"msini{i}_MJ"] = radvel.utils.Msini(chains[f"k{i}"], chains[f"per{i}"], mstar,
                                                     chains[f"e{i}"], Msini_units="jupiter")
    cols = [c for c in chains.columns if re.match(r"(per|k|e|w|msini)\d", c)]
    q = chains[cols].quantile([0.16, 0.5, 0.84]).T
    q["median"], q["-"], q["+"] = q[0.5], q[0.5] - q[0.16], q[0.84] - q[0.5]
    print(q[["median", "-", "+"]].round(4))
    try:
        import corner
        corner.corner(chains[[c for c in cols if not c.startswith("msini")]], show_titles=True,
                      title_fmt=".3f", quantiles=[0.16, 0.5, 0.84])
        plt.show()
    except ImportError:
        print("Install 'corner' for corner plots")

# %%
#######################
# EXOPLANET (PYMC) MODEL BUILDER
#######################
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import exoplanet as xo


def build_xo_model(df, guesses, offsets):
    """
    Same physical model as radvel. Note: exoplanet's omega refers to the planet's
    orbit and differs from radvel's (stellar) omega by pi.
    """
    n = len(guesses)
    idx, names = inst_codes(df)
    tref = np.median(df["t"])
    tt = df["t"].to_numpy() - tref
    y, e = df["rv"].to_numpy(), df["err"].to_numpy()
    P0 = np.array([g["P"] for g in guesses])
    tc0 = np.array([g["tc"] for g in guesses]) - tref
    K0 = np.array([max(g["K"], 2.0) for g in guesses])

    with pm.Model() as model:
        logP = pm.Normal("logP", mu=np.log(P0), sigma=0.02, shape=n)
        P = pm.Deterministic("P", pt.exp(logP))
        tc = pm.Normal("tc", mu=tc0, sigma=P0 / 4, shape=n)
        logK = pm.Uniform("logK", lower=-2, upper=np.log(1000), initval=np.log(K0), shape=n)
        K = pm.Deterministic("K", pt.exp(logK))

        ecs = pm.Uniform("ecs", lower=-1, upper=1, shape=(2, n), initval=0.1 * np.ones((2, n)))
        ecc = pm.Deterministic("ecc", pt.sum(ecs**2, axis=0))
        omega = pm.Deterministic("omega", pt.arctan2(ecs[1], ecs[0]))
        pm.Potential("unit_disk", pt.switch(pt.all(pt.lt(ecc, 1.0-1e-9)), 1e-9, -np.inf))
        # Kipping (2013) prior; equivalently xo.eccentricity.kipping13("ecc_prior", fixed=True, observed=ecc)
        pm.Potential("kipping", pm.logp(pm.Beta.dist(alpha=KIPPING_A, beta=KIPPING_B),
                                        pt.clip(ecc, 1e-8, 1 - 1e-8)).sum())

        gamma = pm.Normal("gamma", mu=np.array([offsets[nm] for nm in names]), sigma=100, shape=len(names))
        log_jit = pm.Normal("log_jit", mu=np.log(2.0), sigma=2.0, shape=len(names))

        orbit = xo.orbits.KeplerianOrbit(period=P, t0=tc, ecc=ecc, omega=omega)
        rv_planets = orbit.get_radial_velocity(tt, K=K)
        pm.Deterministic("rv_planets", rv_planets)
        mu = pt.sum(rv_planets, axis=-1) + gamma[idx]
        sigma = pt.sqrt(e**2 + pt.exp(2 * log_jit)[idx])
        pm.Normal("obs", mu=mu, sigma=sigma, observed=y)
    return model


def sample_xo(model, draws=1500, tune=1500, chains=4):
    with model:
        map_soln = pm.find_MAP(progressbar=False)
        names = [v.name for v in model.value_vars]
        init = {k: v for k, v in map_soln.items() if k in names}
        idata = pm.sample(draws=draws, tune=tune, chains=chains, cores=min(chains, os.cpu_count() or 1),
                          target_accept=0.95, initvals=init, random_seed=42,
                          idata_kwargs={"log_likelihood": True})
    return map_soln, idata


# %%
#######################
# JAXOPLANET FIT AND POSTERIOR SAMPLING
#######################
if RUN_JAXOPLANET:
    jx_args, jx_init, jx_tref = jx_inputs(df, guesses, offsets_guess)
    jx_map_soln = jx_map(jx_args, jx_init)
    print("MAP periods:", jx_map_soln["P"], " K:", jx_map_soln["K"], " e:", jx_map_soln["ecc"])

    jx_mcmc = jx_sample(jx_args, jx_map_soln)
    jx_mcmc.print_summary(exclude_deterministic=False)
    jx_samples = {k: np.asarray(v) for k, v in jx_mcmc.get_samples().items()}

    # Posteriors of the main parameters
    fig, axes = plt.subplots(3, len(guesses), figsize=(4.5 * len(guesses), 7), squeeze=False)
    for j in range(len(guesses)):
        for row, name in enumerate(["P", "K", "ecc"]):
            axes[row, j].hist(jx_samples[name][:, j], bins=50, histtype="step")
            axes[row, j].set_xlabel(f"{name} (planet {j + 1})")
    fig.tight_layout()
    plt.show()

    # Phase-folded plot per planet using the MAP solution
    idx, names = inst_codes(df)
    t_rel = df["t"].to_numpy() - jx_tref
    rv_pl = np.asarray(keplerian_rv(jnp.asarray(t_rel), jx_map_soln["P"], jx_map_soln["tc"],
                                    jx_map_soln["ecc"], jx_map_soln["omega"], jx_map_soln["K"]))
    base = df["rv"].to_numpy() - jx_map_soln["gamma"][idx]
    fig, axes = plt.subplots(1, len(guesses), figsize=(5 * len(guesses), 3.5), squeeze=False)
    for j, ax in enumerate(axes[0]):
        Pj, tcj = float(jx_map_soln["P"][j]), float(jx_map_soln["tc"][j])
        ph = ((t_rel - tcj) / Pj + 0.5) % 1 - 0.5
        order = np.argsort(ph)
        ax.errorbar(ph, base - (rv_pl.sum(axis=0) - rv_pl[j]), df["err"], fmt="o", ms=2.5, alpha=0.6)
        ax.plot(ph[order], rv_pl[j][order], "k-")
        ax.set_title(f"P = {Pj:.3f} d")
        ax.set_xlabel("Phase")
    axes[0, 0].set_ylabel("RV [m/s]")
    fig.tight_layout()
    plt.show()


# %%
#######################
# PERIODOGRAM OF THE RESIDUALS
#######################

# %%
#######################
# RESIDUALS OF THE RADVEL MAP FIT
#######################
def radvel_residuals(post):
    out = []
    for like in post.likelihood.like_list:
        out.append(pd.DataFrame({"t": like.x, "rv": like.residuals(), "err": like.errorbars(),
                                 "inst": like.suffix.lstrip("_")}))
    return pd.concat(out).sort_values("t").reset_index(drop=True)


res_df = radvel_residuals(post)
print("Residual rms per instrument [m/s]:")
print(res_df.groupby("inst")["rv"].std().round(2))

f_r = freq_grid(res_df["t"].to_numpy())
_, p_r, fal_r = ls_periodogram(res_df["t"], res_df["rv"], res_df["err"], f_r)
f_rb, logp_rb = run_bgls(res_df["t"].to_numpy(), res_df["rv"].to_numpy(), res_df["err"].to_numpy(), f_r)

fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
plot_periodogram(ax[0], f_r, p_r, fal=fal_r)
ax[0].set_title(f"GLS of residuals after {len(guesses)}-planet fit")
plot_periodogram(ax[1], f_rb, logp_rb, ylabel="BGLS log P - max")
ax[1].set_ylim(max(logp_rb.min(), -100), 5)
fig.tight_layout()
plt.show()

# %%
#######################
# ITERATIVE SEARCH FOR ADDITIONAL SIGNALS IN THE RESIDUALS
#######################
print("Prewhitening the residuals:")
pw_res = prewhiten(res_df, n_signals=3)
P_candidate = pw_res["periods"][0]
print(f"\nStrongest residual signal: P = {P_candidate:.4f} d (FAP {pw_res['faps'][0]:.1e})")
near = {k: v for k, v in KNOWN_PERIODS.items() if abs(v - P_candidate) / v < 0.02}
in_prot = PROT_RANGE[0] <= P_candidate <= PROT_RANGE[1]
print("Matches known planet:", near or "none", "| inside rotation range:", in_prot)
print("Remember: b/c resonant interactions also leave residual power near 30 and 61 d.")

# %%
#######################
# BAYESIAN MODEL COMPARISON
#######################

# %%
#######################
# CANDIDATE MODELS
#######################
model_periods = {1: pw["periods"][:1],
                 2: pw["periods"][:2],
                 3: pw["periods"][:2] + [P_candidate]}
model_guesses = {}
for n_pl, plist in model_periods.items():
    lf = linear_fit(df, plist)
    model_guesses[n_pl] = (lf["sines"], lf["offsets"])
    print(n_pl, "planet(s):", np.round(plist, 4))

# %%
#######################
# MAXIMUM LIKELIHOOD, BIC AND AIC (RADVEL MAP)
#######################
N_obs = len(df)
comp = {}
posts = {}
for n_pl, (gs, offs) in model_guesses.items():
    p_ = radvel_map(build_radvel_posterior(df, gs, offs))
    posts[n_pl] = p_
    lnL = p_.likelihood.logprob()
    k = len(p_.get_vary_params())
    comp[n_pl] = dict(k=k, lnL_max=lnL, BIC=k * np.log(N_obs) - 2 * lnL, AIC=2 * k - 2 * lnL)
comp_tab = pd.DataFrame(comp).T
comp_tab["dBIC_vs_prev"] = comp_tab["BIC"].diff()
print(comp_tab.round(2))
print("dBIC < -10 is very strong support for the extra planet (Kass & Raftery 1995).")

# %%
#######################
# BAYESIAN EVIDENCE WITH NESTED SAMPLING (DYNESTY)
#######################
# Period priors are restricted to +-3% around each periodogram peak, so this is a
# 'targeted' evidence; kima (below) explores the full period range.
# The prior volume matters: keep the K, gamma and jitter priors identical across models.


def make_nested_problem(post, guesses, offsets):
    names = post.likelihood.name_vary_params()
    lo, hi = [], []
    for nm in names:
        if m := re.fullmatch(r"per(\d+)", nm):
            P = guesses[int(m[1]) - 1]["P"]; lo.append(0.97 * P); hi.append(1.03 * P)
        elif m := re.fullmatch(r"tc(\d+)", nm):
            g = guesses[int(m[1]) - 1]; lo.append(g["tc"] - g["P"] / 2); hi.append(g["tc"] + g["P"] / 2)
        elif re.fullmatch(r"se(cos|sin)w\d+", nm):
            lo.append(-1.0); hi.append(1.0)
        elif re.fullmatch(r"k\d+", nm):
            lo.append(0.0); hi.append(500.0)
        elif nm.startswith("gamma_"):
            o = offsets[nm[6:]]; lo.append(o - 300); hi.append(o + 300)
        elif nm.startswith("jit_"):
            lo.append(0.0); hi.append(30.0)
        else:
            raise ValueError(f"No prior defined for {nm}")
    lo, hi = np.array(lo), np.array(hi)
    ecc_pairs = [(names.index(f"secosw{i}"), names.index(f"sesinw{i}")) for i in range(1, len(guesses) + 1)]

    def prior_transform(u):
        return lo + u * (hi - lo)

    def loglike(x):
        lp = 0.0
        for a, b in ecc_pairs:
            e = x[a] ** 2 + x[b] ** 2
            if e >= 1:
                return -1e300
            lp += stats.beta.logpdf(max(e, 1e-8), KIPPING_A, KIPPING_B)
        val = post.likelihood.logprob_array(x) + lp
        return val if np.isfinite(val) else -1e300

    return names, prior_transform, loglike


if RUN_DYNESTY:
    import dynesty
    for n_pl, (gs, offs) in model_guesses.items():
        names, ptf, llf = make_nested_problem(posts[n_pl], gs, offs)
        sampler = dynesty.NestedSampler(llf, ptf, len(names), nlive=500, bound="multi",
                                        sample="rwalk", rstate=np.random.default_rng(1))
        sampler.run_nested(dlogz=0.5, print_progress=False)
        r = sampler.results
        # The sampler draws (secosw, sesinw) in the square [-1,1]^2 and weights by Beta(e);
        # normalising to the unit-disk Kipping prior adds ln(4/pi) per planet.
        comp_tab.loc[n_pl, "lnZ"] = r.logz[-1] + n_pl * np.log(4 / np.pi)
        comp_tab.loc[n_pl, "lnZ_err"] = r.logzerr[-1]
        print(f"{n_pl} planet(s): lnZ = {comp_tab.loc[n_pl, 'lnZ']:.2f} +- {r.logzerr[-1]:.2f}")
    comp_tab["dlnZ_vs_prev"] = comp_tab["lnZ"].diff()
    print("dlnZ > 5 (Bayes factor > ~150) is strong evidence (Trotta 2008).")


# %%
#######################
# WAIC WITH JAXOPLANET POSTERIORS
#######################
# WAIC estimates how well each model would predict new data (elpd: expected log
# predictive density), penalising models that overfit. It is computed from the
# pointwise log-likelihood of every posterior sample.
if RUN_JX_COMPARISON:
    for n_pl, (gs, offs) in model_guesses.items():
        print(f"Sampling {n_pl}-planet jaxoplanet model...")
        cmp_args, cmp_init, _ = jx_inputs(df, gs, offs)
        cmp_mcmc = jx_sample(cmp_args, jx_map(cmp_args, cmp_init), draws=1000, tune=1000)
        elpd, se = waic(cmp_mcmc, cmp_args)
        comp_tab.loc[n_pl, "elpd_waic"] = elpd
        comp_tab.loc[n_pl, "waic_se"] = se
    comp_tab["d_elpd_vs_prev"] = comp_tab["elpd_waic"].diff()
    print(comp_tab[["elpd_waic", "waic_se", "d_elpd_vs_prev"]].round(2))
    print("Higher elpd_waic is better; a gain of several times the SE favours the larger model.")


# %%
#######################
# MODEL COMPARISON SUMMARY
#######################
print(comp_tab.round(2).to_string())
comp_tab.to_csv(OUTDIR / "model_comparison.csv")

# %%
#######################
# FULL ANALYSIS WITH KIMA (NUMBER OF PLANETS AS A FREE PARAMETER)
#######################

# %%
#######################
# WRITE PER-INSTRUMENT DATA FILES AND SET UP THE KIMA MODEL
#######################
if RUN_KIMA:
    import kima
    from kima import RVData, RVmodel, distributions

    kdir = (OUTDIR / "kima").resolve()
    kdir.mkdir(exist_ok=True)
    files = []
    for inst, g in df.groupby("inst"):
        fn = kdir / f"{inst}.txt"
        np.savetxt(fn, np.c_[g["t"], g["rv"], g["err"]], fmt="%.6f")
        files.append(str(fn))

    data = RVData(files, units="ms", skip=0)
    NPMAX = 6
    kmodel = RVmodel(fix=False, npmax=NPMAX, data=data)
    kmodel.conditional.Pprior = distributions.LogUniform(1.0, 2000.0)
    kmodel.conditional.Kprior = distributions.ModifiedLogUniform(1.0, 500.0)
    kmodel.conditional.eprior = distributions.Kumaraswamy(KIPPING_A, KIPPING_B)  # close match to Kipping Beta
    kmodel.Jprior = distributions.ModifiedLogUniform(0.1, 100.0)

# %%
#######################
# RUN KIMA (DIFFUSIVE NESTED SAMPLING)
#######################
# Use >= 100000 steps for a publication-grade run; this is a classroom setting.
KIMA_STEPS = 50000
if RUN_KIMA:
    cwd = os.getcwd()
    os.chdir(kdir)
    try:
        kima.run(kmodel, steps=KIMA_STEPS, num_threads=4)
        kres = kima.load_results(kmodel)
    finally:
        os.chdir(cwd)

# %%
#######################
# POSTERIOR ON THE NUMBER OF PLANETS AND THE DETECTED PERIODS
#######################
if RUN_KIMA:
    try:
        Np = kres.posterior_sample[:, kres.index_component].astype(int)
        counts = np.bincount(Np, minlength=NPMAX + 1)
        print("Posterior samples per Np:", dict(enumerate(counts)))
        n_detected = 0
        for k in range(1, NPMAX + 1):
            # Uniform prior on Np: posterior ratio = Bayes factor
            bf = counts[k] / counts[k - 1] if counts[k - 1] > 0 else np.inf
            print(f"  BF(Np={k} vs {k - 1}) = {bf:.1f}" if np.isfinite(bf) else
                  f"  BF(Np={k} vs {k - 1}) > {counts[k]} (no samples with {k - 1})")
            if bf > 150:
                n_detected = k
        print(f"Number of planets supported (BF > 150 threshold): {n_detected}")
    except Exception as exc:
        print("Could not compute Np Bayes factors directly:", exc)

    for fn in ["plot_posterior_np", "plot1", "plot_posterior_period", "plot2", "plot_random_samples"]:
        try:
            getattr(kres, fn)()
            plt.show()
        except Exception:
            pass

# %%
#######################
# DISCUSSION PROMPTS FOR STUDENTS
#######################
"""
Compare the number of planets favoured by kima with the radvel/exoplanet model comparison.
Do lnZ, BIC and WAIC agree? Why might they differ (prior volume, asymptotic assumptions)?
Which residual peaks coincide with the rotation range or with activity-indicator peaks?
How do the residuals change when you restrict TIME_WINDOW? What does this say about the
b/c resonance and the need for an N-body fit?
"""
