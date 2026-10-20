# %%
#######################
# INTRODUCTION TO RV PLANET ANALYSIS: 51 PEGASI
#######################
"""
A first look at radial-velocity (RV) planet hunting, using 51 Pegasi, the first
exoplanet found around a Sun-like star (Mayor & Queloz 1995).

Steps
    Download the RVs from DACE
    Remove obvious anomalies
    Find the strongest periodic signal with a BGLS periodogram
    Fit a 1-planet model with radvel (maximum a posteriori, MAP)
    Look for a second signal in the residuals
    Fit a 2-planet model and compare the two with the BIC
    Compare 1 and 2 planets independently with kima (Bayesian evidence)

Install (with uv; use Python 3.11 or 3.12)
    uv add dace-query numpy pandas matplotlib radvel kima
"""

# %%
#######################
# IMPORTS AND SETTINGS
#######################
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize

from RVs_functions import *

STAR = "HD217014"        # 51 Peg; try "51Peg" if DACE does not recognise this name
HOME_DIR = Path.cwd()   # kima changes folder while running; we come back here
OUTDIR = Path("51peg_outputs")
OUTDIR.mkdir(exist_ok=True)

# %%
#######################
# DOWNLOAD THE DATA FROM DACE
#######################
from dace_query.spectroscopy import Spectroscopy

raw = pd.DataFrame(Spectroscopy.get_timeseries(STAR, sorted_by_instrument=False,
                                               output_format="pandas"))
print("Number of rows:", len(raw))
print(raw.groupby(["instrument_name", "drs_id"]).size())

# %%
#######################
# KEEP ONE DATA REDUCTION PER INSTRUMENT
#######################
# DACE can hold the same observations reduced with several pipeline (DRS) versions.
# Keep the version with the most points for each instrument so nothing is counted twice.
raw["drs_id"] = raw["drs_id"].fillna("unknown").astype(str)
counts = (raw.groupby(["instrument_name", "drs_id"]).size().reset_index(name="n")
          .sort_values(["instrument_name", "n", "drs_id"]))
best = counts.groupby("instrument_name").tail(1)[["instrument_name", "drs_id"]]
raw = raw.merge(best, on=["instrument_name", "drs_id"])

# A simple table with the four columns we need: time, RV, RV error, instrument
df = pd.DataFrame({
    "t": raw["rjd"].astype(float),            # BJD - 2400000 [days]
    "rv": raw["rv"].astype(float),            # [m/s]
    "err": raw["rv_err"].astype(float),       # [m/s]
    "inst": raw["instrument_name"].astype(str).str.replace(r"[^A-Za-z0-9]", "", regex=True),
}).sort_values("t").reset_index(drop=True)
print(df.groupby("inst").size())

# %%
#######################
# REMOVE ANOMALIES
#######################
# Anomalies: missing values, non-positive errors, and points far from the rest of their
# instrument's data. We use the median absolute deviation (MAD), which, unlike the
# standard deviation, is not inflated by the outliers themselves. 5 MAD is a generous
# cut: the planet's ~56 m/s signal is well inside it.
n_before = len(df)
df = df[np.isfinite(df["t"]) & np.isfinite(df["rv"]) & np.isfinite(df["err"]) & (df["err"] > 0)]

keep = np.ones(len(df), dtype=bool)
for inst in df["inst"].unique():
    m = (df["inst"] == inst).to_numpy()
    rv, err = df["rv"].to_numpy()[m], df["err"].to_numpy()[m]
    mad_rv = 1.4826 * np.median(np.abs(rv - np.median(rv)))
    mad_err = 1.4826 * np.median(np.abs(err - np.median(err)))
    bad = (np.abs(rv - np.median(rv)) > 5 * mad_rv) | (err > np.median(err) + 5 * mad_err)
    keep[np.where(m)[0][bad]] = False
    print(f"{inst:12s}: {m.sum():4d} points, {bad.sum():3d} removed")
df = df[keep]

# Instruments with very few points cannot constrain their own offset
n_per_inst = df["inst"].value_counts()
df = df[df["inst"].isin(n_per_inst[n_per_inst >= 5].index)].reset_index(drop=True)
print(f"\nKept {len(df)} of {n_before} points")

insts = sorted(df["inst"].unique())
t, rv, err = df["t"].to_numpy(), df["rv"].to_numpy(), df["err"].to_numpy()

# %%
#######################
# PLOT THE DATA
#######################
# Each instrument has its own zero point, so subtract each instrument's median RV.
rv_centered = rv - df.groupby("inst")["rv"].transform("median").to_numpy()

fig, ax = plt.subplots(figsize=(10, 4))
for inst in insts:
    m = (df["inst"] == inst).to_numpy()
    ax.errorbar(t[m], rv_centered[m], err[m], fmt="o", ms=3, label=inst)
ax.set_xlabel("Time [BJD - 2400000]")
ax.set_ylabel("RV - median [m/s]")
ax.legend()
plt.show()

# %%

# %%
#######################
# BGLS PERIODOGRAM OF THE RVS
#######################
# The highest peak of the BGLS periodogram is our first planet candidate.
PMIN, PMAX = 1.1, np.ptp(t)
periods, log_prob = bgls(t, rv_centered, err, PMIN, PMAX)

P1 = periods[np.argmax(log_prob)]
print(f"Highest BGLS peak: P = {P1:.4f} days")

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(periods, log_prob / np.log(10), lw=0.8)        # convert ln to log10
ax.axvline(P1, color="r", ls="--", alpha=0.5)
ax.set_xscale("log")
ax.set_ylim(-50, 1)
ax.set_xlabel("Period [days]")
ax.set_ylabel("log10 relative probability")
ax.set_title("BGLS periodogram of the RVs")
plt.show()

# %%
#######################
# INITIAL GUESS FOR THE ORBIT
#######################
# radvel needs a starting point. Fitting a sine wave at P1 (a linear least-squares fit)
# gives the amplitude K and the time of conjunction Tc. The extra columns let each
# instrument have its own offset.
t0 = np.median(t)
phase = 2 * np.pi * (t - t0) / P1
X = np.column_stack([np.sin(phase), np.cos(phase)] +
                    [(df["inst"] == inst).to_numpy().astype(float) for inst in insts])
coef = np.linalg.lstsq(X / err[:, None], rv / err, rcond=None)[0]
s, c = coef[0], coef[1]
K1_guess = np.hypot(s, c)
tc1_guess = t0 + np.arctan2(c, -s) * P1 / (2 * np.pi)   # v(t) = -K sin(2 pi (t - Tc) / P)
gamma_guess = dict(zip(insts, coef[2:]))
print(f"Guess: P = {P1:.4f} d, K = {K1_guess:.1f} m/s, Tc = {tc1_guess:.3f}")

# %%
#######################
# RADVEL: SET UP THE 1-PLANET MODEL
#######################
# Fitting basis: period, time of conjunction, sqrt(e)cos(w), sqrt(e)sin(w), semi-amplitude.
# Each instrument gets an offset (gamma) and a jitter term (extra noise added to the errors).
import radvel

params1 = radvel.Parameters(1, basis="per tc secosw sesinw k")
params1["per1"] = radvel.Parameter(value=P1)
params1["tc1"] = radvel.Parameter(value=tc1_guess)
params1["secosw1"] = radvel.Parameter(value=0.01)
params1["sesinw1"] = radvel.Parameter(value=0.01)
params1["k1"] = radvel.Parameter(value=K1_guess)
params1["dvdt"] = radvel.Parameter(value=0.0, vary=False)
params1["curv"] = radvel.Parameter(value=0.0, vary=False)
model1 = radvel.RVModel(params1, time_base=t0)

likes = []
for inst in insts:
    m = (df["inst"] == inst).to_numpy()
    like = radvel.likelihood.RVLikelihood(model1, t[m], rv[m], err[m], suffix="_" + inst)
    like.params["gamma_" + inst] = radvel.Parameter(value=gamma_guess[inst])
    like.params["jit_" + inst] = radvel.Parameter(value=2.0)
    likes.append(like)
post1 = radvel.posterior.Posterior(radvel.likelihood.CompositeLikelihood(likes))

post1.priors += [radvel.prior.EccentricityPrior(1)]          # keeps 0 <= e < 1
post1.priors += [radvel.prior.PositiveKPrior(1)]
for inst in insts:
    post1.priors += [radvel.prior.HardBounds("jit_" + inst, 0.0, 50.0)]

# %%
#######################
# RADVEL: MAP FIT OF THE 1-PLANET MODEL
#######################
# Maximise the posterior probability (minimise its negative). Running two optimisers
# in turn helps the fit settle.
for method in ["Powell", "Nelder-Mead", "Powell"]:
    result = optimize.minimize(post1.neglogprob_array, post1.get_vary_params(),
                               method=method, options={"maxiter": 100000, "maxfev": 100000})
    post1.set_vary_params(result.x)
print(post1)

e1 = post1.params["secosw1"].value**2 + post1.params["sesinw1"].value**2
print(f"\nP = {post1.params['per1'].value:.5f} d, K = {post1.params['k1'].value:.2f} m/s, e = {e1:.3f}")

# %%
#######################
# PLOT THE 1-PLANET FIT (PHASE-FOLDED)
#######################
P_fit, tc_fit = post1.params["per1"].value, post1.params["tc1"].value
rv_minus_gamma = rv - df["inst"].map({i: post1.params["gamma_" + i].value for i in insts}).to_numpy()
phase = ((t - tc_fit) / P_fit) % 1

t_model = tc_fit + np.linspace(0, P_fit, 500)
fig, ax = plt.subplots(figsize=(8, 4))
for inst in insts:
    m = (df["inst"] == inst).to_numpy()
    ax.errorbar(phase[m], rv_minus_gamma[m], err[m], fmt="o", ms=3, alpha=0.6, label=inst)
ax.plot(np.linspace(0, 1, 500), model1(t_model), "k-", lw=2)
ax.set_xlabel("Orbital phase")
ax.set_ylabel("RV [m/s]")
ax.set_title(f"51 Peg b: P = {P_fit:.4f} d")
ax.legend()
plt.show()

# %%
#######################
# RESIDUALS OF THE 1-PLANET FIT
#######################
# Residuals = data - (planet + instrument offset). Any remaining periodic signal could be
# another planet (or stellar activity, or an instrumental effect).
res_t, res_rv, res_err, res_inst = [], [], [], []
for like in post1.likelihood.like_list:
    res_t.append(like.x)
    res_rv.append(like.residuals())
    res_err.append(like.errorbars())         # error bars including the fitted jitter
    res_inst += [like.suffix[1:]] * len(like.x)
res_t, res_rv, res_err = np.concatenate(res_t), np.concatenate(res_rv), np.concatenate(res_err)
print(f"RMS of the residuals: {np.std(res_rv):.2f} m/s (median error {np.median(err):.2f} m/s)")

fig, ax = plt.subplots(figsize=(10, 3.5))
ax.errorbar(res_t, res_rv, res_err, fmt="o", ms=3)
ax.axhline(0, color="k")
ax.set_xlabel("Time [BJD - 2400000]")
ax.set_ylabel("Residual RV [m/s]")
plt.show()

# %%
#######################
# BGLS PERIODOGRAM OF THE RESIDUALS
#######################
periods_r, log_prob_r = bgls(res_t, res_rv, res_err, 10, 6000)

P2 = periods_r[np.argmax(log_prob_r)]
print(f"Highest peak in the residuals: P = {P2:.4f} days")

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(periods_r, log_prob_r / np.log(10), lw=0.8)
ax.axvline(P2, color="r", ls="--", alpha=0.5, label=f"{P2:.2f} d")
ax.axvline(P_fit, color="gray", ls=":", label="planet b")
ax.set_xscale("log")
ax.set_ylim(-20, 1)
ax.set_xlabel("Period [days]")
ax.set_ylabel("log10 relative probability")
ax.set_title("BGLS periodogram of the residuals")
ax.legend()
plt.show()

# Questions: is the peak as dominant as in the first periodogram? Is it near 1 day, a year,
# or the length of the data set? Those suggest aliases, seasonal or instrumental effects.

# %%
#######################
# INITIAL GUESS FOR THE SECOND SIGNAL
#######################
t0r = np.median(res_t)
phase = 2 * np.pi * (res_t - t0r) / P2
X = np.column_stack([np.sin(phase), np.cos(phase), np.ones_like(res_t)])
coef = np.linalg.lstsq(X / res_err[:, None], res_rv / res_err, rcond=None)[0]
K2_guess = np.hypot(coef[0], coef[1])
tc2_guess = t0r + np.arctan2(coef[1], -coef[0]) * P2 / (2 * np.pi)
print(f"Guess: P = {P2:.4f} d, K = {K2_guess:.2f} m/s, Tc = {tc2_guess:.3f}")

# %%
#######################
# RADVEL: SET UP THE 2-PLANET MODEL
#######################
# Planet 1 starts from the 1-planet MAP solution; planet 2 from the guess above.
params2 = radvel.Parameters(2, basis="per tc secosw sesinw k")
for name in ["per1", "tc1", "secosw1", "sesinw1", "k1"]:
    params2[name] = radvel.Parameter(value=post1.params[name].value)
params2["per2"] = radvel.Parameter(value=P2)
params2["tc2"] = radvel.Parameter(value=tc2_guess)
params2["secosw2"] = radvel.Parameter(value=0.01)
params2["sesinw2"] = radvel.Parameter(value=0.01)
params2["k2"] = radvel.Parameter(value=max(K2_guess, 1.0))
params2["dvdt"] = radvel.Parameter(value=0.0, vary=False)
params2["curv"] = radvel.Parameter(value=0.0, vary=False)
model2 = radvel.RVModel(params2, time_base=t0)

likes = []
for inst in insts:
    m = (df["inst"] == inst).to_numpy()
    like = radvel.likelihood.RVLikelihood(model2, t[m], rv[m], err[m], suffix="_" + inst)
    like.params["gamma_" + inst] = radvel.Parameter(value=post1.params["gamma_" + inst].value)
    like.params["jit_" + inst] = radvel.Parameter(value=post1.params["jit_" + inst].value)
    likes.append(like)
post2 = radvel.posterior.Posterior(radvel.likelihood.CompositeLikelihood(likes))

post2.priors += [radvel.prior.EccentricityPrior(2)]
post2.priors += [radvel.prior.PositiveKPrior(2)]
for inst in insts:
    post2.priors += [radvel.prior.HardBounds("jit_" + inst, 0.0, 50.0)]

# %%
#######################
# RADVEL: MAP FIT OF THE 2-PLANET MODEL
#######################
for method in ["Powell", "Nelder-Mead", "Powell"]:
    result = optimize.minimize(post2.neglogprob_array, post2.get_vary_params(),
                               method=method, options={"maxiter": 100000, "maxfev": 100000})
    post2.set_vary_params(result.x)
print(post2)

# %%
#######################
# PLOT THE 2-PLANET FIT (EACH SIGNAL PHASE-FOLDED)
#######################

rv_minus_gamma = rv - df["inst"].map({i: post2.params["gamma_" + i].value for i in insts}).to_numpy()
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for n, ax in zip([1, 2], axes):
    other = 2 if n == 1 else 1
    P_n, tc_n = post2.params[f"per{n}"].value, post2.params[f"tc{n}"].value
    rv_this_planet = rv_minus_gamma - model2(t, planet_num=other)   # remove the other planet
    ax.errorbar(((t - tc_n) / P_n) % 1, rv_this_planet, err, fmt="o", ms=3, alpha=0.5)
    ax.plot(np.linspace(0, 1, 500), model2(tc_n + np.linspace(0, P_n, 500), planet_num=n), "k-", lw=2)
    ax.set_title(f"Signal {n}: P = {P_n:.4f} d, K = {post2.params[f'k{n}'].value:.2f} m/s")
    ax.set_xlabel("Orbital phase")
axes[0].set_ylabel("RV [m/s]")
fig.tight_layout()
plt.show()

# %%
#######################
# MODEL COMPARISON WITH THE BIC
#######################
# BIC = k ln(N) - 2 ln(L_max). It rewards a better fit (higher likelihood) but penalises
# every extra free parameter k. Lower BIC is better.
# Delta BIC = BIC(2 planets) - BIC(1 planet):
#     below -10  -> very strong support for the second planet
#     -10 to -6  -> strong;  -6 to -2 -> positive;  above -2 -> not worth the extra parameters
N = len(t)
lnL1, k1 = post1.likelihood.logprob(), len(post1.get_vary_params())
lnL2, k2 = post2.likelihood.logprob(), len(post2.get_vary_params())
BIC1 = k1 * np.log(N) - 2 * lnL1
BIC2 = k2 * np.log(N) - 2 * lnL2

print(f"1 planet : ln L = {lnL1:9.2f}, k = {k1:2d}, BIC = {BIC1:9.2f}")
print(f"2 planets: ln L = {lnL2:9.2f}, k = {k2:2d}, BIC = {BIC2:9.2f}")
print(f"Delta BIC (2 - 1) = {BIC2 - BIC1:.2f}")

# %%
#######################
# KIMA: PREPARE THE DATA FILES
#######################
# kima works on its own and does not use anything from radvel: it explores all
# parameters (periods included) over wide priors and computes the Bayesian evidence Z.
#
# Each instrument is shifted to a median of zero before it goes to kima. DACE RVs can be
# absolute velocities (tens of km/s) with zero points differing by kilometres per second
# between instruments. kima's default offset priors then span that whole range, and a
# 56 m/s planet is easily missed while offsets, or long-period 'planets' imitating them,
# absorb the differences. After centring, every offset is within a few tens of m/s of 0.
import kima
from kima import RVData, RVmodel, distributions

kima_dir = (OUTDIR / "kima").resolve()
kima_dir.mkdir(exist_ok=True)
data_files = []
for inst in insts:
    m = (df["inst"] == inst).to_numpy()
    filename = kima_dir / f"{inst}.txt"
    rv_shifted = rv[m] - np.median(rv[m])
    np.savetxt(filename, np.column_stack([t[m], rv_shifted, err[m]]), fmt="%.6f")
    print(f"{inst:12s}: subtracted median {np.median(rv[m]):12.2f} m/s")
    data_files.append(str(filename))
data = RVData(data_files, units="ms", skip=0)

# %%
#######################
# KIMA: PRIORS
#######################
# Because the data are centred, the systemic velocity and the instrument offsets only need
# to cover +-OFFSET_RANGE m/s. Keep the range comfortably wider than the planet's K so the
# prior does not force offsets that hide part of the signal. kima measures the offsets
# relative to one reference instrument, so there is one fewer offset than instruments.
# The period prior stops at the time span of the data: longer 'periods' are unconstrained
# and are an easy way for the sampler to imitate offsets.
OFFSET_RANGE = 100.0   # m/s
KIMA_C_PRIOR = distributions.Uniform(-OFFSET_RANGE, OFFSET_RANGE)
KIMA_OFFSET_PRIORS = [distributions.Uniform(-OFFSET_RANGE, OFFSET_RANGE)] * (len(insts) - 1)
KIMA_P_PRIOR = distributions.LogUniform(1.0, np.ptp(t))
KIMA_K_PRIOR = distributions.ModifiedLogUniform(1.0, 500.0)

# %%
#######################
# KIMA: RUN WITH EXACTLY 1 PLANET
#######################
# fix=True with npmax=1 means 'exactly one planet'. Each run writes its output files
# into the current folder, so each model gets its own folder.
KIMA_STEPS = 50000   # too few steps is another reason kima can miss a planet

kima_model_1 = RVmodel(fix=True, npmax=1, data=data)
kima_model_1.Cprior = KIMA_C_PRIOR
if hasattr(kima_model_1, "individual_offset_prior"):
    kima_model_1.individual_offset_prior = KIMA_OFFSET_PRIORS
else:
    kima_model_1.offsets_prior = KIMA_OFFSET_PRIORS[0]
kima_model_1.conditional.Pprior = KIMA_P_PRIOR
kima_model_1.conditional.Kprior = KIMA_K_PRIOR
run_dir_1 = kima_dir / "one_planet"
run_dir_1.mkdir(exist_ok=True)
os.chdir(run_dir_1)
kima.run(kima_model_1, steps=KIMA_STEPS, num_threads=4)
kima_results_1 = kima.load_results(kima_model_1)
os.chdir(HOME_DIR)

# %%
#######################
# KIMA: RUN WITH EXACTLY 2 PLANETS
#######################
kima_model_2 = RVmodel(fix=True, npmax=2, data=data)
kima_model_2.Cprior = KIMA_C_PRIOR
if hasattr(kima_model_2, "individual_offset_prior"):
    kima_model_2.individual_offset_prior = KIMA_OFFSET_PRIORS
else:
    kima_model_2.offsets_prior = KIMA_OFFSET_PRIORS[0]
kima_model_2.conditional.Pprior = KIMA_P_PRIOR
kima_model_2.conditional.Kprior = KIMA_K_PRIOR
run_dir_2 = kima_dir / "two_planets"
run_dir_2.mkdir(exist_ok=True)
os.chdir(run_dir_2)
kima.run(kima_model_2, steps=KIMA_STEPS, num_threads=4)
kima_results_2 = kima.load_results(kima_model_2)
os.chdir(HOME_DIR)

# %%
#######################
# KIMA: COMPARE THE EVIDENCES
#######################
# ln Z is the probability of the data averaged over the whole prior. The Bayes factor
# B = Z2 / Z1, i.e. ln B = ln Z2 - ln Z1:
#     ln B > 5   -> strong evidence for 2 planets (B > ~150)
#     ln B < 0   -> the data prefer 1 planet
lnZ1 = kima_results_1.evidence
lnZ2 = kima_results_2.evidence
print(f"kima ln Z (1 planet)  = {lnZ1:.2f}")
print(f"kima ln Z (2 planets) = {lnZ2:.2f}")
print(f"ln Bayes factor (2 vs 1) = {lnZ2 - lnZ1:.2f}")

# Look at what kima found for the periods in the 2-planet run.
# The name of this plotting method differs between kima versions, so try the usual ones
# and otherwise list the plotting methods your version has.
for method_name in ["plot_posterior_period", "plot2"]:
    if hasattr(kima_results_2, method_name):
        getattr(kima_results_2, method_name)()
        plt.show()
        break
else:
    print("Plotting methods in this kima version:",
          [m for m in dir(kima_results_2) if m.startswith("plot")])

# %%
#######################
# KIMA: PHASE-FOLDED SIGNALS OF THE BEST SOLUTION
#######################
# Take the posterior sample with the highest likelihood from each kima run and let kima
# draw each planet's signal folded on its period (the other planet removed). Compare
# these with the radvel plots above: kima found its solution without any starting guess.
# Method names differ between kima versions, so the code tries the common ones.
for label, results in [("1 planet", kima_results_1), ("2 planets", kima_results_2)]:
    # Best (maximum-likelihood) posterior sample
    if hasattr(results, "maximum_likelihood_sample"):
        try:
            best_sample = results.maximum_likelihood_sample(printit=False)
        except TypeError:
            best_sample = results.maximum_likelihood_sample()
    else:
        best_sample = results.posterior_sample[np.argmax(results.posterior_lnlike[:, 1])]

    # Phase-folded plot of that sample
    if hasattr(results, "phase_plot"):
        results.phase_plot(best_sample)
        plt.suptitle(f"kima, {label}: maximum-likelihood solution")
        plt.show()
    else:
        print("No phase_plot in this kima version. Plotting methods available:",
              [m for m in dir(results) if "plot" in m])
        break

# %%
#######################
# QUESTIONS FOR STUDENTS
#######################
"""
Do the BIC and the kima evidence agree about a second planet?
In the 2-planet kima run, is the second period well defined, or spread over many values?
Could the second signal be stellar activity (the star's rotation period is several weeks),
an alias of planet b, or a long-term offset between instruments?
"""
