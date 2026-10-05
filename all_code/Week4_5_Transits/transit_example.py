# %% Importing packages

import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
import astropy.units as u
from astropy.timeseries import BoxLeastSquares
from scipy.signal import savgol_filter
from pytransit import QuadraticModel
import emcee
import corner

# %%
# ==============================================================================
# Search MAST and Download TESS Lightcurves for LHS 1140
# ==============================================================================
search_result = lk.search_lightcurve("LHS 1140", mission="TESS", author="SPOC", exptime=120)
print(search_result)

lc_collection = search_result.download_all()
print(f"Downloaded {len(lc_collection)} sector(s).")

# %%
# ==============================================================================
# Simple lightkurve plot - not so easy to use with gaps!
# ==============================================================================

lc_collection.stitch().plot()

# %%

# ==============================================================================
# Plot Light Curve per Sector (Stacked, No Time Gaps)
# ==============================================================================

# Filter and clean individual sectors
lc_clean_list = [lc.remove_outliers(sigma=5).normalize() for lc in lc_collection]

num_sectors = len(lc_clean_list)
fig, axes = plt.subplots(num_sectors, 1, figsize=(10, 2.2 * num_sectors), sharey=True)
if num_sectors == 1:
    axes = [axes]

for ax, lc in zip(axes, lc_clean_list):
    sec_num = lc.sector
    lc_binned = lc.bin(time_bin_size=1800 * u.s)

    # Raw 120s points (faint black dots)
    lc.plot(ax=ax, ls='none', marker='.', ms=1.5, color='k', alpha=0.2, label=f'Sector {sec_num} (120s)')
    # Binned 1800s points
    lc_binned.plot(ax=ax, ls='none', marker='o', ms=3.5, color='tab:red', label='Binned (1800s)')

    ax.set_ylabel("Norm. Flux")
    ax.set_title(f"LHS 1140 - TESS Sector {sec_num}")

axes[-1].set_xlabel("Time (BJD - 2457000)")
plt.tight_layout()
plt.show()

# %%

# ==============================================================================
# Initial Savitzky-Golay Detrending
# ==============================================================================

window_length = int(16*3600/120)+1  # Using ~16hr window
poly_order = 3

lc_flat_pass1_list = []
for lc in lc_clean_list:
    smoothed = savgol_filter(lc.flux.value, window_length, poly_order)
    lc_flat = lc.copy()
    lc_flat.flux = lc.flux / smoothed
    lc_flat_pass1_list.append(lc_flat)

lc_flat_pass1 = lk.LightCurveCollection(lc_flat_pass1_list).stitch().remove_nans()

# %%
# ==============================================================================
# Initial BLS Transit Search, Periodogram & Zoomed Plot
# ==============================================================================

t1 = lc_flat_pass1.time.value
f1 = lc_flat_pass1.flux.value
ferr1 = lc_flat_pass1.flux_err.value

bls1 = BoxLeastSquares(t1, f1, dy=ferr1)
period_grid = np.linspace(1.0, 30.0, 20000)
durations = np.linspace(0.04, 0.2, 12)

results1 = bls1.power(period_grid, durations)
best_idx1 = np.argmax(results1.power)
p1 = results1.period[best_idx1]
t0_1 = results1.transit_time[best_idx1]
dur1 = results1.duration[best_idx1]

print(f"Pass 1 BLS Peak: Period = {p1:.5f} d | T0 = {t0_1:.5f} BJD | Duration = {dur1*24:.2f} h")
# %%

# Plot 1: BLS Periodogram Pass 1
plt.figure(figsize=(10, 3.5))
plt.plot(results1.period, results1.power, 'k-', lw=1)
plt.axvline(p1, color='r', linestyle='--', label=f'P1: {p1:.4f} d')
plt.xlabel("Period [days]")
plt.ylabel("BLS Power")
plt.title("BLS Periodogram - Pass 1")
plt.legend()
plt.tight_layout()
plt.show()

# %%
# Plot 2: Zoom on Phase-Folded Transit Pass 1
phase1 = ((t1 - t0_1 + 0.5 * p1) % p1) - 0.5 * p1
mask_zoom1 = np.abs(phase1) < 0.2

plt.figure(figsize=(8, 4))
plt.plot(phase1[mask_zoom1], f1[mask_zoom1], 'k.', alpha=0.2, label='120s Data')
lc_pass1_binned = lc_flat_pass1.fold(period=p1, epoch_time=t0_1).bin(time_bin_size=900*u.s)
plt.plot(lc_pass1_binned.time.value, lc_pass1_binned.flux.value, 'ro', ms=4, label='Binned')
plt.xlabel("Time from Mid-Transit [days]")
plt.ylabel("Flattened Flux")
plt.title("Pass 1: Phase-Folded Transit Zoom")
plt.xlim(-0.2, 0.2)
plt.legend()
plt.tight_layout()
plt.show()

# %%
# ==============================================================================
# Mask Pass 1 Transits, Re-flatten & Rerun BLS (Pass 2)
# ==============================================================================

lc_flat_pass2_list = []

for lc in lc_clean_list:
    # 1. Identify in-transit points of Planet 1 using Astropy's transit_mask()
    in_tr1 = bls1.transit_mask(lc.time.value, period=p1, duration=dur1, transit_time=t0_1)

    # 2. Build smooth trend by interpolating over transits for the detrending step
    flux_trend = lc.flux.value.copy()
    if np.any(in_tr1):
        v_idx = np.where(~in_tr1)[0]
        t_idx = np.where(in_tr1)[0]
        flux_trend[in_tr1] = np.interp(t_idx, v_idx, lc.flux.value[v_idx])

    smoothed = savgol_filter(flux_trend, window_length, poly_order)

    # 3. Create flattened light curve
    lc_flat2 = lc.copy()
    lc_flat2.flux = lc.flux / smoothed

    # 4. REMOVE in-transit points of Planet 1 from this sector
    lc_flat2_clipped = lc_flat2[~in_tr1]
    lc_flat_pass2_list.append(lc_flat2_clipped)

lc_flat_pass2 = lk.LightCurveCollection(lc_flat_pass2_list).stitch().remove_nans()

t2 = lc_flat_pass2.time.value
f2 = lc_flat_pass2.flux.value
ferr2 = lc_flat_pass2.flux_err.value

bls2 = BoxLeastSquares(t2, f2, dy=ferr2)
results2 = bls2.power(period_grid, durations)
best_idx2 = np.argmax(results2.power)
p2 = results2.period[best_idx2]
t0_2 = results2.transit_time[best_idx2]
dur2 = results2.duration[best_idx2]

print(f"Pass 2 BLS Peak: Period = {p2:.5f} d | T0 = {t0_2:.5f} BJD | Duration = {dur2*24:.2f} h")
# %%

# Plot Pass 2 Periodogram
plt.figure(figsize=(10, 3.5))
plt.plot(results2.period, results2.power, 'tab:blue', lw=1)
plt.axvline(p2, color='r', linestyle='--', label=f'P2: {p2:.4f} d')
plt.xlabel("Period [days]")
plt.ylabel("BLS Power")
plt.title("BLS Periodogram - Pass 2 (Masked SavGol)")
plt.legend()
plt.tight_layout()
plt.show()

# %%

# Plot Pass 2 Transit Zoom
phase2 = ((t2 - t0_2 + 0.5 * p2) % p2) - 0.5 * p2
mask_zoom2 = np.abs(phase2) < 0.2

plt.figure(figsize=(8, 4))
plt.plot(phase2[mask_zoom2], f2[mask_zoom2], 'k.', alpha=0.2, label='120s Data')
lc_pass2_binned = lc_flat_pass2.fold(period=p2, epoch_time=t0_2).bin(time_bin_size=600*u.s)
plt.plot(lc_pass2_binned.time.value, lc_pass2_binned.flux.value, 'bo', ms=4, label='Binned')
plt.xlabel("Time from Mid-Transit [days]")
plt.ylabel("Flattened Flux")
plt.title("Pass 2: Phase-Folded Transit Zoom")
plt.xlim(-0.2, 0.2)
plt.legend()
plt.tight_layout()
plt.show()

# %%
# ==============================================================================
# Spline Detrending Masking All Identified Transits
# ==============================================================================

lc_spline_list = []
for lc in lc_clean_list:
    bls_sec = BoxLeastSquares(lc.time.value, lc.flux.value)
    mask_tr1 = bls_sec.transit_mask(lc.time.value, p1, dur1, t0_1)
    mask_tr2 = bls_sec.transit_mask(lc.time.value, p2, dur2, t0_2)
    combined_mask = mask_tr1 | mask_tr2

    # Lightkurve flatten() uses spline-fitting on masked data
    lc_spline = lc.flatten(window_length=501, mask=combined_mask)
    lc_spline_list.append(lc_spline)

lc_final = lk.LightCurveCollection(lc_spline_list).stitch().remove_nans()

# %%
# ==============================================================================
# PyTransit Model Fit
# ==============================================================================
from scipy.optimize import minimize

# Instantiate model and bind time data once
tm = QuadraticModel()

# If using exposure integration (e.g., 2-min cadence = 2/1440 days), add exptime and nsamples:
# tm.set_data(time_data, exptime=2.0/1440.0, nsamples=5)
tm.set_data(lc_final.time.value)

# ------------------------------------------------------------------------------
# Define Log-Likelihood & Negative Log-Likelihood Functions
# ------------------------------------------------------------------------------

def log_likelihood(params, time, flux, yerr):
    """
    Computes Gaussian log-likelihood for the PyTransit model.
    """
    k, t0, p, a, b, u1, u2 = params
    # Convert inclination from degrees to radians for PyTransit
    inc_rad = inc_rad = np.arccos(b / a)

    # Evaluate PyTransit model
    try:
        model_flux = tm.evaluate(
            k=k,
            ldc=[u1, u2],
            t0=t0,
            p=p,
            a=a,
            i=inc_rad,
            e=0.0,  # Circular orbit assumption
            w=0.0
        )
    except Exception:
        return -np.inf

    # Standard Gaussian log-likelihood equation
    sigma2 = yerr**2
    return -0.5 * np.sum(((flux - model_flux)**2 / sigma2) + np.log(2 * np.pi * sigma2))


def nll(params, time, flux, yerr):
    """
    Negative Log-Likelihood function for scipy.optimize.minimize.
    """
    ll = log_likelihood(params, time, flux, yerr)
    if not np.isfinite(ll):
        return 1e12  # Return high penalty value if parameters violate priors
    return -ll


# ------------------------------------------------------------------------------
# Optimize Model with Scipy
# ------------------------------------------------------------------------------
# Initial guesses: [k, t0, p, a, inc_deg, u1, u2]
initial_guess = [
    0.075,     # k (Rp/Rs)
    t0_1,     # t0
    p1,       # period [days]
    96.0,     # a/Rs
    0.5,      # b
    0.1,      # u1
    0.3       # u2
]

# Explicit parameter bounds for Scipy
bounds = [
    (0.001, 0.5),      # k
    (t0_1-0.2, t0_1+0.2),# t0
    (p1*0.95,p1*1.05),   # period
    (20.0, 150.0),     # a/Rs
    (0,1.1),           # b
    (0.0, 1.0),        # u1
    (0.0, 1.0)         # u2
]

# Run SciPy Nelder-Mead or L-BFGS-B optimizer
result = minimize(
    nll,
    x0=initial_guess,
    args=(lc_final.time.value, lc_final.flux.value, lc_final.flux_err.value),
    method='L-BFGS-B',
    bounds=bounds
)

# ------------------------------------------------------------------------------
# Extract Results & Plot Fit
# ------------------------------------------------------------------------------
best_params = result.x
print("Optimization Success:", result.success)
print("\nBest-Fit Parameters:")
print(f"  k (Rp/Rs)   : {best_params[0]:.5f}")
print(f"  t0          : {best_params[1]:.5f}")
print(f"  Period      : {best_params[2]:.5f}")
print(f"  a/Rs        : {best_params[3]:.3f}")
print(f"  Inclination : {best_params[4]:.3f} deg")
print(f"  u1, u2      : {best_params[5]:.3f}, {best_params[6]:.3f}")

# %%
# Generate best-fit light curve
best_model = tm.evaluate(
    k=best_params[0],
    ldc=[best_params[5], best_params[6]],
    t0=best_params[1],
    p=best_params[2],
    a=best_params[3],
    i=np.radians(best_params[4])
)

# Phase fold around primary planet transit (P2 / T0_2)
phase_final = ((lc_final.time.value - best_params[1] + 0.5 * best_params[2]) % best_params[2]) - 0.5 * best_params[2]
transit_window = np.abs(phase_final) < 0.2

t_data = phase_final[transit_window]
f_data = lc_final.flux.value[transit_window]
ferr_data = lc_final.flux_err.value[transit_window]

sort_idx = np.argsort(t_data)
t_data, f_data, ferr_data = t_data[sort_idx], f_data[sort_idx], ferr_data[sort_idx]

# Helper for binning raw phase data
def bin_phase_data(x, y, bin_width=0.005):
    bins = np.arange(np.min(x), np.max(x) + bin_width, bin_width)
    idx = np.digitize(x, bins)
    bx = [np.mean(x[idx == i]) for i in range(1, len(bins)) if np.sum(idx == i) > 0]
    by = [np.mean(y[idx == i]) for i in range(1, len(bins)) if np.sum(idx == i) > 0]
    return np.array(bx), np.array(by)

bin_t, bin_f = bin_phase_data(t_data, f_data, bin_width=0.006)
import matplotlib.pyplot as plt
plt.figure(figsize=(8, 4))
plt.errorbar(t_data, f_data, yerr=ferr_data, fmt='.k', alpha=0.2, label='Data')
plt.plot(bin_t, bin_f, '.', alpha=0.8, label='Binned Data')
plt.plot(np.sort(phase_final), best_model[np.argsort(phase_final)], 'r-', lw=2, label='PyTransit Best Fit')
plt.xlim(-0.2,0.2)
plt.xlabel("Time / Phase [days]")
plt.ylabel("Normalized Flux")
plt.legend()
plt.show()


# %%
# ==============================================================================
# MCMC Parameter Sampling with EMCEE & Corner Plot
# ==============================================================================


def log_prior(theta):
    k, t0, p, a, b, u1, u2 = theta

    # Flat, uninformative uniform priors within physical bounds
    if not (0.001 < k < 0.5):         return -np.inf  # Radius ratio
    if not (0.1 < p < 100.0):         return -np.inf  # Orbital period
    if not (2 < a < 250.0):           return -np.inf  # Semi-major axis in Rs
    if not (0.0 <= b <= 1+k):         return -np.inf  # Inclination
    if not (0.0 <= u1 <= 1.0):        return -np.inf  # Limb darkening u1
    if not (0.0 <= u2 <= 1.0):        return -np.inf  # Limb darkening u2
    if (u1 + u2 > 1.0):               return -np.inf  # Physical LDC constraint

    return 0.0

def log_probability(theta, time, flux, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, time, flux, yerr)

ndim = len(best_params)
nwalkers = 32
n_burnin = 500
n_steps = 2000

pos = best_params + 1e-4 * np.random.randn(nwalkers, ndim)

# Ensure initial walker positions respect priors
for i in range(nwalkers):
    while not np.isfinite(log_prior(pos[i])):
        pos[i] = initial_guess + 1e-4 * np.random.randn(ndim)

# %%

# Instantiate EnsembleSampler
sampler = emcee.EnsembleSampler(
    nwalkers,
    ndim,
    log_probability,
    args=(lc_final.time.value, lc_final.flux.value, lc_final.flux_err.value)
)

# Run Burn-in
print(f"--> Running burn-in phase ({n_burnin} steps)...")
state = sampler.run_mcmc(pos, n_burnin, progress=True)
sampler.reset()

# Run Production Chain
print(f"--> Running production MCMC sampling ({n_steps} steps)...")
sampler.run_mcmc(state, n_steps, progress=True)

# %%
# ------------------------------------------------------------------------------
# Diagnostics & Parameter Extraction
# ------------------------------------------------------------------------------
# Extract posterior flat samples
flat_samples = sampler.get_chain(flat=True)

print(f"\nSuccessfully collected {flat_samples.shape[0]} posterior samples.")
print("\n--- Parameter Posterior Estimates (Median +/- 1-Sigma) ---")

labels = ["k (Rp/Rs)", "t0", "Period [d]", "a/Rs", "Inc [deg]", "u1", "u2"]
results_summary = {}
for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    results_summary[labels[i]] = mcmc[1]
    print(f"{labels[i]:12s}: {mcmc[1]:.6f} (+{q[1]:.6f} / -{q[0]:.6f})")


# ------------------------------------------------------------------------------
# Visualizing Corner Plot & Best-Fit Model
# ------------------------------------------------------------------------------
# Corner plot of posteriors
fig = corner.corner(
    flat_samples,
    labels=labels,
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={"fontsize": 10}
)
plt.show()
# %%
# Generate model comparison using posterior samples
plt.figure(figsize=(9, 4.5))
plt.errorbar(time_data, flux_data, yerr=err_data, fmt=".k", alpha=0.15, label="Data")

# Draw 100 random posterior light curves to plot uncertainty band
inds = np.random.randint(len(flat_samples), size=100)
for ind in inds:
    sample = flat_samples[ind]
    sample_model = tm.evaluate(
        k=sample[0],
        ldc=[sample[5], sample[6]],
        t0=sample[1],
        p=sample[2],
        a=sample[3],
        i=np.radians(sample[4])
    )
    plt.plot(time_data, sample_model, "r-", alpha=0.03)

# Best-fit median curve
best_k, best_t0, best_p, best_a, best_inc, best_u1, best_u2 = [results_summary[lbl] for lbl in labels]
median_model = tm.evaluate(
    k=best_k, ldc=[best_u1, best_u2], t0=best_t0, p=best_p, a=best_a, i=np.radians(best_inc)
)

plt.plot(time_data, median_model, "r-", lw=2, label="MCMC Median Fit")
plt.xlabel("Time / Phase [days]")
plt.ylabel("Normalized Flux")
plt.legend()
plt.tight_layout()
plt.show()
