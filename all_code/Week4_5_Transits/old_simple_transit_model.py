# %% Importing packages

import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
from astropy.timeseries import BoxLeastSquares
from scipy.signal import savgol_filter
import batman
import emcee
import corner
import astropy.units as u

# %% ACCESSING DATA FROM LIGHTKURVE

search_result = lk.search_lightcurve("LHS 1140", mission="TESS", author="SPOC")
print(search_result)

# %% Downloading

# - We dont need 20s data (it is very heavy)
# - NB - you may need to modify this for FFI-only lightcurves (where exptime could be 200, 600, 1200, 1800)
lc_collection = search_result[(search_result.exptime.value==120)&(search_result.author=="SPOC")].download_all()

# %% Plotting

# Filter/clean individual light curves in the collection
lc_cleaned_list = [lc.remove_outliers(sigma=5).normalize() for lc in lc_collection]

# 2. Create stacked subplots (one for each sector)
num_sectors = len(lc_cleaned_list)
fig, axes = plt.subplots(num_sectors, 1, figsize=(10, 2.5 * num_sectors), sharey=True)

# Ensure axes is iterable even if there's only 1 sector
if num_sectors == 1:
    axes = [axes]

# 3. Plot each sector in its own subplot
for i, lc in enumerate(lc_cleaned_list):
    ax = axes[i]
    sector_num = lc.sector

    # Create binned data for this sector
    lc_binned = lc.bin(time_bin_size=1800*u.s)

    # Plot raw 120s dots
    lc.plot(
        ax=ax,
        ls='none',
        marker='.',
        ms=1.5,
        color='k',
        alpha=0.25,
        label=f'Sector {sector_num} (120s)'
    )

    # Plot binned 1800s dots
    lc_binned.plot(
        ax=ax,
        ls='none',
        marker='o',
        ms=4,
        color='tab:red',
        label=f'Binned (1800s)'
    )

    ax.set_ylabel("Norm. Flux")
    ax.set_title(f"Gliese 12 - TESS Sector {sector_num}")
    ax.set_ylim(1+2.5*(np.nanpercentile(lc_binned.flux,[1,99])-1))
# Common x-label for the bottom subplot
axes[-1].set_xlabel("Time (BJD - 2457000)")

plt.tight_layout()
plt.show()

# %% Cleaning, Flattening & Normalising

# Process and flatten each sector individually
window_length = int((18*3600/120)+1) # Window length must be an odd integer
poly_order = 3

lc_flattened_list = []
for lc in lc_collection:
    lc_clean = lc.remove_outliers(sigma=5).normalize()

    # Perform smooth detrending using Savitzky-Golay filter
    smoothed_trend = savgol_filter(lc_clean.flux.value, window_length, poly_order)

    # Create flattened copy
    lc_flat = lc_clean.copy()
    lc_flat.flux = lc_clean.flux / smoothed_trend

    lc_flattened_list.append(lc_flat)

# 2. Plot sector by sector from the list
fig, axes = plt.subplots(len(lc_flattened_list), 1, figsize=(10, 2.5 * len(lc_flattened_list)), sharey=True)

if len(lc_flattened_list) == 1:
    axes = [axes]

for ax, lc_sec in zip(axes, lc_flattened_list):
    sec_num = lc_sec.sector
    lc_sec_binned = lc_sec.bin(time_bin_size=1800*u.s)

    lc_sec.plot(ax=ax, ls='none', marker='.', ms=1.5, color='tab:blue', alpha=0.3, label='Flattened (120s)')
    lc_sec_binned.plot(ax=ax, ls='none', marker='o', ms=4, color='tab:orange', label='Binned (1800s)')
    ax.set_ylim(1+4*(np.nanpercentile(lc_sec_binned.flux,[1,99])-1))
    ax.axhline(1.0, color='r', linestyle='--', alpha=0.5)
    ax.set_ylabel("Flux")
    ax.set_title(f"LHS 1140 - Sector {sec_num}")

plt.tight_layout()
plt.show()

# %%

# ==============================================================================
# STEP 4 & 5: Box Least Squares (BLS) Transit Search and Periodogram
# ==============================================================================
print("\n--> Running Box Least Squares (BLS) periodogram...")

time = np.hstack([lcf.time.value for lcf in lc_flattened_list])
flux = np.hstack([lcf.flux.value for lcf in lc_flattened_list])
flux_err = np.hstack([lcf.flux_err.value for lcf in lc_flattened_list])

# Create Astropy BLS model
model = BoxLeastSquares(time, flux, dy=flux_err)

# Define search grid based on Gliese 12 b's known ~12.76-day orbital period
period_grid = np.linspace(1.0, 30.0, 10000)
durations = np.linspace(0.05, 0.2, 10)  # ~1.2 to 4.8 hours duration

bls_results = model.power(period_grid, durations)

# Extract best-fit parameters
best_idx = np.argmax(bls_results.power)
best_period = bls_results.period[best_idx]
best_t0 = bls_results.transit_time[best_idx]
best_duration = bls_results.duration[best_idx]
best_depth = bls_results.depth[best_idx]

print(f"BLS Results:")
print(f"  Period:   {best_period:.5f} days")
print(f"  Epoch T0: {best_t0:.5f} BJD")
print(f"  Duration: {best_duration*24:.2f} hours")
print(f"  Depth:    {best_depth*1e3:.2f} ppt")

# Plot BLS Periodogram
plt.figure(figsize=(10, 4))
plt.plot(bls_results.period, bls_results.power, 'k-', lw=1)
plt.axvline(best_period, color='r', linestyle='--', label=f'Peak Period: {best_period:.4f} d')
plt.xlabel("Period [days]")
plt.ylabel("BLS Power")
plt.title("BLS Periodogram")
plt.legend()
plt.tight_layout()
plt.show()
# %%

# Create in-transit mask based on BLS parameters
transit_mask = model.in_transit(time, best_period, best_duration, best_t0)

# Use Lightkurve's spline/Savitzky-Golay detrending on out-of-transit data
lc_out_of_transit = lc_clean[~transit_mask]
lc_final = lc_clean.flatten(window_length=501, mask=transit_mask)

# Isolate phased/trimmed region around transit (± 0.25 days)
phase = ((lc_final.time.value - best_t0 + 0.5 * best_period) % best_period) - 0.5 * best_period
transit_window_mask = np.abs(phase) < 0.25

t_data = phase[transit_window_mask]
f_data = lc_final.flux.value[transit_window_mask]
ferr_data = lc_final.flux_err.value[transit_window_mask]

# Sort by phase/time
sort_idx = np.argsort(t_data)
t_data, f_data, ferr_data = t_data[sort_idx], f_data[sort_idx], ferr_data[sort_idx]
