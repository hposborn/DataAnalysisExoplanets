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

search_result = lk.search_lightcurve("Gliese 12", mission="TESS", author="SPOC")
print(search_result)

# %% Downloading

#We dont need 20s data (it is very heavy)
# NB - you may need to modify this for FFI-only lightcurves (where exptime could be 200, 600, 1200, 1800)
lc_collection = search_result[(search_result.exptime.value==120)&(search_result.author.value=="SPOC")].download_all()


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

# %% Cleaning & Normalising

lc_nanfree = lc_collection.stitch().remove_nans()
lc_clean = lc_nanfree.remove_outliers(sigma=5)
lc_norm = lc_clean.normalize()


# Create a flattened light curve
lc_flattened = lc_clean.copy()
lc_flattened.flux = lc_clean.flux / smoothed_flux #As we are normalised at 1.0, we can _divide_. At 0.0 we would need to _subtract_.

# 1. Process and flatten each sector individually
window_length = int((12*3600/120)+1) # Window length must be an odd integer
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
    lc_sec_binned = lc_sec.bin(time_bin_size=Quantity(1800, 's'))

    lc_sec.plot(ax=ax, ls='none', marker='.', ms=1.5, color='tab:blue', alpha=0.3, label='Flattened (120s)')
    lc_sec_binned.plot(ax=ax, ls='none', marker='o', ms=4, color='tab:orange', label='Binned (1800s)')

    ax.axhline(1.0, color='r', linestyle='--', alpha=0.5)
    ax.set_ylabel("Flux")
    ax.set_title(f"Gliese 12 - Sector {sec_num}")

plt.tight_layout()
plt.show()

# %%

lc_collection.exptime
# %%
