# %% Initial imports
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord

# %% PART 1 - NEARBY STARS
# Using Astroquery to access gaia data (see https://astroquery.readthedocs.io/en/latest/gaia/gaia.html)

#Let's try to get everything within 20parsec (e.g. parallax >50mas)
job = Gaia.launch_job("select top 5000 "
                      "* from gaiadr3.gaia_source_lite "
                      "WHERE parallax > 50")
results = job.get_results()
print(results)

# %% Check what columns are available
print(results.columns)

# %% Plotting colour vs magnitude

plt.scatter(results['bp_rp'],results['phot_g_mean_mag'],c=results['parallax'],s=5,alpha=0.65)
plt.xlim(-0.5,5.5)
plt.ylim(22,2.5)
plt.xlabel("B-R colour")
plt.ylabel("G magnitude")
plt.colorbar(label="parallax")
plt.show()

# %% Deriving luminosity from apparent brightness & distance

def derive_lum():
    lum = #Finish this function!
    return lum

results['lum'] = derive_lum(results[''],results[''])

plt.scatter(results['bp_rp'],results['lum'],s=3+6*results['parallax']/50,c=results['parallax'],alpha=0.65)
plt.xlim(-0.5,5.5)
plt.yscale('log')
plt.ylim(np.min(results['lum'])*0.9,np.max(results['lum'])*1.1)
plt.xlabel("B-R colour")
plt.ylabel("Stellar luminosity")
plt.colorbar(label="parallax")
plt.show()

# %%
# Deriving stellar radius from both

def derive_rad():
    rad = #Finish this function!
    return rad

results['rad'] = derive_rad(results[''],results[''])

plt.scatter(results['rad'],results['lum'],s=3+6*results['parallax']/50,c=results['parallax'],alpha=0.65)
plt.xlim(np.min(results['rad'])*0.9,np.max(results['rad'])*1.1)
plt.yscale('log')
plt.ylim(np.min(results['lum'])*0.9,np.max(results['lum'])*1.1)
plt.xlabel("B-R colour")
plt.ylabel("Stellar luminosity")
plt.colorbar(label="parallax")
plt.show()


# %% Estimating _stellar temperatures_ from Gaia B & R, and 2MASS colours

# Step 1 - accessing cross-matched 2MASS & Gaia catalogue (T)

# %%

# Deriving luminosity from apparent brightness & distance

# Deriving stellar radius from both

# Deriving stellar mass from logg

# Using the astroquery API (and ADQL queries)

# Finding stars which cluster (young)
