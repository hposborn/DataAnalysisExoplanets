# %% Initial imports
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord

# %% PART 1 - NEARBY STARS
# Using Astroquery to access gaia data (see https://astroquery.readthedocs.io/en/latest/gaia/gaia.html)

col_list="source_id, ra, dec, l, b, parallax, parallax_over_error, pmra, pmdec, radial_velocity, radial_velocity_error, phot_g_mean_mag, bp_rp, bp_g,g_rp, parallax, ruwe, teff_gspphot, logg_gspphot, mh_gspphot, ag_gspphot, distance_gspphot, rv_template_teff, rv_template_logg, rv_template_fe_h, vbroad"

#Let's try to get everything within 20parsec (e.g. parallax >50mas)
job = Gaia.launch_job(f"SELECT TOP 5000 {col_list} FROM gaiadr3.gaia_source WHERE parallax > 50")
result = job.get_results().to_pandas()
print(result)

# %% ADQL SEARCH WITH A DIFFERENT SERVER LOCATION:
import pyvo as vo

# 1. Define the TAP service URL for ARI Heidelberg
tap_url = "https://gaia.ari.uni-heidelberg.de/tap"

# 2. Create the TAP service instance
tap_service = vo.dal.TAPService(tap_url)

query = f"""
SELECT TOP 5000
    {col_list}
FROM gaiadr3.gaia_source
WHERE parallax > 20.0
  AND radial_velocity IS NOT NULL
  AND parallax_over_error > 10.0
"""

# 4. Run synchronously (or use submit_job for async)
result = tap_service.search(query).to_table().to_pandas()
print(result)

# %% Check what columns are available
print(result.columns)

# %% Plotting B-R colour vs G magnitude (and parallax as colour)

plt.scatter(result['bp_rp'],result['phot_g_mean_mag'],c=result['parallax'],s=5,alpha=0.65)
plt.xlim(-0.5,5.5)
plt.ylim(22,2.5)
plt.xlabel("B-R colour")
plt.ylabel("G magnitude")
plt.colorbar(label="parallax")
plt.show()



# %% Deriving luminosity from apparent brightness & distance

def derive_lum(mag,dist):
    absmag = mag - 5*np.log10(dist)+5
    #assuming absmag==bolometric mag...
    lum = 10**(0.4*(4.74 - absmag)) #Finish this function!
    return lum

result['lum']=derive_lum(result['phot_g_mean_mag'],1000/result['parallax'])

plt.scatter(result['bp_rp'],result['lum'],s=3+6*result['parallax']/50,c=result['parallax'],alpha=0.65)
plt.xlim(-0.5,5.5)
plt.yscale('log')
plt.ylim(np.min(result['lum'])*0.9,np.max(result['lum'])*1.1)
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
plt.xlabel("Stellar Radius")
plt.ylabel("Stellar luminosity")
plt.colorbar(label="parallax")
plt.show()


# %% PART 2) LOOKING AT _MOTION_ OF STARS NEAR OUR SUN

col_list="source_id, ra, dec, l, b, parallax, parallax_over_error, pmra, pmdec, radial_velocity, radial_velocity_error, phot_g_mean_mag, bp_rp, bp_g,g_rp, parallax, ruwe, teff_gspphot, logg_gspphot, mh_gspphot, ag_gspphot, distance_gspphot, rv_template_teff, rv_template_logg, rv_template_fe_h, vbroad"

# Let's create a selection function for sunlike stars!
fgk_query = f"""
SELECT TOP 100000 {col_list}
FROM gaiadr3.gaia_source
WHERE parallax > 2.0                      -- Distance <= 500 pc (parallax >= 2.0 mas)
  AND parallax_over_error > 10.0          -- Quality cut for reliable astrometry
  AND bp_rp BETWEEN 0.4 AND 1.8           -- Typical BP-RP colour range for FGK-type (sunlike) stars
  AND (phot_g_mean_mag + 5.0 * LOG10(parallax) - 10.0) BETWEEN 2.6 AND 8.5          -- Absolute magnitude range
  AND ruwe < 1.4                          -- Single star cut
  AND radial_velocity IS NOT NULL         -- Must have a measured radial velocity
  AND radial_velocity_error < 7.5         -- Clean radial velocity measurements
"""

fgk_job = Gaia.launch_job_async(fgk_query)
fgk_results = fgk_job.get_results()

# %% ALTERNATIVE:
# fgk_results = tap_service.search(fgk_query).to_table().to_pandas()

# Now getting real-world 3D space motion
# %% Build a SkyCoord object with all 6 phase-space dimensions
coords = SkyCoord(
    ra=fgk_results['ra'].values* u.deg,
    dec=fgk_results['dec'].values* u.deg,
    distance=(fgk_results['parallax'].values*u.mas).to(u.pc, equivalencies=u.parallax()),
    pm_ra_cosdec=fgk_results['pmra'].values* u.mas / u.yr,
    pm_dec=fgk_results['pmdec'].values* u.mas / u.yr,
    radial_velocity=fgk_results['radial_velocity'].values* u.km / u.s)

from astropy.coordinates import SkyCoord, Galactocentric, CylindricalRepresentation, CylindricalDifferential
# Transform to the Galactocentric frame - 3D Cartesian velocities (v_x, v_y, v_z) relative to the Galactic Center
galactocentric = coords.transform_to(Galactocentric)

# 1) Toomre Diagram
# Get Cartesian Galactocentric coordinates (positions)
x = galactocentric.x.to(u.kpc).value
y = galactocentric.y.to(u.kpc).value
z = galactocentric.z.to(u.kpc).value

# Get Cartesian Galactocentric velocities
v_x = galactocentric.v_x.to(u.km/u.s).value
v_y = galactocentric.v_y.value
v_z = galactocentric.v_z.value

# Calculate Galactocentric radius in the plane (R or rho)
R = np.sqrt(x**2 + y**2)

# Project Cartesian velocities into Cylindrical components cleanly:
# Radial velocity (v_R): motion outward from the Galactic center
v_r = (x * v_x + y * v_y) / R

# Azimuthal / Rotation velocity (v_phi): prograde rotation component
# (Standard convention: positive in the direction of Galactic rotation)
v_phi = (y * v_x - x * v_y) / R

# Total peculiar horizontal velocity for the Toomre diagram
v_tot_pec = np.sqrt(v_r**2 + v_z**2)

mh = fgk_results['mh_gspphot']

# Plotting the Toomre Diagram
plt.figure(figsize=(9, 6))
sc = plt.scatter(v_phi, v_tot_pec, c=mh, cmap='coolwarm', s=4, alpha=0.7, vmin=-1.0, vmax=0.5)
plt.colorbar(sc, label='Metallicity [M/H]')

# 2. Define the local circular velocity (LSR)
v_circ = 232.0  # Standard value (km/s)

# 3. Add the LSR reference line (vertical line at v_circ)
plt.axvline(v_circ, color='royalblue', linestyle='--', linewidth=1.5,
            label=f'LSR ($V_{{\\text{{circ}}}} = {v_circ}$ km/s)')

velocity_levels = [50, 100, 150, 200]
for v_val in velocity_levels:
    circle = plt.Circle((v_circ, 0), v_val, color='crimson', fill=False,
                        linestyle=':', linewidth=1, alpha=0.7)
    plt.gca().add_patch(circle)

    # Place text label neatly along the circle (at a 45-degree angle)
    label_x = v_circ + v_val * np.cos(np.pi / 4)
    label_y = v_val * np.sin(np.pi / 4)
    plt.text(label_x, label_y, f'{v_val} km/s', color='crimson', fontsize=8,
             rotation=-45, backgroundcolor='white', alpha=0.8)

plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.xlabel(r'Galactic Rotation Velocity $V_\phi$ (km/s)', fontsize=12)
plt.ylabel(r'Peculiar Velocity $\sqrt{U^2 + W^2}$ (km/s)', fontsize=12)
plt.title('Toomre Diagram of Solar Neighbourhood Stars', fontsize=14)
plt.xlim(-100, 400)
plt.ylim(0, 250)
plt.grid(True, linestyle=':', alpha=0.5)
plt.show()

# %% Plotting - 2) Overdensities in solar neighbourhood
from scipy import stats
vel = galactocentric.cartesian.differentials['s']

u_vel = -vel.d_x.to(u.km/u.s).value  # Positive toward Galactic Center
v_vel = vel.d_y.to(u.km/u.s).value-230   # Positive in direction of Galactic rotation, subtracted Solar speed
xy = np.vstack([u_vel, v_vel])
kde = stats.gaussian_kde(xy)
density = kde(xy) #This may be slow as it requires a large grid...

idx = density.argsort()
u_sorted, v_sorted, density_sorted = u_vel[idx], v_vel[idx], density[idx]

# %%
# PLOTTING:
plt.figure(figsize=(10, 8))
sc = plt.scatter(u_sorted, v_sorted, c=density_sorted, cmap='plasma', s=12, alpha=0.8)
plt.colorbar(sc, label='Local Density Estimation')

# Adding contour lines to isolate distinct overdensities
xmin, xmax = -80, 80
ymin, ymin_upper = -80, 80
# Create grid for contours
xx, yy = np.mgrid[-80:80:100j, -50:50:100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
z = np.reshape(kde(positions).T, xx.shape)
plt.contour(xx, yy, z, levels=9, colors='cyan', alpha=0.4, linewidths=1)

plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
plt.xlabel('U Velocity (km/s) [Toward Galactic Center]', fontsize=12)
plt.ylabel('V Velocity (km/s) [Toward Galactic Rotation]', fontsize=12)
plt.title('Solar Neighbourhood U-V Velocity Plane & Substructures', fontsize=14)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymin_upper)
plt.grid(True, linestyle=':', alpha=0.5)
plt.tight_layout()
plt.show()

# %% PART 3) LOOKING AT SPECIFIC STAR CLUSTER
from astropy.coordinates import SkyCoord
from astropy import units as u

hyades_location = SkyCoord("04h27m00s +15d50m00s", unit=(u.hourangle,u.deg))

cluster_ra = hyades_location.ra.deg
cluster_dec = hyades_location.dec.deg
cluster_size = 6.0#deg

cluster_query = f"""
SELECT TOP 5000 *
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {cluster_ra}, {cluster_dec}, {cluster_size}))
  AND parallax BETWEEN 18.0 and 28.0      -- Distance 45±10 pc
  AND parallax_over_error > 10.0          -- Quality cut for reliable astrometry
  AND pmra BETWEEN 90.0 AND 115.0         -- Expect roughly similar proper motion
  AND pmdec BETWEEN -30.0 AND -10.0"""

cluster_job = Gaia.launch_job(cluster_query)
cluster_results = cluster_job.get_results()

# %%

# Try plotting the MR diagram (e.g. B-R vs absolute mag)

# %%

# YOUR TURN!
