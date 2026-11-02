# %%
#######################
# GIANT PLANET INTERIORS: BULK METALLICITY
#######################
"""
How much heavy material (anything other than H/He) is inside a giant planet?

For giants we cannot separate iron, rock and water: we infer the bulk heavy-element
mass fraction Z. This needs evolution models, because giants shrink as they cool,
so the age and the stellar irradiation matter as well as mass and radius.

Planets: HD 148797 b and c (low density) and TOI-4504 c (dense)

Steps
    Load masses, radii, insolations and ages from the NASA Exoplanet Archive
    Check each planet is inside the models' validity range
    Mass-radius curves from planetsynth (Müller & Helled 2021) for several Z
    Simple inversion: Z at the planet's mass, age and irradiation
    Posterior on Z with emcee + planetsynth
    (Optional) GASTLI (Acuña et al. 2024): coupled interior-atmosphere retrieval

Setup (uv)
    git clone https://github.com/tiny-hippo/planetsynth.git
    uv add numpy scipy pandas matplotlib emcee corner
    uv add ./planetsynth
  Optional GASTLI (needs a Fortran compiler, e.g. 'brew install gcc' or 'apt install gfortran'):
    uv add gastli h5py
    Download the GASTLI input data (link in the GASTLI installation docs:
    https://gastli.readthedocs.io), unzip it, and set the data path as the docs describe.

Not suitable: hot Jupiters (Teq above ~1000 K), whose radii are inflated by processes
these models do not include, and which are outside GASTLI's atmosphere grid.
"""

# %%
#######################
# IMPORTS AND SETTINGS
#######################
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
rng = np.random.default_rng(42)

PLANET_NAMES = ["HD 148797 b", "HD 148797 c", "TOI-4504 c"]
OUTDIR = Path("interiors_outputs") / "giants"
OUTDIR.mkdir(parents=True, exist_ok=True)

S_EARTH = 1.361e6          # solar constant [erg/s/cm2]
M_EARTH_PER_MJ = 317.83
RUN_EMCEE = True
RUN_GASTLI = False

# Values to use if the archive query fails or a value is missing. Fill these in from
# the discovery papers. Units: Jupiter masses and radii, K, Earth insolation, Gyr.
MANUAL = {
    # "HD 148797 b": dict(M=..., M_err=..., R=..., R_err=..., Teq=..., S=..., age=..., age_err=...),
}

# %%
#######################
# PLANET PARAMETERS FROM THE NASA EXOPLANET ARCHIVE
#######################
names_sql = ", ".join(f"'{n}'" for n in PLANET_NAMES)
query = ("select pl_name, pl_bmassj, pl_bmassjerr1, pl_bmassjerr2, pl_radj, pl_radjerr1, "
         "pl_radjerr2, pl_eqt, pl_insol, st_age, st_ageerr1, st_ageerr2 from pscomppars "
         f"where pl_name in ({names_sql})")
url = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
       + urllib.parse.quote(query) + "&format=csv")

try:
    raw = pd.read_csv(url)
    planets = pd.DataFrame({
        "name": raw["pl_name"],
        "M": raw["pl_bmassj"], "M_err": (raw["pl_bmassjerr1"] - raw["pl_bmassjerr2"]) / 2,
        "R": raw["pl_radj"], "R_err": (raw["pl_radjerr1"] - raw["pl_radjerr2"]) / 2,
        "Teq": raw["pl_eqt"], "S": raw["pl_insol"],
        "age": raw["st_age"], "age_err": (raw["st_ageerr1"] - raw["st_ageerr2"]) / 2,
    })
    print("Loaded from the NASA Exoplanet Archive")
except Exception as exc:
    print(f"Archive query failed ({exc}); using the MANUAL table")
    planets = pd.DataFrame(columns=["name", "M", "M_err", "R", "R_err", "Teq", "S", "age", "age_err"])

# Add planets missing from the archive and fill missing values from MANUAL
for name in PLANET_NAMES:
    if name not in planets["name"].values:
        planets.loc[len(planets), "name"] = name
for name, values in MANUAL.items():
    for key, v in values.items():
        planets.loc[planets["name"] == name, key] = planets.loc[planets["name"] == name, key].fillna(v)
planets = planets.astype({c: float for c in planets.columns if c != "name"}).reset_index(drop=True)

# Derived quantities and sensible defaults
planets["S"] = planets["S"].fillna((planets["Teq"] / 278.6) ** 4)            # zero-albedo estimate
planets["Teq"] = planets["Teq"].fillna(278.6 * planets["S"] ** 0.25)
planets["logF"] = np.log10(planets["S"] * S_EARTH)                          # planetsynth's flux unit
if planets["age"].isna().any():
    print("Missing ages set to 5 +- 3 Gyr: check the literature!")
planets["age"] = planets["age"].fillna(5.0)
planets["age_err"] = planets["age_err"].fillna(3.0)
planets["short"] = planets["name"]
print(planets.round(3).to_string(index=False))

missing = planets[planets[["M", "R"]].isna().any(axis=1)]["name"].tolist()
if missing:
    raise ValueError(f"No mass or radius for {missing}: add them to MANUAL above.")

# %%
#######################
# ARE THE PLANETS INSIDE THE MODEL GRIDS?
#######################
# planetsynth: 0.1-30 Jupiter masses, log10 flux 1-9 (erg/s/cm2), ages 0.01-10 Gyr, and a
# maximum Z that depends on mass. GASTLI's default atmospheres: Teq 100-1000 K.


def z_max(M):
    """Largest bulk Z in the planetsynth grid for mass M [M_jup]."""
    return np.select([M <= 1, M <= 3, M <= 5], [0.8, 0.5, 0.2], 0.04)


for _, p in planets.iterrows():
    notes = []
    if not 0.1 <= p["M"] <= 30:
        notes.append("mass outside planetsynth range (0.1-30 M_J)")
    if not 1 <= p["logF"] <= 9:
        notes.append(f"log F = {p['logF']:.2f} outside planetsynth range (1-9)")
    if p["Teq"] > 1000:
        notes.append(f"Teq = {p['Teq']:.0f} K: above GASTLI's grid; radius may be inflated")
    status = "OK" if not notes else "; ".join(notes)
    print(f"{p['name']:14s} M = {p['M']:.3f} M_J ({p['M'] * M_EARTH_PER_MJ:.0f} M_E), "
          f"R = {p['R']:.3f} R_J, Teq = {p['Teq']:.0f} K, Zmax = {z_max(p['M']):.2f}: {status}")

# %%
#######################
# LOAD PLANETSYNTH
#######################
# planetsynth stores its model grid as a pickled SciPy interpolator. Pickles from older
# SciPy versions do not work with current SciPy, so we rebuild the interpolator from the
# stored grid and values (same data, same method).
from scipy.interpolate import RegularGridInterpolator
from planetsynth import PlanetSynth

pls = PlanetSynth(verbose=False)
_stored = pls.reg.__dict__
pls.reg = RegularGridInterpolator(_stored["grid"], _stored["values"], method=_stored["method"],
                                  bounds_error=_stored["bounds_error"], fill_value=_stored["fill_value"])


def radius_model(M, Z, logF, age_gyr, Z_atm=0.02):
    """
    planetsynth radius [R_jup] for arrays of mass [M_jup], bulk Z, log10 flux and age [Gyr].
    The cooling track is computed on planetsynth's time grid and interpolated (linearly in
    log age) to each planet's age. Returns NaN outside the model range.
    """
    M, Z, logF, age_gyr = np.broadcast_arrays(*map(np.atleast_1d, (M, Z, logF, age_gyr)))
    params = np.column_stack([M, Z, np.minimum(Z_atm, Z), logF]).astype(float)
    tracks = pls.synthesize(params)[:, :, 0]                       # (n, n_times) radius
    logt = np.log10(age_gyr * 1e9)
    x = pls.log_time
    j = np.clip(np.searchsorted(x, logt) - 1, 0, len(x) - 2)
    w = (logt - x[j]) / (x[j + 1] - x[j])
    R = tracks[np.arange(len(M)), j] * (1 - w) + tracks[np.arange(len(M)), j + 1] * w
    R[(logt < x[0]) | (logt > x[-1])] = np.nan
    return R


# %%
#######################
# MASS-RADIUS DIAGRAM
#######################
# Each planet gets its own panel, because the curves depend on its age and irradiation.
mass_grid = np.logspace(-1, np.log10(10), 120)
z_values = [0.0, 0.05, 0.1, 0.2, 0.4, 0.6]
fig, axes = plt.subplots(1, len(planets), figsize=(5 * len(planets), 4.5), squeeze=False)
for ax, (_, p) in zip(axes[0], planets.iterrows()):
    for Z, color in zip(z_values, plt.cm.viridis(np.linspace(0, 0.9, len(z_values)))):
        Mg = mass_grid[z_max(mass_grid) >= Z]
        ax.plot(Mg, radius_model(Mg, Z, min(p["logF"], 9), p["age"]), color=color, label=f"Z = {Z}")
    ax.errorbar(p["M"], p["R"], xerr=p["M_err"], yerr=p["R_err"], fmt="o", color="k", ms=6, zorder=5)
    ax.plot([1.0, 0.299], [1.0, 0.843], "s", color="tab:red", ms=5)        # Jupiter, Saturn
    ax.annotate("J", (1.0, 1.0), xytext=(4, 4), textcoords="offset points", color="tab:red")
    ax.annotate("S", (0.299, 0.843), xytext=(4, 4), textcoords="offset points", color="tab:red")
    ax.set_xscale("log")
    ax.set_xlabel(r"Mass [M$_J$]")
    ax.set_title(f"{p['name']}: {p['age']:.1f} Gyr, Teq {p['Teq']:.0f} K", fontsize=10)
axes[0, 0].set_ylabel(r"Radius [R$_J$]")
axes[0, 0].legend(fontsize=7)
fig.tight_layout()
plt.show()
# Jupiter and Saturn are shown for orientation only: they are older and far less irradiated.

# %%
#######################
# SIMPLE INVERSION: BULK METALLICITY Z
#######################
# At fixed mass, age and irradiation, more heavy elements make a smaller planet. For many
# random draws of mass and radius, compute R(Z) and read off the Z matching the radius.
n_draws = 1500
simple = []
for _, p in planets.iterrows():
    z_grid = np.linspace(0, z_max(p["M"]), 81)
    M_draw = np.clip(rng.normal(p["M"], p["M_err"], n_draws), 0.1, 30)
    R_draw = rng.normal(p["R"], p["R_err"], n_draws)
    R_model = radius_model(np.repeat(M_draw, z_grid.size), np.tile(z_grid, n_draws),
                           min(p["logF"], 9), p["age"]).reshape(n_draws, z_grid.size)
    too_big = R_draw > R_model[:, 0]        # larger than pure H/He: inflated or very young
    too_small = R_draw < R_model[:, -1]     # denser than the grid allows
    ok = ~too_big & ~too_small & np.isfinite(R_model).all(axis=1)
    Z = np.array([np.interp(R_draw[i], R_model[i, ::-1], z_grid[::-1]) for i in np.where(ok)[0]])
    q = np.percentile(Z, [16, 50, 84]) if len(Z) else [np.nan] * 3
    simple.append(dict(planet=p["name"], Z=q[1], minus=q[1] - q[0], plus=q[2] - q[1],
                       M_Z_earth=q[1] * p["M"] * M_EARTH_PER_MJ,
                       frac_too_large=too_big.mean(), frac_too_dense=too_small.mean()))
simple = pd.DataFrame(simple)
print(simple.round(3).to_string(index=False))
print("\nM_Z_earth: heavy-element mass in Earth masses. For comparison, Jupiter's Z is ~0.05-0.1.")
print("Z is computed only from the draws inside the grid. If 'frac_too_large' is high, most")
print("draws are bigger than even a pure H/He planet: Z is then close to zero, or the planet is")
print("inflated by something the models leave out (e.g. extra heating, or a younger age).")

# %%
#######################
# POSTERIOR ON Z WITH EMCEE + PLANETSYNTH
#######################
# Parameters: Z, mass, age and log flux. Priors: uniform Z (0 to the grid maximum),
# Gaussian mass, age and flux (flux error 0.05 dex). Likelihood: Gaussian on the radius.
# The log-probability is vectorised: all walkers are evaluated in one planetsynth call.
import emcee

z_samples = {}
if RUN_EMCEE:
    for _, p in planets.iterrows():
        logF0 = min(p["logF"], 9)

        def log_prob(theta, p=p, logF0=logF0):
            Z, M, age, logF = theta.T
            lp = (-0.5 * ((M - p["M"]) / p["M_err"]) ** 2
                  - 0.5 * ((age - p["age"]) / p["age_err"]) ** 2
                  - 0.5 * ((logF - logF0) / 0.05) ** 2)
            inside = (Z >= 0) & (Z <= z_max(M)) & (M >= 0.1) & (M <= 30) & \
                     (age > 0.05) & (age < 10) & (logF >= 1) & (logF <= 9)
            out = np.full(len(Z), -np.inf)
            if inside.any():
                R = radius_model(M[inside], Z[inside], logF[inside], age[inside])
                ll = -0.5 * ((R - p["R"]) / p["R_err"]) ** 2
                out[inside] = np.where(np.isfinite(ll), lp[inside] + ll, -np.inf)
            return out

        nwalkers = 48
        start_Z = simple.loc[simple["planet"] == p["name"], "Z"].fillna(0.1).values[0]
        p0 = np.column_stack([
            np.clip(start_Z + 0.02 * rng.normal(size=nwalkers), 0.001, z_max(p["M"]) - 0.001),
            p["M"] + 0.1 * p["M_err"] * rng.normal(size=nwalkers),
            np.clip(p["age"] + 0.1 * p["age_err"] * rng.normal(size=nwalkers), 0.1, 9.9),
            logF0 - 0.01 * np.abs(rng.normal(size=nwalkers)),
        ])
        sampler = emcee.EnsembleSampler(nwalkers, 4, log_prob, vectorize=True)
        sampler.run_mcmc(p0, 4000, progress=True)
        z_samples[p["name"]] = sampler.get_chain(discard=1500, flat=True)

    rows = []
    for name, s in z_samples.items():
        Z, M = s[:, 0], s[:, 1]
        MZ = Z * M * M_EARTH_PER_MJ
        rows.append(dict(planet=name, Z_median=np.median(Z), Z_16=np.percentile(Z, 16),
                         Z_84=np.percentile(Z, 84), M_Z_earth=np.median(MZ),
                         M_Z_16=np.percentile(MZ, 16), M_Z_84=np.percentile(MZ, 84)))
    print(pd.DataFrame(rows).round(3).to_string(index=False))

    fig, ax = plt.subplots(figsize=(7, 4))
    for name, s in z_samples.items():
        ax.hist(s[:, 0], bins=60, histtype="step", density=True, label=name)
    ax.set_xlabel("Bulk heavy-element mass fraction Z")
    ax.set_ylabel("Posterior density")
    ax.legend()
    plt.show()

# %%
#######################
# GASTLI: GRID OF COUPLED INTERIOR-ATMOSPHERE EVOLUTION MODELS (OPTIONAL)
#######################
# GASTLI computes thermal evolution with a real atmosphere model. Following the GASTLI
# retrieval tutorial, we first compute a small grid of evolution sequences for each planet
# (its Teq; a range of core mass fractions and masses; solar envelope metallicity), save it,
# and then interpolate it in an MCMC. Each sequence takes a few seconds to a minute, so a
# grid takes ~10-30 min per planet. It is saved and reused on later runs.
# In GASTLI 'CMF' is the core mass fraction; the bulk metal fraction is
# Z = CMF + (1 - CMF) * Z_envelope.
CMF_GRID = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
TINT_GRID = np.array([50.0, 100.0, 200.0, 300.0, 400.0])

gastli_grids = {}
if RUN_GASTLI:
    import gastli.Thermal_evolution as therm

    for _, p in planets.iterrows():
        if p["Teq"] > 1000:
            print(f"{p['name']}: Teq above GASTLI's atmosphere grid, skipped")
            continue
        grid_file = OUTDIR / f"gastli_grid_{p['name'].replace(' ', '_')}.npz"
        M_E = p["M"] * M_EARTH_PER_MJ
        mass_grid_E = M_E + np.array([-3, -1.5, 0, 1.5, 3]) * p["M_err"] * M_EARTH_PER_MJ
        if grid_file.exists():
            gastli_grids[p["name"]] = dict(np.load(grid_file))
            continue
        shape = (len(CMF_GRID), len(mass_grid_E), len(TINT_GRID))
        out = {k: np.full(shape, np.nan) for k in ["Rtot", "Mtot", "age", "Zenv"]}
        for i, cmf in enumerate(CMF_GRID):
            for j, m in enumerate(mass_grid_E):
                print(f"{p['name']}: CMF = {cmf:.2f}, M = {m:.0f} M_E")
                try:
                    te = therm.thermal_evolution(pow_law_formass=0.315)
                    te.main(m, cmf, p["Teq"], TINT_GRID, CO=0.55, log_FeH=0.0)
                    te.solve_thermal_evol_eq()
                    out["Rtot"][i, j] = te.Rtot_TE          # Jupiter radii
                    out["Mtot"][i, j] = te.Mtot_TE          # Earth masses
                    out["age"][i, j] = te.age_points        # Gyr
                    out["Zenv"][i, j] = te.Zenv_TE
                except Exception as exc:
                    print("   model failed:", exc)
        np.savez(grid_file, CMF=CMF_GRID, mass=mass_grid_E, Tint=TINT_GRID, **out)
        gastli_grids[p["name"]] = dict(np.load(grid_file))

# %%
#######################
# GASTLI: RETRIEVAL OF THE BULK METALLICITY (OPTIONAL)
#######################
# Parameters: CMF, input mass and internal temperature (Tint sets the age). Likelihood:
# Gaussian on the total mass, radius and age.
if RUN_GASTLI:
    for _, p in planets.iterrows():
        if p["name"] not in gastli_grids:
            continue
        g = gastli_grids[p["name"]]
        axes_ = (g["CMF"], g["mass"], g["Tint"])
        interp = {k: RegularGridInterpolator(axes_, g[k], bounds_error=False, fill_value=np.nan)
                  for k in ["Rtot", "Mtot", "age", "Zenv"]}
        M_E, M_E_err = p["M"] * M_EARTH_PER_MJ, p["M_err"] * M_EARTH_PER_MJ

        def log_prob_gastli(theta, p=p, interp=interp, M_E=M_E, M_E_err=M_E_err):
            pt = np.atleast_2d(theta)
            R, M, age = (interp[k](pt)[0] for k in ["Rtot", "Mtot", "age"])
            if not np.all(np.isfinite([R, M, age])):
                return -np.inf
            return (-0.5 * ((R - p["R"]) / p["R_err"]) ** 2 - 0.5 * ((M - M_E) / M_E_err) ** 2
                    - 0.5 * ((age - p["age"]) / p["age_err"]) ** 2)

        nwalkers = 32
        p0 = np.column_stack([rng.uniform(0.05, 0.3, nwalkers),
                              M_E + 0.3 * M_E_err * rng.normal(size=nwalkers),
                              rng.uniform(80, 250, nwalkers)])
        sampler = emcee.EnsembleSampler(nwalkers, 3, log_prob_gastli)
        sampler.run_mcmc(p0, 4000, progress=True)
        s = sampler.get_chain(discard=1500, flat=True)
        Zenv = interp["Zenv"](s)
        Z = s[:, 0] + (1 - s[:, 0]) * Zenv
        print(f"{p['name']}: GASTLI Z = {np.nanmedian(Z):.3f} "
              f"(+{np.nanpercentile(Z, 84) - np.nanmedian(Z):.3f} -{np.nanmedian(Z) - np.nanpercentile(Z, 16):.3f})")

# %%
#######################
# QUESTIONS FOR STUDENTS
#######################
"""
How much heavy material does each planet contain, in Earth masses?
Why do the low-density planets need so little heavy material, and the dense one so much?
How sensitive is Z to the age? Try changing the age prior.
Compare planetsynth and GASTLI: what physics differs, and how much does Z change?
"""
