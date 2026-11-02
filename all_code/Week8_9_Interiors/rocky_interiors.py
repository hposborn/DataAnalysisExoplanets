# %%
#######################
# ROCKY PLANET INTERIORS: THE TRAPPIST-1 SYSTEM
#######################
"""
What are the seven TRAPPIST-1 planets made of?

Steps
    Load masses, radii and temperatures from the NASA Exoplanet Archive
    Compare them with mass-radius curves for rocky and water-rich compositions
    Simple inversion: core mass fraction (CMF) from a grid of rocky models
    Posterior sampling with smint: CMF of a dry rocky planet, and the minimum water fraction
    (Optional) ExoMDN: core, mantle, water and gas mass fractions in one go

Setup (uv)
    git clone https://github.com/cpiaulet/smint.git all_code/ (code + model grids)
    uv add numpy scipy pandas matplotlib emcee corner astropy
    uv add all_code/smint
  Optional ExoMDN (needs TensorFlow 2.11, so Python 3.10):
    git clone https://github.com/philippbaumeister/ExoMDN.git
    uv python pin 3.10
    uv add "tensorflow==2.11.*" "tensorflow-probability==0.15.*" "scikit-learn==1.1.1" "numpy==1.26.*" seaborn joblib
    uv add ./ExoMDN
    (on Apple-silicon Macs TensorFlow 2.11 may need "tensorflow-macos==2.11.*" instead)

Key idea for students
    Mass and radius give the mean density. Many different mixtures of iron, rock and
    water share the same density, so the answer is a probability distribution of
    compositions, not a single composition.
"""

# %%
#######################
# IMPORTS AND SETTINGS
#######################
import os
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
rng = np.random.default_rng(42)

HOST = "TRAPPIST-1"
SMINT_DIR = Path(os.environ.get("SMINT_DIR", "smint"))          # the cloned smint repository
MODELS = str(SMINT_DIR / "example" / "smint_models") + "/"      # smint expects a trailing slash
EXOMDN_DIR = Path(os.environ.get("EXOMDN_DIR", "ExoMDN"))       # the cloned ExoMDN repository
OUTDIR = Path("interiors_outputs") / "rocky"
OUTDIR.mkdir(parents=True, exist_ok=True)

RUN_SMINT_CMF = True       # dry rocky planets: iron + rock
RUN_SMINT_WATER = True     # iron-free rock + water
RUN_EXOMDN = False         # needs the ExoMDN environment described above

# %%
#######################
# PLANET PARAMETERS FROM THE NASA EXOPLANET ARCHIVE
#######################
# 'pscomppars' gives one row per planet, combining the most complete set of parameters.
# Values may come from different papers; check the references for serious work.
query = ("select pl_name, pl_bmasse, pl_bmasseerr1, pl_bmasseerr2, pl_rade, pl_radeerr1, "
         "pl_radeerr2, pl_eqt, pl_insol, st_age from pscomppars "
         f"where hostname = '{HOST}' order by pl_orbper")
url = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
       + urllib.parse.quote(query) + "&format=csv")

try:
    raw = pd.read_csv(url)
    planets = pd.DataFrame({
        "name": raw["pl_name"],
        "M": raw["pl_bmasse"], "M_err": (raw["pl_bmasseerr1"] - raw["pl_bmasseerr2"]) / 2,
        "R": raw["pl_rade"], "R_err": (raw["pl_radeerr1"] - raw["pl_radeerr2"]) / 2,
        "Teq": raw["pl_eqt"], "S": raw["pl_insol"],
    })
    print("Loaded from the NASA Exoplanet Archive")
except Exception as exc:
    # Fallback: Agol et al. (2021). Check against the archive before relying on these.
    print(f"Archive query failed ({exc}); using built-in values (Agol et al. 2021)")
    planets = pd.DataFrame({
        "name": [f"TRAPPIST-1 {p}" for p in "bcdefgh"],
        "M": [1.374, 1.308, 0.388, 0.692, 1.039, 1.321, 0.326],
        "M_err": [0.069, 0.056, 0.012, 0.022, 0.031, 0.038, 0.020],
        "R": [1.116, 1.097, 0.788, 0.920, 1.045, 1.129, 0.755],
        "R_err": [0.013, 0.013, 0.011, 0.013, 0.013, 0.014, 0.014],
        "Teq": [398.0, 340.0, 286.0, 250.0, 218.0, 197.0, 172.0],
        "S": [np.nan] * 7,
    })

# Insolation (Earth = 1) from Teq if missing: Teq = 278.6 K * S**0.25 (zero albedo)
planets["S"] = planets["S"].fillna((planets["Teq"] / 278.6) ** 4)
planets["short"] = planets["name"].str.split().str[-1]
print(planets.round(3).to_string(index=False))

# %%
#######################
# MODEL GRIDS FROM SMINT
#######################
# smint ships the rocky models of Zeng et al. (2016) merged into the Aguichine et al. (2021)
# grid. Axes: f'core (iron fraction of the iron + rock part), irradiation temperature,
# water mass fraction and log10(mass). With no water, f'core is the core mass fraction (CMF).
import os
os.environ['SMINT_DIR']="../smint"
from smint import fit_cmf

interp_r = fit_cmf.make_interpolator_A21(MODELS + "Aguichine2021/withZengRocky/")
interp_valid = fit_cmf.make_interpolator_A21(MODELS + "Aguichine2021/", which_quantity="validity")

# Zeng et al. (2016) iron-free rock + condensed water models (water mass fraction in %)
zeng_water = pd.read_csv(MODELS + "t_rock_h2o_Zeng2016.csv", sep=r"\s+")

# %%
#######################
# MASS-RADIUS DIAGRAM
#######################
mass_grid = np.logspace(np.log10(0.2), np.log10(20), 200)
fig, ax = plt.subplots(figsize=(8, 6))

for cmf, ls, label in [(0.9, "-", "90% iron"), (0.33, "-", "Earth-like (33% iron)"),
                       (0.0, "-", "pure rock (no iron)")]:
    pts = np.column_stack([np.full_like(mass_grid, cmf), np.full_like(mass_grid, 400.0),
                           np.zeros_like(mass_grid), np.log10(mass_grid)])
    ax.plot(mass_grid, interp_r(pts), ls, lw=2, label=label)

for wf, color in [("20%h2o", "tab:cyan"), ("50%h2o", "tab:blue"), ("100%h2o", "navy")]:
    ax.plot(zeng_water["Mass"], zeng_water[wf], "--", color=color,
            label=f"rock (no iron) + {wf.replace('h2o', '')} water")

ax.errorbar(planets["M"], planets["R"], xerr=planets["M_err"], yerr=planets["R_err"],
            fmt="o", color="k", ms=5, zorder=5)
for _, p in planets.iterrows():
    ax.annotate(p["short"], (p["M"], p["R"]), xytext=(5, 5), textcoords="offset points")
for name, m, r in [("Earth", 1.0, 1.0), ("Venus", 0.815, 0.949), ("Mars", 0.107, 0.532)]:
    ax.plot(m, r, "s", color="tab:green", ms=6)
    ax.annotate(name, (m, r), xytext=(5, -10), textcoords="offset points", color="tab:green")

ax.set_xscale("log")
ax.set_xlim(0.1, 3)
ax.set_ylim(0.5, 1.5)
ax.set_xlabel(r"Mass [M$_\oplus$]")
ax.set_ylabel(r"Radius [R$_\oplus$]")
ax.set_title("TRAPPIST-1 planets and rocky / water-rich models")
ax.legend(fontsize=8, loc="upper left")
plt.show()

# %%
#######################
# SIMPLE INVERSION: CORE MASS FRACTION
#######################
# For a dry rocky planet, radius at fixed mass shrinks as the iron fraction grows. So for
# each planet: compute R(CMF) at its mass, and read off the CMF that gives the observed R.
# Repeating this for many random draws of mass and radius (within their errors) gives an
# uncertainty. Draws outside the grid are counted separately:
#   too large for pure rock -> the planet needs something lighter (water, atmosphere)
#   too small for 90% iron  -> denser than the models allow
cmf_grid = np.linspace(0.0, 0.9, 91)
n_draws = 4000
simple = []
for _, p in planets.iterrows():
    M_draw = rng.normal(p["M"], p["M_err"], n_draws)
    R_draw = rng.normal(p["R"], p["R_err"], n_draws)
    pts = np.column_stack([np.tile(cmf_grid, n_draws), np.full(n_draws * cmf_grid.size, 400.0),
                           np.zeros(n_draws * cmf_grid.size), np.repeat(np.log10(M_draw), cmf_grid.size)])
    R_model = interp_r(pts).reshape(n_draws, cmf_grid.size)        # radius vs CMF, per draw

    too_big = R_draw > R_model[:, 0]       # larger than pure rock
    too_small = R_draw < R_model[:, -1]    # smaller than 90% iron
    ok = ~too_big & ~too_small
    cmf = np.array([np.interp(R_draw[i], R_model[i, ::-1], cmf_grid[::-1]) for i in np.where(ok)[0]])
    q = np.percentile(cmf, [16, 50, 84]) if len(cmf) else [np.nan] * 3
    simple.append(dict(planet=p["short"], CMF=q[1], minus=q[1] - q[0], plus=q[2] - q[1],
                       frac_too_large=too_big.mean(), frac_too_dense=too_small.mean()))
simple = pd.DataFrame(simple)
print(simple.round(3).to_string(index=False))
print("\nEarth's CMF is about 0.33. 'frac_too_large' > 0 means some of the allowed (M, R)")
print("cannot be explained by iron and rock alone.")

# %%
#######################
# SMINT: POSTERIOR ON THE CORE MASS FRACTION (DRY ROCKY PLANETS)
#######################
# smint samples (CMF, mass) with emcee: a Gaussian prior on the mass, a uniform prior on
# CMF between 0 and 0.9, and a Gaussian likelihood on the radius.
cmf_samples = {}
if RUN_SMINT_CMF:
    for _, p in planets.iterrows():
        params = dict(Mp_earth=p["M"], err_Mp_earth=p["M_err"], Rp_earth=p["R"], err_Rp_earth=p["R_err"],
                      use_KDE_for_Mp_prior=False, nwalkers=40, ndim=2, nsteps=2500, frac_burnin=0.5,
                      save=False, outputdir=str(OUTDIR / "smint"), fname="cmf_" + p["short"])
        params = fit_cmf.setup_priors(params)
        params = fit_cmf.ini_fit(params)
        sampler = fit_cmf.run_fit(params, interp_r, interp_valid)
        chain = sampler.chain[:, int(params["frac_burnin"] * params["nsteps"]):, :]
        cmf_samples[p["short"]] = chain.reshape(-1, 2)[:, 0]

# %%
#######################
# COMPARE THE SIMPLE INVERSION AND THE SMINT POSTERIORS
#######################
if cmf_samples:
    fig, ax = plt.subplots(figsize=(9, 4))
    names = list(cmf_samples)
    ax.violinplot([cmf_samples[n] for n in names], positions=range(len(names)), showmedians=True)
    ax.errorbar(np.arange(len(names)) + 0.25, simple["CMF"], yerr=[simple["minus"], simple["plus"]],
                fmt="s", color="k", label="simple inversion")
    ax.axhline(0.33, color="tab:green", ls="--", label="Earth")
    ax.set_xticks(range(len(names)), names)
    ax.set_ylabel("Core mass fraction")
    ax.set_title("TRAPPIST-1: CMF posteriors (smint) vs simple inversion")
    ax.legend()
    plt.show()
    # Note: smint's prior stops at CMF = 0 and 0.9. A posterior piled up against CMF = 0
    # means the planet is too large for a dry rock: the data want a lighter component.

# %%
#######################
# SMINT: ROCK + WATER
#######################
# The other way to lower a planet's density is to add water. smint's Zeng et al. (2016)
# rock + water models use iron-free silicate rock and condensed (cold) water, which suits
# the temperate TRAPPIST-1 planets. (smint's Aguichine et al. models instead assume a hot,
# irradiated steam/supercritical water layer, valid above ~400 K.)
# Because the rock here contains no iron, this fit gives the MINIMUM water needed:
# an Earth-like iron/rock ratio would need more water to reach the same density.
# Parameters: water mass fraction (in %, uniform prior 0-100) and mass.
from smint import fit_fh2o

interp_h2o = fit_fh2o.make_interpolator_fh2o(zeng_water)
water_samples = {}
if RUN_SMINT_WATER:
    for _, p in planets.iterrows():
        params = dict(Mp_earth=p["M"], err_Mp_earth=p["M_err"], Rp_earth=p["R"], err_Rp_earth=p["R_err"],
                      nwalkers=40, ndim=2, nsteps=2500, frac_burnin=0.5, save=False)
        params["pos0"] = [np.array([abs(rng.normal(10, 5)), rng.normal(p["M"], p["M_err"])])
                          for _ in range(params["nwalkers"])]
        sampler = fit_fh2o.run_fit(params, interp_h2o)
        chain = sampler.chain[:, int(params["frac_burnin"] * params["nsteps"]):, :]
        water_samples[p["short"]] = chain.reshape(-1, 2)[:, 0] / 100     # % -> fraction

    rows = [dict(planet=n, water_median=np.median(w), water_16=np.percentile(w, 16),
                 water_84=np.percentile(w, 84)) for n, w in water_samples.items()]
    print(pd.DataFrame(rows).round(3).to_string(index=False))

    fig, ax = plt.subplots(figsize=(9, 4))
    names = list(water_samples)
    ax.violinplot([water_samples[n] for n in names], positions=range(len(names)), showmedians=True)
    ax.set_xticks(range(len(names)), names)
    ax.set_ylabel("Water mass fraction")
    ax.set_title("TRAPPIST-1: minimum water (iron-free rock + condensed water)")
    plt.show()
    # A water fraction near zero means iron-poor dry rock already explains the planet.
    # The truth could be anywhere between 'less iron than Earth, no water' and 'Earth-like
    # iron/rock plus some water': mass and radius alone cannot tell these apart.

# %%
#######################
# EXOMDN: CORE, MANTLE, WATER AND GAS IN ONE GO (OPTIONAL)
#######################
# ExoMDN is a neural network trained on ~5.6 million interior models. Given mass, radius
# and Teq (with errors) it returns samples of the mass fractions of the iron core, the
# silicate mantle, a water layer and an H/He atmosphere, in about a second.
if RUN_EXOMDN:
    from exomdn import ExoMDN

    exo = ExoMDN(model_path=str(EXOMDN_DIR / "models"), data_path=str(EXOMDN_DIR / "data"))
    exo.load_model(str(EXOMDN_DIR / "models" / "mass_radius_Teq"))

    mdn = {}
    for _, p in planets.iterrows():
        prediction, _, _ = exo.predict_with_error(x=[[p["M"], p["R"], p["Teq"]]],
                                                  errors=[[p["M_err"], p["R_err"], 10.0]],
                                                  samples=(200, 100))
        mdn[p["short"]] = prediction[exo.mf]

    layers = ["core_mf", "mantle_mf", "water_mf", "atmosphere_mf"]
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.5), sharey=False)
    for ax, layer in zip(axes, layers):
        ax.boxplot([mdn[n][layer] for n in mdn], labels=list(mdn), showfliers=False)
        ax.set_title(layer.replace("_mf", ""))
    axes[0].set_ylabel("Mass fraction")
    fig.suptitle("TRAPPIST-1: ExoMDN interior mass fractions")
    fig.tight_layout()
    plt.show()

# %%
#######################
# QUESTIONS FOR STUDENTS
#######################
"""
Are the TRAPPIST-1 planets more or less dense than Earth? Does that mean less iron, or water?
Which planets cannot be explained by iron and rock alone?
The dry fit and the water fit both match the data. What measurement could tell them apart?
How would a prior from the star's Fe/Mg/Si abundances help?
"""
