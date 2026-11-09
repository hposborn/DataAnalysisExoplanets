# %%
#######################
# SUB-NEPTUNE INTERIORS: THE TOI-270 SYSTEM
#######################
"""
Are the TOI-270 planets rocky, water-rich, or rocky with a hydrogen-helium envelope?

Steps
    Load masses, radii, temperatures and insolations from the NASA Exoplanet Archive
    Compare them with mass-radius curves for rock, water worlds and H/He envelopes
    'Gas-rich' fit with smint: Earth-like core + H/He envelope (Lopez & Fortney 2014)
    'Water-rich' fit with smint: iron + rock + irradiated water (Aguichine et al. 2021)
    (Optional) ExoMDN: core, mantle, water and gas fractions together

Setup (uv)
    git clone https://github.com/cpiaulet/smint.git all_code/smint (code + model grids)
    uv add numpy scipy pandas matplotlib emcee corner astropy
    uv add ./all_code/smint
  For ExoMDN, see the setup notes in interiors_rocky_trappist1.py.

Key idea for students
    About 1% of the mass in H/He inflates a planet as much as tens of % of water.
    Both models can fit the same mass and radius: this is the sub-Neptune degeneracy.
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
from astropy.io import ascii as aioascii

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
rng = np.random.default_rng(42)

HOST = "TOI-270"
SMINT_DIR = os.path.join(os.path.dirname(os.path.abspath(os.getcwd())),"smint")
MODELS = str(SMINT_DIR / "example" / "smint_models") + "/"
EXOMDN_DIR = Path(os.environ.get("EXOMDN_DIR", "ExoMDN"))
OUTDIR = Path("interiors_outputs") / "subneptunes"
OUTDIR.mkdir(parents=True, exist_ok=True)

RUN_GAS_FIT = True
RUN_WATER_FIT = True
RUN_EXOMDN = False

# %%
#######################
# PLANET PARAMETERS FROM THE NASA EXOPLANET ARCHIVE
#######################
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
    # Fallback: approximate values after Van Eylen et al. (2021). Check them before use.
    print(f"Archive query failed ({exc}); using approximate built-in values")
    planets = pd.DataFrame({
        "name": ["TOI-270 b", "TOI-270 c", "TOI-270 d"],
        "M": [1.58, 6.15, 4.78], "M_err": [0.26, 0.37, 0.43],
        "R": [1.21, 2.36, 2.13], "R_err": [0.04, 0.06, 0.06],
        "Teq": [528.0, 410.0, 340.0], "S": [np.nan] * 3,
    })

planets["S"] = planets["S"].fillna((planets["Teq"] / 278.6) ** 4)   # zero-albedo estimate
planets["S_err"] = 0.1 * planets["S"]                                 # assume 10% if unknown
planets["short"] = planets["name"].str.split().str[-1]
print(planets.round(3).to_string(index=False))

# %%
#######################
# MODEL GRIDS FROM SMINT
#######################
from smint import fit_fhhe, fit_irrow

# Iron + rock + water (Aguichine et al. 2021, rocky end from Zeng et al. 2016)
interp_r = fit_irrow.make_interpolator_A21(MODELS + "Aguichine2021/withZengRocky/")
interp_valid = fit_irrow.make_interpolator_A21(MODELS + "Aguichine2021/", which_quantity="validity")

# Earth-like core + H/He envelope (Lopez & Fortney 2014); the envelope axis is log10(f_env in %)
lf14_table = aioascii.read(MODELS + "master_table_LF14_20201014.csv")
lf14_radii = np.load(MODELS + "LF14_20201014.npy")
interp_lf14 = fit_fhhe.make_interpolator_LF14(lf14_table, lf14_radii, log_fenv_prior=True)

# %%
#######################
# MASS-RADIUS DIAGRAM
#######################
# Water and H/He curves depend on temperature, so we draw them for TOI-270 c (the middle
# planet); the other planets are hotter (b) or cooler (d).
ref = planets.iloc[len(planets) // 2]
Tirr_ref = float(np.clip(ref["Teq"], 400, 1300))
mass_grid = np.logspace(0, np.log10(20), 150)
fig, ax = plt.subplots(figsize=(8, 6))

for cmf, label in [(0.9, "90% iron"), (0.33, "Earth-like")]:
    pts = np.column_stack([np.full_like(mass_grid, cmf), np.full_like(mass_grid, 400.0),
                           np.zeros_like(mass_grid), np.log10(mass_grid)])
    ax.plot(mass_grid, interp_r(pts), "-", lw=2, color="k" if cmf > 0.5 else "tab:brown", label=label)

for fh2o, color in [(0.1, "tab:cyan"), (0.3, "tab:blue"), (0.5, "navy")]:
    pts = np.column_stack([np.full_like(mass_grid, 0.33), np.full_like(mass_grid, Tirr_ref),
                           np.full_like(mass_grid, fh2o), np.log10(mass_grid)])
    ax.plot(mass_grid, interp_r(pts), "--", color=color, label=f"{fh2o:.0%} water ({Tirr_ref:.0f} K)")

for fenv, color in [(0.5, "gold"), (1.0, "orange"), (2.0, "tab:red"), (5.0, "darkred")]:
    pts = np.column_stack([np.full_like(mass_grid, 1.0), np.full_like(mass_grid, 5.0),
                           np.full_like(mass_grid, np.log10(ref["S"])), np.log10(mass_grid),
                           np.full_like(mass_grid, np.log10(fenv))])
    ax.plot(mass_grid, interp_lf14(pts), ":", lw=2, color=color, label=f"{fenv}% H/He (5 Gyr)")

ax.errorbar(planets["M"], planets["R"], xerr=planets["M_err"], yerr=planets["R_err"],
            fmt="o", color="k", ms=6, zorder=5)
for _, p in planets.iterrows():
    ax.annotate(p["short"], (p["M"], p["R"]), xytext=(6, 6), textcoords="offset points")
ax.set_xscale("log")
ax.set_xlim(1, 20)
ax.set_ylim(0.8, 3.5)
ax.set_xlabel(r"Mass [M$_\oplus$]")
ax.set_ylabel(r"Radius [R$_\oplus$]")
ax.set_title("TOI-270 planets: rock, water worlds and H/He envelopes")
ax.legend(fontsize=7, loc="upper left", ncol=2)
plt.show()

# %%
#######################
# GAS-RICH FIT: EARTH-LIKE CORE + H/He ENVELOPE (SMINT)
#######################
# Parameters: log10 envelope mass fraction (in %), mass, age and insolation. The envelope
# metallicity is fixed, so we run it twice: solar and 50x solar. The age prior is flat
# from 1 to 10 Gyr. The models' envelope fraction stops at 0.01%: a posterior piled up
# there means 'no envelope needed'.
gas = {}
if RUN_GAS_FIT:
    for _, p in planets.iterrows():
        params = dict(Mp_earth=p["M"], err_Mp_earth=p["M_err"], Rp_earth=p["R"], err_Rp_earth=p["R_err"],
                      Sinc_earth=p["S"], err_Sinc_earth=p["S_err"], flat_age=True,
                      age_Gyr_inf=1.0, age_Gyr_sup=10.0, age_Gyr=5.0, err_age_Gyr=3.0,
                      log_fenv_prior=True, nwalkers=50, ndim=4, nsteps=3000, frac_burnin=0.5,
                      save=False, outputdir=str(OUTDIR / "smint"), fname="hhe_" + p["short"])
        params = fit_fhhe.setup_priors(params)
        params = fit_fhhe.ini_fit(params)
        for met in [1, 50]:
            sampler = fit_fhhe.run_fit(params, interp_lf14, met=met)
            chain = sampler.chain[:, int(params["frac_burnin"] * params["nsteps"]):, :]
            gas[(p["short"], met)] = 10 ** chain.reshape(-1, 4)[:, 0] / 100   # % -> fraction

    fig, axes = plt.subplots(1, len(planets), figsize=(4.5 * len(planets), 3.5), squeeze=False)
    for ax, (_, p) in zip(axes[0], planets.iterrows()):
        for met, color in [(1, "tab:blue"), (50, "tab:green")]:
            ax.hist(np.log10(gas[(p["short"], met)]), bins=50, histtype="step", color=color,
                    density=True, label=f"{met}x solar")
        ax.set_xlabel("log10 H/He mass fraction")
        ax.set_title(p["name"])
    axes[0, 0].legend()
    fig.tight_layout()
    plt.show()

# %%
#######################
# WATER-RICH FIT: IRON + ROCK + IRRADIATED WATER (SMINT)
#######################
# Parameters: f'core (iron fraction of the iron + rock part), water mass fraction,
# irradiation temperature and mass. No H/He. The models assume a hot (steam/supercritical)
# water layer and start at 400 K; cooler planets are treated as 400 K.
water = {}
if RUN_WATER_FIT:
    for _, p in planets.iterrows():
        Tirr = float(np.clip(p["Teq"], 400, 1300))
        params = dict(Mp_earth=p["M"], err_Mp_earth=p["M_err"], Rp_earth=p["R"], err_Rp_earth=p["R_err"],
                      Tirr=Tirr, err_Tirr=20.0, Tirr_min=400.0, use_KDE_for_Mp_prior=False,
                      nwalkers=50, ndim=4, nsteps=3000, frac_burnin=0.5, save=False,
                      outputdir=str(OUTDIR / "smint"), fname="water_" + p["short"])
        params = fit_irrow.setup_priors(params)
        params = fit_irrow.ini_fit(params)
        sampler = fit_irrow.run_fit(params, interp_r, interp_valid)
        chain = sampler.chain[:, int(params["frac_burnin"] * params["nsteps"]):, :]
        water[p["short"]] = chain.reshape(-1, 4)

    fig, axes = plt.subplots(1, len(planets), figsize=(4.5 * len(planets), 3.5), squeeze=False)
    for ax, (_, p) in zip(axes[0], planets.iterrows()):
        ax.hist(water[p["short"]][:, 1], bins=50, histtype="step", density=True, color="tab:blue")
        ax.set_xlabel("Water mass fraction")
        ax.set_title(p["name"])
    fig.tight_layout()
    plt.show()

# %%
#######################
# SUMMARY: TWO EXPLANATIONS FOR EACH PLANET
#######################
rows = []
for _, p in planets.iterrows():
    row = dict(planet=p["short"])
    for met in [1, 50]:
        if (p["short"], met) in gas:
            f = gas[(p["short"], met)]
            row[f"H/He % ({met}x sol)"] = f"{100 * np.median(f):.2f} (+{100 * (np.percentile(f, 84) - np.median(f)):.2f} -{100 * (np.median(f) - np.percentile(f, 16)):.2f})"
    if p["short"] in water:
        w = water[p["short"]][:, 1]
        row["water %"] = f"{100 * np.median(w):.0f} (+{100 * (np.percentile(w, 84) - np.median(w)):.0f} -{100 * (np.median(w) - np.percentile(w, 16)):.0f})"
    rows.append(row)
print(pd.DataFrame(rows).to_string(index=False))
print("\nThe same planet needs either ~1% H/He or tens of % water. Which is more likely")
print("depends on formation (inside or outside the snow line) and on atmospheric spectra.")

# %%
#######################
# EXOMDN: ALL FOUR LAYERS TOGETHER (OPTIONAL)
#######################
# ExoMDN lets core, mantle, water and H/He all vary at once, so it shows the full
# degeneracy rather than two separate scenarios.
if RUN_EXOMDN:
    from exomdn import ExoMDN

    exo = ExoMDN(model_path=str(EXOMDN_DIR / "models"), data_path=str(EXOMDN_DIR / "data"))
    exo.load_model(str(EXOMDN_DIR / "models" / "mass_radius_Teq"))
    fig, axes = plt.subplots(1, len(planets), figsize=(4.5 * len(planets), 3.5), squeeze=False)
    for ax, (_, p) in zip(axes[0], planets.iterrows()):
        prediction, _, _ = exo.predict_with_error(x=[[p["M"], p["R"], p["Teq"]]],
                                                  errors=[[p["M_err"], p["R_err"], 10.0]],
                                                  samples=(200, 100))
        ax.boxplot([prediction[c] for c in exo.mf], labels=["core", "mantle", "water", "H/He"],
                   showfliers=False)
        ax.set_title(p["name"])
    axes[0, 0].set_ylabel("Mass fraction")
    fig.tight_layout()
    plt.show()

# %%
#######################
# QUESTIONS FOR STUDENTS
#######################
"""
Which TOI-270 planet looks rocky, and which need either H/He or water?
Why does the H/He fraction change between solar and 50x solar envelope metallicity?
What observations could break the water-vs-gas degeneracy (hint: JWST spectra)?
"""
