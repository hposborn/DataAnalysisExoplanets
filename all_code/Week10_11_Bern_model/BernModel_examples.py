# %%
#######################
# THE BERN MODEL: SYNTHETIC PLANET POPULATIONS (NGPPS)
#######################
"""
Explore the New Generation Planetary Population Synthesis (NGPPS) of the Bern model on
DACE, and turn the synthetic planets into 'observed' populations to compare with reality.

Steps
    List the populations and their columns (dace-query)
    Load the 5 Gyr snapshots of NG74, NG75 and NG76
    Mass-period diagram (radius as symbol size, core mass fraction as colour)
    Formation and evolution tracks of selected planets
    Observational biases: Kepler-like transits, HARPS-like RVs, Roman-like microlensing,
    and a HabEx-like direct-imaging survey
    Compare the 'observed' synthetic planets with real ones (NASA Exoplanet Archive)

Setup (uv)
    uv add dace-query astroquery numpy scipy pandas matplotlib seaborn

Notes for students
    Each population has ~1000 systems around a 1 solar-mass star. NG76 starts each disc
    with 100 planetary embryos; NG74 and NG75 with fewer (see the descriptions printed
    below). The model runs in two stages: formation with N-body interactions (first
    ~20 Myr), then the evolution of each planet in isolation up to 10 Gyr. Orbital
    elements such as eccentricity and inclination are only updated during formation.
"""

# %%
#######################
# IMPORTS AND SETTINGS
#######################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
import urllib.parse

from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.table import Table, Column, MaskedColumn

ARCHIVE_CACHE = Path("archive_cache")


def query_archive(where, columns=None, order=None, table="pscomppars"):
    """
    Query the NASA Exoplanet Archive with astroquery and return a pandas DataFrame.

    where   : ONE simple condition, e.g. "hostname='TRAPPIST-1'" or "discoverymethod='Imaging'"
    columns : list of columns to keep (all columns are downloaded)
    order   : column to sort by

    The archive's firewall blocks some query patterns (e.g. 'and ... is not null'), so the
    request is always 'select * ... where <one condition>'. Do any further filtering
    (ranges, missing values) in pandas. Results are cached in ARCHIVE_CACHE; if the archive
    cannot be reached, the error gives a browser link and the file name to save it as.
    """
    import hashlib

    ARCHIVE_CACHE.mkdir(exist_ok=True)
    cache_file = ARCHIVE_CACHE / (hashlib.md5(f"{table}|{where}".encode()).hexdigest()[:12] + ".csv")
    if cache_file.exists():
        df = pd.read_csv(cache_file, low_memory=False)
    else:
        try:
            result = NasaExoplanetArchive.query_criteria(table=table, select="*", where=where)
        except Exception as err:
            browser_url = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
                           + urllib.parse.urlencode({"query": f"select * from {table} where {where}",
                                                     "format": "csv"}))
            raise RuntimeError(f"Archive query failed: {err}\n"
                               f"Workaround: open this link in a browser and save the CSV as {cache_file}:\n"
                               f"{browser_url}") from err
        result = Table(result)
        # drop non-tabular extras that astroquery adds (e.g. a SkyCoord 'sky_coord' column)
        result.remove_columns([c for c in result.colnames if not isinstance(result[c], (Column, MaskedColumn))])
        df = result.to_pandas()               # plain values (units dropped), masked -> NaN
        df.to_csv(cache_file, index=False)
    if order:
        df = df.sort_values(order)
    if columns:
        df = df[columns]
    return df.reset_index(drop=True)
from dace_query.population import Population

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})
rng = np.random.default_rng(42)

POPULATIONS = ["ng74", "ng75", "ng76"]
MAIN_POP = "ng76"
AGE = "5000000000"          # snapshot age in years (see Population.get_snapshot_ages())
M_STAR = 1.0                # host star mass [M_sun] for the nominal NGPPS populations
M_MIN = 0.015               # [M_earth] drop embryos that barely grew (they start at 0.01)
A_MAX = 300.0               # [AU] drop bodies scattered onto very wide orbits

R_EARTH_PER_RJUP = 11.209
M_EARTH_PER_MJUP = 317.83
M_SUN_PER_M_EARTH = 3.003e-6
R_SUN_CM, AU_CM, R_EARTH_CM = 6.957e10, 1.496e13, 6.371e8

# %%
#######################
# THE POPULATIONS ON DACE
#######################
pops_info = Population.query_database(output_format="pandas")
mask = pops_info.astype(str).apply(lambda col: col.str.lower().isin(POPULATIONS)).any(axis=1)
pd.set_option("display.max_colwidth", 120)
print(pops_info[mask].T if mask.any() else pops_info.head(20).T)

# %%
#######################
# COLUMN NAMES AND LABELS
#######################
# get_columns returns a table with, for each quantity, the name used in requests and a
# human-readable label (and usually a unit). We look quantities up by their label.
cols = Population.get_columns(MAIN_POP, output_format="pandas")
print("Columns of the column table:", list(cols.columns))
print(cols.head(10).to_string())

label_col = next(c for c in cols.columns if cols[c].astype(str).eq("Eccentricity").any())
key_candidates = [c for c in cols.columns if c != label_col and cols[c].astype(str).nunique() == len(cols)
                  and not cols[c].astype(str).str.contains(" ").any()]
key_col = key_candidates[0] if key_candidates else label_col
unit_col = next((c for c in cols.columns if "unit" in c.lower()), None)
print(f"\nlabels in '{label_col}', request names in '{key_col}'"
      + (f", units in '{unit_col}'" if unit_col else ""))

LABEL_TO_KEY = dict(zip(cols[label_col].astype(str), cols[key_col].astype(str)))
LABEL_TO_UNIT = dict(zip(cols[label_col].astype(str), cols[unit_col].astype(str))) if unit_col else {}

# Short names for the quantities we use (label on DACE -> name in this script)
WANTED = {
    "System ID": "system_id",
    "Planet ID": "planet_id",
    "Time": "time",
    "Core mass": "m_core",                                    # M_earth
    "Envelope mass": "m_env",                                 # M_earth
    "Total mass": "m_tot",                                    # M_earth
    "Total radius": "r_tot",                                  # R_jup
    "Transit radius": "r_transit",                            # R_jup (evolution phase only)
    "Semi-major axis": "a",                                   # AU
    "Eccentricity": "ecc",
    "Inclination": "inc",                                     # rad, relative to the initial disc
    "Longitude of the ascending node": "long_node",           # rad
    "Longitude of the pericentre": "long_peri",               # rad
    "Mean anomaly": "mean_anom",                              # rad
    "Status of planet": "status",                             # 0 = still exists
    "Mass fraction of volatiles in the core": "f_ice",
    "Mass fraction of iron in the core's refractory contents": "f_iron",
    "Initial position of the embryo": "a_start",              # AU
    "Temperature at tau=2/3": "t_photo",                      # K
    "Stellar radius": "r_star",
}
missing = [lab for lab in WANTED if lab not in LABEL_TO_KEY]
if missing:
    print("Labels not found on DACE (skipped):", missing)
KEYS = {short: LABEL_TO_KEY[lab] for lab, short in WANTED.items() if lab in LABEL_TO_KEY}
for lab, short in WANTED.items():
    if lab in LABEL_TO_KEY:
        print(f"  {short:10s} <- {KEYS[short]:30s} [{LABEL_TO_UNIT.get(lab, '?')}]  {lab}")


def rename_to_short(df):
    """Rename DACE columns (request names or labels) to the short names used here."""
    mapping = {v: k for k, v in KEYS.items()}
    mapping.update({lab: short for lab, short in WANTED.items()})
    return df.rename(columns=mapping)


# %%
#######################
# LOAD THE 5 GYR SNAPSHOTS
#######################
snap_fields = ["system_id", "planet_id", "m_core", "m_env", "m_tot", "r_tot", "r_transit", "a",
               "ecc", "inc", "long_node", "status", "f_ice", "f_iron", "a_start", "r_star"]
snaps = {}
for pop in POPULATIONS:
    df = Population.get_snapshots(pop, AGE, columns=[KEYS[f] for f in snap_fields if f in KEYS],
                                  output_format="pandas")
    snaps[pop] = rename_to_short(pd.DataFrame(df))
    print(f"{pop}: {len(snaps[pop])} bodies")

# Quick look at the value ranges: check the units before using them
print(snaps[MAIN_POP].describe().T[["min", "50%", "max"]])

# %%
#######################
# CLEAN THE SNAPSHOTS AND DERIVE QUANTITIES
#######################
# Count the simulated stars BEFORE cleaning: systems whose bodies are all removed still
# count as stars when computing detections per star.
N_SYSTEMS = {}
for pop, df in snaps.items():
    N_SYSTEMS[pop] = df["system_id"].nunique()
    df = df[(df["status"] == 0) & (df["m_tot"] > M_MIN) & (df["a"] > 0) & (df["a"] < A_MAX)].copy()
    df["P_days"] = 365.25 * np.sqrt(df["a"] ** 3 / M_STAR)
    r = df["r_transit"] if "r_transit" in df and (df["r_transit"] > 0).any() else df["r_tot"]
    df["R_earth"] = np.where(r > 0, r, df["r_tot"]) * R_EARTH_PER_RJUP
    df["cmf"] = df["m_core"] / df["m_tot"]                   # solid core / total mass
    df["m_water"] = df["m_core"] * df["f_ice"]
    df["m_iron"] = df["m_core"] * (1 - df["f_ice"]) * df["f_iron"]
    df["m_silicate"] = df["m_core"] * (1 - df["f_ice"]) * (1 - df["f_iron"])
    df["m_gas"] = df["m_env"]
    # Stellar radius: convert to solar radii whatever the unit (cm or R_sun)
    if "r_star" in df and df["r_star"].median() > 1e6:
        df["r_star"] = df["r_star"] / R_SUN_CM
    df["r_star"] = df.get("r_star", pd.Series(1.0, index=df.index)).where(lambda x: x > 0, 1.0)
    snaps[pop] = df
    print(f"{pop}: {len(df)} surviving planets in {df['system_id'].nunique()} systems")

# %%
#######################
# MASS-PERIOD DIAGRAM AT 5 GYR
#######################
fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(6 * len(POPULATIONS), 5), sharey=True, squeeze=False)
for ax, pop in zip(axes[0], POPULATIONS):
    df = snaps[pop].sort_values("R_earth")
    sc = ax.scatter(df["P_days"], df["m_tot"], s=3 * np.clip(df["R_earth"], 0.3, 15) ** 1.5,
                    c=df["cmf"], cmap="viridis", vmin=0, vmax=1, alpha=0.7, edgecolors="none")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Period [days]")
    ax.set_title(f"{pop.upper()} at {int(AGE) / 1e9:.0f} Gyr ({len(df)} planets)")
axes[0, 0].set_ylabel(r"Mass [M$_\oplus$]")
fig.colorbar(sc, ax=axes[0].tolist(), label="Core mass fraction (solid core / total mass)")
for R in [1, 4, 11]:
    axes[0, -1].scatter([], [], s=3 * R ** 1.5, color="gray", label=f"{R} R$_\\oplus$")
axes[0, -1].legend(loc="lower right", fontsize=8, title="Radius")
plt.show()

# %%
#######################
# SELECT PLANETS FOR THE TRACKS
#######################
# Planets near 5 Earth radii and 50 Earth masses, inside P = 50 days (edit to explore).
R_TARGET, R_TOL = 5.0, 0.5        # +-25%
M_TARGET, M_TOL = 50.0, 5.0       # +-30%
P_MAX = 50.0
N_TRACKS_MAX = 10 #

df = snaps[MAIN_POP]
sel = df[(np.abs(df["R_earth"] / R_TARGET - 1) < R_TOL) & (np.abs(df["m_tot"] / M_TARGET - 1) < M_TOL)
         & (df["P_days"] < P_MAX)]
print(f"{len(sel)} planets selected in {MAIN_POP}")
if len(sel) > N_TRACKS_MAX:
    sel = sel.sample(N_TRACKS_MAX, random_state=1)
print(sel[["system_id", "planet_id", "m_tot", "m_core", "m_env", "R_earth", "P_days", "a_start", "f_ice"]]
      .round(3).to_string(index=False))
if len(sel) == 0:
    print("No planets selected: widen R_TOL / M_TOL or change the targets.")

# %%
#######################
# DOWNLOAD THE TRACKS
#######################
track_fields = ["time", "m_core", "m_env", "m_tot", "r_tot", "r_transit", "a", "ecc", "inc",
                "status", "f_ice", "f_iron", "t_photo"]
tracks = {}
for _, p in sel.iterrows():
    sid, pid = int(p["system_id"]), int(p["planet_id"])
    tr = Population.get_track(MAIN_POP, sid, pid, columns=[KEYS[f] for f in track_fields if f in KEYS],
                              output_format="pandas")
    tracks[(sid, pid)] = rename_to_short(pd.DataFrame(tr)).sort_values("time")
print(f"Downloaded {len(tracks)} tracks")

# %%
#######################
# PLOT THE TRACKS
#######################
# The vertical band marks the approximate end of the formation stage (~20 Myr): after it,
# eccentricity and inclination are no longer updated.
panels = [("m_tot", r"Total mass [M$_\oplus$]", "log"), ("m_env", r"Envelope mass [M$_\oplus$]", "log"),
          ("r_tot", r"Radius [R$_\oplus$]", "log"), ("a", "Semi-major axis [AU]", "log"),
          ("ecc", "Eccentricity", "linear"), ("f_ice", "Ice fraction of the core", "linear")]
fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for (sid, pid), tr in tracks.items():
    t = tr["time"].clip(lower=1e4)
    for ax, (col, label, scale) in zip(axes.flat, panels):
        if col not in tr:
            continue
        y = tr[col] * (R_EARTH_PER_RJUP if col == "r_tot" else 1)
        ax.plot(t, y, lw=1, label=f"{sid}-{pid}")
for ax, (col, label, scale) in zip(axes.flat, panels):
    ax.set_xscale("log")
    ax.set_yscale(scale)
    ax.set_ylabel(label)
    ax.axvspan(1e7, 3e7, color="gray", alpha=0.15)
for ax in axes[1]:
    ax.set_xlabel("Time [yr]")
axes[0, 0].legend(fontsize=7, ncol=2, title="system-planet")
fig.suptitle(f"{MAIN_POP.upper()}: planets near {R_TARGET} R$_\\oplus$, {M_TARGET} M$_\\oplus$, P < {P_MAX} d")
fig.tight_layout()
plt.show()

# %%
#######################
# OBSERVATIONAL BIASES: TARGETS AND VIEWING GEOMETRY
#######################
# A real survey observes many different stars. For every synthetic system and every one of
# N_VIEWS 'viewings', we draw (1) a random orientation of the system on the sky and (2) a
# target star from a survey-like sample (magnitude, size, noise, observing strategy).
# Averaging over the viewings gives each planet a detection PROBABILITY, which reflects the
# whole range of targets rather than one idealised star.
#
# The Bern planets formed around a 1 M_sun star, so orbits keep that stellar mass; the
# targets vary in brightness, radius and noise as Sun-like survey targets do.
from scipy.special import ndtr          # standard normal CDF, used for smooth detection thresholds

N_VIEWS = 50


def draw_power_law_magnitudes(size, mag_range, slope, rng):
    """Magnitudes with dN/dm ~ 10**(slope*m) (slope = 0.6 for stars uniform in space)."""
    lo, hi = 10 ** (slope * np.array(mag_range))
    return np.log10(lo + rng.uniform(0, 1, size) * (hi - lo)) / slope


def random_orientations(df, n_views, rng):
    """
    Sky inclination of every planet for n_views random viewing directions.
    Each system gets an isotropic orientation; planets are tilted from the system plane by
    their model inclination: i_sky ~ i_sys + inc * cos(node - node_observer).
    Also returns a random argument of periastron (the model's is measured in the disc frame).
    Returns cos_i, sin_i, omega with shape (n_views, n_planets), and the system index per planet.
    """
    sys_codes = pd.factorize(df["system_id"])[0]
    n_sys = sys_codes.max() + 1
    cos_sys = rng.uniform(0, 1, (n_views, n_sys))[:, sys_codes]
    node_obs = rng.uniform(0, 2 * np.pi, (n_views, n_sys))[:, sys_codes]
    inc = df["inc"].fillna(0).to_numpy() if "inc" in df else 0.0
    node = df["long_node"].fillna(0).to_numpy() if "long_node" in df else 0.0
    i_sky = np.arccos(cos_sys) + inc * np.cos(node - node_obs)
    omega = rng.uniform(0, 2 * np.pi, (n_views, len(df)))
    return np.abs(np.cos(i_sky)), np.abs(np.sin(i_sky)), omega, sys_codes


# %%
#######################
# A KEPLER-LIKE TRANSIT SURVEY
#######################
# Targets: Sun-like dwarfs with Kepler magnitudes 11-16 (more faint stars than bright ones),
# radii scattered around 1 R_sun, and a noise level (CDPP over 6.5 h) combining photon noise
# with a star-dependent variability term. Detection uses the Multiple Event Statistic (MES,
# the S/N of all observed transits) with a smooth recovery curve: 50% at MES = 7.1 and ~95%
# by MES ~ 9, similar in shape to Kepler's injection-recovery tests (Christiansen et al.
# 2016). At least 3 transits must fall in the data, which has gaps (duty cycle).
KEPLER = dict(T_obs_yr=4.0, duty_cycle=0.92, kepmag_range=(11.0, 16.0),
              photon_ppm_at_kp12=20.0, star_noise_median_ppm=20.0, star_noise_dex=0.25,
              r_star_median=1.0, r_star_dex=0.08, mes_50=7.1, mes_width=1.0, min_transits=3)


def draw_kepler_targets(size, rng, cfg=KEPLER):
    kepmag = draw_power_law_magnitudes(size, cfg["kepmag_range"], 0.6, rng)
    r_star = cfg["r_star_median"] * 10 ** (cfg["r_star_dex"] * rng.normal(size=size))
    star_noise = cfg["star_noise_median_ppm"] * 10 ** (cfg["star_noise_dex"] * rng.normal(size=size))
    photon = cfg["photon_ppm_at_kp12"] * 10 ** (0.2 * (kepmag - 12))
    return dict(kepmag=kepmag, r_star=r_star, cdpp=np.sqrt(photon**2 + star_noise**2))


def kepler_detection_probability(R_p_earth, P_days, a_AU, ecc, cos_i, sin_i, omega, star, rng, cfg=KEPLER):
    """
    Probability of detection for one target and orientation. Includes the transit geometry
    (eccentric orbits: impact parameter and duration scale with (1-e^2)/(1+e sin w)), grazing
    transits, the number of transits observed, and the MES recovery curve.
    Returns (probability, MES).
    """
    k = R_p_earth * R_EARTH_CM / (star["r_star"] * R_SUN_CM)
    a_over_R = a_AU * AU_CM / (star["r_star"] * R_SUN_CM)
    ecc_factor = (1 - ecc**2) / (1 + ecc * np.sin(omega))
    b = a_over_R * cos_i * ecc_factor
    transits = b < 1 + k
    chord = np.sqrt(np.clip((1 + k) ** 2 - b**2, 0, None))
    T14_hr = (24 * P_days / np.pi * np.arcsin(np.clip(chord / (a_over_R * np.clip(sin_i, 1e-6, None)), 0, 1))
              * np.sqrt(1 - ecc**2) / (1 + ecc * np.sin(omega)))
    depth_ppm = 1e6 * k**2 * np.clip((1 + k - b) / (2 * k), 0, 1)      # reduced for grazing transits
    n_possible = np.floor(cfg["T_obs_yr"] * 365.25 / P_days + rng.uniform(0, 1, np.shape(P_days)))
    n_obs = rng.binomial(n_possible.astype(int), cfg["duty_cycle"])
    mes = depth_ppm / star["cdpp"] * np.sqrt(n_obs * T14_hr / 6.5)
    p = ndtr((mes - cfg["mes_50"]) / cfg["mes_width"]) * (n_obs >= cfg["min_transits"]) * transits
    return p, mes


# %%
#######################
# A HARPS-LIKE RADIAL-VELOCITY SURVEY
#######################
# Targets: bright Sun-like stars (V = 6-10). Each has its own photon noise (scaling with
# brightness), stellar jitter (activity, granulation: ~1-3 m/s for quiet stars), number of
# RVs and time baseline. Detection: the S/N of the semi-amplitude, K / sigma_K with
# sigma_K = sigma_RV * sqrt(2/N), passes a smooth threshold around 5, and the period is
# shorter than the baseline.
HARPS = dict(V_range=(6.0, 10.0), sigma_ph_at_V8=0.5, t_exp_min=15.0, jitter_median=1.5, jitter_dex=0.2,
             n_rv_median=80, n_rv_dex=0.2, baseline_yr=(5.0, 15.0), snr_50=5.0, snr_width=0.7)


def rv_semi_amplitude(m_p_earth, P_days, sin_i, ecc=0.0, M_star_sun=1.0):
    """RV semi-amplitude K [m/s]."""
    m_jup = m_p_earth / M_EARTH_PER_MJUP
    M_tot = M_star_sun + m_p_earth * M_SUN_PER_M_EARTH
    return 28.4329 / np.sqrt(1 - ecc**2) * m_jup * sin_i * M_tot ** (-2 / 3) * (P_days / 365.25) ** (-1 / 3)


def draw_harps_targets(size, rng, cfg=HARPS):
    V = draw_power_law_magnitudes(size, cfg["V_range"], 0.6, rng)
    sigma_ph = cfg["sigma_ph_at_V8"] * 10 ** (0.2 * (V - 8)) * np.sqrt(15.0 / cfg["t_exp_min"])
    jitter = cfg["jitter_median"] * 10 ** (cfg["jitter_dex"] * rng.normal(size=size))
    n_rv = np.maximum(10, np.round(cfg["n_rv_median"] * 10 ** (cfg["n_rv_dex"] * rng.normal(size=size))))
    baseline = rng.uniform(*cfg["baseline_yr"], size)
    return dict(V=V, sigma_rv=np.sqrt(sigma_ph**2 + jitter**2), jitter=jitter, n_rv=n_rv, baseline_yr=baseline)


def harps_detection_probability(m_p_earth, P_days, ecc, sin_i, star, cfg=HARPS):
    """Returns (probability, S/N of K)."""
    K = rv_semi_amplitude(m_p_earth, P_days, sin_i, ecc, M_STAR)
    snr = K / (star["sigma_rv"] * np.sqrt(2.0 / star["n_rv"]))
    p = ndtr((snr - cfg["snr_50"]) / cfg["snr_width"]) * (P_days < star["baseline_yr"] * 365.25)
    return p, snr


# %%
#######################
# 'OBSERVE' THE SYNTHETIC SYSTEMS WITH KEPLER AND HARPS
#######################
for pop, df in snaps.items():
    cos_i, sin_i, omega, sys_codes = random_orientations(df, N_VIEWS, rng)
    n_sys_obs = sys_codes.max() + 1
    # one target star per system and viewing, shared by all its planets
    kep_targets = {k: v[:, sys_codes] for k, v in draw_kepler_targets((N_VIEWS, n_sys_obs), rng).items()}
    rv_targets = {k: v[:, sys_codes] for k, v in draw_harps_targets((N_VIEWS, n_sys_obs), rng).items()}

    ecc = df["ecc"].fillna(0).clip(0, 0.95).to_numpy() if "ecc" in df else 0.0
    p_kep, _ = kepler_detection_probability(df["R_earth"].to_numpy(), df["P_days"].to_numpy(), df["a"].to_numpy(),
                                            ecc, cos_i, sin_i, omega, kep_targets, rng)
    p_rv, _ = harps_detection_probability(df["m_tot"].to_numpy(), df["P_days"].to_numpy(), ecc, sin_i, rv_targets)

    df["p_kepler"] = p_kep.mean(axis=0)
    df["p_harps"] = p_rv.mean(axis=0)
    df["msini_example"] = df["m_tot"] * sin_i[0]            # one viewing direction, for plots
    n_sys = N_SYSTEMS[pop]                                  # all simulated stars, for yields per star
    print(f"{pop}: {df['p_kepler'].sum() / n_sys:.3f} Kepler and {df['p_harps'].sum() / n_sys:.2f} "
          f"HARPS detections per star (expected)")

# %%
#######################
# THE RANGE OF SENSITIVITY ACROSS EACH SURVEY
#######################
# Each target has its own sensitivity. Here: the smallest planet each Kepler target could
# detect (MES = 7.1, central transit, P = 10 d), and the smallest m sin i each HARPS target
# could detect (K/sigma_K = 5, P = 10 d).
kt = draw_kepler_targets(5000, rng)
P_ref = 10.0
a_ref = (M_STAR * (P_ref / 365.25) ** 2) ** (1 / 3)
T14_ref = 24 * P_ref / np.pi * kt["r_star"] * R_SUN_CM / (a_ref * AU_CM)
n_tr_ref = KEPLER["T_obs_yr"] * 365.25 * KEPLER["duty_cycle"] / P_ref
R_min = kt["r_star"] * R_SUN_CM / R_EARTH_CM * np.sqrt(
    KEPLER["mes_50"] * kt["cdpp"] * 1e-6 / np.sqrt(n_tr_ref * T14_ref / 6.5))

ht = draw_harps_targets(5000, rng)
K_min = HARPS["snr_50"] * ht["sigma_rv"] * np.sqrt(2.0 / ht["n_rv"])
msini_min = K_min / rv_semi_amplitude(1.0, P_ref, 1.0)                # K scales linearly with mass

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
sc = axes[0].scatter(kt["kepmag"], R_min, s=3, c=kt["r_star"], cmap="coolwarm")
fig.colorbar(sc, ax=axes[0], label=r"Stellar radius [R$_\odot$]")
axes[0].set_xlabel("Kepler magnitude")
axes[0].set_ylabel(r"Smallest detectable radius at P = 10 d [R$_\oplus$]")
axes[0].set_title("Kepler-like targets")
sc = axes[1].scatter(ht["V"], msini_min, s=3, c=ht["jitter"], cmap="viridis")
fig.colorbar(sc, ax=axes[1], label="Stellar jitter [m/s]")
axes[1].set_xlabel("V magnitude")
axes[1].set_ylabel(r"Smallest detectable m sin i at P = 10 d [M$_\oplus$]")
axes[1].set_title("HARPS-like targets")
fig.tight_layout()
plt.show()

# %%
#######################
# WHAT THE SURVEYS WOULD SEE
#######################
df = snaps[MAIN_POP]
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, prob, ycol, ylab, title in [
        (axes[0], "p_kepler", "R_earth", r"Radius [R$_\oplus$]", "Kepler-like transit survey"),
        (axes[1], "p_harps", "m_tot", r"Mass [M$_\oplus$]", "HARPS-like RV survey")]:
    ax.scatter(df["P_days"], df[ycol], s=2, color="lightgray", label="all planets")
    d = df[df[prob] > 1e-4]
    sc = ax.scatter(d["P_days"], d[ycol], s=8, c=d[prob], cmap="plasma", norm=plt.matplotlib.colors.LogNorm(),
                    label="detectable")
    fig.colorbar(sc, ax=ax, label="Detection probability")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Period [days]")
    ax.set_ylabel(ylab)
    ax.set_title(f"{MAIN_POP.upper()}: {title}")
    ax.legend(loc="upper left", fontsize=8)
fig.tight_layout()
plt.show()

# %%
#######################
# MICROLENSING: A ROMAN-LIKE SURVEY
#######################
# A planet is found when it perturbs the light curve of a background (source) star
# magnified by its host (the lens). Here the light-curve physics is approximated:
#   - lens distances follow the Galactic disc + bulge density, weighted by the event rate;
#     each event has its own lens-source proper motion and source brightness;
#   - the planet's perturbation lies at u_c = |s - 1/s| from the lens (s: projected
#     separation in Einstein radii) and must fall within the region the source crosses
#     (u_c < u_max), which gives the classic 'lensing zone' around s = 1;
#   - the perturbation has a size w ~ 4 sqrt(q) Einstein radii (q: planet/lens mass ratio),
#     so the source crosses it with probability ~ w / (2 u_max);
#   - the crossing lasts w * t_E and is sampled every 15 min; a detection needs
#     Delta chi^2 > 160 (as in Penny et al. 2019), with photometric noise set by the source
#     magnitude and reduced by the magnification, and the signal diluted when the source
#     star is larger than the perturbation (finite-source effect).
# Most lenses are M dwarfs: set MICROLENS_POP to an M-dwarf NGPPS population (Burn et al.
# 2021; see the population list at the top) and MICROLENS_STAR_MASS to its stellar mass.
MICROLENS_POP = "ngm11"
MICROLENS_STAR_MASS = 0.5      # CHECK: set to NGM11's host-star mass, from its description printed below

if MICROLENS_POP is not None:
    MICROLENS_POP = MICROLENS_POP.lower()            # DACE population IDs are lower case, e.g. 'ng76'
    ml_info = pops_info[pops_info.astype(str).apply(lambda col: col.str.lower().eq(MICROLENS_POP)).any(axis=1)]
    print(ml_info.T if len(ml_info) else f"'{MICROLENS_POP}' not found in the DACE population list")
    print(f"Using a lens (host) mass of {MICROLENS_STAR_MASS} M_sun: check it matches the description above.")
ROMAN = dict(D_S_kpc=8.0, cadence_min=15.0, sigma_at_H21=0.01, H_range=(18.0, 23.0), u_max=1.0,
             mu_rel_median=6.0, mu_rel_dex=0.2, delta=0.3, dchi2_min=160.0, R_source_sun=1.0)


def einstein_radius_mas(M_lens_sun, D_L_kpc, D_S_kpc=ROMAN["D_S_kpc"]):
    """Angular Einstein radius [mas]; multiplied by D_L [kpc] it gives the physical radius in AU."""
    return np.sqrt(8.144 * M_lens_sun * (1.0 / D_L_kpc - 1.0 / D_S_kpc))


def draw_lens_distances(size, rng, D_S_kpc=ROMAN["D_S_kpc"]):
    """Lens distances towards the bulge: (disc + bulge density) x D^2 x theta_E, sampled on a grid."""
    D = np.linspace(0.2, D_S_kpc - 0.05, 2000)
    disc = np.exp(-(D_S_kpc - D) / 2.6)                         # exponential disc, scale length 2.6 kpc
    bulge = 5.0 * np.exp(-0.5 * ((D - D_S_kpc) / 0.8) ** 2)     # bulge centred on the source distance
    weight = (disc + bulge) * D**2 * np.sqrt(1.0 / D - 1.0 / D_S_kpc)
    cdf = np.cumsum(weight) / weight.sum()
    return np.interp(rng.uniform(0, 1, size), cdf, D)


def draw_roman_events(size, rng, cfg=ROMAN):
    return dict(D_L=draw_lens_distances(size, rng, cfg["D_S_kpc"]),
                mu_rel=cfg["mu_rel_median"] * 10 ** (cfg["mu_rel_dex"] * rng.normal(size=size)),
                H=draw_power_law_magnitudes(size, cfg["H_range"], 0.3, rng))


def roman_detection_probability(q, s, theta_E_mas, event, cfg=ROMAN):
    """Probability that a planet (mass ratio q, separation s) is detected in one event."""
    u_c = np.abs(s - 1.0 / np.clip(s, 1e-6, None))
    w = 4 * np.sqrt(q) / (np.clip(s, 1e-6, None) * np.sqrt(np.abs(s**2 - 1) + q ** (1 / 3)))
    p_cross = np.clip(w / (2 * cfg["u_max"]), 0, 1) * (u_c < cfg["u_max"])
    u = np.clip(u_c, 1e-3, None)
    magnification = (u**2 + 2) / (u * np.sqrt(u**2 + 4))
    t_E_days = theta_E_mas / event["mu_rel"] * 365.25
    n_points = w * t_E_days * 24 * 60 / cfg["cadence_min"]
    theta_star_mas = cfg["R_source_sun"] * 0.00465 / cfg["D_S_kpc"]          # source angular radius
    rho = theta_star_mas / theta_E_mas
    delta = cfg["delta"] * np.minimum(1.0, (w / (2 * rho)) ** 2)            # finite-source dilution
    sigma = cfg["sigma_at_H21"] * 10 ** (0.2 * (event["H"] - 21)) / np.sqrt(magnification)
    dchi2 = n_points * (delta / sigma) ** 2
    return p_cross * (dchi2 > cfg["dchi2_min"])


# Load the microlensing population (same cuts as the main snapshots)
if MICROLENS_POP is not None:
    ml_raw = rename_to_short(pd.DataFrame(Population.get_snapshots(
        MICROLENS_POP, AGE, columns=[KEYS[f] for f in snap_fields if f in KEYS], output_format="pandas")))
    n_sys_ml = ml_raw["system_id"].nunique()             # all simulated stars, before cleaning
    ml = ml_raw[(ml_raw["status"] == 0) & (ml_raw["m_tot"] > M_MIN) & (ml_raw["a"] > 0)
                & (ml_raw["a"] < A_MAX)].copy()
    M_LENS, ml_name = MICROLENS_STAR_MASS, MICROLENS_POP.upper()
else:
    ml, M_LENS, ml_name = snaps[MAIN_POP].copy(), M_STAR, f"{MAIN_POP.upper()} (solar-mass hosts)"
    n_sys_ml = N_SYSTEMS[MAIN_POP]

# Each system acts as the lens of N_VIEWS events, each with its own lens distance, proper
# motion, source brightness, orientation and orbital phase
cos_i, sin_i, omega, sys_codes = random_orientations(ml, N_VIEWS, rng)
events = {k: v[:, sys_codes] for k, v in draw_roman_events((N_VIEWS, sys_codes.max() + 1), rng).items()}
phase = rng.uniform(0, 2 * np.pi, cos_i.shape)
r_perp = ml["a"].to_numpy() * np.sqrt(np.cos(phase) ** 2 + (cos_i * np.sin(phase)) ** 2)
theta_E = einstein_radius_mas(M_LENS, events["D_L"])
q = ml["m_tot"].to_numpy() * M_SUN_PER_M_EARTH / M_LENS
ml["q"] = q
ml["p_roman"] = roman_detection_probability(q, r_perp / (theta_E * events["D_L"]), theta_E, events).mean(axis=0)
print(f"{ml_name}: median Einstein radius {np.median(theta_E * events['D_L']):.2f} AU, "
      f"median t_E {np.median(theta_E / events['mu_rel'] * 365.25):.0f} d")
print(f"Expected planets detected per microlensing event: {ml['p_roman'].sum() / n_sys_ml:.4f}")

# Real microlensing planets (masses and projected separations as listed in the archive)
real_ml = query_archive("discoverymethod='Microlensing'", columns=["pl_name", "pl_bmasse", "pl_orbsmax"])
real_ml = real_ml.dropna(subset=["pl_bmasse", "pl_orbsmax"])

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(ml["a"], ml["m_tot"], s=2, color="lightgray", label="all synthetic planets")
d = ml[ml["p_roman"] > 1e-4]
sc = ax.scatter(d["a"], d["m_tot"], s=8, c=d["p_roman"], cmap="plasma",
                norm=plt.matplotlib.colors.LogNorm(), label="detectable")
ax.scatter(real_ml["pl_orbsmax"], real_ml["pl_bmasse"], marker="*", s=40, color="k", label="real microlensing")
fig.colorbar(sc, ax=ax, label="Detection probability per event")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Semi-major axis (synthetic) / projected separation (real) [AU]")
ax.set_ylabel(r"Mass [M$_\oplus$]")
ax.set_title(f"Roman-like microlensing: {ml_name}")
ax.legend(fontsize=8)
plt.show()

# %%
#######################
# DIRECT IMAGING: A HABEX-LIKE SURVEY
#######################
# Mature planets shine by reflected starlight. Targets: the n_stars nearest Sun-like stars
# (uniform in volume between d_min and d_max), with apparent magnitudes from their distance.
# A planet is detected if it lies outside the inner working angle (IWA), inside the outer
# working angle (OWA), above the contrast curve (which degrades towards the IWA), and is
# brighter than the faintest magnitude reachable in a realistic exposure. The albedo varies
# from planet to planet. Values are close to the HabEx concept; treat them as adjustable.
HABEX = dict(n_stars=150, d_min_pc=3.0, d_max_pc=15.0, IWA_arcsec=0.062, OWA_arcsec=1.0,
             contrast_floor=1e-10, V_lim_planet=31.0, albedo_range=(0.2, 0.5), M_V_star=4.83, M_V_scatter=0.3)


def lambert_phase(beta):
    """Lambert-sphere phase function for phase angle beta [rad] (1 at full phase)."""
    return (np.sin(beta) + (np.pi - beta) * np.cos(beta)) / np.pi


def reflected_contrast(R_p_earth, r_AU, beta, albedo=0.3):
    """Planet/star flux ratio in reflected light."""
    return albedo * lambert_phase(beta) * (R_p_earth * R_EARTH_CM / (r_AU * AU_CM)) ** 2


def contrast_curve(sep_arcsec, cfg=HABEX):
    """Faintest detectable contrast vs separation: the floor, degrading near the IWA."""
    return np.where(sep_arcsec > cfg["IWA_arcsec"],
                    cfg["contrast_floor"] * (1 + (cfg["IWA_arcsec"] / np.clip(sep_arcsec, 1e-6, None)) ** 3),
                    np.inf)


def draw_habex_targets(size, rng, cfg=HABEX):
    u = rng.uniform(0, 1, size)
    dist = (cfg["d_min_pc"] ** 3 + u * (cfg["d_max_pc"] ** 3 - cfg["d_min_pc"] ** 3)) ** (1 / 3)
    V = cfg["M_V_star"] + 5 * np.log10(dist / 10) + cfg["M_V_scatter"] * rng.normal(size=size)
    return dict(dist_pc=dist, V=V)


def habex_detection(R_p_earth, a_AU, cos_i, sin_i, phase, albedo, star, cfg=HABEX):
    """Returns (detected, separation [arcsec], contrast) for circular orbits."""
    r_proj = a_AU * np.sqrt(np.cos(phase) ** 2 + (cos_i * np.sin(phase)) ** 2)
    beta = np.arccos(np.clip(sin_i * np.sin(phase), -1, 1))          # star-planet-observer angle
    contrast = reflected_contrast(R_p_earth, a_AU, beta, albedo)
    sep = r_proj / star["dist_pc"]
    V_planet = star["V"] - 2.5 * np.log10(np.clip(contrast, 1e-30, None))
    detected = (sep < cfg["OWA_arcsec"]) & (contrast > contrast_curve(sep, cfg)) & (V_planet < cfg["V_lim_planet"])
    return detected, sep, contrast


di = snaps[MAIN_POP]
cos_i, sin_i, omega, sys_codes = random_orientations(di, N_VIEWS, rng)
targets = {k: v[:, sys_codes] for k, v in draw_habex_targets((N_VIEWS, sys_codes.max() + 1), rng).items()}
phase = rng.uniform(0, 2 * np.pi, cos_i.shape)
albedo = rng.uniform(*HABEX["albedo_range"], cos_i.shape)
det, sep, contrast = habex_detection(di["R_earth"].to_numpy(), di["a"].to_numpy(), cos_i, sin_i, phase, albedo, targets)
di["p_habex"] = det.mean(axis=0)
di["sep_example"], di["contrast_example"] = sep[0], contrast[0]
n_sys = N_SYSTEMS[MAIN_POP]

# 'Exo-Earth candidates': rocky-size planets in the habitable zone of a Sun-like star
eec = di["R_earth"].between(0.5, 1.5) & di["a"].between(0.95, 1.67)
print(f"Expected detections for {HABEX['n_stars']} stars: {HABEX['n_stars'] * di['p_habex'].sum() / n_sys:.1f} planets, "
      f"of which {HABEX['n_stars'] * di.loc[eec, 'p_habex'].sum() / n_sys:.2f} exo-Earth candidates")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
ax.scatter(di["sep_example"], di["contrast_example"], s=2, color="lightgray", label="one random view")
ax.scatter(di.loc[eec, "sep_example"], di.loc[eec, "contrast_example"], s=8, color="tab:green",
           label="exo-Earth candidates")
sep_grid = np.logspace(np.log10(HABEX["IWA_arcsec"] * 1.01), np.log10(HABEX["OWA_arcsec"]), 100)
ax.plot(sep_grid, contrast_curve(sep_grid), color="tab:red", label="contrast curve")
ax.axvline(HABEX["IWA_arcsec"], color="k", ls="--", label="IWA / OWA")
ax.axvline(HABEX["OWA_arcsec"], color="k", ls="--")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(1e-13, 1e-6)
ax.set_xlabel("Angular separation [arcsec]")
ax.set_ylabel("Planet / star contrast (reflected light)")
ax.legend(fontsize=8)

ax = axes[1]
ax.scatter(di["a"], di["m_tot"], s=2, color="lightgray", label="all synthetic planets")
d = di[di["p_habex"] > 0]
sc = ax.scatter(d["a"], d["m_tot"], s=8, c=d["p_habex"], cmap="plasma",
                norm=plt.matplotlib.colors.LogNorm(), label="detectable")
fig.colorbar(sc, ax=ax, label="Detection probability")
real_di = query_archive("discoverymethod='Imaging'", columns=["pl_name", "pl_bmasse", "pl_orbsmax"])
real_di = real_di.dropna(subset=["pl_bmasse", "pl_orbsmax"])
ax.scatter(real_di["pl_orbsmax"], real_di["pl_bmasse"], marker="*", s=40, color="k",
           label="real imaged planets (young, thermal)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Semi-major axis [AU]")
ax.set_ylabel(r"Mass [M$_\oplus$]")
ax.legend(fontsize=8)
fig.suptitle(f"HabEx-like reflected-light survey of {HABEX['n_stars']} nearby Sun-like stars ({MAIN_POP.upper()})")
fig.tight_layout()
plt.show()
# Note: today's imaged planets are young, hot giants seen in thermal emission at wide
# separations, a very different regime from a reflected-light survey of mature systems.

# %%
#######################
# REAL PLANETS FROM THE NASA EXOPLANET ARCHIVE
#######################
# Planets around roughly Sun-like stars. pl_bmasse is the true mass where known, else m sin i.
real_kepler = query_archive("disc_facility='Kepler'", columns=["pl_name", "pl_orbper", "pl_rade", "st_mass"])
real_kepler = real_kepler[real_kepler["st_mass"].between(0.8, 1.2)].dropna(subset=["pl_orbper", "pl_rade"])
real_rv = query_archive("discoverymethod='Radial Velocity'", columns=["pl_name", "pl_orbper", "pl_bmasse", "st_mass"])
real_rv = real_rv[real_rv["st_mass"].between(0.8, 1.2)].dropna(subset=["pl_orbper", "pl_bmasse"])
print(f"{len(real_kepler)} Kepler planets, {len(real_rv)} RV planets")

# %%
#######################
# SYNTHETIC VS REAL: TRANSIT SURVEY
#######################
# Compare shapes (normalised histograms): the real Kepler sample covers ~150 000 stars of
# many magnitudes, the synthetic one ~1000 stars at Kp = 13.
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
for pop, d in snaps.items():
    axes[0].hist(np.log10(d["P_days"]), bins=40, range=(-0.5, 3), weights=d["p_kepler"],
                 density=True, histtype="step", lw=1.5, label=pop.upper())
    axes[1].hist(np.log10(d["R_earth"]), bins=40, range=(-0.5, 1.5), weights=d["p_kepler"],
                 density=True, histtype="step", lw=1.5, label=pop.upper())
axes[0].hist(np.log10(real_kepler["pl_orbper"]), bins=40, range=(-0.5, 3), density=True,
             color="k", alpha=0.25, label="Kepler (real)")
axes[1].hist(np.log10(real_kepler["pl_rade"]), bins=40, range=(-0.5, 1.5), density=True,
             color="k", alpha=0.25, label="Kepler (real)")
axes[0].set_xlabel("log10 Period [d]")
axes[1].set_xlabel(r"log10 Radius [R$_\oplus$]")
axes[0].legend(fontsize=8)
d = snaps[MAIN_POP]
axes[2].scatter(real_kepler["pl_orbper"], real_kepler["pl_rade"], s=3, color="k", alpha=0.3, label="Kepler (real)")
d = d[d["p_kepler"] > 0.05]
axes[2].scatter(d["P_days"], d["R_earth"], s=6, color="tab:red", alpha=0.6, label=f"{MAIN_POP.upper()} detectable")
axes[2].set_xscale("log")
axes[2].set_yscale("log")
axes[2].set_xlabel("Period [d]")
axes[2].set_ylabel(r"Radius [R$_\oplus$]")
axes[2].legend(fontsize=8)
fig.tight_layout()
plt.show()

# %%
#######################
# KDE COMPARISON IN LOG PERIOD - LOG RADIUS
#######################
# Kernel density estimates (smoothed 2D histograms) of three distributions:
#   intrinsic : every synthetic planet (what the model makes)
#   observed  : synthetic planets weighted by their Kepler-like detection probability
#   Kepler    : real planets found by Kepler around roughly Sun-like stars
# Comparing intrinsic vs observed shows the survey's bias; observed vs Kepler tests the model.
# Contours enclose 20%, 40%, 60% and 80% of each distribution.
import seaborn as sns

d = snaps[MAIN_POP]
d = d[(d["P_days"] > 0.5) & (d["P_days"] < 1000)]
synth = pd.DataFrame({"log10 P [d]": np.log10(d["P_days"]), "log10 R [R_E]": np.log10(d["R_earth"]),
                      "weight": d["p_kepler"]})
real = pd.DataFrame({"log10 P [d]": np.log10(real_kepler["pl_orbper"]),
                     "log10 R [R_E]": np.log10(real_kepler["pl_rade"])})
xlim, ylim = (-0.3, 3.0), (-0.4, 1.4)
levels = [0.2, 0.4, 0.6, 0.8]          # seaborn levels are probability-mass thresholds from the outside in

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharex=True, sharey=True)
ax = axes[0]
sns.kdeplot(data=synth, x="log10 P [d]", y="log10 R [R_E]", ax=ax, levels=levels, fill=True,
            cmap="Greys", alpha=0.6, label="intrinsic")
sns.kdeplot(data=synth[synth["weight"] > 0], x="log10 P [d]", y="log10 R [R_E]", weights="weight",
            ax=ax, levels=levels, color="tab:red", linewidths=1.5)
ax.set_title(f"{MAIN_POP.upper()}: intrinsic (grey) vs observed by Kepler (red)")

ax = axes[1]
sns.kdeplot(data=real, x="log10 P [d]", y="log10 R [R_E]", ax=ax, levels=levels, fill=True,
            cmap="Blues", alpha=0.6)
sns.kdeplot(data=synth[synth["weight"] > 0], x="log10 P [d]", y="log10 R [R_E]", weights="weight",
            ax=ax, levels=levels, color="tab:red", linewidths=1.5)
ax.scatter(real["log10 P [d]"], real["log10 R [R_E]"], s=1, color="tab:blue", alpha=0.3)
ax.set_title(f"Real Kepler planets (blue) vs {MAIN_POP.upper()} observed (red)")

for ax in axes:
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    for R in [1.5, 4.0]:                                 # rough super-Earth / sub-Neptune / giant borders
        ax.axhline(np.log10(R), color="k", lw=0.5, ls=":")
fig.tight_layout()
plt.show()

# 1D marginals of the same three distributions
fig, axes = plt.subplots(1, 2, figsize=(13, 4))
for ax, col in zip(axes, ["log10 P [d]", "log10 R [R_E]"]):
    sns.kdeplot(synth[col], ax=ax, color="grey", fill=True, alpha=0.3, label=f"{MAIN_POP.upper()} intrinsic")
    sns.kdeplot(x=synth.loc[synth["weight"] > 0, col], weights=synth.loc[synth["weight"] > 0, "weight"],
                ax=ax, color="tab:red", lw=2, label=f"{MAIN_POP.upper()} observed")
    sns.kdeplot(real[col], ax=ax, color="tab:blue", lw=2, label="Kepler (real)")
    ax.set_xlabel(col)
axes[0].legend(fontsize=8)
fig.tight_layout()
plt.show()
# Questions: where does the bias move the peak of the radius distribution? Does the model
# reproduce the radius valley near 1.5-2 R_E, or the drop in small planets at short periods?

# %%
#######################
# SYNTHETIC VS REAL: RV SURVEY
#######################
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for pop, d in snaps.items():
    axes[0].hist(np.log10(d["m_tot"]), bins=40, range=(-0.5, 4), weights=d["p_harps"],
                 density=True, histtype="step", lw=1.5, label=pop.upper())
axes[0].hist(np.log10(real_rv["pl_bmasse"]), bins=40, range=(-0.5, 4), density=True,
             color="k", alpha=0.25, label="RV planets (real)")
axes[0].set_xlabel(r"log10 Mass [M$_\oplus$]")
axes[0].legend(fontsize=8)
d = snaps[MAIN_POP]
axes[1].scatter(real_rv["pl_orbper"], real_rv["pl_bmasse"], s=4, color="k", alpha=0.3, label="RV planets (real)")
d = d[d["p_harps"] > 0.05]
axes[1].scatter(d["P_days"], d["msini_example"], s=6, color="tab:blue", alpha=0.6,
                label=f"{MAIN_POP.upper()} detectable (m sin i)")
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_xlabel("Period [d]")
axes[1].set_ylabel(r"M or M sin i [M$_\oplus$]")
axes[1].legend(fontsize=8)
fig.tight_layout()
plt.show()

# %%
#######################
# QUESTIONS FOR STUDENTS
#######################
"""
Where do the water-rich planets end up, and where did they start (a_start)?
Which features of the real Kepler period and radius distributions does the model reproduce?
Which does it miss (e.g. the radius valley, the number of hot Jupiters)?
How do the results change from NG74 to NG76, i.e. with the number of embryos per disc?
How do the detection yields change for a fainter star, fewer RVs, or a shorter mission?
Which planets does each technique find, and which does no technique find?
Why do microlensing and direct imaging probe planets that transits and RVs miss?
"""
