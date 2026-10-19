# %%
#######################
# L1 PERIODOGRAM OF GJ 876 (OPTIONAL EXTRA)
#######################
"""
The l1 periodogram (Hara, Boue, Laskar & Correia 2017, MNRAS 464, 1220) fits many
sinusoids at once with a sparsity (l1) penalty. Because it does not search for one
signal at a time, it suffers much less from aliasing than Lomb-Scargle.

The code is not pip-installable, so it needs a one-off manual setup.

Setup steps (run in a terminal)
    1. Run the main notebook (gj876_rv_analysis.py) at least up to the cleaning
       section. It writes gj876_outputs/gj876_rv_clean.csv, which this script reads.

    2. Clone the repository next to this script:
           git clone https://github.com/nathanchara/l1periodogram.git

    3. Install its dependencies (numpy, scipy, matplotlib; check the repository
       README for any others, e.g. a requirements file):
           pip install numpy scipy matplotlib

    4. Look inside the cloned folder. The Python modules (l1periodogram_v1.py and
       helpers) sit in a subfolder, usually l1periodogram_codes/. If the folder name
       differs, edit L1_CODE_DIR below.

    5. Run this script from the folder that contains both it and the clone:
           python l1_periodogram_gj876.py
       or run it cell by cell in an editor that understands # %% cells.

The repository README and its example notebook are the reference for the API.
The calls below follow that example; if a newer version renames something, the
error message will point to the line to adapt.
"""

# %%
#######################
# PATHS AND IMPORTS
#######################
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO_DIR = Path("l1periodogram")                    # the git clone
L1_CODE_DIR = REPO_DIR / "l1periodogram_codes"      # folder holding l1periodogram_v1.py
DATA_FILE = Path("gj876_outputs") / "gj876_rv_clean.csv"
PMIN = 1.2                                          # shortest period searched [d]
KNOWN_PERIODS = {"d": 1.9378, "c": 30.08, "b": 61.08, "e": 123.8}

if not L1_CODE_DIR.exists():
    raise FileNotFoundError(f"{L1_CODE_DIR} not found: clone the repository first (see the header).")
sys.path.insert(0, str(L1_CODE_DIR.resolve()))
sys.path.insert(0, str(REPO_DIR.resolve()))

import l1periodogram_v1

# %%
#######################
# LOAD THE CLEANED DATA
#######################
df = pd.read_csv(DATA_FILE).sort_values("t").reset_index(drop=True)
t = df["t"].to_numpy()
y = df["rv"].to_numpy()
err = df["err"].to_numpy()
names = sorted(df["inst"].unique())
idx = df["inst"].map({n: i for i, n in enumerate(names)}).to_numpy()
print(f"{len(df)} points from {len(names)} instruments: {names}")

# %%
#######################
# NOISE MODEL AND INSTRUMENT OFFSETS
#######################
# Each instrument gets its own unpenalised offset: one column of ones per instrument (MH0).
MH0 = np.column_stack([(idx == j).astype(float) for j in range(len(names))])

# Diagonal covariance: error bars plus a per-instrument jitter. The default 3 m/s is a
# placeholder; the 'scatter versus error bars' table in the main notebook gives an
# estimate for each instrument, so paste those values in.
JITTER = {n: 3.0 for n in names}                     # m/s, edit as needed
jit = df["inst"].map(JITTER).to_numpy()
V = np.diag(err**2 + jit**2)

# %%
#######################
# RUN THE L1 PERIODOGRAM
#######################
c = l1periodogram_v1.l1p_class(t, y)
c.set_model(omegamax=2 * np.pi / PMIN, V=V, MH0=MH0, Nphi=1)
c.l1_perio(numerical_method="lars",
           significance_evaluation_methods=["fap", "evidence_laplace"],
           max_n_significance_tests=8, verbose=0, plot_output=True)
plt.show()

# %%
#######################
# CUSTOM PLOT WITH KNOWN PLANETS MARKED
#######################
try:
    periods = 2 * np.pi / c.omegas[c.omegas > 0]
    amp = c.smoothed_solution[c.omegas > 0]
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(periods, amp, lw=0.8)
    for name, P in KNOWN_PERIODS.items():
        ax.axvline(P, color="tab:red", alpha=0.35)
        ax.text(P, 0.97, name, color="tab:red", ha="center", va="top",
                transform=ax.get_xaxis_transform())
    ax.set_xscale("log")
    ax.set_xlabel("Period [d]")
    ax.set_ylabel("l1 amplitude [m/s]")
    ax.set_title("l1 periodogram of GJ 876")
    plt.show()
except AttributeError as exc:
    print("Attribute names differ in this version of the code; use the built-in plot above.", exc)

# %%
#######################
# QUESTIONS FOR STUDENTS
#######################
"""
Which peaks does the l1 periodogram recover, and how do they compare with the GLS
prewhitening sequence in the main notebook?
Are there fewer alias peaks than in the GLS periodogram? Why?
How do the results change if you increase the jitter values?
"""
