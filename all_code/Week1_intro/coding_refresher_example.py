# %% INITIAL IMPORTS

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import corner
from scipy.optimize import minimize
import emcee
import arviz as az
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import warnings
warnings.filterwarnings("ignore")

# Set global publication visual defaults
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.autolayout": True
})

# %% NASA EXOPLANET ARCHIVE ACCESS
# Using Astroquery to access nasa exoplanet archive data (see https://astroquery.readthedocs.io)

def fetch_sub_neptunes() -> pd.DataFrame:
    """Fetch transiting planets with precise radii and mass measurements."""

    cols="pl_name,pl_rade,pl_radeerr1,pl_bmasse,pl_bmasseerr1,pl_eqt,sy_vmag,sy_gaiamag,st_mass,st_teff,st_rad"
    wh="default_flag=1 AND tran_flag=1 AND pl_rade IS NOT NULL AND pl_rade > 0 AND pl_bmasse IS NOT NULL AND pl_bmasse > 0 AND pl_bmasseerr1 IS NOT NULL AND pl_bmasseerr1 > 0 AND pl_radeerr1 IS NOT NULL"

    table = NasaExoplanetArchive.query_criteria(table="ps",where=wh,select=cols)
    df = table.to_pandas()

    # Filter for Sub-Neptunes with precise measurements (<20% mass error)
    mask = (
        (df['pl_rade'] >= 1.8) &
        (df['pl_rade'] <= 6.0) &
        ((df['pl_bmasseerr1'] / df['pl_bmasse']) <= 0.20)
    )
    return df[mask].dropna().reset_index(drop=True)

df = fetch_sub_neptunes()

print(df)

# %% Plotting
fig, ax = plt.subplots(figsize=(10, 6))

# __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
# \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
#  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
#   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
#   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|



# %% MODIFY THIS PLOT

## Choose different variables

# __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
# \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
#  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
#   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
#   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|


## Include other planet types (giant planets, rocky planets)

# __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
# \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
#  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
#   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
#   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|


## Include error bars

# __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
# \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
#  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
#   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
#   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|


## Add seaborn kde (e.g.):

ax=sns.kdeplot(x=df['pl_bmasse'], y=df['pl_rade'],bw_adjust=0.9,fill=True)
scatter = ax.errorbar(
            df['pl_bmasse'],
            df['pl_rade'],fmt='.k',yerr=df['pl_radeerr1'],xerr=df['pl_bmasseerr1'],
            ecolor='#999',alpha=0.5)
plt.show()

# %% DEFINING MODELS

x = np.log(df['pl_rade'].values)
y = np.log(df['pl_bmasse'].values)
yerr = df['pl_bmasseerr1'].values / df['pl_bmasse'].values #log errors ~ rel error

def linear_model(theta, x):
    # Define a linear model from parameter array theta
    # __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
    # \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
    #  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
    #   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
    #   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|
    return None

def quadratic_model(theta, x):
    # Define a quadratic model from parameter array theta
    # __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
    # \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
    #  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
    #   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
    #   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|
    return None


# %% LOG LIKELIHOOD

def log_likelihood(theta, x, y, yerr, model=linear_model, **kwargs):
    # __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
    # \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
    #  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
    #   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
    #   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|
    return None

# %%

def fit_frequentist(x, y, yerr, **kwargs):
    """Scipy Least-Squares / MLE fit."""
    nll = lambda theta: log_likelihood(theta, x, y, yerr, **kwargs)*-1
    if 'model' in kwargs and kwargs['model']==quadratic_model:
        init=np.tile(1,3)
    else:
        init=[1.0,0.5]
    res = minimize(nll, init, method='Nelder-Mead')
    return res.x

bfit = fit_frequentist(x,y,yerr)
loglik = log_likelihood(bfit,x,y,yerr)
print(bfit, loglik)

# %% MODEL COMPARISON USING BIC FOR 2 and 3-PARAM MODELS

b3fit = fit_frequentist(x,y,yerr,model=quadratic_model)
loglik3 = log_likelihood(b3fit,x,y,yerr,model=quadratic_model)


def BIC(loglik,n_params,nsamps):
    # __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
    # \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
    #  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
    #   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
    #   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|
    return None

print("2-parameter:",calc_bic(loglik,2,len(x)),"| 3-parameter:",calc_bic(loglik3,3,len(x))," (Lowest is preferred)")

# %% PLOTTING BEST-FIT

fig, ax = plt.subplots(figsize=(10, 6))

sizes = (15 - df['sy_vmag'].clip(upper=15)) * 15 + 20

scatter = ax.errorbar(np.exp(x),np.exp(y),
            yerr=yerr*y,
            fmt='.',
            alpha=0.8,
            ecolor='#ddd',
            linewidth=0.5
        )
ax.set_ylabel(r'Planetary Mass $M_{\text{Earth}}$ [K]')
ax.set_xlabel(r'Planetary Radius $R_\oplus$ [$R_\text{Earth}$]')
ax.set_title('Best-fit power law model')
ax.grid(True, linestyle='--', alpha=0.5)

ax.plot(np.arange(1.5,6,0.1),
        np.exp(linear_model(bfit,
                            np.log(np.arange(1.5,6,0.1)))),
        '--')

plt.savefig("best_MR_fit.png", dpi=300)
plt.show()

# %% ADD A PRIOR FUNCTION

def log_prior_linear(theta):
    """
    Calculates log prior using a combination of:
    1. Uniform prior on slope 'a' (must be physical: 0 < a < 5)
    2. Gaussian prior on intercept 'b' centered on 0.5 with sigma=2.5
    """

    a, b = theta

    # __   __   ___    _   _   ____         _   _   _   ____    _   _   _____
    # \ \ / /  / _ \  | | | | |  _ \       | | | \ | | |  _ \  | | | | |_   _|
    #  \ V /  | | | | | | | | | |_) |      | | |  \| | | |_) | | | | |   | |
    #   | |   | |_| | | |_| | |  _ <       | | | |\  | |  __/  | |_| |   | |
    #   |_|    \___/   \___/  |_| \_\      |_| |_| \_| |_|      \___/    |_|

    return log_prior_b

def log_probability_lin(theta, x, y, yerr):
    """Log posterior = Log Prior + Log Likelihood."""
    lp = log_prior_linear(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

# %% PERFORMING MCMC SAMPLING WITH EMCEE

def run_mcmc_linear(x,y,yerr, p0_guess=[1.5, 0.5], nwalkers=32, nsteps=20000):
    """Sample posterior for linear model using emcee."""
    ndim = len(p0_guess)
    pos = p0_guess + 1e-4 * np.random.randn(nwalkers, ndim)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability_lin,args=[x,y,yerr])
    sampler.run_mcmc(pos, nsteps, progress=True)

    samples = sampler.get_chain(discard=5000, thin=15, flat=True)

    # Save per-observation log likelihoods for ArviZ model evaluation
    log_per_obs = np.array([
        [-0.5 * (((y[i] - linear_model(s, x[i])) / yerr[i])**2 + np.log(2 * np.pi * yerr[i]**2))
            for i in range(len(x))]
        for s in samples
    ])

    return samples, log_per_obs

trace, log = run_mcmc_linear(x,y,yerr, bfit)

# %% PLOTTING THE SAMPLES AS A "CORNER"

fig = corner.corner(
        trace,
        labels=[r"$\gamma$ (Slope)", r"$\ln(M_0)$ (Intercept)"],
        show_titles=True
    )


# %% PLOTTING PLANETS WITH A BEST-FIT REGION

# First getting the 1- and 2-sigma regions as a function of a new x (using all trace outputs)
xplot=np.arange(1.5,6,0.1)
ally=trace[:,0][None,:]+np.log(xplot)[:,None]*trace[:,1][None,:]
yprcnts = np.exp((np.nanpercentile(ally,[5,16,50,84,95],axis=1)))

fig, ax = plt.subplots(figsize=(10, 6))

sizes = (15 - df['sy_vmag'].clip(upper=15)) * 15 + 20

scatter = ax.errorbar(np.exp(x),np.exp(y),
            yerr=yerr,
            fmt='.',
            alpha=0.8,
            ecolor='#888',
            linewidth=0.5
        )
ax.set_ylabel(r'Planetary Mass $M_{\text{Earth}}$ [K]')
ax.set_xlabel(r'Planetary Radius $R_\oplus$ [$R_\text{Earth}$]')
ax.set_title('Best-fit power law model')
ax.grid(True, linestyle='--', alpha=0.5)

plt.fill_between(xplot,yprcnts[0],yprcnts[4],
        color='C5',alpha=0.15)
plt.fill_between(xplot,yprcnts[1],yprcnts[3],
        color='C5',alpha=0.15)

plt.savefig("best_MR_fit_region.png", dpi=300)
plt.show()

# %%
