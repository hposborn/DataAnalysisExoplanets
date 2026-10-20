import numpy as np
import pandas as pd

def mad_std(x):
    """Robust standard deviation from the median absolute deviation."""
    x = np.asarray(x, float)
    return 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))


def find_col(df, candidates):
    """Return the first column name from candidates that exists in df."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def inst_codes(df):
    """Integer instrument index for each row and the ordered instrument names."""
    names = sorted(df["inst"].unique())
    idx = df["inst"].map({n: i for i, n in enumerate(names)}).to_numpy()
    return idx, names


def design_matrix(t, idx, n_inst, periods, t0):
    """Columns: one offset per instrument, then a sin/cos pair per period."""
    cols = [(idx == j).astype(float) for j in range(n_inst)]
    for P in periods:
        ph = 2 * np.pi * (t - t0) / P
        cols += [np.sin(ph), np.cos(ph)]
    return np.column_stack(cols)


def linear_fit(df, periods=(), jitter=0.0, y=None):
    """
    Weighted linear least squares: per-instrument offsets + sinusoids at fixed periods.
    Returns the model, residuals, offsets and, per period, (P, K, tc) where tc is the
    time of conjunction for a circular orbit, v(t) = -K sin(2 pi (t - tc) / P),
    which matches the radvel/exoplanet convention.
    """
    t = df["t"].to_numpy()
    y = df["rv"].to_numpy() if y is None else np.asarray(y)
    e = df["err"].to_numpy()
    idx, names = inst_codes(df)
    t0 = np.median(t)
    X = design_matrix(t, idx, len(names), periods, t0)
    w = 1.0 / np.sqrt(e**2 + jitter**2)
    beta, *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)
    model = X @ beta
    offsets = dict(zip(names, beta[: len(names)]))
    sines = []
    for k, P in enumerate(periods):
        s, c = beta[len(names) + 2 * k], beta[len(names) + 2 * k + 1]
        phic = np.arctan2(c, -s)
        sines.append(dict(P=P, K=np.hypot(s, c), tc=t0 + phic * P / (2 * np.pi)))
    return dict(model=model, resid=y - model, offsets=offsets, sines=sines, t0=t0)



def jx_inputs(df, guesses, offsets):
    """Model arguments (times relative to the median) and starting values from the guesses."""
    idx, names = inst_codes(df)
    tref = np.median(df["t"])
    P0 = np.array([g["P"] for g in guesses])
    tc0 = np.array([g["tc"] for g in guesses]) - tref
    K0 = np.array([max(g["K"], 2.0) for g in guesses])
    gamma0 = np.array([offsets[nm] for nm in names])
    args = (df["t"].to_numpy() - tref, df["rv"].to_numpy(), df["err"].to_numpy(),
            idx, len(names), P0, tc0, K0, gamma0)
    init = dict(logP=np.log(P0), tc=tc0, logK=np.log(K0), ecc=np.full(len(P0), 0.05),
                w_vec=np.vstack([np.ones(len(P0)), np.zeros(len(P0))]),
                gamma=gamma0, log_jit=np.full(len(names), np.log(2.0)))
    return args, init, tref


def keplerian_rv(t, P, tc, ecc, omega, K):
    """RV of each planet, shape (n_planets, n_times). Same parameterisation as radvel."""
    from jaxoplanet.core.kepler import kepler
    import jax.numpy as jnp
    P, tc, ecc, omega, K = (jnp.atleast_1d(x)[:, None] for x in (P, tc, ecc, omega, K))
    # time of periastron from the time of conjunction (true anomaly pi/2 - omega)
    f_c = jnp.pi / 2 - omega
    E_c = 2 * jnp.arctan2(jnp.sqrt(1 - ecc) * jnp.sin(f_c / 2), jnp.sqrt(1 + ecc) * jnp.cos(f_c / 2))
    tp = tc - P / (2 * jnp.pi) * (E_c - ecc * jnp.sin(E_c))
    M = jnp.mod(2 * jnp.pi * (t[None, :] - tp) / P, 2 * jnp.pi)
    sinf, cosf = kepler(M, jnp.broadcast_to(ecc, M.shape))
    return K * (cosf * jnp.cos(omega) - sinf * jnp.sin(omega) + ecc * jnp.cos(omega))

def rv_model(t, y, err, idx, n_inst, P0, tc0, K0, gamma0):
    """
    Priors: narrow log-normal on P around the periodogram peak, wide on tc, log-uniform K.
    Eccentricity: the Kipping (2013) Beta(0.867, 3.03) prior, sampled directly.
    omega: from an unconstrained 2D vector (avoids the 0/2pi wrap-around that hurts NUTS).
    """
    import numpyro
    import numpyro.distributions as dist
    import jax.numpy as jnp
    from jaxoplanet.orbits.keplerian import Central, System
    KIPPING_A, KIPPING_B = 0.867, 3.03
    n = len(P0)
    # Priors: narrow around the periodogram period, broad on tc and K
    logP = numpyro.sample("logP", dist.Normal(jnp.log(P0), 0.02))
    P = numpyro.deterministic("P", jnp.exp(logP))
    tc = numpyro.sample("tc", dist.Normal(tc0, P0 / 4))
    logK = numpyro.sample("logK", dist.Uniform(-2.0, jnp.log(1000.0)).expand([n]))
    K = numpyro.deterministic("K", jnp.exp(logK))
    # Kipping (2013) eccentricity prior
    ecc = numpyro.sample("ecc", dist.Beta(KIPPING_A, KIPPING_B).expand([n]))
    # omega from a 2D vector, which avoids the wrap-around at +-pi that troubles NUTS
    w_vec = numpyro.sample("w_vec", dist.Normal(0.0, 1.0).expand([2, n]))
    omega = numpyro.deterministic("omega", jnp.arctan2(w_vec[1], w_vec[0]))
    # One offset and one jitter per instrument
    gamma = numpyro.sample("gamma", dist.Normal(gamma0, 100.0))
    log_jit = numpyro.sample("log_jit", dist.Normal(jnp.log(2.0), 2.0).expand([n_inst]))

    # The planetary system
    system = System(Central())
    for i in range(n):
        system = system.add_body(period=P[i], time_transit=tc[i], eccentricity=ecc[i],
                                 omega_peri=omega[i], radial_velocity_semiamplitude=K[i])
    rv_planets = numpyro.deterministic("rv_planets", system.radial_velocity(t))  # (n_planets, n_times)

    mu = rv_planets.sum(axis=0) + gamma[idx]
    sigma = jnp.sqrt(err**2 + jnp.exp(2 * log_jit)[idx])
    numpyro.sample("obs", dist.Normal(mu, sigma), obs=y)


def jx_map(args, init, seed=0):
    """MAP: minimise the negative log posterior with L-BFGS, using JAX gradients."""
    from numpyro.infer.util import initialize_model
    from jax.flatten_util import ravel_pytree
    import jax
    import jax.numpy as jnp
    from scipy import optimize
    info = initialize_model(jax.random.PRNGKey(seed), rv_model, model_args=args,
                            init_strategy=init_to_value(values=init))
    x0, unravel = ravel_pytree(info.param_info.z)
    nlp = jax.jit(lambda x: info.potential_fn(unravel(x)))
    grad = jax.jit(jax.grad(lambda x: info.potential_fn(unravel(x))))
    res = optimize.minimize(lambda x: float(nlp(x)), np.asarray(x0),
                            jac=lambda x: np.asarray(grad(x)), method="L-BFGS-B")
    map_soln = info.postprocess_fn(unravel(jnp.asarray(res.x)))
    return {k: np.asarray(v) for k, v in map_soln.items()}


def jx_sample(args, init_values, draws=1500, tune=1500, chains=4, seed=42):
    import jax
    from numpyro.infer import MCMC, NUTS, init_to_value
    sample_sites = ["logP", "tc", "logK", "ecc", "w_vec", "gamma", "log_jit"]
    kernel = NUTS(rv_model, target_accept_prob=0.95,
                  init_strategy=init_to_value(values={k: init_values[k] for k in sample_sites}))
    mcmc = MCMC(kernel, num_warmup=tune, num_samples=draws, num_chains=chains,
                chain_method="parallel", progress_bar=True)
    mcmc.run(jax.random.PRNGKey(seed), *args)
    return mcmc


def run_jaxoplanet(df, guesses, offsets, draws=1500, tune=1500, chains=4):
    """Sample the model with NUTS, starting from the periodogram guesses."""
    from numpyro.infer import MCMC, NUTS, init_to_value
    import jax
    idx, names = inst_codes(df)
    tref = np.median(df["t"])                 # times relative to the middle of the data
    P0 = np.array([g["P"] for g in guesses])
    tc0 = np.array([g["tc"] for g in guesses]) - tref
    K0 = np.array([max(g["K"], 2.0) for g in guesses])
    gamma0 = np.array([offsets[nm] for nm in names])
    args = (df["t"].to_numpy() - tref, df["rv"].to_numpy(), df["err"].to_numpy(),
            idx, len(names), P0, tc0, K0, gamma0)
    start = dict(logP=np.log(P0), tc=tc0, logK=np.log(K0), ecc=np.full(len(P0), 0.05),
                 w_vec=np.vstack([np.ones(len(P0)), np.zeros(len(P0))]),
                 gamma=gamma0, log_jit=np.full(len(names), np.log(2.0)))
    mcmc = MCMC(NUTS(rv_model, target_accept_prob=0.95, init_strategy=init_to_value(values=start)),
                num_warmup=tune, num_samples=draws, num_chains=chains, chain_method="parallel")
    mcmc.run(jax.random.PRNGKey(42), *args)
    return mcmc, args, tref

def waic(mcmc, args):
    """WAIC (elpd and its standard error) from the pointwise log-likelihood of the samples."""
    import numpyro
    ll = np.asarray(numpyro.infer.log_likelihood(rv_model, mcmc.get_samples(), *args)["obs"])
    lppd_i = np.logaddexp.reduce(ll, axis=0) - np.log(ll.shape[0])
    elpd_i = lppd_i - ll.var(axis=0, ddof=1)
    return elpd_i.sum(), np.sqrt(len(elpd_i) * elpd_i.var())

def freq_grid(t, pmin=12, pmax=3000, ofac=10):
    T = np.ptp(t)
    return np.arange(1.0 / pmax, 1.0 / pmin, 1.0 / (ofac * T))


def ls_periodogram(t, y, err, freq):
    from astropy.timeseries import LombScargle
    ls = LombScargle(t, y, err)
    power = ls.power(freq)
    fal = ls.false_alarm_level([0.1, 0.01, 0.001],
                               minimum_frequency=freq.min(), maximum_frequency=freq.max())
    return ls, power, fal


def plot_periodogram(ax, freq, power, label=None, fal=None, marks=[],
                     ylabel="LS power", prot=True, prot_range=(80.0, 100.0), **kw):
    ax.plot(1.0 / freq, power, lw=0.8, label=label, **kw)
    if fal is not None:
        for lvl, ls_, lab in zip(fal, [":", "--", "-."], ["10%", "1%", "0.1%"]):
            ax.axhline(lvl, color="k", ls=ls_, lw=0.7, label=f"FAP {lab}")
    for name, P in (marks or {}).items():
        ax.axvline(P, color="tab:red", alpha=0.35, lw=1)
        ax.text(P, 0.97, name, color="tab:red", ha="center", va="top",
                transform=ax.get_xaxis_transform(), fontsize=8)
    if prot:
        ax.axvspan(*prot_range, color="tab:orange", alpha=0.12)
    ax.set_xscale("log")
    ax.set_xlabel("Period [d]")
    ax.set_ylabel(ylabel)


def prewhiten(df, n_signals=4, pmin=1.2, pmax=3000, y=None, verbose=True):
    """
    Iterative prewhitening: find the highest LS peak, add a sinusoid at that period
    to a joint linear fit (with instrument offsets), repeat on the residuals.
    """
    from astropy.timeseries import LombScargle

    t, e = df["t"].to_numpy(), df["err"].to_numpy()
    freq = freq_grid(t, pmin, pmax)
    periods, faps, powers = [], [], []
    for _ in range(n_signals):
        fit = linear_fit(df, periods, y=y)
        ls = LombScargle(t, fit["resid"], e)
        p = ls.power(freq)
        k = np.argmax(p)
        fap = ls.false_alarm_probability(p[k], minimum_frequency=freq.min(),
                                         maximum_frequency=freq.max())
        periods.append(1.0 / freq[k])
        faps.append(fap)
        powers.append(p)
        if verbose:
            print(f"  signal {len(periods)}: P = {periods[-1]:10.4f} d   FAP = {fap:.2e}")
    final = linear_fit(df, periods, y=y)
    return dict(periods=periods, faps=faps, powers=powers, freq=freq, fit=final)


def offset_subtracted(df, offsets):
    return df["rv"].to_numpy() - df["inst"].map(offsets).to_numpy()



def bgls(t, y, err, pmin, pmax, oversample=5):
    """
    BGLS periodogram (Mortier et al. 2015).

    t, y, err  : times [d], RVs and RV errors [m/s] (numpy arrays)
    pmin, pmax : shortest and longest period to search [d]
    oversample : frequency grid points per peak width (1 / time span)

    Returns (periods, log_prob), where log_prob is the natural log of the probability
    relative to the best period (so the highest peak is 0).
    """
    freqs = np.arange(1 / pmax, 1 / pmin, 1 / (oversample * np.ptp(t)))
    tt = t - t.min()                     # small numbers keep the trigonometry precise
    w = 1 / err**2                       # weights
    W, Y = w.sum(), (w * y).sum()
    log_prob = np.zeros(freqs.size)
    for i in range(0, freqs.size, 2000):           # 2000 frequencies at a time to save memory
        omega = 2 * np.pi * freqs[i:i + 2000, None]
        theta = 0.5 * np.arctan2((w * np.sin(2 * omega * tt)).sum(1),
                                 (w * np.cos(2 * omega * tt)).sum(1))
        x = omega * tt - theta[:, None]
        C, S = (w * np.cos(x)).sum(1), (w * np.sin(x)).sum(1)
        YC, YS = (w * y * np.cos(x)).sum(1), (w * y * np.sin(x)).sum(1)
        CC, SS = (w * np.cos(x)**2).sum(1), (w * np.sin(x)**2).sum(1)
        K = (C**2 * SS + S**2 * CC - W * CC * SS) / (2 * CC * SS)
        L = (Y * CC * SS - C * YC * SS - S * YS * CC) / (CC * SS)
        M = (YC**2 * SS + YS**2 * CC) / (2 * CC * SS)
        log_prob[i:i + 2000] = -0.5 * np.log(CC * SS * np.abs(K)) + M - L**2 / (4 * K)
    return 1 / freqs, log_prob - log_prob.max()
