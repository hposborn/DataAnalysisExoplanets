scatter = ax.scatter(df['pl_bmasse'],df['pl_rade'],
            c=np.log10(df['pl_eqt']),
            s=(15 - df['sy_vmag'].clip(upper=15)) * 15 + 20, cmap='viridis', alpha=0.8, edgecolors='k',linewidth=0.5
        )
cbar = fig.colorbar(scatter, ax=ax)
cbar.set_label(r'$\log_{10}(\text{Teff } [K])$')
ax.set_xlabel(r'Mass $M_{\text{Earth}}$ [$M_\text{Earth}$]')
ax.set_ylabel(r'Planetary Radius $R_\oplus$ [$R_\text{Earth}$]')
ax.set_title('Sub-Neptune Population Overview')
ax.grid(True, linestyle='--', alpha=0.5)
plt.savefig("population_4d.png", dpi=300)
plt.show()

def linear_model(theta, x):
    # Define a linear model from parameter array theta
    a, b = theta
    return a * x + b

def quadratic_model(theta, x):
    # Define a quadratic model from parameter array theta
    a, b, c = theta
    return a * x**2 + b * x + c

def log_likelihood(theta, x, y, yerr, model=linear_model, **kwargs):
    """Gaussian log likelihood."""
    ymodel = model(theta, x)
    log_lik = -0.5 * np.sum(((y - ymodel) / yerr)**2 + np.log(2 * np.pi * yerr**2))
    return log_lik

def BIC(loglik,n_params,nsamps):
    return 2 * loglik + n_params * np.log(nsamps)

def log_prior_linear(theta):
    """
    Calculates log prior using a combination of:
    1. Uniform prior on slope 'a' (must be physical: 0 < a < 5)
    2. Gaussian prior on intercept 'b' centered on 0.5 with sigma=2.5
    """
    a, b = theta

    # 1. Uniform Prior component
    if not (0.0 < a < 5.0):
        return -np.inf  # Reject parameters outside range

    # 2. Gaussian Prior component: log N(b | mu=0.5, sigma=1.0)
    mu_b, sigma_b = 0.5, 5.0
    log_prior_b = -0.5 * ((b - mu_b) / sigma_b)**2 - np.log(sigma_b * np.sqrt(2 * np.pi))

    return log_prior_b
