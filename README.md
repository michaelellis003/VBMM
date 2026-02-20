# VBMM

An R package for Bayesian inference on Gaussian mixture models using Coordinate Ascent Variational Inference (CAVI). Supports both single-level Dirichlet process mixtures (DPGMM) and hierarchical mixtures of Dirichlet processes (MDPGMM).

## Installation

```r
# install.packages("devtools")
devtools::install_github("michaelellis003/VBMM")
```

## Models

**DPGMM** (B = 1) -- Dirichlet process Gaussian mixture model with stick-breaking prior. Performs nonparametric density estimation and clustering.

**MDPGMM** (B > 1) -- Mixture of Dirichlet process Gaussian mixture models. Supports hierarchical clustering with B known groups, each containing an unknown number of sub-clusters sharing the same Gaussian kernels.

## Usage

```r
library(VBMM)

# Simulate mixture data
mix <- simulate_mixture(
    N = 1000, B = 1, K = 3,
    mu = c(-20, 5, 25),
    sigmasq = c(0.5, 0.5, 0.5),
    phi = 1,
    p = list(c(1/3, 1/3, 1/3))
)

# Fit DPGMM with variational inference
output <- vbmm(y = mix$y, B = 1, K = 20, max_iter = 50)

# Diagnostics
plot_elbo(output)           # ELBO convergence
plot_q_mu(output)           # Posterior component means
plot_q_sigmasq(output)      # Posterior component variances
```

### Hierarchical model (B > 1)

```r
mix <- simulate_mixture(
    N = 1000, B = 2, K = 3,
    mu = c(-20, 5, 25),
    sigmasq = c(0.5, 0.5, 0.5),
    phi = c(1/2, 1/2),
    p = list(c(1/3, 1/3, 1/3), c(1/2, 1/4, 1/4))
)

output <- vbmm(y = mix$y, w = mix$w, B = 2, K = 20, max_iter = 50)
```

## Mathematical Details

The full derivation of the variational updates and ELBO is in [`VBMM.pdf`](VBMM.pdf).

## License

MIT
