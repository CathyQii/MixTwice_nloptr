# MixTwice v3.0

## Major Changes

This version introduces an updated constrained optimization
framework and expanded posterior inference outputs.

- replaced the previous `"EM-pava"` and `"AugLag"`
  optimization options with a unified optimization
  framework using EM + PAVA initialization followed by
  sequential least-squares programming (SLSQP) via
  `nloptr::slsqp()`.

- simplified the `mixtwice()` interface by reducing
  user-specified optimization choices under a unified
  estimation framework for stability.

- added posterior inference outputs.

  - the full posterior probability array is returned as
    `post_prob`, where `post_prob[i, k, j]` represents
    \(P(\theta_i = a_k, \sigma_i^2 = b_j^2 \mid x_i, s_i^2)\).

  - the marginal posterior probability matrix over the
    theta grid is returned as `posterior_theta`, where
    `posterior_theta[i, k]` represents
    \(P(\theta_i = a_k \mid x_i, s_i^2)\).

  - approximate posterior credible intervals are returned
    as `q025` and `q0975`.

  - posterior means and posterior standard deviations are
    returned as the numeric vectors `post_mean` and
    `post_sd`.

## Minor Changes

- improved numerical stability for sparse and weak-signal
  data settings.

- improved robustness of variance-grid construction when
  the variance estimates have limited spread.

- updated package documentation.