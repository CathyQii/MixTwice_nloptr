mixtwice <- function(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df) {
  
  theta0 <- thetaHat
  s20 <- s2
  p0 <- length(theta0)
  
  theta <- theta0
  s2 <- s20
  p <- length(theta)
  
  cc <- max(abs(theta)) * 1.1
  if (!is.finite(cc) || cc <= 0) cc <- 1
  
  gridtheta <- (cc / Btheta) * seq(-Btheta, Btheta, by = 1)
  
  sd_min <- sqrt(min(s2))
  sd_max <- sqrt(max(s2))
  if (sd_max == sd_min) {
    gridsigma <- rep(sd_min, Bsigma2)
  } else {
    gridsigma <- seq(sd_min, sd_max, length.out = Bsigma2)
  }
  
  ltheta <- length(gridtheta)
  lsigma <- length(gridsigma)
  
  grid1 <- rep(gridtheta, each = lsigma)
  grid2 <- rep(gridsigma, ltheta)
  
  likelihood.theta <- matrix(NA_real_, nrow = p, ncol = length(grid1))
  likelihood.s2 <- matrix(NA_real_, nrow = p, ncol = length(grid2))
  
  for (i in seq_len(p)) {
    likelihood.theta[i, ] <- stats::dnorm(theta[i], mean = grid1, sd = grid2)
    likelihood.s2[i, ] <- stats::dchisq(df * s2[i] / grid2^2, df = df) * df / grid2^2
  }
  
  likelihood.conditional <- likelihood.theta * likelihood.s2
  
  L <- function(x) {
    xtheta <- x[1:ltheta]
    xsigma <- x[(ltheta + 1):(ltheta + lsigma)]
    
    w <- as.vector(xsigma %o% xtheta)
    dens <- pmax(as.vector(likelihood.conditional %*% w), 1e-300)
    -sum(log(dens))
  }
  
  G <- function(x) {
    g <- numeric(ltheta)
    h <- numeric(lsigma)
    
    xtheta <- x[1:ltheta]
    xsigma <- x[(ltheta + 1):(ltheta + lsigma)]
    
    w <- as.vector(xsigma %o% xtheta)
    d <- pmax(as.vector(likelihood.conditional %*% w), 1e-300)
    
    for (i in seq_len(ltheta)) {
      cols <- (lsigma * (i - 1) + 1):(lsigma * i)
      g[i] <- -sum((likelihood.conditional[, cols, drop = FALSE] %*% xsigma) / d)
    }
    
    for (j in seq_len(lsigma)) {
      cols <- seq(j, j + (ltheta - 1) * lsigma, by = lsigma)
      h[j] <- -sum((likelihood.conditional[, cols, drop = FALSE] %*% xtheta) / d)
    }
    
    c(g, h)
  }
  
  idx0 <- Btheta + 1
  
  hin <- function(x) {
    g <- x[1:ltheta]
    c(
      g[1:Btheta] - g[2:(Btheta + 1)],
      -g[(Btheta + 1):(ltheta - 1)] + g[(Btheta + 2):ltheta]
    )
  }
  
  hinjac <- function(x) {
    J <- matrix(0, nrow = ltheta - 1, ncol = ltheta + lsigma)
    
    for (i in seq_len(Btheta)) {
      J[i, i] <- 1
      J[i, i + 1] <- -1
    }
    
    for (i in (Btheta + 1):(ltheta - 1)) {
      J[i, i] <- -1
      J[i, i + 1] <- 1
    }
    
    J
  }
  
  heq <- function(x) {
    c(
      sum(x[1:ltheta]) - 1,
      sum(x[(ltheta + 1):(ltheta + lsigma)]) - 1
    )
  }
  
  heqjac <- function(x) {
    J <- matrix(0, nrow = 2, ncol = ltheta + lsigma)
    J[1, 1:ltheta] <- 1
    J[2, (ltheta + 1):(ltheta + lsigma)] <- 1
    J
  }
  
  loglik <- log(pmax(likelihood.conditional, 1e-300))
  t1 <- rep(seq_len(ltheta), each = lsigma)
  t2 <- rep(seq_len(lsigma), ltheta)
  
  est.theta <- rep(1 / ltheta, ltheta)
  est.sigma <- rep(1 / lsigma, lsigma)
  
  for (iter in 1:100) {
    est <- as.vector(est.sigma %o% est.theta)
    vv <- t(t(loglik) + log(pmax(est, 1e-300)))
    vv <- exp(vv - apply(vv, 1, max))
    vv <- vv / rowSums(vv)
    
    est <- colMeans(vv)
    est.theta <- as.numeric(tapply(est, t1, sum))
    est.sigma <- as.numeric(tapply(est, t2, sum))
    
    if (requireNamespace("Iso", quietly = TRUE)) {
      tmp <- Iso::ufit(y = est.theta, x = seq_along(est.theta), lmode = idx0)
      est.theta <- tmp$y
    }
    
    est.theta[est.theta < 0] <- 0
    est.sigma[est.sigma < 0] <- 0
    est.theta <- est.theta / sum(est.theta)
    est.sigma <- est.sigma / sum(est.sigma)
  }
  
  x0 <- c(est.theta, est.sigma)
  
  fit <- nloptr::slsqp(
    x0 = x0,
    fn = L,
    gr = G,
    lower = rep(0, ltheta + lsigma),
    upper = rep(1, ltheta + lsigma),
    hin = hin,
    hinjac = hinjac,
    heq = heq,
    heqjac = heqjac,
    control = list(
      xtol_rel = 1e-8,
      ftol_rel = 1e-10,
      maxeval = 3000
    ),
    deprecatedBehavior = FALSE,
    nl.info = FALSE
  )
  
  est.theta <- fit$par[1:ltheta]
  est.theta[est.theta < 0] <- 0
  est.theta <- est.theta / sum(est.theta)
  
  est.sigma <- fit$par[(ltheta + 1):(ltheta + lsigma)]
  est.sigma[est.sigma < 0] <- 0
  est.sigma <- est.sigma / sum(est.sigma)
  
  est.matrix <- outer(est.theta, est.sigma)
  est.array <- as.vector(t(est.matrix))
  
  likelihood.theta <- matrix(NA_real_, nrow = p0, ncol = length(grid1))
  likelihood.s2 <- matrix(NA_real_, nrow = p0, ncol = length(grid2))
  
  for (i in seq_len(p0)) {
    likelihood.theta[i, ] <- stats::dnorm(theta0[i], mean = grid1, sd = grid2)
    likelihood.s2[i, ] <- stats::dchisq(df * s20[i] / grid2^2, df = df) * df / grid2^2
  }
  
  likelihood.conditional <- likelihood.theta * likelihood.s2
  
  post_prob <- array(NA_real_, dim = c(p0, ltheta, lsigma))
  posterior_theta <- matrix(NA_real_, nrow = p0, ncol = ltheta)
  post_mean <- numeric(p0)
  post_sd <- numeric(p0)
  q025 <- numeric(p0)
  q0975 <- numeric(p0)
  
  for (i in seq_len(p0)) {
    ddd <- likelihood.conditional[i, ] * est.array
    ddd <- ddd / sum(ddd)
    
    post_prob[i, , ] <- matrix(ddd, nrow = ltheta, ncol = lsigma, byrow = TRUE)
    
    UUU <- rowSums(post_prob[i, , ])
    posterior_theta[i, ] <- UUU
    
    post_mean[i] <- sum(gridtheta * UUU)
    post_sd[i] <- sqrt(sum((gridtheta^2) * UUU) - post_mean[i]^2)
    
    cdf_theta <- cumsum(UUU)
    q025[i] <- gridtheta[which(cdf_theta >= 0.025)[1]]
    q0975[i] <- gridtheta[which(cdf_theta >= 0.975)[1]]
  }
  
  lfdr <- posterior_theta[, idx0]
  
  lfsr <- numeric(p0)
  for (i in seq_len(p0)) {
    left  <- sum(posterior_theta[i, 1:idx0])
    right <- sum(posterior_theta[i, idx0:ltheta])
    lfsr[i] <- min(left, right)
  }
  
  list(
    grid.theta = gridtheta,
    grid.sigma2 = gridsigma^2,
    mix.theta = est.theta,
    mix.sigma2 = est.sigma,
    post_prob = post_prob,
    posterior_theta = posterior_theta,
    lfdr = lfdr,
    lfsr = lfsr,
    post_mean = post_mean,
    post_sd = post_sd,
    q025 = q025,
    q0975 = q0975
  )
}
