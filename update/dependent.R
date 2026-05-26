library(ashr)
library(limma)
library(MixTwice)
library(ggplot2)
library(patchwork)
l   <- 10000
pi0 <- 0.3
n   <- 10
df  <- 2 * (n - 1)
Btheta  <- 15
Bsigma2 <- 10
cuts <- seq(0.01, 0.5, by = 0.01)

# sigma_i^2
z <- rbinom(l, 1, 0.88)
sigma2 <- numeric(l)
sigma2[z == 1] <- 0.11 * rbeta(sum(z == 1), 0.85, 6.5)
sigma2[z == 0] <- 0.11 * rbeta(sum(z == 0), 2.0, 4.2)

# gamma_i = 1: non-null
gamma <- rbinom(l, 1, 1 - pi0)

# theta_i depends on sigma_i^2
theta <- numeric(l)
idx <- which(gamma == 1)
theta[idx] <- rnorm(length(idx), mean = 0, sd = sqrt(0.5 * sigma2[idx]))

thetaHat <- rnorm(l, mean = theta, sd = sqrt(sigma2))

# s2_i | sigma_i^2
s2 <- (2 * sigma2 / n) * rchisq(l, df = df) / df
se_hat <- sqrt(s2)

### fit 6 methods
fit_ash <- ash(
  thetaHat, se_hat,
  lik    = lik_t(df = df),
  method = "fdr"
)

sq <- squeezeVar(var = s2, df = df)
se_post <- sqrt(sq$var.post)
df_post <- df + sq$df.prior

fit_two <- ash(
  thetaHat, se_post,
  lik    = lik_t(df = df_post),
  method = "fdr"
)

fit_mix <- mixtwice_v3(
  thetaHat, s2,
  Btheta  = Btheta,
  Bsigma2 = Bsigma2,
  df      = df
)

fit_hu <- ash(
  thetaHat, se_hat,
  method      = "fdr",
  lik         = lik_t(df = df),
  mixcompdist = "halfuniform"
)

fit_two_hu <- ash(
  thetaHat, se_post,
  lik         = lik_t(df = df_post),
  method      = "fdr",
  mixcompdist = "halfuniform"
)

model <- mixtwice_v2(
  thetaHat = thetaHat,
  s2       = s2,
  Btheta   = Btheta,
  Bsigma2  = Bsigma2,
  df       = df,
  prop     = 0.1,
  "AugLag"
)

lfdr_list <- list(
  "MixTwice (3.0)"      = fit_mix$lfdr,
  "MixTwice (2.0)"     = model$lfdr,
  "ASH-t"               = get_lfdr(fit_ash),
  "Two-t"               = get_lfdr(fit_two),
  "ASH-t(asymmetry)"    = get_lfdr(fit_hu),
  "Two-t (asymmetry)"   = get_lfdr(fit_two_hu)
)

# true null indicator
true_null <- (theta == 0)

# compute FP, TP, total discoveries, FDP
fdp_tab <- do.call(rbind, lapply(names(lfdr_list), function(m) {
  
  lfdr <- lfdr_list[[m]]
  
  FP <- sapply(cuts, function(cut) {
    sum(true_null & (lfdr <= cut), na.rm = TRUE)
  })
  
  TP <- sapply(cuts, function(cut) {
    sum(!true_null & (lfdr <= cut), na.rm = TRUE)
  })
  
  R <- sapply(cuts, function(cut) {
    sum(lfdr <= cut, na.rm = TRUE)
  })
  
  data.frame(
    setup = "dependent theta-sigma setup",
    pi0 = pi0,
    n = n,
    method = m,
    cutoff = cuts,
    FP = FP,
    TP = TP,
    total_discoveries = R,
    FDP = ifelse(R == 0, NA, FP / R)
  )
}))

# FDP plot
p_fdp <- ggplot(fdp_tab,
                aes(x = cutoff,
                    y = FDP,
                    color = method)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("Dependent setup: pi0 = ", pi0, ", n = ", n),
    x = "fdr",
    y = "fdp"
  )

p_fdp
