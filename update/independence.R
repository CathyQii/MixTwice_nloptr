library(ashr)
library(limma)
library(MixTwice)
library(ggplot2)
library(patchwork)

nsamp   <- 10000
pi0     <- 0.3
n       <- 10
df      <- 2 * n - 2
Btheta  <- 15
Bsigma2 <- 10
cuts    <- seq(0.1, 0.7, by = 0.1)

g_alt_list <- list(
  "spiky" = normalmix(
    pi   = c(0.20, 0.30, 0.30, 0.20),
    mean = c(0, 0, 0, 0),
    sd   = c(0.30, 0.80, 1.50, 3.00)
  ),
  
  "near-normal" = normalmix(
    pi   = c(0.55, 0.30, 0.15),
    mean = c(0, 0, 0),
    sd   = c(0.50, 0.90, 1.40)
  ),
  
  "flattop" = normalmix(
    pi   = rep(1/7, 7),
    mean = c(-1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2),
    sd   = c(0.30, 0.45, 0.60, 0.70, 0.60, 0.45, 0.30)
  ),
  
  "skew" = normalmix(
    pi   = c(0.30, 0.30, 0.25, 0.15),
    mean = c(0, 0.8, 1.8, 3.2),
    sd   = c(0.40, 0.80, 1.40, 2.20)
  ),
  
  "big-normal" = normalmix(
    pi   = c(0.35, 0.40, 0.25),
    mean = c(0, 0, 0),
    sd   = c(1.00, 2.20, 3.50)
  ),
  
  "bimodal" = normalmix(
    pi   = c(0.25, 0.25, 0.25, 0.25),
    mean = c(-2.5, -2.5, 2.5, 2.5),
    sd   = c(0.50, 1.50, 0.50, 1.50)
  )
)

all_results <- list()
plot_list <- list()

for (nm in names(g_alt_list)) {
  
  g_alt <- g_alt_list[[nm]]
  
  is_null <- rbinom(nsamp, 1, pi0)
  comp <- sample(seq_along(mixprop(g_alt)),
                 nsamp,
                 replace = TRUE,
                 prob = mixprop(g_alt))
  
  theta <- ifelse(
    is_null == 1,
    0,
    rnorm(nsamp, comp_mean(g_alt)[comp], comp_sd(g_alt)[comp])
  )
  
  sigma2 <- 1 / rgamma(nsamp, shape = 5, rate = 5)
  
  thetaHat <- rnorm(nsamp, mean = theta, sd = sqrt(sigma2))
  s2 <- sigma2 * rchisq(nsamp, df = df) / df
  se_hat <- sqrt(s2)
  
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
    "ASH-t"               = fit_ash$result$lfdr,
    "Two-t"               = fit_two$result$lfdr,
    "ASH-t(asymmetry)"    = fit_hu$result$lfdr,
    "Two-t (asymmetry)"   = fit_two_hu$result$lfdr
  )
  
  true_null <- (theta == 0)
  
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
      g_alt = nm,
      n = n,
      method = m,
      cutoff = cuts,
      FP = FP,
      TP = TP,
      total_discoveries = R,
      FDP = ifelse(R == 0, NA, FP / R)
    )
  }))
  
  all_results[[nm]] <- fdp_tab
  
  plot_list[[nm]] <- ggplot(fdp_tab,
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
      title = paste0(nm, " (n = ", n, ")"),
      x = "fdr",
      y = "fdp"
    )
}

all_results_df <- do.call(rbind, all_results)

wrap_plots(plot_list, ncol = 2)

t_stat <- thetaHat / sqrt(s2)
z <- qnorm(pt(t_stat, df = df))
