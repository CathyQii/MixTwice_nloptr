library(ashr)
library(limma)
library(MixTwice)
library(ggplot2)
library(patchwork)

# LOAD DATA

obj <- load("peptideRAStats.RData")

thetaHat <- array1_weak$thetaHat
s2       <- array1_weak$s2
n        <- 12
df       <- 2 * n - 2
se_hat   <- sqrt(s2)

Btheta  <- 15
Bsigma2 <- 10

#FIT ALL 6 METHODS
fit_ash <- ash(thetaHat,
               lik = lik_t(df = df),
               se_hat)


sq      <- squeezeVar(var = s2, df = df)
se_post <- sqrt(sq$var.post)
df_post <- df + sq$df.prior
fit_two <- ash(thetaHat, se_post,
               lik    = lik_t(df = df_post),
               method = "fdr")


fit_mix <- mixtwice_v3(thetaHat, s2,
                           Btheta  = Btheta,
                           Bsigma2 = Bsigma2,
                           df      = df)

fit_hu <- ash(thetaHat, se_hat,
              method      = "fdr",
              lik         = lik_t(df = df),
              mixcompdist = "halfuniform"
              )


fit_two_hu <- ash(thetaHat, se_post,
                  lik         = lik_t(df = df_post),
                  method      = "fdr",
                  mixcompdist = "halfuniform")


model <- mixtwice_v2(
  thetaHat = thetaHat,
  s2       = s2,
  Btheta   = Btheta,
  Bsigma2  = Bsigma2,
  df       = df,
  prop     = 0.1,
  "AugLag"
)


# Effect-size CDF
plot(
  model$grid.theta,
  cumsum(model$mix.theta),
  type = "s",
  xlab = "grid of effect size",
  ylab = "cdf of effect size"
)

plot(
  fit_mix$grid.theta,
  cumsum(model$mix.theta),
  type = "s",
  xlab = "grid of effect size",
  ylab = "cdf of effect size"
)


# ASH-t
curve(mixcdf(fit_ash$fitted_g, x),
      from = min(thetaHat),
      to   = max(thetaHat),
      xlab = "effect size",
      ylab = "cdf of effect distribution",
      main = "ASH-t")

# Two-step ASH-t
curve(mixcdf(fit_two$fitted_g, x),
      from = min(thetaHat),
      to   = max(thetaHat),
      main = "Two-step ASH-t")

# ASH-HU
curve(mixcdf(fit_hu$fitted_g, x),
      from = min(thetaHat),
      to   = max(thetaHat),
      main = "ASH-HU")

# Two-step ASH-HU
curve(mixcdf(fit_two_hu$fitted_g, x),
      from = min(thetaHat),
      to   = max(thetaHat),
      main = "Two-step ASH-HU")





#############################

# zoomed-in x-range
xmin <- -0.02
xmax <-  0.08
xgrid <- seq(xmin, xmax, length.out = 500)

# empty plot
plot(NULL,
     xlim = c(xmin, xmax),
     ylim = c(0, 1),
     xlab = "effect size",
     ylab = "CDF of g(theta)",
     main = "Estimated Effect Size Distributions")

# MixTwice (2.0)
lines(model$grid.theta,
      cumsum(model$mix.theta),
      type = "s",
      col = "black",
      lwd = 2)

# MixTwice (3.0)
lines(fit_mix$grid.theta,
      cumsum(fit_mix$mix.theta),
      type = "s",
      col = "gray40",
      lwd = 2,
      lty = 2)

# ASH-t
lines(xgrid,
      mixcdf(fit_ash$fitted_g, xgrid),
      col = "blue",
      lwd = 2)

# Two-step ASH-t
lines(xgrid,
      mixcdf(fit_two$fitted_g, xgrid),
      col = "red",
      lwd = 2)

# ASH-HU
lines(xgrid,
      mixcdf(fit_hu$fitted_g, xgrid),
      col = "darkgreen",
      lwd = 2)

# Two-step ASH-HU
lines(xgrid,
      mixcdf(fit_two_hu$fitted_g, xgrid),
      col = "purple",
      lwd = 2)

# legend
legend("bottomright",
       legend = c("MixTwice (2.0)",
                  "MixTwice (3.0)",
                  "ASH-t",
                  "Two-t",
                  "ASH-t (asymmetry)",
                  "Two-t (asymmetry)"),
       col = c("black", "gray40", "blue", "red", "darkgreen", "purple"),
       lwd = 2,
       lty = c(1, 2, 1, 1, 1, 1),
       cex = 0.9,
       bty = "n")








# Variance CDF
plot(
  model$grid.sigma2,
  cumsum(model$mix.sigma2),
  type = "s",
  xlab = "grid of squared standard error",
  ylab = "cdf of squared standard error"
)








lfdr_list <- list(
  "MixTwice (3.0)" = fit_mix$lfdr,
  "MixTwice (2.0)" = model$lfdr,
  "ASH-t"             = get_lfdr(fit_ash),
  "Two-t"      = get_lfdr(fit_two),
  "ASH-t(asymmetry)"          = get_lfdr(fit_hu),
  "Two-t (asymmetry)"       = get_lfdr(fit_two_hu)
)


clamp01 <- function(x) pmin(pmax(x, 0), 1)


plot_lfdr_real <- function(thetaHat, s2, lfdr_vec, title_txt) {
  ok <- is.finite(thetaHat) & is.finite(s2) &
    is.finite(lfdr_vec) & s2 > 0
  dat <- data.frame(
    s2       = s2[ok],
    thetaHat = thetaHat[ok],
    lfdr     = clamp01(lfdr_vec[ok])
  )
  dat <- dat[order(dat$lfdr, decreasing = TRUE), ]
  ggplot(dat, aes(s2, thetaHat, color = lfdr)) +
    geom_point(size = 0.45, alpha = 0.7) +
    scale_color_gradientn(
      colours = c("#b2182b", "#ef8a62", "#fddbc7",
                  "#d1e5f0", "#67a9cf", "#2166ac"),
      limits  = c(0, 1),
      name    = "lfdr"
    ) +
    theme_classic(base_size = 12) +
    labs(title = title_txt,
         x     = expression(hat(s)^2),
         y     = expression(hat(theta)))
}

panels <- lapply(names(lfdr_list), function(nm) {
  plot_lfdr_real(thetaHat, s2, lfdr_list[[nm]], nm)
})

print(
  wrap_plots(panels, ncol = 3, guides = "collect") &
    theme(legend.position = "right")
)


plot_lfdr_sig <- function(thetaHat, s2, lfdr_vec, title_txt, threshold = 0.1) {
  ok <- is.finite(thetaHat) & is.finite(s2) &
    is.finite(lfdr_vec) & s2 > 0
  dat <- data.frame(
    s2       = s2[ok],
    thetaHat = thetaHat[ok],
    lfdr     = clamp01(lfdr_vec[ok])
  )
  dat_bg  <- dat
  dat_sig <- dat[dat$lfdr <= threshold, ]
  ggplot() +
    geom_point(data = dat_bg,
               aes(x = s2, y = thetaHat),
               color = "grey80", size = 0.4, alpha = 0.5) +
    geom_point(data = dat_sig,
               aes(x = s2, y = thetaHat),
               color = "#2ca25f", size = 1.2, alpha = 0.9) +
    theme_classic(base_size = 12) +
    labs(title    = title_txt,
         subtitle = paste0(nrow(dat_sig), " peptides found"),
         x        = expression(hat(s)^2),
         y        = expression(hat(theta)))
}

panels_sig <- lapply(names(lfdr_list), function(nm) {
  plot_lfdr_sig(thetaHat, s2, lfdr_list[[nm]], nm)
})

print(
  wrap_plots(panels_sig, ncol = 3) +
    plot_annotation(
      title   = "Significant peptides at lfdr \u2264 0.1",
      caption = "Green = lfdr \u2264 0.1   |   Grey = not significant"
    )
)



t_stat <- thetaHat / sqrt(s2)
z <- qnorm(pt(t_stat, df = df))
hist(z, breaks = 100, freq = FALSE)
lines(density(z), col = "red", lwd = 2)
abline(v = 0, col = "blue", lty = 2)
