######## Ash-t not sensitive to the structure
library(ashr)
library(limma)
library(MixTwice)
library(ggplot2)
library(patchwork)

set.seed(1)

l   <- 10000
pi0 <- 0.8       # 20% non-null
n   <- 10
df  <- 2 * (n - 1)
Btheta  <- 15
Bsigma2 <- 10
cut <- 0.1

# underlying variance
z <- rbinom(l, 1, 0.88)
sigma2_raw <- numeric(l)
sigma2_raw[z == 1] <- 0.11 * rbeta(sum(z == 1), 0.85, 6.5)
sigma2_raw[z == 0] <- 0.11 * rbeta(sum(z == 0), 2.0, 4.2)
sigma2_raw <- pmax(sigma2_raw, 1e-8)

# variance of thetaHat
r2 <- 2 * sigma2_raw / n

# non-nulls concentrated among low variance peptides
gamma <- rep(0, l)
low_var_idx <- order(r2)[1:round((1 - pi0) * l)]
gamma[low_var_idx] <- 1

# weak positive effects in low-variance region
theta <- numeric(l)
theta[gamma == 1] <- rnorm(
  sum(gamma == 1),
  mean = 0.055,
  sd   = sqrt(0.10 * r2[gamma == 1])
)

# observed effect and estimated SE
thetaHat <- rnorm(l, mean = theta, sd = sqrt(r2))
s2 <- r2 * rchisq(l, df = df) / df
se_hat <- sqrt(s2)

# fit 6 methods
fit_ash <- ash(thetaHat, se_hat, lik = lik_t(df = df), method = "fdr")

sq <- squeezeVar(var = s2, df = df)
se_post <- sqrt(sq$var.post)
df_post <- df + sq$df.prior

fit_two <- ash(thetaHat, se_post, lik = lik_t(df = df_post), method = "fdr")

fit_mix <- mixtwice_v3(
  thetaHat, s2,
  Btheta = Btheta,
  Bsigma2 = Bsigma2,
  df = df
)

fit_hu <- ash(
  thetaHat, se_hat,
  method = "fdr",
  lik = lik_t(df = df),
  mixcompdist = "halfuniform"
)

fit_two_hu <- ash(
  thetaHat, se_post,
  lik = lik_t(df = df_post),
  method = "fdr",
  mixcompdist = "halfuniform"
)

model <- mixtwice_v2(
  thetaHat = thetaHat,
  s2 = s2,
  Btheta = Btheta,
  Bsigma2 = Bsigma2,
  df = df,
  prop = 0.1,
  "AugLag"
)

lfdr_list <- list(
  "MixTwice (3.0)"      = fit_mix$lfdr,
  "MixTwice (2.0)"     = model$lfdr,
  "ASH-t"               = get_lfdr(fit_ash),
  "Two-step ASH"        = get_lfdr(fit_two),
  "ASH-t HU"            = get_lfdr(fit_hu),
  "Two-step HU"         = get_lfdr(fit_two_hu)
)

# selected-peptide scatterplot
plot_df <- do.call(rbind, lapply(names(lfdr_list), function(m) {
  lfdr <- lfdr_list[[m]]
  data.frame(
    method = m,
    thetaHat = thetaHat,
    s2 = s2,
    lfdr = lfdr,
    selected = lfdr <= cut
  )
}))

counts <- aggregate(selected ~ method, plot_df, sum)
names(counts)[2] <- "n_found"
plot_df <- merge(plot_df, counts, by = "method")
plot_df$label <- paste0(plot_df$method, "\n", plot_df$n_found, " peptides found")

ggplot(plot_df, aes(x = s2, y = thetaHat)) +
  geom_point(color = "grey80", size = 0.35, alpha = 0.5) +
  geom_point(
    data = subset(plot_df, selected),
    color = "#2ca25f", size = 0.7, alpha = 0.9
  ) +
  facet_wrap(~ label, ncol = 3) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("Significant peptides at lfdr <= ", cut),
    x = expression(hat(s)^2),
    y = expression(hat(theta))
  )

# FDP plot
cuts <- seq(0.01, 0.5, by = 0.01)
true_null <- (gamma == 0)

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
    setup = "dependent low-variance weak positive effects",
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

p_fdp <- ggplot(fdp_tab,
                aes(x = cutoff,
                    y = FDP,
                    color = method)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("Dependent setup: pi0 = ", pi0, ", n = ", n),
    x = "fdr",
    y = "fdp"
  )

p_fdp



###### 

######## mimic weak-signal real data:
######## global z-shift + subtle index shift
######## asymmetry ASH not sensitive to the structure

library(ashr)
library(limma)
library(MixTwice)
library(ggplot2)
library(patchwork)

set.seed(1)

l  <- 150000
n  <- 10
df <- 2 * n - 2

Btheta  <- 15
Bsigma2 <- 10
cut <- 0.1

group <- rep(0, l)
shift_start <- 100001
group[shift_start:l] <- 1

gamma <- as.integer(group == 1)

sigma2_raw <- numeric(l)

sigma2_raw[group == 0] <- 0.08 * rbeta(sum(group == 0), 0.9, 6.5)
sigma2_raw[group == 1] <- 0.08 * rbeta(sum(group == 1), 1.2, 5.8)

sigma2_raw <- pmax(sigma2_raw, 1e-8)
r2 <- 2 * sigma2_raw / n


theta <- numeric(l)

# shift
theta[group == 0] <- rnorm(
  sum(group == 0),
  mean = 0.60 * sqrt(r2[group == 0]),
  sd   = 0.03 * sqrt(r2[group == 0])
)

# shifted block
theta[group == 1] <- rnorm(
  sum(group == 1),
  mean = 0.90 * sqrt(r2[group == 1]),
  sd   = 0.04 * sqrt(r2[group == 1])
)

# Observed effect and estimated SE
thetaHat <- rnorm(l, mean = theta, sd = sqrt(r2))

s2 <- r2 * rchisq(l, df = df) / df
s2 <- pmax(s2, 1e-10)
se_hat <- sqrt(s2)

zscore <- thetaHat / se_hat

# plot z 
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

plot(
  zscore,
  pch = 1,
  cex = 0.45,
  xlab = "Index",
  ylab = "z",
  ylim = c(-4, 5),
  main = ""
)

abline(v = shift_start, col = "red", lty = 2, lwd = 2)

hist(
  zscore,
  breaks = 90,
  probability = TRUE,
  col = "grey80",
  border = "grey40",
  xlim = c(-4, 4),
  ylim = c(0, 0.43),
  main = "Histogram of z",
  xlab = "z",
  ylab = "Density"
)

curve(
  dnorm(x, mean = mean(zscore), sd = sd(zscore)),
  add = TRUE,
  col = "red",
  lwd = 2
)

abline(v = 0, col = "royalblue", lty = 2, lwd = 2)

par(mfrow = c(1, 1))

cat("\nMean z by block:\n")
print(tapply(zscore, group, mean))

fit_ash <- ash(
  betahat = thetaHat,
  sebetahat = se_hat,
  lik = lik_t(df = df),
  method = "fdr"
)

sq <- limma::squeezeVar(var = s2, df = df)

se_post <- sqrt(as.numeric(sq$var.post))

df_post <- df + sq$df.prior
df_post <- median(as.numeric(df_post), na.rm = TRUE)

fit_two <- ash(
  betahat = thetaHat,
  sebetahat = se_post,
  lik = lik_t(df = df_post),
  method = "fdr"
)


fit_hu <- ash(
  betahat = thetaHat,
  sebetahat = se_hat,
  method = "fdr",
  lik = lik_t(df = df),
  mixcompdist = "halfuniform"
)


fit_two_hu <- ash(
  betahat = thetaHat,
  sebetahat = se_post,
  lik = lik_t(df = df_post),
  method = "fdr",
  mixcompdist = "halfuniform"
)


fit_mix <- mixtwice_v3(
  thetaHat,
  s2,
  Btheta = Btheta,
  Bsigma2 = Bsigma2,
  df = df
)

model <- mixtwice_v2(
  thetaHat = thetaHat,
  s2 = s2,
  Btheta = Btheta,
  Bsigma2 = Bsigma2,
  df = df,
  prop = 0.1,
  "AugLag"
)


lfdr_list <- list(
  "ASH-t"           = get_lfdr(fit_ash),
  "ASH-t HU"        = get_lfdr(fit_hu),
  "MixTwice (3.0)"  = fit_mix$lfdr,
  "MixTwice (2.0)" = model$lfdr,
  "Two-step ASH"    = get_lfdr(fit_two),
  "Two-step HU"     = get_lfdr(fit_two_hu)
)


result <- data.frame(
  method = names(lfdr_list),
  n_found = sapply(lfdr_list, function(x) sum(x <= cut, na.rm = TRUE)),
  true_positive = sapply(lfdr_list, function(x) sum(x <= cut & gamma == 1, na.rm = TRUE)),
  false_positive = sapply(lfdr_list, function(x) sum(x <= cut & gamma == 0, na.rm = TRUE))
)

print(result)



plot_df <- do.call(rbind, lapply(names(lfdr_list), function(m) {
  lfdr <- lfdr_list[[m]]
  data.frame(
    method   = m,
    thetaHat = thetaHat,
    s2       = s2,
    lfdr     = lfdr,
    selected = lfdr <= cut
  )
}))

# count discoveries
counts <- aggregate(selected ~ method, plot_df, sum)
names(counts)[2] <- "n_found"

plot_df <- merge(plot_df, counts, by = "method")

plot_df$label <- paste0(plot_df$method, "\n", plot_df$n_found, " peptides found")


ggplot(plot_df, aes(x = s2, y = thetaHat)) +
  geom_point(color = "grey80", size = 0.35, alpha = 0.5) +
  geom_point(
    data = subset(plot_df, selected),
    color = "#2ca25f",
    size = 0.7,
    alpha = 0.9
  ) +
  facet_wrap(~ label, ncol = 3) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("Significant peptides at lfdr <= ", cut),
    x = expression(hat(s)^2),
    y = expression(hat(theta))
  )


# FDP

cuts <- seq(0.01, 0.5, by = 0.01)

# true null indicator
true_null <- (gamma == 0)

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
    setup = "dependent low-variance weak positive effects",
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

p_fdp <- ggplot(fdp_tab,
                aes(x = cutoff,
                    y = FDP,
                    color = method)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("Dependent setup: pi0 = ", pi0, ", n = ", n),
    x = "fdr",
    y = "fdp"
  )

p_fdp


