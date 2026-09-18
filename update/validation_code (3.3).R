library(ashr)
library(MixTwice)
library(dplyr)
library(ggplot2)

set.seed(123)


#set up

p <- 100
pi0 <- 0.8
n_values <- c(2, 4, 5, 10, 15, 20, 50)

# non-null effects magnitude across flat-top, near-normal, spiky
target_nonnull_rms <- 0.7


# Original effect-distribution shapes
g_scenarios_raw <- list(
  
  spiky = normalmix(
    pi   = c(0.4, 0.2, 0.2, 0.2),
    mean = c(0, 0, 0, 0),
    sd   = c(0.25, 0.5, 1, 2)
  ),
  
  near_normal = normalmix(
    pi   = c(2/3, 1/3),
    mean = c(0, 0),
    sd   = c(1, 2)
  ),
  
  flattop = normalmix(
    pi   = rep(1/7, 7),
    mean = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
    sd   = rep(0.5, 7)
  )
)


# Rescale every non-null distribution to the same magnitude
rescale_normalmix <- function(g, target_rms) {
  
  weights <- mixprop(g)
  means   <- comp_mean(g)
  sds     <- comp_sd(g)
  
  current_rms <- sqrt(
    sum(
      weights * (means^2 + sds^2)
    )
  )
  
  scale_factor <- target_rms / current_rms
  
  normalmix(
    pi   = weights,
    mean = means * scale_factor,
    sd   = sds * scale_factor
  )
}


g_scenarios <- lapply(
  g_scenarios_raw,
  rescale_normalmix,
  target_rms = target_nonnull_rms
)


# Verify
scenario_rms <- sapply(
  g_scenarios,
  function(g) {
    
    sqrt(
      sum(
        mixprop(g) *
          (comp_mean(g)^2 + comp_sd(g)^2)
      )
    )
  }
)

print(scenario_rms)


# generation

rtheta <- function(g) {
  
  z <- sample(
    seq_len(ncomp(g)),
    size = p,
    replace = TRUE,
    prob = mixprop(g)
  )
  
  theta <- rnorm(
    p,
    mean = comp_mean(g)[z],
    sd   = comp_sd(g)[z]
  )
  theta[runif(p) < pi0] <- 0
  theta
}



# RMSE = sqrt(mean((posterior mean - true theta)^2))
get_rmse <- function(theta, postmean) {
  
  sqrt(
    mean(
      (postmean - theta)^2
    )
  )
}


#Simulation
results <- list()
k <- 1


for (scenario in names(g_scenarios)) {
  
  theta <- rtheta(
    g_scenarios[[scenario]]
  )
  
  variance_component <- sample(
    1:3,
    size = p,
    replace = TRUE,
    prob = c(0.4, 0.4, 0.2)
  )
  
  sigma2_true <- c(
    0.25,
    1.00,
    4.00
  )[variance_component]
  
  
  z_effect <- rnorm(p)
  
  variance_percentile <- runif(
    p,
    min = 1e-6,
    max = 1 - 1e-6
  )
  
  #loop over n
  for (n in n_values) {
    df <- 2 * n - 2
    se2_true <- 2 * sigma2_true / n
    thetaHat <- theta +
      sqrt(se2_true) * z_effect
    s2 <- se2_true *
      qchisq(
        variance_percentile,
        df = df
      ) / df
    
    s2 <- pmax(
      s2,
      1e-12
    )
    # ASH
    fit_ash <- ash(
      betahat = thetaHat,
      sebetahat = sqrt(s2),
      lik = lik_normal(),
      mixcompdist = "uniform",
      method = "fdr"
    )
    
    rmse_ash <- get_rmse(
      theta = theta,
      postmean = fit_ash$result$PosteriorMean
    )
    
    results[[k]] <- data.frame(
      scenario = scenario,
      n = n,
      df = df,
      method = "ASH",
      RMSE = rmse_ash
    )
    
    k <- k + 1
    
    # Two-step ASH
    sq <- limma::squeezeVar(
      var = s2,
      df = df
    )
    
    fit_two <- ash(
      betahat = thetaHat,
      sebetahat = sqrt(sq$var.post),
      lik = lik_t(
        df = df + sq$df.prior
      ),
      mixcompdist = "uniform",
      method = "fdr"
    )
    
    rmse_two <- get_rmse(
      theta = theta,
      postmean = fit_two$result$PosteriorMean
    )
    
    results[[k]] <- data.frame(
      scenario = scenario,
      n = n,
      df = df,
      method = "Two-step ASH",
      RMSE = rmse_two
    )
    
    k <- k + 1
    
    # MixTwice
    
    theta_boundary <- 1.1 *
      max(abs(thetaHat))
    
    target_grid_step <- 0.25 *
      median(sqrt(s2))
    
    Btheta_use <- max(
      60L,
      ceiling(
        theta_boundary /
          target_grid_step
      )
    )
    
    fit_mt <- mixtwice_v3(
      thetaHat = thetaHat,
      s2 = s2,
      Btheta = Btheta_use,
      Bsigma2 = 50,
      df = df
    )
    
    rmse_mt <- get_rmse(
      theta = theta,
      postmean = fit_mt$post_mean
    )
    
    results[[k]] <- data.frame(
      scenario = scenario,
      n = n,
      df = df,
      method = "MixTwice",
      RMSE = rmse_mt
    )
    
    k <- k + 1
  }
}


results <- bind_rows(results)



results <- results |>
  mutate(
    
    scenario = factor(
      scenario,
      levels = c(
        "flattop",
        "near_normal",
        "spiky"
      ),
      labels = c(
        "Flat-top",
        "Near-normal",
        "Spiky"
      )
    ),
    
    method = factor(
      method,
      levels = c(
        "ASH",
        "Two-step ASH",
        "MixTwice"
      )
    )
  )


print(
  results |>
    arrange(
      scenario,
      n,
      method
    )
)

# plot 1
p_rmse <- ggplot(
  results,
  aes(
    x = df,
    y = RMSE,
    color = method,
    linetype = method,
    shape = method,
    group = method
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.7
  ) +
  facet_wrap(
    ~ scenario,
    nrow = 1
  ) +
  scale_x_log10(
    breaks = c(2, 6, 18, 38, 98)
  ) +
  labs(
    x = "Residual degrees of freedom",
    y = "Posterior mean RMSE",
    color = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "top",
    panel.spacing = unit(
      1,
      "lines"
    )
  )

p_rmse

# plot 2
rmse_agreement <- results |>
  group_by(
    scenario,
    df
  ) |>
  summarise(
    RMSE_range =
      max(RMSE) - min(RMSE),
    .groups = "drop"
  )


print(rmse_agreement)


p_agreement <- ggplot(
  rmse_agreement,
  aes(
    x = df,
    y = RMSE_range
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.7
  ) +
  facet_wrap(
    ~ scenario,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_x_log10(
    breaks = c(2, 6, 18, 38, 98)
  ) +
  labs(
    x = "Residual degrees of freedom",
    y = "Range of RMSE across methods"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    panel.spacing = unit(
      1,
      "lines"
    )
  )

p_agreement