\name{mixtwice}
\alias{mixtwice}
\title{Large-scale hypothesis testing by variance mixing}
\description{
\code{mixtwice} performs large-scale hypothesis testing when each testing unit
provides an estimated effect and its associated squared standard error. It fits
a nonparametric empirical-Bayes mixture model on both the underlying effect sizes
and the underlying variances, and returns posterior summaries including local
false discovery rates and local false sign rates.
}
\usage{
mixtwice(thetaHat, s2, Btheta, Bsigma2, df)
}
\arguments{
  \item{thetaHat}{
    Numeric vector of estimated effect sizes across testing units.
  }
  \item{s2}{
    Numeric vector of estimated squared standard errors corresponding to
    \code{thetaHat}.
  }
  \item{Btheta}{
    Integer controlling the effect-size grid. The procedure constructs a
    symmetric grid with \code{2 * Btheta + 1} support points centered at 0.
  }
  \item{Bsigma2}{
    Integer controlling the variance grid. The procedure constructs a positive
    grid with \code{Bsigma2} support points for the underlying variances.
  }
  \item{df}{
    Degrees of freedom in the chi-squared model for the observed squared
    standard errors.
  }
}
\details{
\code{mixtwice} takes estimated effects and squared standard errors over a set
of testing units and fits a two-layer nonparametric empirical-Bayes model.

For the effect-size distribution, the method constructs a symmetric finite grid
centered at 0, with boundary given by \code{1.1 * max(abs(thetaHat))}. For the
variance distribution, it constructs a finite grid spanning the observed range
of standard deviations and then squares the grid values to obtain variances.

The mixing distribution on effect sizes is constrained to be unimodal with mode
at 0, while the mixing distribution on variances is constrained only to be a
valid probability distribution. Estimation is carried out using a short
EM-PAVA warm start followed by constrained likelihood optimization using SLSQP
via the \code{nloptr} package.

Posterior inference is then computed for each testing unit, including posterior
probabilities on the effect-size and variance grids, posterior means, posterior
standard deviations, credible interval endpoints, local false discovery rates
(\code{lfdr}), and local false sign rates (\code{lfsr}).
}
\value{
A list with components:
\item{grid.theta}{
  Support grid for the effect-size mixing distribution.
}
\item{grid.sigma2}{
  Support grid for the variance mixing distribution.
}
\item{mix.theta}{
  Estimated mixing distribution on \code{grid.theta}.
}
\item{mix.sigma2}{
  Estimated mixing distribution on \code{grid.sigma2}.
}
\item{post_prob}{
  Posterior probabilities over the joint effect-size by variance grid for each
  testing unit.
}
\item{LFDR}{
  Posterior probabilities over the effect-size grid for each testing unit.
}
\item{lfdr}{
  Local false discovery rate for each testing unit, defined as the posterior
  probability of the null grid point.
}
\item{lfsr}{
  Local false sign rate for each testing unit.
}
\item{post_mean}{
  Posterior mean of the effect size for each testing unit.
}
\item{post_sd}{
  Posterior standard deviation of the effect size for each testing unit.
}
\item{q025}{
  Lower endpoint of the 95\% posterior credible interval for each testing unit.
}
\item{q0975}{
  Upper endpoint of the 95\% posterior credible interval for each testing unit.
}
\item{fit}{
  Optimization output returned by \code{nloptr::slsqp()}.
}
}
\references{
Zheng et al. \emph{MixTwice: Large scale hypothesis testing for peptide arrays
by variance mixing}. \bold{Bioinformatics}, 2021.
}
\author{
Zihao Zheng, Michael A. Newton
}
\note{
This updated version uses a grid-based empirical-Bayes model with constrained
optimization and posterior inference for effect sizes and variances.
}
\seealso{
\code{\link{peptide_data}}
}
\examples{
set.seed(1)

l <- 1000
neach <- 20
pi0 <- 0.8

signal1 <- rep(0, l)
signal2 <- signal1
signal2[1:round((1 - pi0) * l)] <- rnorm(round((1 - pi0) * l), mean = 0, sd = 3)

sigma2 <- rep(1, l)
sigma <- sqrt(sigma2)

data1 <- data2 <- matrix(NA, nrow = l, ncol = neach)

for (i in 1:l) {
  data1[i, ] <- rnorm(neach, mean = signal1[i], sd = sigma[i])
  data2[i, ] <- rnorm(neach, mean = signal2[i], sd = sigma[i])
}

thetaHat <- rowMeans(data2) - rowMeans(data1)

sd1 <- apply(data1, 1, sd)
sd2 <- apply(data2, 1, sd)

s2 <- sd1^2 / neach + sd2^2 / neach

fit <- mixtwice(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df = 2 * neach - 2)

plot(fit$grid.theta, cumsum(fit$mix.theta), type = "s",
     xlab = "grid.theta", ylab = "ecdf of theta", lwd = 2)

lines(ecdf(signal2 - signal1), lwd = 0.5, lty = 2, col = "red")

legend("topleft", legend = c("fit.mix", "true.mix"),
       col = c("black", "red"), lwd = 1, pch = 19)

oo <- order(fit$lfdr)

x1 <- sum(cumsum(fit$lfdr[oo]) / (1:l) <= 0.05)
x2 <- sum(cumsum(fit$lfdr[oo]) / (1:l) <= 0.05 & (signal2 != 0)[oo])
x3 <- sum(signal2 != 0)
x4 <- sum(cumsum(fit$lfdr[oo]) / (1:l) <= 0.05 & (signal2 == 0)[oo])

x4 / x1
x2 / x3

max(fit$mix.theta)
}
\keyword{High density peptide array}
\keyword{Large-scale hypothesis testing}
\keyword{Mixing distribution}
\keyword{Local false discovery rate}