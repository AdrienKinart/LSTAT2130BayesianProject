# ==================
#       SETUP
# ==================

# Require the necessary packages

if (!require(mvtnorm)) {install.packages("mvtnorm");require(mvtnorm)}
if (!require(EnvStats)) {install.packages("EnvStats");require(EnvStats)}
if (!require(R2WinBUGS)) {install.packages("R2WinBUGS");require(R2WinBUGS)}
if (!require(coda)) {install.packages("coda");require(coda)}
if (!require(rjags)) {install.packages("rjags");require(rjags)}

# Setup of the table

data <- matrix(
  data = c(
    25, 69, 65, 106, 80, 106, 136, 94, 76, 46,
    17, 36, 47, 58, 47, 53, 59, 54, 33, 21
  ),
  nrow = 2, byrow = T
)

rownames(data) <- c("Flanders", "Wallonia")
colnames(data) <- c(
  "<1200", "[1200-1500)", "[1500-1800)", "[1800-2300)", "[2300-2700)",
  "[2700-3300)", "[3300-4000)", "[4000-4900)", "[4900-6000)", ">=6000"
)

flanders <- data[1,]; wallonia <- data[2,]
flanders.n <- sum(flanders) ; wallonia.n <- sum(wallonia)
flanders.prop <- prop.table(flanders) ; wallonia.prop <- prop.table(wallonia)

# Intervals of the frequency table
intervals = c(0, 1200,1500,1800,2300,2700,3300,4000,4900,6000, Inf)
intervals.mean = c(600,1350,1650,2050,2500,3000,3650,4450,5450,6000)

# Utility function to compute the prob of being between high and low
pgammadiff = function(low, high, kappa, lambda){
  pgamma(high, kappa, lambda) - pgamma(low, kappa, lambda)
}

# =====

alpha = 0.05

# Utility fonctions
kappa = function(phi){1/phi}
lambda = function(phi, mu) {1/(phi*mu)}

# Estimate mu and phi

obtain_the_estimated_parameters <- function(row, boot = F, simulations = 10000) {
  if (boot) {
    sim_prop <- rowSums(rmultinom(simulations, sum(data), prop.table(data)[row, ]))
  } else {
    sim_prop <- data[row, ]
  }
  
  sim_data <- unlist(sapply(1:ncol(data), function(i) {
    rep(intervals.mean[i], sim_prop[i])
  }))
  est_gam <- egamma(sim_data)
  
  est_kappa <- est_gam$parameters["shape"]
  est_lambda <- 1 / est_gam$parameters["scale"]
  
  est_mu <- est_kappa / est_lambda
  est_phi <- 1 / est_kappa
  
  return(list(mu = est_mu, phi = est_phi, kappa = est_kappa, lambda = est_lambda, sim = sim_data))
}

estim_params_fl <- obtain_the_estimated_parameters(1, T)
estim_params_wal <- obtain_the_estimated_parameters(2, T)

estim_mu_fl <- estim_params_fl$mu
estim_phi_fl <- estim_params_fl$phi
estim_fl <- estim_params_fl$sim
estim_mu_wal <- estim_params_wal$mu
estim_phi_wal <- estim_params_wal$phi
estim_wal <- estim_params_wal$sim

# ===================
# QUESTION 2 : PRIORS
# ===================

# Set the priors
mu_prior <- 3000
sigma_prior <- 306.12

# ===============================
# QUESTION 3 : lpost(theta, freq)
# ===============================

lpost <- function(theta, freq) {
  # Transform mu and phi -> kappa and lambda
  kappa <- kappa(theta[2])
  lambda <- lambda(theta[2], theta[1])
  
  # Likelihood : sum_j^n[x_j * ln(probability_of_being_in_interval_j)]
  LL <- sum(sapply(1:length(freq), function(j) {
    freq[j] * log(pgammadiff(low = intervals[j], high = intervals[j + 1], kappa, lambda))
  }))
  
  # Log posterior: LL + priors
  lpi <- dnorm(theta[1], mu_prior, sigma_prior, log = T) + dunif(theta[2], 0, 10, log = T)
  lpost <- LL + lpi
  
  names(lpost) <- "lpost"
  return(lpost)
}

# ====================
# QUESTION 4 : Laplace
# ====================

# Starting values
inits <- c(mu = mu_prior, phi = 0.01)

# Laplace approximation
laplace_approx <- function(inits, frequencies) {
  # Try to maximize lpost and retrieve the hessian
  fit <- optim(inits, lpost, control = list(fnscale = -1), hessian = T, freq = frequencies)
  params <- fit$par
  param_cov_mat <- solve(-fit$hessian)
  # Multivariate normal
  samples <- rmvnorm(50000, params, param_cov_mat)
  
  list(samples = samples, parameters = params, cov = param_cov_mat)
}

laplace_fl <- laplace_approx(inits, flanders)
laplace_fl.mcmc <- mcmc(laplace_fl$samples)

# Credible intervals of Laplace (HPD vs QB) for the mean
marginal_post_Fl <- rnorm(1e6, laplace_fl$parameters["mu"], sqrt(laplace_fl$cov[1, 1]))

QB_LP <- quantile(marginal_post_Fl, probs = c(alpha / 2, 1 - alpha / 2))
HPD_LP <- HPDinterval(as.mcmc(marginal_post_Fl), prob = 1 - alpha)

# =======================
# QUESTION 5 : Metropolis
# =======================

metropolis_algorithm <- function(M, theta, sd.propositions, factors, frequencies) {
  # Make sure args are vectors
  theta <- as.vector(theta)
  sd.propositions <- as.vector(sd.propositions)
  factors <- as.vector(factors)
  frequencies <- as.vector(frequencies)
  
  # Init an array of dim (M, 2) with theta as first elem
  thetas <- array(dim = c(M + 1, 2))
  thetas[1, ] <- theta
  
  accepted <- c(0, 0)
  
  # variance adjustment via a factor
  sigma_mu <- factors[1] * sd.propositions[1]
  sigma_phi <- factors[2] * sd.propositions[2]
  
  for (i in 2:(M + 1)) {
    # Retrieve the previous value
    theta_mu <- theta_phi <- this_theta <- thetas[i - 1, ]
    
    # Proposition of mu and phi, "independently"
    theta_mu[1] <- this_theta[1] + rnorm(1, 0, sigma_mu)
    theta_phi[2] <- this_theta[2] + rnorm(1, 0, sigma_phi)
    
    thresholds <- c(
      min(exp(lpost(theta_mu, frequencies) - lpost(this_theta, frequencies)), 1),
      min(exp(lpost(theta_phi, frequencies) - lpost(this_theta, frequencies)), 1)
    )
    
    accepts <- runif(2) <= thresholds
    
    # Update if accepted
    # otherwise use the previous values
    thetas[i, ] <- this_theta
    if (accepts[1]) thetas[i, 1] <- theta_mu[1]
    if (accepts[2]) thetas[i, 2] <- theta_phi[2]
    accepted <- accepted + as.integer(accepts)
  }
  
  # We decided to do the burnin outside
  
  colnames(thetas) <- c("mu", "phi")
  names(accepted) <- c("mu", "phi")
  return(list(theta = thetas, accept_rate = accepted / M))
}

burnin <- function(array, percent) {
  tail(array, -percent * nrow(array))
}

# Initial propositions using laplace
sd_hat_mu <- sqrt(laplace_fl$cov[1, 1])
sd_hat_phi <- sqrt(laplace_fl$cov[2, 2])

# =====
# Run the algorithm (take some time)
metro_fl1 <- metropolis_algorithm(5e4,
  theta = c(mu = estim_mu_fl, phi = estim_phi_fl),
  sd.propositions = c(
    mu = sd_hat_mu,
    phi = sd_hat_phi
  ),
  factors = c(mu = 2.75, phi = 2.75),
  frequencies = flanders
)
# metro_fl1$accept_rate

# Remove burnin samples (10%)
metro_fl1.theta <- burnin(metro_fl1$theta, 0.1)
metro_fl1.theta.mcmc <- as.mcmc(metro_fl1.theta)

# =====
# Convergence analysis
traceplot(metro_fl1.theta.mcmc)
acf(metro_fl1.theta[, 1]); pacf(metro_fl1.theta[, 1])
acf(metro_fl1.theta[, 2]); pacf(metro_fl1.theta[, 2])
effectiveSize(metro_fl1.theta.mcmc)

# Gelman rubin
# Run another metropolis with different initial theta
metro_fl2 <- metropolis_algorithm(50000,
  theta = c(mu = 4000, phi = 0.1),
  sd.propositions = c(
    mu = sd_hat_mu,
    phi = sd_hat_phi
  ),
  factors = c(mu = 2.75, phi = 2.75),
  frequencies = flanders
)

# Don't apply burnin as we want to see the convergence
metro_fl2.theta <- metro_fl2$theta 
metro_fl2.theta.mcmc <- as.mcmc(metro_fl2.theta)

traceplot(list(
  head(as.mcmc(metro_fl1$theta), 1000),
  head(as.mcmc(metro_fl2$theta), 1000)
))

gelman.diag(list(head(metro_fl1.theta.mcmc, 1000), head(metro_fl2.theta.mcmc, 1000)))$psrf
gelman.plot(list(
  head(as.mcmc(metro_fl1$theta), 2000),
  head(as.mcmc(metro_fl2$theta), 2000)
))

# Geweke
geweke.diag(mcmc((metro_fl1.theta)))
# Z-Score plot
geweke.plot(metro_fl1.theta.mcmc, nbins = 50)

# =====
# Credible interval for mu1 and comparison with Laplace approx

# Credible intervals of the Metropolis for the mean
QB_ME <- quantile(metro_fl1.theta[, 1], probs = c(alpha / 2, 1 - alpha / 2))
HPD_ME <- HPDinterval(as.mcmc(metro_fl1.theta[, 1]), prob = 1 - alpha)

# Comparison of Laplace vs Metropolis
credible_intervals_comparison <- function(alpha) {
  colors <- rainbow(n = 4, alpha = .6)
  
  densplot(laplace_fl.mcmc[, 1], main = "Intervals comparison for mu in Flanders", xlab="Income", sub = expression(paste(alpha, " = ", 0.05)))
  legend("topright",
         legend = c("Quantiles Laplace", "Quantiles Metropolis", "HPD Laplace", "HPD Metropolis"), col = colors, lty = c(2, 2, 1, 1), cex = .7
  )
  
  abline(v = c(QB_LP[1], QB_LP[2]), col = colors[1], lty = 2)
  abline(v = c(QB_ME[1], QB_ME[2]), col = colors[2], lty = 2)
  abline(v = c(HPD_LP[1], HPD_LP[2]), col = colors[3], lty = 1)
  abline(v = c(HPD_ME[1], HPD_ME[2]), col = colors[4], lty = 1)
  
  # return(list(QB=QB,HPD=HPD))
}

credible_intervals_comparison(.05)

# ===============================
# QUESTION 6 : JAGS with Flanders
# ===============================

# JAGS Metropolis Model 
metro_model <- function() {
  # Utils
  kappa <- 1 / phi
  lambda <- 1 / (phi * mu)
  
  # Probabilities
  for (i in 1:10) {
    pi[i] <- pgamma(intervals[i + 1], kappa, lambda) - pgamma(intervals[i], kappa, lambda)
  }
  
  # Likelihood
  y ~ dmulti(pi, n)
  
  # Priors
  mu ~ dnorm(3000, pow(306.12, -2)) # 2nd term is precision, which is 1/sd^2
  phi ~ dunif(0, 10)
}

model.file <- "resources/metropolis.bug"
write.model(metro_model, model.file)

# =====
# Test with Flanders
jags_fl <- jags.model(
  file = model.file,
  inits = list(list(mu = estim_mu_fl, phi = estim_phi_fl)),
  data = list(n = flanders.n, intervals = intervals, y = flanders),
  n.chains = 1,
  quiet = T
)

update(jags_fl, 1000)

out_fl <- coda.samples(model = jags_fl, c("mu", "phi"), n.iter = 50000)
out_fl.matrix <- as.matrix(out_fl)

# =====
# Convergence analysis : same methods as previously
traceplot(out_fl)
acf(out_fl.matrix[, 1]); pacf(out_fl.matrix[, 1])
acf(out_fl.matrix[, 2]); pacf(out_fl.matrix[, 2])
effectiveSize(out_fl)

# Gelman Rubin
# Same as previous : define another model
jags_fl2 <- jags.model(
  file = model.file,
  inits = list(list(mu = 4000, phi = 0.1)),
  data = list(n = flanders.n, intervals = intervals, y = flanders),
  n.chains = 1,
  quiet = T
)

update(jags_fl2, 1000)
out_fl2 <- coda.samples(model = jags_fl2, c("mu", "phi"), n.iter = 45000) 
out_fl2.matrix <- as.matrix(out_fl2)

# Traceplot of the two metropolis, comparison
traceplot(list(
  head(as.mcmc(out_fl.matrix), 1000),
  head(as.mcmc(out_fl2.matrix), 1000)
))

gelman.diag(list(head(as.mcmc(out_fl.matrix), 1000), head(as.mcmc(out_fl2.matrix), 1000)))$psrf
gelman.plot(list(
  head(as.mcmc(out_fl.matrix), 2000),
  head(as.mcmc(out_fl2.matrix), 2000)
))

geweke.diag(out_fl.matrix)
geweke.plot(out_fl, nbins = 50)

# =====
# Credible intervals of JAGS model
QB_JA_FL <- quantile(out_fl.matrix[, 1], probs = c(alpha / 2, 1 - alpha / 2))
HPD_JA_FL <- HPDinterval(mcmc(out_fl.matrix[, 1]), prob = 1 - alpha)

# ===============================
# QUESTION 7 : JAGS with Wallonia
# ===============================

jags_wal <- jags.model(
  file = model.file,
  inits = list(list(mu = estim_mu_wal, phi = estim_phi_wal)),
  data = list(n = wallonia.n, intervals = intervals, y = wallonia),
  n.chains = 1,
  quiet = T
)

update(jags_wal, 1000)

out_wal <- coda.samples(model = jags_wal, c("mu", "phi"), n.iter = 50000)
out_wal.matrix <- as.matrix(out_wal)

# =====
# Convergence analysis : same methods as previously
traceplot(out_wal)
acf(out_wal.matrix[, 1]); pacf(out_wal.matrix[, 1])
acf(out_wal.matrix[, 2]); pacf(out_wal.matrix[, 2])
effectiveSize(out_wal)

# Gelman Rubin
# Same as previous : define another model
jags_wal2 <- jags.model(
  file = model.file,
  inits = list(list(mu = 4000, phi = 0.1)),
  data = list(n = 425, intervals = intervals, y = wallonia),
  n.chains = 1,
  quiet = T
)

update(jags_wal2, 1000)

out_wal2 <- coda.samples(model = jags_wal2, c("mu", "phi"), n.iter = 45000)
out_wal2.matrix <- as.matrix(out_wal2)

# Traceplot of the two metropolis, comparison
traceplot(list(
  head(as.mcmc(out_wal.matrix), 1000),
  head(as.mcmc(out_wal2.matrix), 1000)
))

gelman.diag(list(head(as.mcmc(out_wal.matrix), 1000), head(as.mcmc(out_wal2.matrix), 1000)))$psrf
gelman.plot(list(
  head(as.mcmc(out_wal.matrix), 2000),
  head(as.mcmc(out_wal2.matrix), 2000)
))

geweke.diag(out_wal.matrix)
geweke.plot(out_wal, nbins = 50)

# =====
# Credible intervals of JAGS model with wallonia
QB_JA_WAL <- quantile(out_wal.matrix[, 1], probs = c(alpha / 2, 1 - alpha / 2))
HPD_JA_WAL <- HPDinterval(out_wal[,1], prob = 1 - alpha)

# =============================
# QUESTION 8 : Difference of mu
# =============================

# Credible intervals of the difference
diff_of_mu <- out_fl.matrix[, 1] - out_wal.matrix[, 1]
QB_DIFF <- quantile(diff_of_mu, probs = c(alpha / 2, 1 - alpha / 2))
HPD_DIFF <- HPDinterval(as.mcmc(diff_of_mu), prob = 1 - alpha)




