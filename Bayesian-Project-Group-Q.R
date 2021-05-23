
# Require the necessary packages

if(!require(mnormt)){ install.packages("mnormt"); require(mnormt)}
if(!require(EnvStats)){ install.packages("EnvStats"); require(EnvStats)}
if(!require(R2WinBUGS)){ install.packages("R2WinBUGS"); require(R2WinBUGS)}
if(!require(coda)){ install.packages("coda"); require(coda)}
if(!require(rjags)){ install.packages("rjags"); require(rjags)}

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

# Utility fonctions
kappa = function(phi){1/phi}
lambda = function(phi, mu) {1/(phi*mu)}

# Estimate mu and phi by "bootstrapping" the rmultinom

obtain_the_estimated_parameters = function(simulations, row){
  sim_prop = rowSums(rmultinom(simulations, sum(data), prop.table(data)[row,]))
  sim_data = unlist(sapply(1:ncol(data), function(i){rep(intervals.mean[i],sim_prop[i])}))
  est_gam = egamma(sim_data)
  
  est_kappa = as.double(est_gam$parameters["shape"])
  est_lambda = 1/as.double(est_gam$parameters["scale"])
  
  est_mu = est_kappa/est_lambda
  est_phi = 1 / est_kappa
  
  return(list(mu=est_mu, phi=est_phi, sim=sim_data))
}

estim_params_fl <- obtain_the_estimated_parameters(10000, 1)
estim_params_wal <- obtain_the_estimated_parameters(10000, 2)

estim_mu_fl <- estim_params_fl$mu ; estim_phi_fl <- estim_params_fl$phi; estim_fl <- estim_params_fl$sim
estim_mu_wal <- estim_params_wal$mu ; estim_phi_wal <- estim_params_wal$phi; estim_wal <- estim_params_wal$sim