install.packages("EnvStats")
library(EnvStats)

mu = 50
phi=100
gammad = dgamma(x,1/phi,1/(phi*mu))
gammap = pgamma(x,1/phi,1/(phi*mu))
plot(gammad)
plot(gammap)


mu = 2500
phi = 1.0
kappa = 1/phi
lambda = 1/(phi*mu) 
pgamma(2700,kappa,lambda)-pgamma(2300,kappa,lambda)

plot(pgamma(2700,kappa,lambda))

kappaFct <- function(phi){1/phi}
lambdaFct <- function(phi,mu){1/(phi*mu)}

muExample <- 2500
phiExample <- 1
kappaExample <- kappaFct(phiExample)
lambdaExample <- lambdaFct(muExample,phiExample)

pgamma(2700,kappaExample,lambdaExample)
pgamma(2700,kappaExample,lambdaExample)-pgamma(2300,kappaExample,lambdaExample)


x <- rgamma(10000,shape = 1/phi,scale =1/(phi*mu))

xtest = rep(1200, 25)

muflander = 3089


# Flander 
x1 <- rep(1200, 25)
x1 <- append(x1, runif(69, 1200, 1500))
x1 <- append(x1, runif(65, 1500, 1800))
x1 <- append(x1, runif(106, 1800, 2300))
x1 <- append(x1, runif(80, 2300, 2700))
x1 <- append(x1, runif(106, 2700, 3300))
x1 <- append(x1, runif(136, 3300, 4000))
x1 <- append(x1, runif(94, 4000, 4900))
x1 <- append(x1, runif(76, 4900, 6000))
x1 <- append(x1, rep(6000, 46))
length(x1)

shape_test = egamma(x1)$parameters[1]
scale_test = egamma(x1)$parameters[2]

shape_test*scale_test


plot(density(x1), ylim=c(0,0.0003))
lines(density(rgamma(1000000,shape = shape_test, 1/scale_test)), col="red")


# Walonia 

x2 <- rep(1200, 17)
x2 <- append(x2, runif(36, 1200, 1500))
x2 <- append(x2, runif(47, 1500, 1800))
x2 <- append(x2, runif(58, 1800, 2300))
x2 <- append(x2, runif(47, 2300, 2700))
x2 <- append(x2, runif(53, 2700, 3300))
x2 <- append(x2, runif(59, 3300, 4000))
x2 <- append(x2, runif(54, 4000, 4900))
x2 <- append(x2, runif(33, 4900, 6000))
x2 <- append(x2, rep(6000, 21))
length(x2)

shape_test2 = egamma(x2)$parameters[1]
scale_test2 = egamma(x2)$parameters[2]

shape_test2*scale_test2


plot(density(x2), ylim=c(0,0.0003))
lines(density(rgamma(1000000,shape = shape_test2, 1/scale_test2)), col="red")





