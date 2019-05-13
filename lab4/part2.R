library(mvtnorm)
library(readr)
library(matlib)
library(LaplacesDemon)
library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

bidderData<-read.table("ebayNumberOfBidderData.dat",header=TRUE)


ARfunction <- function(n, mu, omega, sigma2) {
  x = matrix(0, n, 1)
  x[1, 1] = mu
  noise <- rnorm(n-1,0,sigma2)
  for(i in 1:(n-1)) {
    x[i+1, 1] = mu + omega*(x[i, 1] - mu) + noise[i]
  }
  return(x)
}

simulatedX <- ARfunction(100, 10, -1, 2)
plot(simulatedX, type='l')

for(i in 1:21) {
  simulatedX <- ARfunction(100, 10, -11/10 + i/10, 2)
  plot(simulatedX, type='l')
}


### Part b

simulatedX <- ARfunction(1000, 10, 0.3, 1)
simulatedY <- ARfunction(1000, 10, 0.95, 1)

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit <- stan(file = '8schools.stan', data = schools_dat)

print(fit)
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = TR UE) # return a list of arrays 
mu <- la$mu 

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit, permuted = FALSE) 

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)
