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

for(i in 1:21) {
  simulatedX <- ARfunction(200, 10, -11/10 + i/10, 2)
  plot(simulatedX, type='l')
  Sys.sleep(0.5)
}
###
simulatedX <- ARfunction(100000, 10, 1, 1)
plot(simulatedX, type='l')


### Part b

yI = c(4,5,6,4,0,2,5,3,8,6,10,8)
Ny = length(yI)

normal_dat <- list(N = Ny, 
                    y = yI)
fit <- stan(file = 'lec1stan.stan', data = normal_dat)

print(fit)

###

simulatedX <- c(ARfunction(1000, 10, 0.3, 1))
simulatedY <- c(ARfunction(1000, 10, 0.95, 1))

Nx = length(simulatedX)
Ny = length(simulatedY)

x_dat <- list(N = Nx,
              X = simulatedX,
              initialMuMu = 0,
              muSigma2 = 1,
              initialOmegaMu = 0,
              omegaSigma2 = 1,
              sigma2Nu = 0,
              sigma2Sigma = 1)

y_dat <- list(N = Ny,
              X = simulatedY,
              initialMuMu = 0,
              muSigma2 = 100000,
              initialOmegaMu = 0,
              omegaSigma2 = 100000,
              sigma2Nu = 0,
              sigma2Sigma = 100000)

fitX <- stan(file = 'mystan.stan', data = x_dat)
print(fitX)

fitY <- stan(file = 'mystan.stan', data = y_dat)
print(fitY)

plot(fitX)
pairs(fitX, pars = c("muRandom", "sigma2Random", "omegaRandom", "lp__"))
###

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

### part c
campyData<-c(t(read.table("campy.dat",header=TRUE)))

Ny = length(campyData)

campy_dat <- list(N = Ny,
              y = campyData,
              initialMuMu = 0,
              muSigma2 = 100000,
              initialOmegaMu = 0,
              omegaSigma2 = 100000,
              sigma2Nu = 0,
              sigma2Sigma = 100000)
fitP <- stan(file = 'poisson.stan', data = campy_dat)

print(fitP)

### part d
campy_dat <- list(N = Ny,
                  y = campyData,
                  initialMuMu = 0,
                  muSigma2 = 100000,
                  initialOmegaMu = 0,
                  omegaSigma2 = 100000,
                  sigma2Nu = 0,
                  sigma2Sigma = 1)
fitP <- stan(file = 'poisson.stan', data = campy_dat)

print(fitP)

