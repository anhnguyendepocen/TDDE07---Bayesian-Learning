library(mvtnorm)
library(readr)
library(matlib)
library(LaplacesDemon)
library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



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

#What omega does is how much the previous value affect the current value and in which direction
### Part b

simulatedX <- c(ARfunction(1000, 10, 0.3, 1))
simulatedY <- c(ARfunction(1000, 10, 0.95, 1))

Nx = length(simulatedX)
Ny = length(simulatedY)

x_dat <- list(N = Nx,
              X = simulatedX)

y_dat <- list(N = Ny,
              X = simulatedY)

fitX <- stan(file = 'mystan.stan', data = x_dat)
print(fitX)
summaryFitX = summary(fitX)$summary
dataX = extract(fitX)

meanX = summaryFitX[,1]
highX = summaryFitX[,8]
lowX = summaryFitX[,4]
efficientX = summaryFitX[, 9]
print(meanX)
print(highX)
print(lowX)
print(efficientX)

plot(dataX$muRandom, dataX$omegaRandom)

library(MASS)
den3d <- kde2d(dataX$muRandom, dataX$omegaRandom)
library(plotly)
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z/length(dataX$muRandom)) %>% add_surface()

plot(cumsum(dataX$muRandom)/seq(1,length(dataX$muRandom)), type='l')
plot(cumsum(dataX$omegaRandom)/seq(1,length(dataX$omegaRandom)), type='l')

fitY <- stan(file = 'mystan.stan', data = y_dat)
print(fitY)
summaryFitY = summary(fitY)$summary
dataY = extract(fitY)

meanY = summaryFitY[,1]
highY = summaryFitY[,8]
lowY = summaryFitY[,4]
efficientY = summaryFitY[, 9]
print(meanY)
print(highY)
print(lowY)
print(efficientY)

plot(dataY$muRandom, dataY$omegaRandom)

library(MASS)
den3d <- kde2d(dataY$muRandom, dataY$omegaRandom)
library(plotly)
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z/length(dataY$muRandom)) %>% add_surface()

plot(cumsum(dataY$muRandom)/seq(1,length(dataY$muRandom)), type='l')
plot(cumsum(dataY$omegaRandom)/seq(1,length(dataY$omegaRandom)), type='l')
###

#pairs(fit, pars = c("mu", "tau", "lp__"))

#la <- extract(fit, permuted = TR UE) # return a list of arrays 
#mu <- la$mu 

### return an array of three dimensions: iterations, chains, parameters 
#a <- extract(fit, permuted = FALSE) 

### use S3 functions on stanfit objects
#a2 <- as.array(fit)
#m <- as.matrix(fit)
#d <- as.data.frame(fit)

### part c
campyData<-c(t(read.table("campy.dat",header=TRUE)))

Ny = length(campyData)

campy_dat <- list(N = Ny,
              y = campyData)
fitP <- stan(file = 'poisson.stan', data = campy_dat)

summaryFitP = summary(fitP)$summary
rowsCols = dim(summaryFitP)
highTheta = c()
meanTheta = c()
lowTheta = c()
for(i in 1:(rowsCols[1]-4)) {
  highTheta = c(highTheta, exp(summaryFitP[i, 8]))
  meanTheta = c(meanTheta, exp(summaryFitP[i, 1]))
  lowTheta = c(lowTheta, exp(summaryFitP[i, 4]))
}

plot(meanTheta, type='l', ylim=c(min(lowTheta), 50))
lines(highTheta, col='red')
lines(lowTheta, col='red')

### part d
campy_dat <- list(N = Ny,
                  y = campyData,
                  sigma2Nu = 50,
                  sigma2Sigma = 0.01)
fitP <- stan(file = 'poissonPrior.stan', data = campy_dat)

summaryFitP = summary(fitP)$summary
rowsCols = dim(summaryFitP)
rowsCols[1] = rowsCols[1] - 4
highTheta = c()
meanTheta = c()
lowTheta = c()
for(i in 1:rowsCols[1]) {
  highTheta = c(highTheta, exp(summaryFitP[i, 8]))
  meanTheta = c(meanTheta, exp(summaryFitP[i, 1]))
  lowTheta = c(lowTheta, exp(summaryFitP[i, 4]))
}

plot(meanTheta, type='l', ylim=c(min(lowTheta), 50))
lines(highTheta, col='red')
lines(lowTheta, col='red')