library(mvtnorm)
library(readr)
library(matlib)
library(LaplacesDemon)

rainfallData<-read.table("rainfall.dat",header=TRUE)
plotColors = list('red', 'green', 'blue', 'yellow', 'black', 'orange')
# Setup 
mu1 <- 1
mu2 <- -1
rho <- 0.9
mu <- c(mu1,mu2)
Sigma = matrix(c(1,rho,rho,1),2,2)
nDraws <- 1000 # Number of draws

v0 = 2
mu0 = 20
sigma0 = 1
tau0 = 20

currentMu = mu0
currentSigma = sigma0
currentTau = tau0
currentV = v0

n = nrow(rainfallData)
averageRain = sum(rainfallData)/nrow(rainfallData)
sigma = var(rainfallData-averageRain)[1]
gibbsDraws <- matrix(0,nDraws,2)

for (i in 1:nDraws) {
  currentTau = 1/((n/currentSigma) + 1/(tau0)^2)
  w = (n/currentSigma)/(n/currentSigma + 1/(tau0)^2)
  currentMu = w*averageRain + (1 - w)*mu0
  currentMu <- rnorm(1, currentMu, currentTau)
  gibbsDraws[i,1] <- currentMu

  currentV = n + v0
  currentSigma = (v0*sigma0 + sum((rainfallData - currentMu)^2)/(n + v0))
  currentSigma <- rinvchisq(n = 1, df = currentV, scale = currentSigma)
  gibbsDraws[i, 2] <- currentSigma
}

averages <- matrix(0, nDraws, 2)
for (i in (1:((nDraws-500)/10))) {
  averages[i, 1] = sum(gibbsDraws[(i*10):(i*10+500), 1])/500
  averages[i, 2] = sum(gibbsDraws[(i*10):(i*10+500), 2])/500
}
#plot(averages[1:(nDraws-500)/10, 1], type='l')
#lines(1:nDraws, matrix(averageRain, (nDraws-500)/10, 1))

#plot(averages[1:(nDraws-500)/10, 2], type='l')
#lines(1:nDraws, matrix(sigma, (nDraws-500)/10, 1))

plot(1:nDraws, gibbsDraws[,1], type = "l", col = 'red', ylab='y', 
     lwd = 1, axes=FALSE, xlab = 'MCMC iteration', xlim = c(0,nDraws), ylim = c(minY, maxY), 
     main = 'Raw - Gibbs')
axis(side = 1, at = seq(0, nDraws, by = 250))
axis(side = 2, at = seq(minY, maxY, by = 0.5))

hist(gibbsDraws[,1], freq = FALSE, main='Gibbs draws', ylim = c(0,0.5), xlab='x')
lines(seq(-2,4,by=0.01),dnorm(seq(-2,4,by=0.01), mean = 1), col = 'orange', 
      lwd = 1)

cusumData =  cumsum(gibbsDraws[,1])/seq(1,nDraws)
minY = floor(min(cusumData))
maxY = ceiling(max(cusumData))
plot(1:nDraws, cusumData, type = "l", col = 'orange', ylab='Cumulative estimate', 
     lwd = 1, axes=FALSE, xlab = 'MCMC iteration', xlim = c(0,nDraws), 
     ylim = c(minY,maxY), main = 'Cusum - Gibbs')
lines(seq(1,nDraws),1*matrix(1,1,nDraws),col = 'orange', lwd=1)
axis(side = 1, at = seq(0, nDraws, by = 250))
axis(side = 2, at = seq(minY, maxY, by = 0.5))

a = acf(gibbsDraws[,1], main='Gibbs draws', lag.max = 20, plot = F)
barplot(height = a$acf[-1], names.arg=seq(1,20), col = 'orange')

#dev.off()

gibbsDraws1 = c(mu0, gibbsDraws[,1])
gibbsDraws2 = c(sigma0, gibbsDraws[,2])
plot(gibbsDraws1,gibbsDraws2, col ='red')
#plot(gibbsDraws[,1],gibbsDraws[,2])
plot(gibbsDraws1[1:50],gibbsDraws2[1:50], type ='s', col ='orange')

# Plotting the cumulative path of estimates of Pr(theta1>0, theta2>0)
par(mfrow=c(2,1))
plot(cumsum(gibbsDraws[,1]>32 & gibbsDraws[,2]>1600)/seq(1,nDraws),type="l", main='Gibbs draws', xlab='Iteration number', ylab='', ylim = c(0,1))