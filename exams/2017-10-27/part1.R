### Part 1

thetaGrid = seq(0.01, 15, length=1000)
bayesPointsD = c()
bayesPointsX = c()
for(theta in thetaGrid) {
  currentRes = 0
  for(y in yProp){
    currentRes = currentRes + log(dbeta(x = y, theta, theta))
  }
  currentRes = currentRes + log(dexp(theta, rate = 1))
  bayesPointsD = c(bayesPointsD, currentRes)
  bayesPointsX = c(bayesPointsX, theta)
}

totalSum = sum(exp(bayesPointsD))
distribution = exp(bayesPointsD)/(totalSum*(15/1000))

plot(x = bayesPointsX, y = distribution)

#b)

thetaLikelihood <- function(thetaVec, xVals) {
  currentRes = 0
  print(thetaVec)
  for(y in xVals) {
    currentRes = currentRes + log(dbeta(x = y, thetaVec[1], thetaVec[2]))
  }
  currentRes = currentRes + log(dexp(thetaVec[1], rate = 1)) + log(dexp(thetaVec[2], rate = 1))
  return(currentRes)
}
initVal <- as.vector(rep(1,2)); 
logPost = thetaLikelihood;
OptimResults<-optim(initVal,logPost,gr=NULL, yProp,method=c("L-BFGS-B"),control=list(fnscale=-1),hessian=TRUE)

std <- sqrt(diag(-solve(OptimResults$hessian)))
betatilde = OptimResults$par

### Part 2
## a)
library("mvtnorm")
omega0 = solve(100*diag(dim(X)[2]))
randomDraws <- BayesLinReg(X = as.matrix(X), y = as.matrix(y), mu_0 = matrix(rep(0, dim(X)[2])), Omega_0 = omega0, v_0 = 1, sigma2_0 = 36, nIter=5000)
meansBeta = c()
meansSigma = c()
stdsBeta = c()
stdsSigma = c()
for(i in 1:dim(X)[2]) {
  currentMeanBeta = mean(randomDraws$betaSample[,i])
  currentMeanSigma = mean(randomDraws$sigma2Sample[,i])
  currentStdBeta = sqrt(var(randomDraws$betaSample[,i]))
  currentStdSigma = sqrt(var(randomDraws$sigmaSample[,i]))
  
  meansBeta = c(meansBeta, currentMeanBeta)
  meansSigma = c(meansSigma, currentMeanSigma)
  meansBeta = c(meansBeta, currentMeanBeta)
  meansSigma = c(meansSigma, currentMeanSigma)
  
  
}

mean14 = mean(randomDraws$betaSample[,14])
std14 = sqrt(var(randomDraws$betaSample[,14]))
stdMean = sqrt(mean(randomDraws$sigma2Sample))
quantile(randomDraws$betaSample[,14], c(0.025, 0.975))

### part b)
x9 = X[9,]
num = 10000
rDraws = rnorm(num, 0, 1)
#rDrawsBeta = runif(num, 0, 1)

valDiff = rep(mean14, num)*x9[14]*-0.3 + rDraws
hist(valDiff, n=30)

### PART 4

y = c(195, 191, 196, 197, 189)
sigma = 10
stdSigma = sqrt(sigma^2 + sigma^2/2)
rDraws = rnorm(1000, mean = mean(y), sd = stdSigma)
hist(rDraws, n=30)

### b)
rDraws = rnorm(100000, mean = mean(y), sd = stdSigma)
rDraws230 = rDraws[rDraws > 230]

prob = length(rDraws230)/length(rDraws)
1 - (1 - prob)^365

### c)

getExpectedLoss <- function(a, meanY, sigmaY, numberOfDays) {
  prob = pnorm(10*a, meanY, sigmaY)^numberOfDays
  notCrashCost = a*prob
  crashCost = (a+100)*(1 - prob)
  return(notCrashCost + crashCost)
}

aGrid = seq(20, 40, length=1000)

lossA = c()
xVal = c()
minVal = 1000
minX = 0
for(a in aGrid){
  currentLoss = getExpectedLoss(a, mean(y), stdSi, 365)
  lossA = c(lossA, currentLoss)
  xVal = c(xVal, a)
  if(currentLoss < minVal) {
    minVal = currentLoss
    minX = a
  }
}
plot(x = xVal, y = lossA)