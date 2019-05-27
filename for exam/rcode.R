#### LAB 1

### PART 1

library(manipulate)

alpha = 2
beta = 2
s = 14
num = 20
f = num - s
p = s/num

prior = dbeta(alpha, beta, )
posterior = rbeta(alpha + s, beta + f, n = num)

print(paste0("ground truth std: ", sqrt(p*(1-p))))
print(paste0("ground truth Mean: ", p))

BetaRandomPosteriorPriorPlotandMeanStdComparison <- function(a,b,n,p){
  xGrid <- seq(0.001, 0.999, by=0.001)
  normalizedLikelihood = dbeta(xGrid, 20*p+1, 20*(1-p)+1)
  prior = dbeta(xGrid, a, b)
  posterior = dbeta(xGrid, a+20*p, b+20*(1-p))
  posteriorAlpha = a+20*p
  posteriorBeta = b+20*(1-p)
  posteriorRandom <- rbeta(n, posteriorAlpha, posteriorBeta)
  maxDensity <- max(normalizedLikelihood, prior, posterior) # Use to make the y-axis high enough
  
  plot(xGrid, normalizedLikelihood, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "theta", 
       ylab = 'Density', main = 'Bernoulli model - Beta(a,b) prior')
  lines(xGrid, posterior, lwd = 3, col = "red")
  lines(xGrid, prior, lwd = 3, col = "green")
  legend(x = 0.01, y = maxDensity*0.95, legend = c("Likelihood (normalized)", "Prior", "Posterior"), col = c("blue","green","red"), lwd = c(3,3,3), cex = 0.7)
  
  posteriorMean = posteriorAlpha/(posteriorAlpha + posteriorBeta)
  print(paste0('Posterior mean GT: ', posteriorMean))
  print(paste0("ground truth std: ", sqrt(posteriorAlpha*posteriorBeta/(((posteriorAlpha+posteriorBeta)^2)*(posteriorAlpha+posteriorBeta+1)))))
  print(paste0("std: ", sd(posteriorRandom)))
  print(paste0("Mean: ", mean(posteriorRandom)))
  
  #b part
  posteriorAlpha = a+20*p
  posteriorBeta = b+20*(1-p)
  posteriorRandomDraw <- rbeta(10000, posteriorAlpha, posteriorBeta)
  posteriorRandomOverCondition = which(posteriorRandomDraw < 0.4)
  probabilityCondition = length(posteriorRandomOverCondition)/length(posteriorRandomDraw)
  trueProbability = pbeta(0.4, posteriorAlpha, posteriorBeta)
  
  print(paste0("propability condition with random: ", probabilityCondition))
  print(paste0("ground truth probability: ", trueProbability))
  
  #c part
  logOdds = log(posteriorRandomDraw/(1 - posteriorRandomDraw))
  density(logOdds)
  # Kommentera bort hist om du vill köra det bort kommenterade ovan, om båda delarna plottas så kommer manipulate-panelen inte upp.
  #hist(logOdds)
  #hist(posteriorRandomDraw)
  #hist(posteriorRandomDraw/(1 - posteriorRandomDraw))
}


manipulate(
  BetaRandomPosteriorPriorPlotandMeanStdComparison(alpha,beta,num,p),
  num = slider(1, 10000, step=1, initial = 2, label = "Number of trials")
)

### PART 2

library(manipulate)
library(LaplacesDemon)
u = 3.5
y = c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
logY = log(y)
meanY = sum(y)/length(y)
meanLogY = sum(logY)/length(logY)

getDistrubutionFromPoints <- function(x){
  n = length(x)
  sortedX = sort(x)
  cumPercent = (1:n)/n
  densityPercent = (cumPercent[2:n] - cumPercent[1:n-1])/(sortedX[2:n] - sortedX[1:(n-1)])
  
  dataMatrix = matrix(1:(2*n), nrow=n, ncol=2, byrow=TRUE)
  dataMatrix[ ,1] = sortedX
  dataMatrix[ ,2] = cumPercent
  
  return(dataMatrix)
}

getIntervalData <- function(x, lowerPercent, upperPercent){
  lowerFlag = FALSE
  upperFlag = FALSE
  n = length(x)
  sortedX = sort(x)
  i = 0
  for(data in x) {
    i = i + 1
    if(i/n >= lowerPercent && !lowerFlag) {
      lowerFlag = TRUE
      lowerIndex = i
    }
    else if(i/n > upperPercent){
      upperIndex = i - 1
      upperFlag = TRUE
      break;
    }
    
  }
  if(!upperFlag){
    upperIndex = n
  }
  return(sortedX[lowerIndex:upperIndex])
}


dScaledInvChi2 <- function(x, v, tau2) {
  densities = array(0, length(x))
  for(i in 1:length(x)) {
    densities[i] = (((tau2*v/2)^(v/2)/gamma(v/2))*(exp(-v*tau2/(2*x[i]))/(x[i]^(1+v/2))))
  }
  return(densities)
}

drawFromPosterior <- function(mean, n, logy){
  xGrid <- seq(0.001, 0.999, by=0.001)
  tau2 = (sum(logy - mean)^2)/n
  posterior = dScaledInvChi2(xGrid, n, tau2)
  #posterior = dinvchisq(xGrid, n, tau2, log=FALSE)
  print(posterior)
  maxDensity <- max(posterior)
  
  posteriorRandomDraw <- rinvchisq(10000,  n, tau2)
  
  #plotting
  attach(mtcars)
  par(mfrow=c(4,2))
  plot(xGrid, posterior, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "theta", 
       ylab = 'Density', main = 'Chi Posterior of Variance')
  plot(sort(posteriorRandomDraw))
  hist(posteriorRandomDraw)
  
  #part b)
  giniTheory = 2*pnorm(posteriorRandomDraw/sqrt(2), mean = 0, sd = 1) - 1
  
  logGini = log(giniTheory)
  distrubtionGiniLog = getDistrubutionFromPoints(logGini)
  distrubtionGini = getDistrubutionFromPoints(giniTheory)
  plot(distrubtionGiniLog)
  hist(logGini)
  plot(distrubtionGini)
  hist(giniTheory)
  
  #part c)
  lowerData = getIntervalData(giniTheory, 0, 0.025)
  middleData = getIntervalData(giniTheory, 0.025, 0.975)
  upperData = getIntervalData(giniTheory, 0.975, 1)
  
  print(density(middleData))
  print(paste0("Lower end of interval: ", min(middleData)))
  print(paste0("Upper end of interval: ", max(middleData)))
}

num = length(logY)
draws = drawFromPosterior(u, num, logY)

### PART 3

mu = 2.39
obsRad = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

kGrid <- seq(0.001, 10, by=0.001)

priorDist <- dexp(kGrid)
kappa = kGrid

dvonmises <- function(y, mu, kappa){
  res = 1
  for (i in y){
    res = res*(exp(kappa*cos(i-mu))/(2*pi*besselI(kappa,0)))
  }
  return(res)
}

likeliHood <- dvonmises(obsRad, mu, kappa)
posterior = likeliHood*priorDist

posterior = (posterior/sum(posterior))/0.001
likeliHood = (likeliHood/sum(likeliHood))/0.001

plot(kGrid, priorDist, type = 'l', lwd = 3, col = "blue", xlab = "k", 
     ylab = 'Density', main = 'von Mises - Wind direction')
lines(kGrid, likeliHood, lwd = 3, col = "red")
lines(kGrid, posterior, lwd = 3, col = "green")
legend(x = 5, y = 0.8, legend = c("Likelihood (normalized)", "Prior", "Posterior"), 
       col = c("red","blue","green"), lwd = c(3,3,3), cex = 0.7)    


#### LAB 2

### PART 1

library(mvtnorm)
library(readr)
library(matlib)
library(LaplacesDemon)

TempLinkoping <- read_csv("TempLinkoping.csv")
#View(TempLinkoping)

x = TempLinkoping['time']
y = TempLinkoping['temp']
x['time2'] = x^2
x['1'] = x['time']/x['time']
ph = x['time']
#Just to have the betas in right order
x['time'] = x['1']
x['1'] = x['time2']
x['time2'] = ph

matrix_x = data.matrix(x)
matrix_y = data.matrix(y)

mu0 = c(-11,85,-70)
omega0 = matrix(c(0.03, 0, 0, 0, 0.01, 0, 0, 0, 0.03), 3, 3)
v0 = 3
sigmasq0 = 0.03

e = rnorm(1, mean = 0, sd = sigmasq0)


betahat = inv((t(matrix_x)%*%matrix_x))%*%(t(matrix_x)%*%matrix_y) 
randomSigma2 <- rinvchisq(n = 100, df = v0, scale = sigmasq0)
randomBetas <- c()
plot(TempLinkoping, col="black") 
for(singleSigma in randomSigma2) {
  randomBeta <- rmvt(n = 1,mu = t(mu0), S = singleSigma*inv(omega0))
  ys = c()
  xs = cbind(matrix(1, 1000, 1), matrix(1:1000, 1000, 1)/1000, (matrix(1:1000, 1000, 1)/1000)^2)
  ys = xs%*%randomBeta[1, ]
  #for (i in 1:1000) {
  #  y = randomBeta[k, 1] + randomBeta[k, 2]*i/1000 + randomBeta[k, 3]*(i/1000)^2
  #  ys = c(ys, y)
  #  xs = c(xs, i/1000)
  #}
  lines(matrix(1:1000, 1000, 1)/1000, ys, col="red") 
}


### part b)

muN = inv(t(matrix_x)%*%matrix_x + omega0)%*%(t(matrix_x)%*%matrix_x%*%betahat + omega0%*%mu0)
omegaN = t(matrix_x)%*%matrix_x + omega0
vN = v0 + 366
vNsigmaN2 = v0*sigmasq0 + (t(matrix_y)%*%matrix_y + t(mu0)%*%omega0%*%mu0 - t(muN)%*%omegaN%*%muN)

randomSigma2 <- rinvchisq(n = 1000, df = vN, scale = (vNsigmaN2/vN))
randomBetas <- c()
plot(TempLinkoping, col="black")
sigmas = data.frame(randomSigma2)
y_df = data.frame(matrix(1, 1000, 1))
betas0 = c()
betas1 = c()
betas2 = c()
sigma2s = c()
ibeta = 0
for(singleSigma in randomSigma2) {
  randomBeta <- rmvt(n = 1,mu = t(muN), S = singleSigma*inv(omegaN))
  ibeta = ibeta + 1
  ys = c()
  for (i in 1:1000) {
    y = randomBeta[1, 1] + randomBeta[1, 2]*i/1000 + randomBeta[1, 3]*(i/1000)^2
    ys = c(ys, y)
  }
  betas0 = c(betas0, randomBeta[1, 1])
  betas1 = c(betas1, randomBeta[1, 2])
  betas2 = c(betas2, randomBeta[1, 3])
  sigma2s = c(sigma2s, singleSigma)
  y_df[paste0("trial", ibeta)] <- data.frame(ys)
}
xs = c()
for(i in 1:1000) {
  xs = c(xs, i/1000)
}
y_df = subset(y_df, select = -c(1) )

mediany = matrix(1, 1000, 1)
lowery = matrix(1, 1000, 1)
uppery = matrix(1, 1000, 1)
for (row in 1:nrow(y_df)) {
  mediany[row] = median(as.numeric(as.vector(y_df[row, ])))
  lowery[row] = quantile(x = as.numeric(as.vector(y_df[row, ])), probs = 0.025)
  uppery[row] = quantile(x = as.numeric(as.vector(y_df[row, ])), probs = 0.975)
}
plot(TempLinkoping, col="black")
lines(xs, mediany, col="blue")
lines(xs, lowery, col="red")
lines(xs, uppery, col="red")
hist(betas0, nclass=30)
hist(betas1, nclass=30)
hist(betas2, nclass=30)
hist(sigma2s, nclass=30)

### part c

timeOptimum = -betas1/(2*betas2)
hist(timeOptimum, nclass=30)

### part d

## If there is a higher order but we are more certain that higher order parameters are not needed we can set the omega-values high as it creates a stronger prior 
## and the mu-prior values at 0.


### PART 2

womenWork<-read.table("WomenWork.dat",header=TRUE)  # Spam data from Hastie et al.
glmModel <- glm(Work ~ 0 + ., data = womenWork, family = binomial)

### b

chooseCov <- c(1:8) # Here we choose which covariates to include in the model
tau <- 10; # Prior scaling factor such that Prior Covariance = (tau^2)*I
###########     END USER INPUT    ################


# install.packages("mvtnorm") # Loading a package that contains the multivariate normal pdf
library("mvtnorm") # This command reads the mvtnorm package into R's memory. NOW we can use dmvnorm function.

# Loading data from file
Data<-read.table("WomenWork.dat",header=TRUE)  # Spam data from Hastie et al.
y <- as.vector(womenWork[,1]); # Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
X <- as.matrix(womenWork[,2:9]);
covNames <- names(womenWork)[2:length(names(womenWork))];
X <- X[,chooseCov]; # Here we pick out the chosen covariates.
covNames <- covNames[chooseCov];
nPara <- dim(X)[2];

# Setting up the prior
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara);

LogPostLogistic <- function(betaVect,y,X,mu,Sigma){ 
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  
  logLik <- sum( linPred*y -log(1 + exp(linPred)));
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  return(logLik + logPrior)
}

initVal <- as.vector(rep(0,dim(X)[2])); 
logPost = LogPostLogistic;
OptimResults<-optim(initVal,logPost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Printing the results to the screen
names(OptimResults$par) <- covNames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian))) # Computing approximate standard deviations.
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian)))
names(approxPostStd) <- covNames # Naming the coefficient by covariates

betatilde = OptimResults$par
print("Betatilde: ")
print(betatilde)
print("Jacobiany beta: ")
print(approxPostStd)

print("Intervall NSmallChild: ")
upperb = betatilde["NSmallChild"] + 1.96*approxPostStd["NSmallChild"]
lowerb = betatilde["NSmallChild"] - 1.96*approxPostStd["NSmallChild"]
print(upperb)
print(lowerb)

### part c
covarMatrix = -inv(OptimResults$hessian)

ladyInput = c(1, 10, 8, 10, (10/10)^2, 40, 1, 1)
workOrNots = c()
workOrNotsBinary = c()
betas = rmvt(n = 1000,mu = matrix(betatilde), S = covarMatrix)
for(row in 1:nrow(betas)) {
  working = exp(ladyInput %*% betas[row, ])/(1 + exp(ladyInput %*% betas[row, ]))
  workingBinary = round(working)
  workOrNots = c(workOrNots, working)
  workOrNotsBinary = c(workOrNotsBinary, workingBinary)
}
hist(workOrNots, n=30)
hist(workOrNotsBinary, n=30)

#### LAB 3

### PART 1

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

### PART 1b

# Estimating a simple mixture of normals
# Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com

##########    BEGIN USER INPUT #################
# Data options
rawData <- read.table("rainfall.dat",header=TRUE)
x <- as.matrix(rawData)

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(30,nComp) # Prior mean of mu
tau2Prior <- rep(10,nComp) # Prior std of mu
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 1000 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 0 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0.2,0.8,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x, n = 50)$density))
stored_mu = matrix(0, nIter, nComp)
stored_sigma = matrix(0, nIter, nComp)

for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
    stored_mu[k, j] = mu[j]
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
    stored_sigma[k, j] = sigma2[j]
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%50 ==0)){
    effIterCount <- effIterCount + 1
    hist(x, breaks = 50, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  
}

hist(x, breaks = 50, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDens, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)

par(mfrow=c(2,2))
plot(stored_mu[, 1], type='l', main = 'mu 1')
plot(stored_mu[, 2], type='l', main = 'mu 2')
plot(stored_sigma[, 1], type='l', main = 'sigma 1')
plot(stored_sigma[, 2], type='l', main = 'sigma 2')

par(mfrow=c(2,2))
plot(stored_mu[20:nIter, 1], type='l', main = 'mu 1')
plot(stored_mu[20:nIter, 2], type='l', main = 'mu 2')
plot(stored_sigma[20:nIter, 1], type='l', main = 'sigma 1')
plot(stored_sigma[20:nIter, 2], type='l', main = 'sigma 2')

### PART 2

library(mvtnorm)
library(readr)
library(matlib)
library(LaplacesDemon)

rawData <- read.table("ebayNumberOfBidderData.dat",header=TRUE)

y <- as.vector(rawData[,1]); # Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
X <- rawData[,2:10];

model <- glm(y ~ PowerSeller+VerifyID+Sealed+Minblem+MajBlem+LargNeg+LogBook+MinBidShare, data=X, family=poisson())

## Significant variables are VerifyID, Sealed, MajBlem, LogBook and MinBidShare

# b)
xMatrix = as.matrix(X)
sigmaGPrior = 100*solve((t(xMatrix)%*%xMatrix))
tau = 10
chooseCov <- c(1:9)

covNames <- names(rawData)[2:length(names(rawData))];
xMatrix <- xMatrix[,chooseCov]; # Here we pick out the chosen covariates.
covNames <- covNames[chooseCov];
nPara <- dim(X)[2];

# Setting up the prior
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara);

PoiPost <- function(theta,y,X, mu, SigmaGPrior) {
  nPara <- length(theta);
  linPred <- X%*%theta;
  logPoiLik <- sum( linPred*y -exp(linPred) - log(factorial(y)));
  if (abs(logPoiLik) == Inf) logPoiLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logBetaPrior <- dmvnorm(theta, matrix(0,nPara,1), SigmaGPrior, log=TRUE);
  return(logPoiLik + logBetaPrior) 
}

initVal <- as.vector(rep(0,dim(X)[2])); 
logPost = PoiPost;
OptimResults<-optim(initVal,logPost,gr=NULL,y,xMatrix,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

approxPostStd <- sqrt(diag(-solve(OptimResults$hessian)))
names(approxPostStd) <- covNames # Naming the coefficient by covariates

betatilde = OptimResults$par
print("Betatilde: ")
print(betatilde)
print("Jacobiany beta: ")
print(approxPostStd)

### part c)

RWMSampler <- function(logPostFunc, n,  c, covar,  ...) {
  currentTheta = as.vector(rep(0, dim(covar)[1]))
  draws = matrix(0, nrow = n, ncol = dim(covar)[1])
  oldProbability = logPostFunc(currentTheta, ...)
  for(i in 1:n) {
    currentDraw <- rmvt(1, mu = currentTheta, S = c*covar)
    newProbability <- logPostFunc(as.vector(currentDraw), ...)
    
    alpha = min(1, newProbability/oldProbability)
    uniformDraw = runif(1, 0, 1)
    if(uniformDraw >= alpha) {
      oldProbability = newProbability
      currentTheta = currentDraw
    }
    draws[i,] = currentDraw
  }
  return(draws)
}

myDraws <- RWMSampler( PoiPost, 10000, 20, diag(diag(-solve(OptimResults$hessian))), y, xMatrix, mu, Sigma)

hist(myDraws[,9])

for(i in 1:9){
  plot(myDraws[,i], type='s')
  a = c(rep(betatilde[i],length(myDraws)))
  lines(a, col='red')
}

#### LAB 4

### R-code

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

### C-stan code

data {
  int<lower=0> N;
  int  y[N];
  real sigma2Nu;
  real sigma2Sigma;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real x[N];
  real muRandom;
  real omegaRandom;
  real sigma2Random;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  sigma2Random ~ scaled_inv_chi_square(sigma2Nu, sigma2Sigma);
  //sigma2Random ~ normal(sigma2Nu, sigma2Sigma);
  
  
  x[1] ~ normal(muRandom, sigma2Random);
  y[1] ~ poisson(exp(x[1]));
  for(i in 2:N) {
    x[i] ~ normal(muRandom + omegaRandom*(x[i-1] - muRandom), sigma2Random);
    y[i] ~ poisson(exp(x[i]));
  }
}