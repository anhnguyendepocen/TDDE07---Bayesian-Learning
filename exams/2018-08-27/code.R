sigma2 = 0.04
mu0 = 5
sigma0 = 1


xGrid = seq(5.1, 5.4, 0.0001)
dists = c()


logPost <- function(data, mu, sigma2) {
  likelihood = 0
  for(point in data) {
    likelihood = likelihood + dlnorm(x = point, meanlog = mu, sdlog = sqrt(sigma2), log = TRUE)
  }
  prior = dnorm(mu, 5, 1, log = TRUE)
  return(likelihood + prior)
}
for(x in xGrid) {
  dists = c(dists, exp(logPost(lions, x, sigma2)))
}
distsNorm = dists / (sum(dists)*0.0001)
plot(x = xGrid, y = distsNorm)

##

sigma20 = 0.04
v0 = 5

### PART 2

logReg <- function(x, betaVect) {
  upper = exp(t(x)%*%betaVect)
  lower = 1 + exp(t(x)%*%betaVect)
  return(upper/lower)
}

LogPostLogistic <- function(betaVect,y,X,mu,Sigma){ 
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  
  logLik <- sum( linPred*y -log(1 + exp(linPred)));
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, matrix(mu,nPara,1), Sigma, log=TRUE);
  return(logLik + logPrior)
  
}

library('mvtnorm')
tau = 50

X = titanic[,2:6]
y = titanic[,1]

Sigma = tau*diag(dim(X)[2])
mu = 0*rep(0, (dim(X)[2]))
initVal <- as.vector(rep(0,dim(X)[2])); 
logPost = LogPostLogistic;
OptimResults <- optim(par = initVal,logPost,gr=NULL,as.matrix(y),as.matrix(X),mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
OptimResults


names(OptimResults$par) <- covNames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian))) # Computing approximate standard deviations.
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian)))
names(approxPostStd) <- covNames # Naming the coefficient by covariates

betatilde = OptimResults$par
print("Betatilde: ")
print(betatilde)
print("Jacobiany beta: ")
print(approxPostStd)
xVals = c()
for(i in 1:5) {
  xGrid = seq(betatilde[i] - 5*approxPostStd[i], betatilde[i] + 5*approxPostStd[i], length = 10000)
  dist = c()
  xVals = c(xVals, xGrid)
  for(x in xGrid) {
    dist = c(dist, dnorm(x, betatilde[i], approxPostStd[i]))
  }
  plot(x = xGrid, y = dist, main = paste0("beta: ", i, sep = ""))
}

## b)

probs = c()
for(i in 1:5) {
  probs = c(probs, pnorm(0, betatilde[i], approxPostStd[i]))
}
print(pnorm(0, betatilde[2], approxPostStd[2]))

## c)
library("Rlab")
library('mvtnorm')
ladyInput = c(1, 1, 0, 1, 0)
manInput = c(1, 1, 1, 0, 0)
covarMatrix = -solve(OptimResults$hessian)
survivedLady = c()
survivedMan = c()
survivedLadyB = c()
survivedManB = c()
ourScenario = c()
ourScenarioB = c()
betas = rmvnorm(n = 10000, sigma = as.matrix(covarMatrix), mean = matrix(betatilde))
for(row in 1:nrow(betas)) {
  survivedLadyT = exp(ladyInput %*% betas[row, ])/(1 + exp(ladyInput %*% betas[row, ]))
  survivedManT = exp(manInput %*% betas[row, ])/(1 + exp(manInput %*% betas[row, ]))
  scenario = survivedLadyT*(1-survivedManT)
  survivedLady = c(survivedLady, survivedLadyT)
  survivedMan = c(survivedMan, survivedManT)
  ourScenario = c(ourScenario, scenario)
  
  survivedLadyB = c(survivedLadyB, rbern(1, survivedLadyT))
  survivedManB = c(survivedManB, rbern(1, survivedManT))
  ourScenarioB = c(ourScenarioB, rbern(1, scenario))
}
hist(survivedMan, n=30)
hist(survivedLady, n=30)
hist(ourScenario, n=30)

hist(survivedManB, n=30)
hist(survivedLadyB, n=30)
hist(ourScenarioB, n=30)


### PART 3

## c)
x = c(2, 1, 12)
alpha = 0.5
beta = 0.5

theta = 0.5

priorProb1 = 0.1
priorProb2 = 0.9


model1 <- function(theta, alpha, beta, data) {
  n = length(data)
  meanX = mean(data)
  alpha1 = n + alpha
  beta1 = meanX*n + beta
  return(gamma(alpha1 + beta1)/(gamma(alpha1)*gamma(beta1))*theta^(alpha1 - 1)*(1 - theta)^(beta1 - 1))
}

model2 <- function(theta, data) {
  likelihood = sum(log((theta - 1)^(data)*theta))
  return(exp(likelihood))
}

thetaGrid = seq(0, 1, length = 10000)

distVals = c()
maxProb = 0
bestTheta = 0
for(theta2 in thetaGrid) {
  val = model1(theta2, alpha, beta, x)
  distVals = c(distVals, val)
  if(val > maxProb) {
    maxProb = val
    bestTheta = theta2
  }
}
plot(thetaGrid, distVals)

logmarginallikelihood1 = exp(3*5*log(1 - theta) + log(theta)*3)
n = length(x)
meanX = mean(x)
alpha1 = n + alpha
beta1 = n*meanX + beta
logmarginallikelihood2 = lbeta(alpha, beta) - lbeta(alpha1, beta1) 
(model1(bestTheta, alpha, beta, x)*priorProb1)/(likeli*priorProb2)
exp(logmarginallikelihood1)*priorProb1/(exp(logmarginallikelihood2)*priorProb2)

### PART 4
## a)

## b)
library('LaplacesDemon')
prob <- function(x1, x2, x3, prob) {
  return((gamma(x1 + x2 + x3)/gamma(x1)*gamma(x2 + x3))*(prob)^(x1 - 1)*(1 - prob)^(x2+x3 - 1))
}

prob(185, 68, 150, 0.5)
alpha = c(185, 68, 150)
draws = rdirichlet(10000, alpha)

As = draws[,1]
n = length(As)
As50 = As[As > 0.5]
length(As50)/n


###


Bs = draws[,2]
Cs = draws[,3]

AsMaj = As[As >  pmax(Bs, Cs)]

length(AsMaj)/n

### PART 1


