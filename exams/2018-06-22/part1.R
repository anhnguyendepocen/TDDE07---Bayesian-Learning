## PART 1
nXs = 10
beta = 2
aGrid = seq(0.01, 20, length=1000)

diffDist = c()



for(a in aGrid) {
  xGrid = seq(0.0)
  priorMean = a*beta
  xGrid = 
  
  
}


rdraws 

## b)

## Part 2
# a)
library('mvtnorm')
mu0 = c(0, 0, 0)
omega0 = 0.01*diag(3)
#omega0[3, 3] = 1000000
v0 = 1
sigma20 = 10000
fishL = fish$length
fishX = fish[,2:4]


bDraws = BayesLinReg(X = as.matrix(fishX), y = as.matrix(fishL), mu_0 = mu0, Omega_0 = omega0, v_0 = v0, sigma2_0 = sigma20, nIter = 10000)

for(i in 1:3){
  hist(bDraws$betaSample[,i],n=30, probability = TRUE)
}
## b)

for(i in 1:3){
  print(quantile(bDraws$betaSample[,i],c(0.05, 0.95)))
}

## c)

## d)

sigmaM = mean(bDraws$sigma2Sample)
expectedLength <- 
rDrawFish <- function(n, betaVect, xVals, sigma) {
  pred = betaVect%*%xVals
  draws = rnorm(n, pred, sigma)
  return(draws)
}
library('Rlab')
randomB <- function(betaVects, probability, sigmas) {
  draws = c()
  for(i in 1:length(betaVects[,1])) {
    bernD = rbern(1, probability)
    if(bernD == 1) {
      draws = c(draws, rDrawFish(1, betaVects[i,], c(1, 30, 30), sqrt(sigmas[i])))
    }
    else {
      draws = c(draws, rDrawFish(1, betaVects[i,], c(1, 100, 30), sqrt(sigmas[i])))
    }
  }
  return(draws)
}
betas = c(mean(bDraws$betaSample[,1]), mean(bDraws$betaSample[,2]), mean(bDraws$betaSample[,3]))
draws = randomB(bDraws$betaSample, 0.5, bDraws$sigma2Sample)
hist(draws, n=50, probability=TRUE)

### PART 4

## a)

sigma = 100
l = 200
logPostFunc <- function(mu, sigma, data, l) {
  prior = log(1)
  upper = sum(log(dnorm((data - mu)/sigma)) - log(sigma) - log(1 - pnorm((l-mu)/sigma)))
  #lower = n*log(sigma*(1 - pnorm(((l-mu)/sigma))))
  loglikelihood = upper #- lower
  return(loglikelihood + prior)
}

xGrid = seq(100, 400, 1)
logPosts = c()
filterData = sulfur[sulfur > l]
xVals = c()
for(x in xGrid) {
  print(x)
  logPosts = c(logPosts, logPostFunc(x, sigma, filterData, l))
  xVals = rbind(xVals, x)
}
plot(x = xGrid, y = logPosts)

## JESPER/JONTE kod

logPostFunc = function(x, mu, sigma, L) {
  distribution <- sum(log(dnorm((x-mu)/sigma)) - log(sigma) - log(1 - pnorm((L-mu)/sigma))) + 0
  return(distribution)
}

xGrid <- seq(100, 400, 1)
x <- sulfur[sulfur > 200]
sigma <- 100
L <- 200

distribution <- rep(0, length(xGrid))
for (mu in xGrid) {
  print(mu)
  distribution[mu-99] <- logPostFunc(x, mu, sigma, L)
}

plot(xGrid, exp(distribution))

##
library('rstan')
censModel <- '
data {
  int<lower=0> T;       // Total number of time points
  int<lower=0> T_cens;  // Number of censored time points
  real x[T];            // Partly censored data
  real<upper=max(x)> L; // Lower truncation point
}
parameters {
  real mu;
  real<lower=0> sigma;
  real<upper=L> x_cens[T_cens]; // Censored values
}
model {
  int t_cens = 0;
  for (t in 1:T){
    if (x[t] > L) 
      x[t] ~ normal(mu,sigma);
    else {
      t_cens += 1;
      x_cens[t_cens] ~ normal(mu,sigma);
    }
  }
}
'


options(mc.cores = parallel::detectCores())
stanObj <- stan(file = 'mystan.stan', data=censModel, iter = 1000)
traceplot(stanObj)

