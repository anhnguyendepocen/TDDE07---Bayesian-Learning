library(mvtnorm)
library(readr)
library(matlib)
library(LaplacesDemon)

TempLinkoping <- read_csv("TempLinkoping.csv")
#View(TempLinkoping)

### part a), 
## In this part we tested the hyperparameters for our prior so that the lines becomes sort of how we think they will be
## We performed normal regression-method to get the betas and then used a invchisq-distribution to draw sigma values,
## these betas and sigmas was then used to draw from a multivariate normal dist betas for our line. 
## total of 100 lines were drawn to showcase the performance of the priors
##


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


betahat = solve((t(matrix_x)%*%matrix_x))%*%(t(matrix_x)%*%matrix_y) 
randomSigma2 <- rinvchisq(n = 100, df = v0, scale = sigmasq0)
randomBetas <- c()
plot(TempLinkoping, col="black") 
for(singleSigma in randomSigma2) {
  randomBeta <- rmvt(n = 1,mu = t(mu0), S = singleSigma*inv(omega0))
  ys = c()
  xs = cbind(matrix(1, 1000, 1), matrix(1:1000, 1000, 1)/1000, (matrix(1:1000, 1000, 1)/1000)^2)
  ys = xs%*%randomBeta[1, ]
  lines(matrix(1:1000, 1000, 1)/1000, ys, col="red") 
}


### part b)
## Used the data together with the prior made in a) to compute a posterior distribution of the betas, and then plotted the median
## and the 95% credible interval of the line


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
## derivated the regression line by normal derivation and wanted the derivative to be at 0, 
## this gave the x-value for the maximum point
timeOptimum = -betas1/(2*betas2)
hist(timeOptimum, nclass=30)

### part d

## If there is a higher order but we are more certain that higher order parameters are not needed we can set the omega-values high as it creates a stronger prior 
## and the mu-prior values at 0.


