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

mu0 = c(-11,92,-80)
omega0 = matrix(c(0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01), 3, 3)
v0 = 3
sigmasq0 = 0.03

e = rnorm(1, mean = 0, sd = sigmasq0)


betahat = inv((t(matrix_x)%*%matrix_x))%*%(t(matrix_x)%*%matrix_y) 
randomSigma2 <- rinvchisq(n = 10, df = v0, scale = sigmasq0)
randomBetas <- c()
plot(TempLinkoping, col="black") 
for(singleSigma in randomSigma2) {
  randomBeta <- rmvt(n = 10,mu = t(mu0), S = singleSigma*inv(omega0))
  
  for(k in 1:10) {
    ys = c()
    xs = c()
    for (i in 1:1000) {
      y = randomBeta[k, 1] + randomBeta[k, 2]*i/1000 + randomBeta[k, 3]*(i/1000)^2
      ys = c(ys, y)
      xs = c(xs, i/1000)
    }
    lines(xs, ys, col="red") 
  }
}


### part b)




muN = inv(t(matrix_x)%*%matrix_x + omega0)%*%(t(matrix_x)%*%matrix_x%*%betahat + omega0%*%mu0)
omegaN = t(matrix_x)%*%matrix_x + omega0
vN = v0 + 3
vNsigmaN2 = v0*sigmasq0 + (t(matrix_y)%*%matrix_y + t(mu0)%*%omega0%*%mu0 - t(muN)%*%omegaN%*%muN)

randomSigma2 <- rinvchisq(n = 100, df = vN, scale = (vNsigmaN2/vN))
randomBetas <- c()
plot(TempLinkoping, col="black")
sigmas = data.frame(randomSigma2)
y_df = data.frame(matrix(1, 1000, 1))
betas0 = c()
betas1 = c()
betas2 = c()
ibeta = 0
for(singleSigma in randomSigma2) {
  randomBeta <- rmvt(n = 10,mu = muN, S = singleSigma*inv(omegaN))
  for(k in 1:10) {
    ibeta = ibeta + 1
    ys = c()
    xs = c()
    for (i in 1:1000) {
      y = randomBeta[k, 1] + randomBeta[k, 2]*i/1000 + randomBeta[k, 3]*(i/1000)^2
      ys = c(ys, y)
      xs = c(xs, i/1000)
    }
    betas0 = c(betas0, randomBeta[k, 1])
    betas1 = c(betas1, randomBeta[k, 2])
    betas2 = c(betas2, randomBeta[k, 3])
    y_df[paste0("trial", ibeta)] <- data.frame(ys)
  }
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

