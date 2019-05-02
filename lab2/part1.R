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

matrix_x = data.matrix(x)
matrix_y = data.matrix(y)

mu0 = c(-10,100,-100)
omega0 = matrix(0.01, 3, 3)
onesM = one(3, 3)
v0 = 4
sigmasq0 = 1

e = rnorm(1, mean = 0, sd = sigmasq0)

plot(TempLinkoping)

betahat = inv((t(matrix_x)%*%matrix_x))%*%(t(matrix_x)%*%matrix_y) 

muN = inv(t(matrix_x)%*%matrix_x + omega0)%*%(t(matrix_x)%*%matrix_x%*%betahat + omega0%*%mu0)
omegaN = t(matrix_x)%*%matrix_x + omega0
vN = v0 + 3
vNsigmaN2 = v0*sigmasq0 + (t(matrix_y)%*%matrix_y + t(mu0)%*%omega0%*%mu0 - t(muN)%*%omegaN%*%muN)

randomSigma2 <- 
randomBeta <- rnorm(n=10,mean = muN, sigmas)


