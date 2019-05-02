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

mu0 = c(-10,100,-100)
omega0 = matrix(0.01, 3, 3)
v0 = 4
sigmasq0 = 1

e = rnorm(1, mean = 0, sd = sigmasq0)


betahat = inv((t(matrix_x)%*%matrix_x))%*%(t(matrix_x)%*%matrix_y) 

muN = inv(t(matrix_x)%*%matrix_x + omega0)%*%(t(matrix_x)%*%matrix_x%*%betahat + omega0%*%mu0)
omegaN = t(matrix_x)%*%matrix_x + omega0
vN = v0 + 3
vNsigmaN2 = v0*sigmasq0 + (t(matrix_y)%*%matrix_y + t(mu0)%*%omega0%*%mu0 - t(muN)%*%omegaN%*%muN)

randomSigma2 <- rinvchisq(n = 10, df = vN, scale = (vNsigmaN2/vN))
randomBetas <- c()
plot(TempLinkoping, col="black") 
for(singleSigma in randomSigma2) {
  randomBeta <- rmvt(n = 10,mu = muN, S = singleSigma*inv(omegaN))

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