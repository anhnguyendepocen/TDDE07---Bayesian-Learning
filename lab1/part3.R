#Part 3
mu = 2.39
obsRad = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

kGrid <- seq(0.01, 20, by=0.01)

priorDist <- rexp(kGrid, 1)
plot(priorDist)
kappa = kGrid

dvonmises <- function(y, mu, kappa){
  res = 1
  for (i in y){
      res = res*(exp(kappa*cos(i-mu))/(2*pi*besselI(kappa,1)))
  }
  return(res)
}

likeliHood <- dvonmises(obsRad, mu, kappa)
plot(likeliHood)
posterior = likeliHood*priorDist

posterior = (posterior/sum(posterior))/0.01

plot(posterior)
likeliHood
