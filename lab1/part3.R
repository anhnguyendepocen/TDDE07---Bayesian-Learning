#Part 3
mu = 2.39
obsRad = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

kGrid <- seq(0.001, 10, by=0.001)

priorDist <- dexp(kGrid)
plot(kGrid, priorDist, )
kappa = kGrid

dvonmises <- function(y, mu, kappa){
  res = 1
  for (i in y){
    res = res*(exp(kappa*cos(i-mu))/(2*pi*besselI(kappa,0)))
  }
  return(res)
}

likeliHood <- dvonmises(obsRad, mu, kappa)
plot(kGrid, likeliHood/(sum(likeliHood)*0.001))
posterior = likeliHood*priorDist

posterior = (posterior/sum(posterior))/0.001

plot(kGrid, posterior)

md = Mode(posterior)
index = which(posterior == max(posterior))
k = kGrid[index]
print(k)