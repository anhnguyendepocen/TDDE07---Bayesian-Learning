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