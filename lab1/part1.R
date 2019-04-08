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