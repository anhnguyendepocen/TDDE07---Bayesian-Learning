library(manipulate)
library(LaplacesDemon)
u = 3.5
y = c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
logY = log(y)
meanY = sum(y)/length(y)
meanLogY = sum(logY)/length(logY)

getCdfValue <- function(y, target){
  sortedY = sort(y)
  n = length(sortedY)
  i = 0
  for(dataPoint in sortedY){
    i = i + 1
    if(dataPoint >= target){
      return(i/n)
    }
  }
}

getDistrubutionFromPoints <- function(x){
  n = length(x)
  sortedX = sort(x)
  cumPercent = (1:n)/n
  densityPercent = (cumPercent[2:n] - cumPercent[1:n-1])/(sortedX[2:n] - sortedX[1:(n-1)])
  #dataMatrix = matrix(1:(2*n-2), nrow=n-1, ncol=2, byrow=TRUE)
  #dataMatrix[ ,1] = sortedX[2:n]
  #dataMatrix[ ,2] = densityPercent
  
  dataMatrix = matrix(1:(2*n), nrow=n, ncol=2, byrow=TRUE)
  dataMatrix[ ,1] = sortedX
  dataMatrix[ ,2] = cumPercent
  
  return(dataMatrix)
}

getIntervalData <- function(x, lowerPercent, upperPercent){
  lowerFlag = FALSE
  upperFlag = FALSE
  n = length(x)
  sortedX = sort(x)
  i = 0
  for(data in x) {
    i = i + 1
    if(i/n >= lowerPercent && !lowerFlag) {
      lowerFlag = TRUE
      lowerIndex = i
    }
    else if(i/n > upperPercent){
      upperIndex = i - 1
      upperFlag = TRUE
      break;
    }
    
  }
  if(!upperFlag){
    upperIndex = n
  }
  return(sortedX[lowerIndex:upperIndex])
}


drawFromPosterior <- function(mean, n, logy){
  xGrid <- seq(0.001, 0.999, by=0.001)
  tau2 = (sum(logy - mean)^2)/n
  posterior = dinvchisq(xGrid, n, tau2, log=FALSE)
  
  maxDensity <- max(posterior)
  
  posteriorRandomDraw <- rinvchisq(10000,  n, tau2)
  
  #plotting
  attach(mtcars)
  par(mfrow=c(2,2))
  plot(xGrid, posterior, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "theta", 
       ylab = 'Density', main = 'Chi Posterior of Variance')
  plot(sort(posteriorRandomDraw))
  
  #part b)
  giniTheory = 2*pnorm(posteriorRandomDraw/sqrt(2), mean = 0, sd = 1) - 1
  
  logGini = log(giniTheory)
  distrubtionGiniLog = getDistrubutionFromPoints(logGini)
  distrubtionGini = getDistrubutionFromPoints(giniTheory)
  plot(distrubtionGiniLog)
  plot(distrubtionGini)
  
  #part c)
  lowerData = getIntervalData(giniTheory, 0, 0.025)
  middleData = getIntervalData(giniTheory, 0.025, 0.975)
  upperData = getIntervalData(giniTheory, 0.975, 1)
  
  print(density(middleData))
  print(min(middleData))
  print(max(middleData))
}

num = length(logY)
draws = drawFromPosterior(u, num, logY)