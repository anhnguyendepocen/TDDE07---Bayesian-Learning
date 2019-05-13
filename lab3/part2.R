library(mvtnorm)
library(readr)
library(matlib)
library(LaplacesDemon)

rawData <- read.table("ebayNumberOfBidderData.dat",header=TRUE)

y <- as.vector(rawData[,1]); # Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
X <- rawData[,2:10];


model <- glm(y ~ PowerSeller+VerifyID+Sealed+Minblem+MajBlem+LargNeg+LogBook+MinBidShare, data=X, family=poisson())

## Significant variables are VerifyID, Sealed, MajBlem, LogBook and MinBidShare

# b)
xMatrix = as.matrix(X)
sigmaGPrior = 100*solve((t(xMatrix)%*%xMatrix))
tau = 10
chooseCov <- c(1:9)

covNames <- names(rawData)[2:length(names(rawData))];
xMatrix <- xMatrix[,chooseCov]; # Here we pick out the chosen covariates.
covNames <- covNames[chooseCov];
nPara <- dim(X)[2];

# Setting up the prior
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara);

PoiPost <- function(betaVect,y,X, mu, SigmaGPrior) {
  nPara <- length(betaVect);
  
  linPred <- X%*%betaVect;
  
  logPoiLik <- sum( linPred*y -exp(linPred) - log(factorial(y)));
  if (abs(logPoiLik) == Inf) logPoiLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logBetaPrior <- dmvnorm(betaVect, matrix(0,nPara,1), SigmaGPrior, log=TRUE);
  return(logPoiLik + logBetaPrior) 
}

initVal <- as.vector(rep(0,dim(X)[2])); 
logPost = PoiPost;
OptimResults<-optim(initVal,logPost,gr=NULL,y,xMatrix,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

approxPostStd <- sqrt(diag(-solve(OptimResults$hessian)))
names(approxPostStd) <- covNames # Naming the coefficient by covariates

betatilde = OptimResults$par
print("Betatilde: ")
print(betatilde)
print("Jacobiany beta: ")
print(approxPostStd)

### part c)

RWMSampler <- Function(theta, logPostFunc, ...) {
  
}