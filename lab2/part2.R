

womenWork<-read.table("WomenWork.dat",header=TRUE)  # Spam data from Hastie et al.
glmModel <- glm(Work ~ 0 + ., data = womenWork, family = binomial)

### b

chooseCov <- c(1:8) # Here we choose which covariates to include in the model
tau <- 10; # Prior scaling factor such that Prior Covariance = (tau^2)*I
###########     END USER INPUT    ################


# install.packages("mvtnorm") # Loading a package that contains the multivariate normal pdf
library("mvtnorm") # This command reads the mvtnorm package into R's memory. NOW we can use dmvnorm function.

# Loading data from file
Data<-read.table("SpamReduced.dat",header=TRUE)  # Spam data from Hastie et al.
y <- as.vector(womenWork[,1]); # Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
X <- as.matrix(womenWork[,2:9]);
covNames <- names(womenWork)[2:length(names(womenWork))];
X <- X[,chooseCov]; # Here we pick out the chosen covariates.
covNames <- covNames[chooseCov];
nPara <- dim(X)[2];

# Setting up the prior
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara);

LogPostLogistic <- function(betaVect,y,X,mu,Sigma){
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  
  logLik <- sum( linPred*y -log(1 + exp(linPred)));
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  return(logLik + logPrior)
}

initVal <- as.vector(rep(0,dim(X)[2])); 
logPost = LogPostLogistic;
OptimResults<-optim(initVal,logPost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Printing the results to the screen
names(OptimResults$par) <- covNames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian))) # Computing approximate standard deviations.
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian)))
names(approxPostStd) <- covNames # Naming the coefficient by covariates

betatilde = OptimResults$par
print("Betatilde: ")
print(betatilde)
print("Jacobiany beta: ")
print(approxPostStd)

print("Intervall NSmallChild: ")
upperb = betatilde["NSmallChild"] + 1.64*approxPostStd["NSmallChild"]
lowerb = betatilde["NSmallChild"] - 1.64*approxPostStd["NSmallChild"]
print(upperb)
print(lowerb)

### part c
covarMatrix = -inv(OptimResults$hessian)

ladyInput = c(1, 10, 8, 10, (10/10)^2, 40, 1, 1)
workOrNots = c()
betas = rmvt(n = 1000,mu = matrix(betatilde), S = covarMatrix)
for(row in 1:nrow(betas)) {
  working = exp(ladyInput %*% betas[row, ])
  workOrNots = c(workOrNots, working)
}
