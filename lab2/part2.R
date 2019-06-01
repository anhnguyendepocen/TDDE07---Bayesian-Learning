### a

## used glm to fit a model to the data in womenWork, family = binomial sets it so that it fits for a logistic regression
## the 0 added makes sure that it does not add an intercept as the data already has one and the dot adds all the features in our data.

womenWork<-read.table("WomenWork.dat",header=TRUE)  # Spam data from Hastie et al.
glmModel <- glm(Work ~ 0 + ., data = womenWork, family = binomial)

### b

## In this part we defined a likelihood function which optim could optimize for to find the optimal betas
## In this function we had defined a prior and used the likelihood which was for the relevant distribution, in our case the poisson distribution
## We then used the result from optim to get the hessian and mode(where the mode is the optimal betas) and then got the std and estimate.


chooseCov <- c(1:8) # Here we choose which covariates to include in the model
tau <- 10; # Prior scaling factor such that Prior Covariance = (tau^2)*I
###########     END USER INPUT    ################


# install.packages("mvtnorm") # Loading a package that contains the multivariate normal pdf
library("mvtnorm") # This command reads the mvtnorm package into R's memory. NOW we can use dmvnorm function.

# Loading data from file
Data<-read.table("WomenWork.dat",header=TRUE)  # Spam data from Hastie et al.
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
upperb = betatilde["NSmallChild"] + 1.96*approxPostStd["NSmallChild"]
lowerb = betatilde["NSmallChild"] - 1.96*approxPostStd["NSmallChild"]
print(upperb)
print(lowerb)

### part c
## In this part we used the results from before to firstly to random draws from the normal distrubtion of our betas,
## with the estimated betas as mu and std as the one we got from b), from these draws we could make a prediction for the lady
## afterwards we stored the estimated probability that the lady work but also used the probability in a bernoulli distribution to draw
## whether the lady works or not, in order to get the binary distribution on how many cases the lady would work and how many cases she would not work.


library("Rlab")
covarMatrix = -solve(OptimResults$hessian)

ladyInput = c(1, 10, 8, 10, (10/10)^2, 40, 1, 1)
workOrNots = c()
workOrNotsBinary = c()
betas = rmvt(n = 1000, covarMatrix, matrix(betatilde))
for(row in 1:nrow(betas)) {
  working = exp(ladyInput %*% betas[row, ])/(1 + exp(ladyInput %*% betas[row, ]))
  workingBinary = rbern(1, working)
  workOrNots = c(workOrNots, working)
  workOrNotsBinary = c(workOrNotsBinary, workingBinary)
}
hist(workOrNots, n=30)
hist(workOrNotsBinary, n=30)