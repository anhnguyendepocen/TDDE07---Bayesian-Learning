
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] X;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real muRandom;
  real<lower=0> sigma2Random;
  real omegaRandom;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  //muRandom ~ normal(initialMuMu, muSigma2);
  #sigma2Random ~ scaled_inv_chi_square(sigma2Nu, sigma2Sigma);
  //sigma2Random ~ normal(sigma2Nu, sigma2Sigma);
  //omegaRandom ~ normal(initialOmegaMu, omegaSigma2);
  for(i in 2:N) {
    X[i] ~  normal(muRandom + omegaRandom*(X[i-1] - muRandom), sigma2Random);
  }
}


