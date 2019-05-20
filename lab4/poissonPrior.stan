
data {
  int<lower=0> N;
  int  y[N];
  real sigma2Nu;
  real sigma2Sigma;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real x[N];
  real muRandom;
  real omegaRandom;
  real sigma2Random;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  sigma2Random ~ scaled_inv_chi_square(sigma2Nu, sigma2Sigma);
  //sigma2Random ~ normal(sigma2Nu, sigma2Sigma);
  
  
  x[1] ~ normal(muRandom, sigma2Random);
  y[1] ~ poisson(exp(x[1]));
  for(i in 2:N) {
    x[i] ~ normal(muRandom + omegaRandom*(x[i-1] - muRandom), sigma2Random);
    y[i] ~ poisson(exp(x[i]));
  }
}
