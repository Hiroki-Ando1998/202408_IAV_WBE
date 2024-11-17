//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n;
  vector[n] p;
  vector[n] IAV;
  vector[n] nos;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real c1;
  real intercept;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  intercept ~ normal(0, 100);
  c1 ~ normal(0, 100);

  vector[n] x1;
  vector[n] x2;
  vector[n] x3;

  for (i in 1:n){
  x1[i] = exp(intercept + c1*IAV[i]);
  x2[i] = x1[i] / (1 + x1[i]);
  x3[i] = sqrt(x2[i]*(1-(x2[i]))/nos[i]);
  p[i] ~ normal(x2[i], x3[i]);
  }
}

generated quantities{
  vector[n] log_lik;
  vector[n] x4;
  vector[n] x5;
  vector[n] x6;
  
  for (i in 1:n){
  x4[i] = exp(intercept + c1*IAV[i]);
  x5[i] = x4[i] / (1 + x4[i]);
  x6[i] = sqrt(x5[i]*(1-(x5[i]))/n);
  
  log_lik[i] = normal_lpdf(p[i]| x5[i], x6[i]);
}
}



