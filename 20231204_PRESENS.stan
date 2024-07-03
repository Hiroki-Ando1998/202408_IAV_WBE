
data {
  int<lower=0> n;
  vector[n] IAV;
  int<lower=0> na;
  vector[na] WBE;
  int pick_colum [na];
  int<lower=0> nb;
  vector[nb] t;
}


parameters {
  real<lower=0> A;
  real<lower=0> a;
  real<lower=a> b;
  real<lower=0> sigma;
}


transformed parameters {
  vector[nb] feces;
  for (i in 1:nb){
    feces[i] = A * a / (b - a) * exp(-a * t[i]) * (1 - exp((a - b) * t[i]));
    }
}
//21

model {
  A ~ normal(1000, 1000);
  a ~ normal(0.3, 3);
  b ~ normal(20, 20); 
  sigma ~ normal(0, 10);
  
  vector[na] WBE_estimated;
  int x [na];

  for (i in 1:na){
   x[i] = pick_colum[i];
   
   WBE_estimated[i] = IAV[x[i]] * feces[1] + IAV[x[i] - 1] * feces[2] + IAV[x[i] - 2] * feces[3] + IAV[x[i] - 3] * feces[4] + IAV[x[i] - 4] * feces[5] + IAV[x[i] - 5] * feces[6] + IAV[x[i] - 6] * feces[7]
   + IAV[x[i] - 7] * feces[8] +  IAV[x[i] - 8] * feces[9] + IAV[x[i] - 9] * feces[10] + IAV[x[i] - 10] * feces[11] + IAV[x[i] - 11] * feces[12]
   + IAV[x[i] - 12] * feces[13] + IAV[x[i] - 13] * feces[14] + IAV[x[i] - 14] * feces[15] + IAV[x[i] - 15] * feces[16] 
   + IAV[x[i] - 16] * feces[17] + IAV[x[i] - 17] * feces[18] + IAV[x[i] - 18] * feces[19] + IAV[x[i] - 19] * feces[20] + IAV[x[i] - 20] *feces[21]
   + IAV[x[i] - 21] * feces[22] + IAV[x[i] - 22] * feces[23] + IAV[x[i] - 23] * feces[24] + IAV[x[i] - 24] * feces[25] + IAV[x[i] - 25] *feces[26];
   WBE[i] ~ normal(log10(WBE_estimated[i]), sigma);
}
}


generated quantities{
  vector[na] log_lik;
  vector[na] WBE_E;
  int x [na];
  
  vector[nb] F;
  for (i in 1:nb){
    F[i] = A * a / (b - a) * exp(-a * t[i]) * (1 - exp((a - b) * t[i]));
    }
  
  for (i in 1:na){
  x[i] = pick_colum[i];
  
  WBE_E[i] = IAV[x[i]] * F[1] + IAV[x[i] - 1] * F[2] + IAV[x[i] - 2] * F[3] + IAV[x[i] - 3] * F[4] + IAV[x[i] - 4] * F[5] + IAV[x[i] - 5] * F[6] + IAV[x[i] - 6] * F[7]
   + IAV[x[i] - 7] * F[8] +  IAV[x[i] - 8] * F[9] + IAV[x[i] - 9] * F[10] + IAV[x[i] - 10] * F[11] + IAV[x[i] - 11] * F[12]
   + IAV[x[i] - 12] * F[13] + IAV[x[i] - 13] * F[14] + IAV[x[i] - 14] * F[15] + IAV[x[i] - 15] * F[16] 
   + IAV[x[i] - 16] * F[17] + IAV[x[i] - 17] * F[18] + IAV[x[i] - 18] * F[19] + IAV[x[i] - 19] * F[20] + IAV[x[i] - 20] *F[21]
   + IAV[x[i] - 21] * F[22] + IAV[x[i] - 22] * F[23] + IAV[x[i] - 23] * F[24] + IAV[x[i] - 24] * F[25] + IAV[x[i] - 25] *F[26];
  log_lik[i] = normal_lpdf(WBE[i]| log10(WBE_E[i]), sigma);
  
  }
}  

  




