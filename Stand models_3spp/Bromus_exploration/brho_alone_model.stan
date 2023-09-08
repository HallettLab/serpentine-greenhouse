// Beverton-Holt growth model

data{
  int<lower = 1> N;
  int Fecundity[N];
  real intra[N];
  real bg;

}

parameters{
  real<lower = 0> lambda;
  //real alpha_pler;
  //real alpha_brho;
  //real alpha_lapl;
}

model{
  // create a vector of predictions
  vector[N] F_hat;

  // set priors
  //lambda ~ uniform(0, 200);
  lambda ~ gamma(0.001, 0.001);



  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda*intra[i]*bg; 
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat);
}






