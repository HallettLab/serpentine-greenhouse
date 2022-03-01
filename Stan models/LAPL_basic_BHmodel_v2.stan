// Beverton-Holt growth model

data{
  int<lower = 1> N;
  int<lower = 1>P;
  int Fecundity [N];
  int Plot[N];
}

parameters{
  real epsilon[P];
  real<lower = 0> sigma;
  real<lower = 0> lambda;


  //real<lower = 0>  alpha_pler;
  //real<lower = 0>  alpha_brho;
  //real<lower = 0>  alpha_lapl;
  //real<lower = 0>  alpha_femi;


}

model{
  // create a vector of predictions
  vector[N] F_hat;
  vector[N] F_hat2;

  // set priors
  lambda ~ gamma(0.001, 0.001);



  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda; 
    F_hat2[i] = F_hat[i]*epsilon[Plot[i]];
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat2);
}
