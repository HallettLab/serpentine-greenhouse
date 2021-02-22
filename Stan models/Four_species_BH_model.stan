// Beverton-Holt growth model

data{
  int<lower = 1> N;
  int Fecundity[N];
  vector[N] intra;
  vector[N] pler;
  vector[N] brho;
  vector[N] lapl;
  vector[N] femi;
}

parameters{
  real<lower = 0> lambda;
  real alpha_pler;
  real alpha_brho;
  real alpha_lapl;
  real alpha_femi;

  //real<lower = 0>  alpha_pler;
  //real<lower = 0>  alpha_brho;
  //real<lower = 0>  alpha_lapl;
  //real<lower = 0>  alpha_femi;


}

model{
  // create a vector of predictions
  vector[N] F_hat;

  // set priors
  alpha_pler ~ normal(0, 1000);
  alpha_brho ~ normal(0, 1000);
  alpha_lapl ~ normal(0, 1000);
  alpha_femi ~ normal(0, 1000);
  lambda ~ gamma(0.001, 0.001);


  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda*intra[i] / (1 + alpha_pler*pler[i] + alpha_brho*brho[i]+ alpha_lapl*lapl[i]+ alpha_femi*femi[i]);
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat);
}





