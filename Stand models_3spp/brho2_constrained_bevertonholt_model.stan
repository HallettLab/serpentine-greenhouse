// Beverton-Holt growth model

data{
  int<lower = 1> N;
  int Fecundity[N];
  vector[N] intra;
  vector[N] pler;
  vector[N] brho;
  vector[N] lapl;
  int<lower = 1>P;
  int Plot[N];
  real bg;
  real pg;
  real lg;
}

parameters{
  real epsilon[P];
  real<lower = 0> sigma;
  real<lower = 0> lambda;
  real alpha_pler;
  real alpha_brho;
  real alpha_lapl;

  //real<lower = 0>  alpha_pler;
  //real<lower = 0>  alpha_brho;
  //real<lower = 0>  alpha_lapl;



}

model{
  // create a vector of predictions
  vector[N] F_hat;
  vector[N] F_hat2;

  // set priors
  sigma ~ gamma(0.001, 0.001);
  epsilon ~ gamma(sigma, sigma);
  alpha_pler ~ normal(0, 1);
  alpha_brho ~ normal(0, 1);
  alpha_lapl ~ normal(0, 1);
  //lambda ~ normal(0, 1000);
  lambda ~ gamma(0.001, 0.001);



  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda*intra[i]*bg / (1 + alpha_pler*pler[i]*pg + alpha_brho*brho[i]*bg+ alpha_lapl*lapl[i]*lg);
    F_hat2[i] = F_hat[i]*epsilon[Plot[i]];
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat2);
}





