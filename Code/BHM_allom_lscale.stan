data{
  int<lower=0> N;
  real log_x[N];
  real log_y[N];
  
  //Measurment error of x measurements (assumed known for now)
  real<lower=0> tau;
  
  //Vector of 0/1 for control/drought
  vector[N] drought;
  
  //Vector of accession numbers from 1-M
  int accession[N];
  int<lower=0> M;
}


parameters{

  
  real<lower=0> sigma;
  
  //Accession level parameters
  vector[M] beta_0;
  vector[M] beta_drought;
  vector[M] beta_x;
  vector[M] beta_interaction;
  
  //Population level parameters
  real pop_0;
  real pop_drought;
  real pop_x;
  real pop_interaction;
  
  real<lower=0> sigma_0;
  real<lower=0> sigma_drought;
  real<lower=0> sigma_x;
  real<lower=0> sigma_interaction;
  
  //Priors in true x measurement error model
  vector[N] true_logx;
  real x_mean;
  real<lower=0> x_sd;
}

model{
  int indexi;
  
  //Hierarchical priors
  beta_0 ~ normal(pop_0, sigma_0);
  beta_drought ~ normal(pop_drought, sigma_drought);
  beta_x ~ normal(pop_x, sigma_x);
  beta_interaction ~ normal(pop_interaction, sigma_interaction);
  
  //Measurment error model on x measurements (x-axis)
  true_logx ~ normal(x_mean, x_sd); //Model for the true x measurements (on log scale)
  log_x ~ normal(true_logx, tau); //Model for the observed x measurements around the true value
  
  //Likelihood
  for(i in 1:N){
    indexi = accession[i];
    log_y[i] ~ normal(beta_0[indexi] + beta_drought[indexi] * drought[i] + beta_x[indexi] * true_logx[i] + 
                          beta_interaction[indexi] * drought[i] * true_logx[i], sigma); //The actual allometric model
  }
  
  //Variance priors
  sigma ~ cauchy(0, 10);
  sigma_0 ~ cauchy(0, 10);
  sigma_drought ~ cauchy(0, 10);
  sigma_x ~ cauchy(0, 10);
  sigma_interaction ~ cauchy(0, 10);

}

generated quantities{
  
  vector[M] b_drought;
  vector[M] b_control;
  vector[M] k_drought;
  vector[M] k_control;
  real b_control_pop;
  real b_drought_pop;
  real k_control_pop;
  real k_drought_pop;
  
  b_control = exp(beta_0);
  b_drought = exp(beta_drought+beta_0);
  k_control = beta_x;
  k_drought = beta_interaction + k_control;
  
  b_control_pop = exp(pop_0);
  b_drought_pop = exp(pop_drought + pop_0);
  k_control_pop = pop_x;
  k_drought_pop = pop_x + pop_interaction;
  
}
