functions {
  vector ode_meat(real t, vector y, vector theta) {
  vector[2] dydt;
  // State 2: meat temperature in C
  dydt[2] = theta[2] * (y[1] - y[2]);
  // State 1: water temperature in C
  dydt[1] = theta[1] * (theta[3] - y[1]) - theta[4] * dydt[2];
  return dydt;
}

}

data {
  int<lower=1> N;
  array[N] real time;
  array[N] real temp;
  real water;
  real meat;
  real meat0;
  vector[2] prior_k_cooler;
  vector[2] bounds_k_cooler;
  vector[2] prior_k_meat;
  vector[2] bounds_k_meat;
  vector[2] prior_z;
  vector[2] prior_sigma;
  real rel_tol;
  real abs_tol;
  int max_steps;
}
transformed data {
  real meat_to_water_ratio = meat / water;
}

parameters {
  real<lower=bounds_k_cooler[1], upper=bounds_k_cooler[2]> k_cooler;
  real<lower=bounds_k_meat[1], upper=bounds_k_meat[2]> k_meat;
  real<lower=0, upper=100> z;
  real<lower=1e-2> sigma;
  vector<lower=0, upper=100>[2] y0;
}

transformed parameters {
  vector[4] theta;
  theta[1] = k_cooler;
  theta[2] = k_meat;
  theta[3] = z;
  theta[4] = meat_to_water_ratio;
}

model {
  // Integrate in time
  array[N - 1] vector[2] y = ode_rk45_tol(ode_meat, y0, time[1], time[2:N],
                                          rel_tol, abs_tol, max_steps, theta);

  // Priors
  k_cooler ~ cauchy(prior_k_cooler[1], prior_k_cooler[2]);
  k_meat ~ cauchy(prior_k_meat[1], prior_k_meat[2]);
  z ~ cauchy(prior_z[1], prior_z[2]);
  sigma ~ normal(prior_sigma[1], prior_sigma[2]);

  // Likelihood of initial conditions
  y0[1] ~ normal(temp[1], sigma);
  y0[2] ~ normal(meat0, sigma);

  // Likelihood of future states
  for (t in 1:N - 1) {
    temp[t + 1] ~ normal(y[t][1], sigma);
  }
}

generated quantities {
  array[N - 1] vector[2] y_sim = ode_rk45_tol(ode_meat, y0, time[1], time[2:N],
                                              rel_tol, abs_tol, max_steps, theta);
}