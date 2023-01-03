functions {
  real linterp(real t, vector t_bounds, vector bounds) {
    real frac = (t - t_bounds[1]) / (t_bounds[2] - t_bounds[1]);
    return bounds[1] + frac * (bounds[2] - bounds[1]);
  }

  vector ode_ref(real t, vector y, real k_cooler,
                 vector t_bounds_ambient, vector bounds_ambient) {
    vector[1] dydt;
    real ambient = linterp(t, t_bounds_ambient, bounds_ambient);
    dydt[1] = (ambient - y[1]) / k_cooler;
    return dydt;
  }
}

data {
  // Priors and boundaries
  vector[2] prior_k_cooler;
  vector[2] bounds_k_cooler;
  vector[2] prior_sigma;

  // Reference time series
  int<lower=1> T_ref;
  array[T_ref] real time_ref;
  array[T_ref] real ambient_ref;
  array[T_ref] real temp_ref;

  // Integrator config
  real rel_tol;
  real abs_tol;
  int max_steps;
}

transformed data {
  array[T_ref - 1] vector[2] t_bounds_ambient_ref;
  array[T_ref - 1] vector[2] bounds_ambient_ref;
  for (t in 1:T_ref - 1) {
    t_bounds_ambient_ref[t][1] = time_ref[t];
    t_bounds_ambient_ref[t][2] = time_ref[t + 1];

    bounds_ambient_ref[t][1] = ambient_ref[t];
    bounds_ambient_ref[t][2] = ambient_ref[t + 1];
  }
}

parameters {
  real<lower=bounds_k_cooler[1], upper=bounds_k_cooler[2]> k_cooler;
  real<lower=1e-2> sigma;
  real<lower=0, upper=100> y0_ref;
}

transformed parameters {
}

model {
  // Reference integration
  vector[1] ylast_ref;
  ylast_ref[1] = y0_ref;
  for (t in 2:T_ref) {
    array[1] real t1;
    t1[1] = time_ref[t];
    array[1] vector[1] y = ode_rk45_tol(ode_ref,
                                        ylast_ref, time_ref[t - 1], t1,
                                        rel_tol, abs_tol, max_steps,
                                        k_cooler,
                                        t_bounds_ambient_ref[t - 1],
                                        bounds_ambient_ref[t - 1]);
    temp_ref[t] ~ normal(y[1][1], sigma);
  }

  // Priors
  k_cooler ~ cauchy(prior_k_cooler[1], prior_k_cooler[2]);
  sigma ~ inv_gamma(prior_sigma[1], prior_sigma[2]);

  // Likelihood of initial conditions
  y0_ref ~ normal(temp_ref[1], sigma);
}

generated quantities {
  array[T_ref] real y_sim_ref;
  y_sim_ref[1] = y0_ref;
  for (t in 2:T_ref) {
    array[1] real t1;
    t1[1] = time_ref[t];
    vector[1] y_last;
    y_last[1] = y_sim_ref[t - 1];
    array[1] vector[1] yout = ode_rk45_tol(ode_ref, y_last, time_ref[t - 1], t1,
                            rel_tol, abs_tol, max_steps,
                            k_cooler,
                            t_bounds_ambient_ref[t - 1],
                            bounds_ambient_ref[t - 1]);
    y_sim_ref[t] = yout[1][1];
  }
}
