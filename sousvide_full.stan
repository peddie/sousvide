// New model:
//
//  - Multiple meats parameterised by thickness and volume
//
//  - Ambient temperature must be measured concurrently, using the
//    same sensor.  This removes the offset parameter `z`.
//
//  - Linear interpolation between ambient temperatures

functions {
  vector ode_meat(real t, vector y, real k_cooler,
                  vector t_bounds_ambient,
                  vector bounds_ambient,
                  int N_meat,
                  vector k_meat,
                  vector mass_ratio_meat) {
  vector[1 + N_meat] dydt;

  real frac_ambient = (t - t_bounds_ambient[1]) /
    (t_bounds_ambient[2] - t_bounds_ambient[1]);
  real range_ambient = bounds_ambient[2] - bounds_ambient[1];
  real ambient = bounds_ambient[1] + frac_ambient * range_ambient;
  
  // Compute water cooling to ambient sink
  dydt[1] = (ambient - y[1]) / k_cooler; 

  for (m in 1:N_meat) {
    // meat temperature in C
    dydt[m + 1] = (y[1] - y[m + 1]) / k_meat[m];
    // Compute energy balance between water and meat
    dydt[1] = dydt[1] - mass_ratio_meat[m] * dydt[m + 1];
  }
  return dydt;
  }

}

data {
  int<lower=1> T;
  int<lower=1> N_meat;
  array[T] real time;
  array[T] real ambient;
  array[T] vector[1 + N_meat] temps;
  array[T] real add;
  real water_kg;
  vector[N_meat] meat_kg;
  vector[2] prior_k_cooler;
  vector[2] bounds_k_cooler;
  vector[2] prior_k_meat;
  vector[2] bounds_k_meat;
  vector[2] prior_sigma;
  real rel_tol;
  real abs_tol;
  int max_steps;
}
transformed data {
  int N_states = N_meat + 1;
  vector[N_meat] mass_ratio_meat = meat_kg / water_kg;

  // Precompute ambient temperature pairs for interpolation
  //
  // Indices are one off from timestep, because we don't need to
  // integrate to get to t0.
  array[T - 1] vector[2] t_bounds_ambient;
  array[T - 1] vector[2] bounds_ambient;
  for (t in 1:T - 1) {
    t_bounds_ambient[t][1] = time[t];
    t_bounds_ambient[t][2] = time[t + 1];

    bounds_ambient[t][1] = ambient[t];
    bounds_ambient[t][2] = ambient[t + 1];
  }
}

parameters {
  real<lower=bounds_k_cooler[1], upper=bounds_k_cooler[2]> k_cooler;
  vector<lower=bounds_k_meat[1], upper=bounds_k_meat[2]>[N_meat] k_meat;
  real<lower=1e-2> sigma;
  vector<lower=0, upper=100>[1 + N_meat] y0;
}

transformed parameters {
}

model {
  vector[N_states] ylast = y0;
  for (t in 2:T) {
    // Add hot water after the initial measurements
    ylast[1] = ylast[1] + add[t - 1];
    array[1] real t1;
    t1[1] = time[t];

    array[1] vector[N_states] y = ode_rk45_tol(ode_meat,
                                      ylast, time[t - 1], t1,
                                      rel_tol, abs_tol, max_steps,
                                      k_cooler,
                                      t_bounds_ambient[t - 1],
                                      bounds_ambient[t - 1],
                                      N_meat,
                                      k_meat,
                                      mass_ratio_meat);
    temps[t] ~ normal(y[1], sigma);
  }

  // Priors
  k_cooler ~ cauchy(prior_k_cooler[1], prior_k_cooler[2]);
  k_meat ~ cauchy(prior_k_meat[1], prior_k_meat[2]);
  sigma ~ inv_gamma(prior_sigma[1], prior_sigma[2]);

  // Likelihood of initial conditions
  y0 ~ normal(temps[1], sigma);
}

generated quantities {
  array[T] vector[N_states] y_sim;
  y_sim[1] = y0;
  for (t in 2:T) {
    y_sim[t - 1] = y_sim[t - 1] + add[t - 1];
    array[1] real t1;
    t1[1] = time[t];
    array[1] vector[N_states] yout = ode_rk45_tol(ode_meat, y_sim[t - 1], time[t - 1], t1,
                            rel_tol, abs_tol, max_steps,
                            k_cooler,
                            t_bounds_ambient[t - 1],
                            bounds_ambient[t - 1],
                            N_meat,
                            k_meat,
                            mass_ratio_meat);
    y_sim[t] = yout[1];
  }
}
