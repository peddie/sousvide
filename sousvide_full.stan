// New model:
//
//  - Multiple meats parameterised by thickness and volume
//
//  - Ambient temperature must be measured concurrently, using the
//    same sensor.  This removes the offset parameter `z`.
//
//  - Linear interpolation between ambient temperatures
//
//  - Accepts a separate reference time series for jointly inferring
//    the cooler time constant, since the constant becomes hard to
//    observe when you constantly top up the water.

functions {
  real linterp(real t, vector t_bounds, vector bounds) {
    real frac = (t - t_bounds[1]) / (t_bounds[2] - t_bounds[1]);
    return bounds[1] + frac * (bounds[2] - bounds[1]);
  }

  vector ode_meat(real t, vector y, real k_cooler,
                  vector t_bounds_ambient,
                  vector bounds_ambient,
                  int N_meat,
                  vector k_meat,
                  vector mass_ratio_meat) {
    vector[1 + N_meat] dydt;

    real ambient = linterp(t, t_bounds_ambient, bounds_ambient);

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

  vector ode_ref(real t, vector y, real k_cooler,
                 vector t_bounds_ambient, vector bounds_ambient) {
    vector[1] dydt;
    real ambient = linterp(t, t_bounds_ambient, bounds_ambient);
    dydt[1] = (ambient - y[1]) / k_cooler;
    return dydt;
  }
}

data {
  // Constants
  int<lower=1> T;
  int<lower=1> N_meat;
  real water_kg;
  vector[N_meat] meat_kg;

  // Measurement info
  array[T] real time;
  array[T] real ambient;
  array[T] vector[1 + N_meat] temps;
  array[T] real add;

  // Priors and boundaries
  vector[2] prior_k_cooler;
  vector[2] bounds_k_cooler;
  vector[2] prior_k_meat;
  vector[2] bounds_k_meat;
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
  vector<lower=bounds_k_meat[1], upper=bounds_k_meat[2]>[N_meat] k_meat;
  real<lower=1e-2> sigma;
  vector<lower=0, upper=100>[1 + N_meat] y0;
  real<lower=0, upper=100> y0_ref;
}

transformed parameters {
}

model {
  // Meat integration
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

  // Reference integration
  vector[1] ylast_ref;
  ylast_ref[1] = y0_ref;
  for (t in 2:T_ref) {
    array[1] real t1;
    t1[1] = time[t];
    array[1] vector[1] y = ode_rk45_tol(ode_ref,
                                        ylast_ref, time[t - 1], t1,
                                        rel_tol, abs_tol, max_steps,
                                        k_cooler,
                                        t_bounds_ambient_ref[t - 1],
                                        bounds_ambient_ref[t - 1]);
    temp_ref[t] ~ normal(y[1][1], sigma);
  }

  // Priors
  k_cooler ~ cauchy(prior_k_cooler[1], prior_k_cooler[2]);
  k_meat ~ cauchy(prior_k_meat[1], prior_k_meat[2]);
  sigma ~ inv_gamma(prior_sigma[1], prior_sigma[2]);

  // Likelihood of initial conditions
  y0 ~ normal(temps[1], sigma);
  y0_ref ~ normal(temp_ref[1], sigma);
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
