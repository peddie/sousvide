library(tidyverse)
library(brms)
library(loo)
library(lubridate)
library(cmdstanr)
library(bayesplot)
library(invgamma)

times <- c(7 * 60 + 35,
           8 * 60,
           11 * 60 + 35,
           15 * 60 + 44,
           16 * 60 + 36,
           17 * 60 + 42)

measurements <- tibble(
    time = times - times[1],
    temp = c(54.5,
             53.6,
             45.8,
             40.4,
             39.4,
             38.3))

priors <- c(
    brms::prior(cauchy(25, 5), class = b, nlpar = a, lb = 0),
    brms::prior(cauchy(1e3, 1e2), class = b, nlpar = k, lb = 0),
    brms::prior(cauchy(28, 1), class = b, nlpar = z, lb = 0)
)

formula <- brms::bf(temp ~ a * exp(-time / k) + z,
                    a + k + z ~ 1,
                    nl = TRUE)

fit <- brms::brm(formula = formula,
                 data = measurements,
                 prior = priors,
                 chains = 8,
                 cores = 8,
                 control = list(max_treedepth = 15,
                                adapt_delta = 0.99),
                 iter = 5000)

summary(fit)
pairs(fit)
loo(fit, moment_match = TRUE)

draws <- fit$fit %>% as.data.frame %>% tibble

## 350g of meat
## 10L cooler
## Tap temperature ~50C
## Target temperature ~58C
t_kettle <- 98
t_target <- 57
t_tap <- 50
capacity_kg <- 10
meat_kg <- 0.35
water_kg <- capacity_kg - meat_kg

## t_tap * frac + t_kettle * (1 - frac) = t_target
frac <- (t_target - t_kettle) / (t_tap - t_kettle)
tap_fill <- water * frac
kettle_fill <- water * (1 - frac)

print(sprintf(
    "Fill cooler to %.2fL with tap water and boil %.2fL in the kettle.", 
    tap_fill,
    kettle_fill))

work <- NULL
work$temps <- c(62.3,
                60.1,
                58.8,
                57.6,
                56.5)
work$time_strings <- 
    c("15:25",
      "15:51",
      "16:20",
      "16:51",
      "17:23")
work$times_raw <- work$time_strings %>%
    lubridate::fast_strptime("%H:%M") %>%
    (\(x) lubridate::hour(x) * 60 + lubridate::minute(x))
work$times <- work$times_raw - work$times_raw[1]
work$measurements <- tibble(
    time = work$times,
    temp = work$temps)

priors <- c(
    brms::prior(cauchy(35, 1), class = b, nlpar = a, lb = 0),
    brms::prior(cauchy(1e3, 5e1), class = b, nlpar = k, lb = 0),
    brms::prior(cauchy(28, 1e-1), class = b, nlpar = z, lb = 0)
)

formula <- brms::bf(temp ~ z + a * exp(-time / k),
                    a + k + z ~ 1,
                    nl = TRUE)

fit <- brms::brm(formula = formula,
                 data = work$measurements,
                 prior = priors,
                 chains = 8,
                 cores = 8,
                 control = list(max_treedepth = 15,
                                adapt_delta = 0.99),
                 iter = 5000)

summary(fit)
pairs(fit)
loo(fit, moment_match = TRUE)

draws <- fit$fit %>% as.data.frame %>% tibble

## Full ODE model

meat <- NULL
meat$meat0 <- 10
meat$time <- work$times
meat$temp <- work$temps
meat$N <- length(meat$time)
meat$water <- water_kg
meat$meat <- meat_kg
## These values are inverse time constants, so invert to get time
## constants in minutes.
meat$prior_k_cooler <- c(1e-3, 3e-4);
meat$bounds_k_cooler <- c(1e-4, 1e-2);
meat$prior_k_meat <- c(2e-3, 1e-5);
meat$bounds_k_meat <- c(1e-4, 1e-1);
meat$prior_z <- c(28, 1e-1);
meat$prior_sigma <- c(0, 1);
meat$rel_tol <- 1e-8
meat$abs_tol <- 1e-8
meat$max_steps <- 5000

model <- cmdstanr::cmdstan_model("sousvide.stan")

fit <- model$sample(data = meat,
                    seed = 2222,
                    adapt_delta = 0.99,
                    max_treedepth = 16,
                    chains = 4,
                    parallel_chains = 4)

fit$summary()
bayesplot::mcmc_pairs(fit$draws(),
                      pars = vars(-contains("y_sim"),
                                  -contains("lp__")),
                      np = nuts_params(fit)## ,
                      ## condition = pairs_condition(nuts = "divergent__")
                      )

## Simulation

meat_sim <- meat
meat_sim$theta <- c(1e-3, 1e-1, 28, meat_sim$meat / meat_sim$water)

model <- cmdstanr::cmdstan_model("sousvide_sim.stan")
fit <- model$sample(data = meat_sim,
                    fixed_param = TRUE,
                    seed = 2222,
                    chains = 1)



## Raw data from 2-meat run
convert_time_strings <- function(time_strings) {
    times_raw <- time_strings %>%
        lubridate::fast_strptime("%H:%M") %>%
        (\(x) lubridate::hour(x) * 60 + lubridate::minute(x))
    return(times_raw - times_raw[1])
}

plot_invgamma <- function(prior, bounds) {
    dists <- seq(bounds[1], bounds[2], length.out = 200)
    plot(dists,
         invgamma::dinvgamma(dists,
                             prior[1],
                             scale = prior[2]),
         type = "l")
}

plot_cauchy <- function(prior, bounds) {
    distk <- seq(bounds[1], bounds[2], length.out = 200)
    plot(distk,
         dcauchy(distk,
                 prior[1],
                 prior[2]),
         type = "l")
}

reference_run <- tibble(time=convert_time_strings(c("15:10", "17:11", "18:13", "19:23")),
                        temp=c(75.3, 69.0, 67.0, 64.0),
                        ambient=seq(27.5, 24.0, length.out = 4))

full_run <- tibble(time=convert_time_strings(c("14:33", "14:54", "15:14", "16:10", "16:48", "17:27")),
                   ambient=c(28, 27.8, 27.7, 27.5, 27.2, 27.0),
                   ## add=c(0, 0, 4.1, 0, 4.8, 0),
                   ## water=c(58.4, 55.4, 51.0, 53.7, 52.7, 55.3),
                   ## thin=c(6.5, 54.0, 53.0, 53.2, 52.5, 55.0),
                   add=c(0, 0, 2.1, 0, 4.8, 0),
                   water=c(58.4, 55.4, 53.0, 53.7, 52.7, 55.3),
                   thin=c(6.5, 54.0, 51.0, 53.2, 52.5, 55.0),
                   thick=c(0.7, 42.8, 45.6, 52.2, 52.5, 54.3),
                   )

added_to_water <- full_run %>%
    dplyr::filter(add != 0) %>%
    dplyr::mutate(time = time + 0.5, water = water + add, add = NULL)

dplyr::bind_rows(full_run, added_to_water) %>%
    dplyr::select(-c(add)) %>%
    tidyr::pivot_longer(cols = -c(time), names_to = "type", values_to = "temp") %>%
    ggplot(aes(x = time)) +
    geom_line(aes(y = temp, color=type)) +
    labs(x = "Time [minutes]",
         y = "Temperature [C]",
         title = "Sous vide with two meats")

## Constants
full_meat <- NULL
full_meat$capacity_kg <- 10
full_meat$meat_kg <- c(0.350, 0.500)  ## Fairly approximate
full_meat$water_kg <- full_meat$capacity_kg - sum(full_meat$meat_kg)
full_meat$thickness_cm <- c(0.4, 3)

## Reformulate tablular data for Stan
full_meat$time <- full_run$time
full_meat$T <- length(full_meat$time)
full_meat$ambient <- full_run$ambient
full_meat$add <- full_run$add
full_meat$temps <- full_run %>%
    dplyr::select(-c(ambient, time, add)) %>%
    as.matrix
full_meat$N_meat <- ncol(full_meat$temps) - 1

## Reference time series
full_meat$time_ref <- reference_run$time
full_meat$T_ref <- length(full_meat$time_ref)
full_meat$ambient_ref <- reference_run$ambient
full_meat$temp_ref <- reference_run$temp

## time constants in minutes
## Cauchy priors
##
## A mystery remains here: without a relatively tight prior and a high
## lower bound on the cooler time constant, somehow the posterior
## clusters much shorter, around 250 minutes.  The reference model
## alone gives the expected results (closer to 1000 minutes).  Why is
## this happening?
full_meat$prior_k_cooler <- c(1000, 50);
full_meat$bounds_k_cooler <- c(700, 2000);
full_meat$prior_k_meat <- c(50, 20);
full_meat$bounds_k_meat <- c(1, 100);
## Inverse gamma prior
full_meat$prior_sigma <- c(3, 1);
full_meat$rel_tol <- 1e-5
full_meat$abs_tol <- 1e-5
full_meat$max_steps <- 5000

plot_invgamma(full_meat$prior_sigma, c(1e-2, 10))

plot_cauchy(full_meat$prior_k_meat, full_meat$bounds_k_meat)

plot_cauchy(full_meat$prior_k_cooler, full_meat$bounds_k_cooler)

ref_model <- cmdstanr::cmdstan_model("sousvide_ref.stan")
ref_fit <- ref_model$sample(data = full_meat,
                            seed = 2222,
                            chains = 4,
                            parallel_chains = 4)
ref_fit$summary()
bayesplot::mcmc_pairs(ref_fit$draws(),
                      pars = vars(-contains("y_sim"),
                                  -contains("lp__")),
                      np = nuts_params(ref_fit)## ,
                      ## condition = pairs_condition(nuts = "divergent__")
                      )

model <- cmdstanr::cmdstan_model("sousvide_full.stan")
fit <- model$sample(data = full_meat,
                    seed = 2222,
                    chains = 4,
                    parallel_chains = 4)

fit$summary()
bayesplot::mcmc_pairs(fit$draws(),
                      pars = vars(-contains("y_sim"),
                                  -contains("lp__")),
                      np = nuts_params(fit)## ,
                      ## condition = pairs_condition(nuts = "divergent__")
                      )



simulated <- fit$draws("y_sim", format="df") %>% tibble

tidy_posterior <- function(simulated, index, label) {
    sim_draws <- simulated %>%
        dplyr::select(starts_with("y_sim[") & ends_with(paste0(index, "]"))) %>%
        data.frame %>%
        t

    colnames(sim_draws) <- paste0("draw_", seq(1, ncol(sim_draws)))

    draw_table <- as_tibble(sim_draws)

    draw_table$time <- full_meat$time
    return(draw_table %>% 
        tidyr::pivot_longer(cols = -c(time),
                            names_to = "draw",
                            values_to = label))
}

specs <- c(1 = "water", 2 = "thin", 3 = "thick")

posterior <- tidy_posterior(simulated, 1, "water") %>%
    dplyr::inner_join(tidy_posterior(simulated, 2, "thin"), by = c("time", "draw")) %>%
    dplyr::inner_join(tidy_posterior(simulated, 3, "thick"), by = c("time", "draw")) 
    
posterior %>%
    ggplot(aes(x = time, group = draw)) +
    geom_line(aes(y = water),
              alpha = 0.01,
              color = "red") +
    geom_line(aes(y = thin),
              alpha = 0.01,
              color = "green") +
    geom_line(aes(y = thick),
              alpha = 0.01,
              color = "blue") +
    labs(x = "Time [minutes]",
         y = "Temperature [C]",
         title = paste("Sous vide"))
