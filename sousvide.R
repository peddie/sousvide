library(tidyverse)
library(brms)
library(loo)
library(lubridate)
library(cmdstanr)
library(bayesplot)

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
