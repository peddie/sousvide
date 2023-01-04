library(tidyverse)
library(lubridate)
library(cmdstanr)
library(bayesplot)
library(invgamma)

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

tidy_posterior <- function(inputs, simulated, index, label) {
    sim_draws <-
        simulated %>%
        dplyr::select(starts_with("y_sim[") & ends_with(paste0(index, "]"))) %>%
        data.frame %>%
        t

    colnames(sim_draws) <- paste0("draw_", seq(1, ncol(sim_draws)))

    draw_table <- as_tibble(sim_draws)

    full_times = c()
    for (t in 2:length(inputs$time)) {
        new <- seq(inputs$time[t - 1],
                            inputs$time[t],
                            length.out = inputs$T_subsample + 1)[-(inputs$T_subsample + 1)]
        full_times <- c(full_times, new)
    }
    full_times <- c(full_times, inputs$time[length(inputs$time)])

    draw_table$time <- full_times
    return(draw_table %>%
        tidyr::pivot_longer(cols = -c(time),
                            names_to = "draw",
                            values_to = label))
}

tidy_posterior_reference <- function(inputs, simulated, label) {
    sim_draws <-
        simulated %>%
        dplyr::select(starts_with("y_sim_ref")) %>%
        data.frame %>%
        t

    colnames(sim_draws) <- paste0("draw_", seq(1, ncol(sim_draws)))

    draw_table <- as_tibble(sim_draws)

    draw_table$time <- inputs$time_ref
    return(draw_table %>%
        tidyr::pivot_longer(cols = -c(time),
                            names_to = "draw",
                            values_to = label))
}

lengthen <- function(data) {
    return(tidyr::pivot_longer(data,
                               cols = -c(time),
                               names_to = "type",
                               values_to = "temp"))
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

region_colors <-
    c("water" = "red",
      "thin" = "green",
      "thick" = "blue",
      "ambient" = "purple",
      "reference" = "black")

measurement_points <- full_run %>%
    dplyr::select(c(water, time, add)) %>%
    dplyr::filter(add != 0) %>%
    dplyr::mutate(time = time + 0.5, water = water + add, add = NULL) %>% lengthen %>%
    dplyr::bind_rows(full_run %>% dplyr::select(-c(add))  %>% lengthen)

measurement_points %>%
    ggplot(aes(x = time)) +
    geom_line(aes(y = temp, color=type)) +
    labs(x = "Time [minutes]",
         y = "Temperature [C]",
         title = "Sous vide with two meats")

full_meat <- NULL

## Constants
full_meat$capacity_kg <- 10
full_meat$meat_kg <- c(0.350, 0.7)  ## Fairly approximate
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

## Simulation
full_meat$T_subsample <- 50

## time constants in minutes
## Cauchy priors
full_meat$prior_k_cooler <- c(1000, 300);
full_meat$bounds_k_cooler <- c(100, 2000);
full_meat$prior_k_meat <- c(50, 20);
full_meat$bounds_k_meat <- c(1, 100);
## Inverse gamma prior
full_meat$prior_sigma <- c(3, 1);
full_meat$rel_tol <- 1e-5
full_meat$abs_tol <- 1e-5
full_meat$max_steps <- 5000

## pdf("sousvide.pdf", width=30, height=20)

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
mk_tidy <- function(index, label) { return(tidy_posterior(full_meat, simulated, index, label)) }
posterior <- mk_tidy(1, "water") %>%
    dplyr::inner_join(mk_tidy(2, "thin"), by = c("time", "draw")) %>%
    dplyr::inner_join(mk_tidy(3, "thick"), by = c("time", "draw"))

means <- fit$draws(format = "df") %>%
    as_tibble %>%
    dplyr::select(starts_with("k_")) %>%
    colMeans

posterior_line_alpha <- 4e-3
posterior %>%
    ggplot(aes(x = time, group = draw)) +
    geom_line(aes(y = water,
                  color = "water"),
              alpha = posterior_line_alpha) +
    geom_line(aes(y = thin,
                  color = "thin"),
              alpha = posterior_line_alpha) +
    geom_line(aes(y = thick,
                  color = "thick"),
              alpha = posterior_line_alpha) +
    geom_point(data = measurement_points,
               inherit.aes = FALSE,
               aes(x = time,
                   y = temp,
                   color = type,
                   shape = type)) +
    scale_color_manual(values = region_colors) +
    guides(color = guide_legend(
               override.aes = list(alpha = 1),
               title = "Region"),
           shape = guide_legend(
               title = "Region")) +
    labs(x = "Time [minutes]",
         y = "Temperature [C]",
         title = paste("Sous vide with", full_meat$N_meat, "meat cuts"),
         subtitle = "lines are draws from posterior; points are measurements",
         caption = sprintf(
             "time constant means [minutes]\nk_cooler: %.0f\nk_thin: %.2f\nk_thick: %.2f",
             means["k_cooler"], means["k_meat[1]"], means["k_meat[2]"])
         )

ref_means <- ref_fit$draws(format = "df") %>%
    as_tibble %>%
    dplyr::select(starts_with("k_")) %>%
    colMeans
ref <- tidy_posterior_reference(full_meat,
                                ref_fit$draws("y_sim_ref", format="df") %>% tibble,
                                "reference")
ref %>%
    ggplot(aes(x = time, group = draw)) +
    geom_line(aes(y = reference,
                  color = "reference"),
              alpha = posterior_line_alpha) +
    geom_point(data = reference_run,
               inherit.aes = FALSE,
               aes(x = time, y = temp)) +
    guides(color = guide_legend(
               override.aes = list(alpha = 1),
               title = "Region")) +
    scale_color_manual(values = region_colors) +
    labs(x = "Time [minutes]",
         y = "Temperature [C]",
         title = "Reference sous vide experiment",
         caption = paste("k_cooler:", ref_means["k_cooler"], "minutes")
         )


ref_joint <- tidy_posterior_reference(full_meat,
                                      fit$draws("y_sim_ref", format="df") %>% tibble,
                                      "reference")
ref_joint %>%
    ggplot(aes(x = time, group = draw)) +
    geom_line(aes(y = reference,
                  color = "reference"),
              alpha = posterior_line_alpha) +
    geom_point(data = reference_run,
               inherit.aes = FALSE,
               aes(x = time, y = temp)) +
    guides(color = guide_legend(
               override.aes = list(alpha = 1),
               title = "Region")) +
    scale_color_manual(values = region_colors) +
    labs(x = "Time [minutes]",
         y = "Temperature [C]",
         title = "Reference sous vide experiment (joint params)",
         caption = paste("k_cooler:", means["k_cooler"], "minutes")
         )

## dev.off()
