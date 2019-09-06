####################################
#                                  #
# Warning: The script might run a  #
# while depending on your machine! #
#                                  #
####################################

# -----------------
# Loading Libraries
# -----------------

library(magrittr)
library(dplyr)
library(purrr)
library(tibble)
library(reshape2)
library(ggplot2)
library(coda)
library(devtools)
install_github("sbfnk/fitR")
library(fitR)

# -----------------------
# Calculating growth rate
# -----------------------

start <- c(-15000, 0.01)
end <- c(0, 2.8)

prior_data <- data.frame(rbind(start, end))

colnames(prior_data) <- c("BCE", "dens")

new_data <- data.frame(BCE = seq(-15000, 0, by = 1))
model_base <- lm(log(dens) ~ BCE, data = prior_data)
predict_base <- exp(predict(model_base, newdata = new_data))
predict_base_df <- cbind(new_data, predict_base)

base_df <- data.frame(Pt = predict_base,
                      Ptp1 = c(tail(predict_base, -1), NA))

rate <- mean(log(base_df$Pt / base_df$Ptp1), na.rm = T)

# ------------------
# Preparing the data
# ------------------

## Time slots
time_slots <- seq(-1500, -6500, by = -100)

## Expert knowledge
expert_data <- read.csv(
  "data/demography_studies_paper_kruk_neustupny_net_reduction.csv")

## 14C data
data_exp <- expert_data[expert_data$Grossregion != "NE", ]

data_exp_cleaned <- data_exp %>%
  select(Zeit_start_BCE, Zeit_ende_BCE, Wert_global) %>%
  na.omit

data_exp_cleaned <- data_exp_cleaned %>%
  mutate(t_start = Zeit_start_BCE * -1,
         t_end = Zeit_ende_BCE * -1,
         value = Wert_global)

data_14c <- read.csv("data/spd_europe.csv", row.names = 1)

data_14c$calBCE <- 1950 - data_14c$calBP
data_14c$calBCE <- data_14c$calBCE + 0.5

data_14c <- data_14c[data_14c$calBCE %in% time_slots, ]

data_14c$diff <- c(diff(data_14c$PrDens), NA)

diff14c <- data_14c$diff

names(diff14c) <- data_14c$calBCE

likelihood_matrix <- t(data_exp_cleaned %>% apply(1, function(x)
  sapply(time_slots, function(z) {
    if (x[4] <= z & x[5] >= z){return(x[6])} else {return(NA)}
  }) %>% unlist))

colnames(likelihood_matrix) <- time_slots

# ------------------
# Likelihood
# ------------------

likelihood <- function(param){

  likelihood <- log(1)
  x <- param[[1]]

  year <- param[[2]]
  prev_slot <- param[[3]]
  years <- as.numeric(names(prev_slot)) - year

  likelihood_collector <- vector()

  #------ > 0 ---------
  if (x < 0) {return(-Inf)}

  #------ > expert knowledge ---------
  likelihood_collector <- c(likelihood_matrix[, as.character(year)] %>%
                              na.omit %>% sapply(function(z){
                                dnorm(x, z, 1, log = T)
                              })) %>% unlist

  #------ exp_growth ---------
  prog <- prev_slot * exp(rate * years)
  prob <- dnorm(x, prog, 1, log = T)
  likelihood_collector <- c(likelihood_collector, prob)

  #------ 14C ---------
  prob_14c <- dnorm(diff14c[as.character(year)] + x, prev_slot, 1, log = T)

  likelihood_collector <- c(likelihood_collector, prob_14c)

  if (length(likelihood_collector > 0)) {
    likelihood <- sum(likelihood_collector)
  }

  return(likelihood)
}

# -----------------
#       Prior
# -----------------

prior <- function(param){
  x <- param[[1]]
  return(dnorm(x, 0, 10e6))
}

# -----------------
# Proposal function
# -----------------

proposalfunction <- function(param){
  return(rnorm(1, param, 1))
}

# -----------------
#   MCMC function
# -----------------

run_metropolis_MCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations + 1, length(time_slots)))
  colnames(chain) <- time_slots
  chain[1, ] <- startvalue

  for (i in 1:iterations){
    this_run <- chain[i, ]

    for (curr_data in as.character(time_slots)){
      this_slot_idx <- which(names(this_run) == curr_data)
      if (this_slot_idx == 1) {
        chain[i + 1, this_slot_idx] <- chain[i, this_slot_idx]
        next()
      }
      this_slot_idx <- which(names(this_run) == curr_data)

      proposal <- proposalfunction(chain[i, this_slot_idx])

      prev_slot <- this_run[this_slot_idx - 1]

      probab <- exp(likelihood(list(proposal,
                                   as.numeric(curr_data),
                                   prev_slot)) +
                     prior(list(proposal)) -
                     likelihood(list(chain[i, this_slot_idx],
                                     as.numeric(curr_data),
                                     prev_slot)) -
                     prior(list(chain[i, this_slot_idx]))
      )

      if (runif(1) < probab){
        chain[i + 1, this_slot_idx] <- proposal
      }else{
        chain[i + 1, this_slot_idx] <- chain[i, this_slot_idx]
      }
    }
  }
  return(mcmc(chain))
}

# ---------------
# Set up the run
# ---------------

iterations <- 10000
startvalue <- c(2.8, rep(1, length(time_slots) - 1))
curr_data <- as.character(time_slots)[2]

# ---------------
# Run the MCMC
# ---------------

chain <- run_metropolis_MCMC(startvalue, iterations)

# ----------------
# Some diagnostics
# ----------------

acceptanceRate <- 1 - rejectionRate(chain[, -1])
acceptanceRate

effectiveSize(chain[, -1])

raftery.diag(chain[, -1])

# ----------------
# burn in and thin
# ----------------

mcmc.trace.burned <- burnAndThin(chain[, -1], burn = nrow(chain) / 10)

mcmc.trace.burned.thinned <- burnAndThin(mcmc.trace.burned, thin = 5)

# -----------------------
# Summarising the results
# -----------------------

mcmc.trace.burned_mean <- mcmc.trace.burned.thinned %>%
  data.frame() %>%
  summarise_all(mean)

mcmc.trace.burned_q1 <- mcmc.trace.burned.thinned %>%
  data.frame() %>%
  summarise_all(list(Q1 = quantile), probs = 0.025)

mcmc.trace.burned_q2 <- mcmc.trace.burned.thinned %>%
  data.frame() %>%
  summarise_all(list(Q1 = quantile), probs = 0.975)


mcmc.trace.burned_df <- tibble(date = as.numeric(time_slots[-1]),
                               mean = t(mcmc.trace.burned_mean)[, 1],
                               q1 = t(mcmc.trace.burned_q1)[, 1],
                               q2 = t(mcmc.trace.burned_q2)[, 1])

# -----------------------
# Visualising the results
# -----------------------

ggplot(mcmc.trace.burned_df) +
  geom_line(aes(x = date, y = mean), color = "red") +
  geom_line(aes(x = date, y = q1), color = "grey") +
  geom_line(aes(x = date, y = q2), color = "grey")

last_plot() +
  geom_segment(data = data_exp_cleaned,
               aes(x = t_start,
                   xend = t_end,
                   y = value,
                   yend = value),
               col = "blue") +
  coord_cartesian(xlim = c(-6500, -1500),
                  ylim = c(0, 10)) +
  geom_point(data = data_exp_cleaned,
             aes(x = t_start,
                 y = value),
             col = "blue")

johannes_spline <- read.csv(
  "data/johannes_spline_for_engauge.csv", dec = ",")

spline.d <- as.data.frame(spline(johannes_spline$x,
                                 johannes_spline$Curve1))

last_plot() + geom_line(data = spline.d,
                        aes(x = x, y = y))

last_plot() +
  geom_ribbon(data = data_14c,
              aes(x = 1950 - calBP, ymin = 0, ymax = PrDens),
              alpha = .25) +
  geom_line(data = data_14c,
            aes(x = 1950 - calBP, y = PrDens), alpha = .25)
