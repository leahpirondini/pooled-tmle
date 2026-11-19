# ##################################################
# Script: Simulate longitudinal data
# Date: 01/06/2026
# Author: Leah Pirondini
# Notes: 
# ##################################################

library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)

expit <- function(x) 1 / (1 + exp(-x))

#' simulate longitudinal data with irregularly observed time-varying covariate
#'
#' @param n             sample size
#' @param n.obs         number of time points
#' @param N_depends_on  which version of the covariate drives the visit process N: options are "Lstar" or "L"
#' @param Lmod.*        parameters of the model for the true covariate Lstar
#' @param Amod.*        parameters of the treatment model
#' @param Nmod.*        parameters of the visit-process model
#' @param hazmod.*      parameters of the event hazard model
#' @param n.lags        number of lags to create for A, L and N
#' @param verbose       if TRUE, print simple checks of N, A and L
#'
#' @return              a long-format data frame, one row per id-time until the first event

sim_long_dat <- function(n,
                         n.obs,
                         N_depends_on = c("Lstar", "L"),
                         Lmod.beta.0, Lmod.beta.L, Lmod.beta.A, Lmod.beta.U, Lmod.sd,
                         Amod.beta.0, Amod.beta.L, Amod.beta.Lsq, Amod.beta.A, Amod.beta.N,
                         Nmod.beta.0, Nmod.beta.L, Nmod.beta.A, Nmod.beta.N, Nmod.beta.U,
                         hazmod.beta.0, hazmod.beta.A, hazmod.beta.L, hazmod.beta.U,
                         n.lags,
                         verbose = FALSE) {
  
  N_depends_on <- match.arg(N_depends_on)
  stopifnot(n.obs >= n.lags)
  
  # -------------------
  # data generation
  # -------------------
  
  N     <- matrix(nrow = n, ncol = n.obs)
  Lstar <- matrix(nrow = n, ncol = n.obs)
  L     <- matrix(nrow = n, ncol = n.obs)
  A     <- matrix(nrow = n, ncol = n.obs)
  Y     <- matrix(nrow = n, ncol = n.obs)
  
  U <- rnorm(n, 0, 0.1)
  
  N[, 1]     <- 1
  Lstar[, 1] <- rnorm(n, 0 + U, Lmod.sd)
  L[, 1]     <- Lstar[, 1]
  A[, 1]     <- rbinom(n, 1, expit(Amod.beta.0 + Amod.beta.L * L[, 1]))
  
  for (k in 2:n.obs) {
    
    # covariate that the visit process depends on
    L_for_N <- if (N_depends_on == "Lstar") Lstar[, k - 1] else L[, k - 1]
    
    N[, k] <- rbinom(n, 1, expit(Nmod.beta.0 + 
                                   Nmod.beta.L * L_for_N +
                                   Nmod.beta.A * A[, k - 1] +
                                   Nmod.beta.N * N[, k - 1] +
                                   Nmod.beta.U * U))
    Lstar[, k] <- rnorm(n, Lmod.beta.L * Lstar[, k - 1] +
                            Lmod.beta.A * A[, k - 1] +
                            Lmod.beta.0 * (k - 1) +
                            Lmod.beta.U * U,
                            Lmod.sd)
    
    L[, k] <- ifelse(N[, k] == 0, L[, k - 1], Lstar[, k])
    
    A[, k] <- rbinom(n, 1, expit(Amod.beta.0 + 
                                   Amod.beta.L * L[, k] +
                                   Amod.beta.Lsq*(L[, k])^2 +
                                   Amod.beta.A * A[, k - 1] +
                                   Amod.beta.N * N[, k]))
  }
  
  # event indicators (hazard depends on current and previous two time points)
  Y[, 1] <- rbinom(n, 1, expit(hazmod.beta.0 + 
                                 hazmod.beta.A * A[, 1] +
                                 hazmod.beta.L * Lstar[, 1] + 
                                 hazmod.beta.U * U))
  Y[, 2] <- rbinom(n, 1, expit(hazmod.beta.0 +
                                 hazmod.beta.A * (A[, 1] + A[, 2]) +
                                 hazmod.beta.L * (Lstar[, 1] + Lstar[, 2]) +
                                 hazmod.beta.U * U))
  for (k in 3:n.obs) {
    Y[, k] <- rbinom(n, 1, expit(hazmod.beta.0 +
                                   hazmod.beta.A * (A[, k - 2] + A[, k - 1] + A[, k]) +
                                   hazmod.beta.L * (Lstar[, k - 2] + Lstar[, k - 1] + Lstar[, k]) +
                                   hazmod.beta.U * U))
  }
  
  # once an event has occurred, Y stays 1
  Y <- t(apply(Y, 1, cummax))
  
  # ----------------------
  # convert to long format
  # ----------------------
  
  colnames(A)     <- paste0("A",     0:(n.obs - 1))
  colnames(Lstar) <- paste0("Lstar", 0:(n.obs - 1))
  colnames(L)     <- paste0("L",     0:(n.obs - 1))
  colnames(N)     <- paste0("N",     0:(n.obs - 1))
  colnames(Y)     <- paste0("Y",     0:(n.obs - 1))
  
  dat <- data.frame(id = 1:n, A, L, Lstar, N, U, Y)
  
  dat.long <- dat %>%
    pivot_longer(cols = -c(id, U),
                 names_to = c("Variable", "time"),
                 names_pattern = "([A-Za-z]+)([0-9]+)",
                 names_transform = list(time = as.numeric)) %>%
    pivot_wider(id_cols = c(id, time, U),
                names_from = Variable,
                values_from = value) %>%
    arrange(id, time)
  
  # drop rows after the event has occurred (keep the first event row)
  dat.long <- dat.long %>%
    group_by(id) %>%
    mutate(Y_sum = cumsum(Y)) %>%
    ungroup() %>%
    filter(Y_sum < 2) %>%
    dplyr::select(-Y_sum) %>%
    rename(event = Y)
  
  # "observed" Lstar
  dat.long <- dat.long %>%
    mutate(Lstar.obs = ifelse(N == 0, NA, Lstar))
  
  # -----------------
  # add lag variables
  # -----------------
  
  lag_vars <- c("A", "L", "N")
  
  dat.long <- dat.long %>% group_by(id)
  for (var in lag_vars) {
    for (i in seq_len(n.lags)) {
      dat.long <- dat.long %>%
        mutate(!!paste0(var, "lag", i) := lag(.data[[var]], n = i))
    }
  }
  dat.long <- dat.long %>% ungroup()
  
  # missing lags (before time 0) set to 0
  lag_cols <- paste0(rep(lag_vars, each = n.lags), "lag", seq_len(n.lags))
  dat.long <- dat.long %>%
    mutate(across(all_of(lag_cols), ~ ifelse(is.na(.x), 0, .x)))
  
  # ----------
  # checks
  # ----------
  if (verbose) {
    print(table(dat.long$N, useNA = "always"))
    print(table(dat.long$A, useNA = "always"))
    print(summary(dat.long$L))
  }
  
  dat.long
}

# Example:
dat.long <- sim_long_dat(
   n = 5000, 
   n.obs = 20, 
   N_depends_on = "L",
   Lmod.beta.0 = 0.05,
   Lmod.beta.L = 0.5, 
   Lmod.beta.A = -0.5,
   Lmod.beta.U = 1,
   Lmod.sd = 1,
   Amod.beta.0 = -0.5,
   Amod.beta.L = 0.6, 
   Amod.beta.Lsq = 0.05,
   Amod.beta.A = 0.7, 
   Amod.beta.N = 3, 
   Nmod.beta.0 = -1,
   Nmod.beta.L = 0.3,
   Nmod.beta.A = 0.5,
   Nmod.beta.N = 0.5,
   Nmod.beta.U = 0,
   hazmod.beta.0 = -3,
   hazmod.beta.L = 0.3, 
   hazmod.beta.A = -0.3,
   hazmod.beta.U = 0.05,
   n.lags = 3,
   verbose = TRUE
 )
 

 
