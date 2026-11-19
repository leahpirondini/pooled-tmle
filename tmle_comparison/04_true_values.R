# ##################################################
# Script: Simulate true survival under a static treatment strategy
# Notes: companion to sim_long_dat():
#        same data-generating process,
#        but A is set by intervention
# ##################################################

library(dplyr)

expit <- function(x) 1 / (1 + exp(-x))

#' simulate true survival probabilities under a static treatment strategy
#'
#' @param n.true        sample size
#' @param n.obs         number of time points
#' @param A.value       treatment value(s) to set: a single value (1 = always
#'                      treat, 0 = never treat) or a vector of length n.obs
#'                      giving a fixed treatment at each time point
#' @param N_depends_on  which version of the covariate drives the visit process N: options are "Lstar" or "L"
#' @param Lmod.*        parameters of the model for the true covariate Lstar
#' @param Nmod.*        parameters of the visit-process model
#' @param hazmod.*      parameters of the event hazard model
#'
#' @details Uses the same Lmod, Nmod and hazmod parameters as sim_long_dat().
#'          The treatment model (Amod.*) is not needed since A is set.
#'
#' @return              a one-row data frame of survival probabilities at times 1:n.obs

sim_true_surv <- function(n.true,
                          n.obs,
                          A.value,
                          N_depends_on = c("Lstar", "L"),
                          Lmod.beta.0, Lmod.beta.L, Lmod.beta.A, Lmod.beta.U, Lmod.sd,
                          Nmod.beta.0, Nmod.beta.L, Nmod.beta.A, Nmod.beta.N, Nmod.beta.U,
                          hazmod.beta.0, hazmod.beta.A, hazmod.beta.L, hazmod.beta.U) {
  
  N_depends_on <- match.arg(N_depends_on)
  stopifnot(n.obs >= 3, length(A.value) %in% c(1, n.obs))
  
  # -------------------
  # data generation
  # -------------------
  
  # treatment is set by intervention (scalar is recycled across time points)
  A <- matrix(A.value, nrow = n, ncol = n.obs, byrow = TRUE)
  
  N     <- matrix(nrow = n, ncol = n.obs)
  Lstar <- matrix(nrow = n, ncol = n.obs)
  L     <- matrix(nrow = n, ncol = n.obs)
  Y     <- matrix(nrow = n, ncol = n.obs)
  
  U <- rnorm(n, 0, 0.1)
  
  N[, 1]     <- 1
  Lstar[, 1] <- rnorm(n, 0 + U, Lmod.sd)
  L[, 1]     <- Lstar[, 1]
  
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
  
  # -----------------------
  # outcomes
  # -----------------------
  
  surv.prob <- 1 - colMeans(Y)
  
  out <- as.data.frame(t(surv.prob))
  colnames(out) <- 1:n.obs
  out
}


pars <- list(
   n.obs = 20, 
   N_depends_on = "L",
   Lmod.beta.0 = 0.05, 
   Lmod.beta.L = 0.5, 
   Lmod.beta.A = -0.5, 
   Lmod.beta.U = 1, 
   Lmod.sd = 1,
   Nmod.beta.0 = -1, 
   Nmod.beta.L = 0.3, 
   Nmod.beta.A = 0.5, 
   Nmod.beta.N = 0.5, 
   Nmod.beta.U = 0,
   hazmod.beta.0 = -3, 
   hazmod.beta.L = 0.3, 
   hazmod.beta.A = -0.3, 
   hazmod.beta.U = 0.05
 )


summary.true.always <- do.call(sim_true_surv, c(list(n.true = 1e6, A.value = 1), pars))
summary.true.never  <- do.call(sim_true_surv, c(list(n.true = 1e6, A.value = 0), pars))

