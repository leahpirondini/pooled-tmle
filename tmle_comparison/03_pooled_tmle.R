# --------------------------------------------------------------------
# --------------------------------------------------------------------
# Pooled TMLE to estimate parameters of MSM
# --------------------------------------------------------------------
# --------------------------------------------------------------------

# notes:
# - targets parameters of a chosen marginal structural working model
# - then uses the fitted MSM to estimate mean counterfactual outcomes at each time-point
# - in this example the goal is to estimate the effect of static treatment strategies: always treat vs never treat
# - uses the function ltmleMSM within the ltmle package

library(ltmle)

# ------------------
# general set up
# ------------------

n.obs <- 20
n <- 5000

# select only the variables needed
dat.tmle <- dat.long %>%
  select(id,time,L,A,N,event) %>%
  rename(Y=event)

# convert to wide format
dat.tmle <- dat.tmle %>%
  pivot_wider(
    id_cols = id,
    names_from = time,
    values_from = c(L, A, N, Y),
    names_sep = ""
  )

# set Y=1 after event has occurred
dat.tmle <- dat.tmle %>%
  mutate(across(all_of(paste0("Y", 0:(n.obs-1))), ~ifelse(is.na(.),1,.))
  )

# change name of N vars
dat.tmle <- dat.tmle %>% select(-N0)
names(dat.tmle)[names(dat.tmle) %in% paste0("N", 1:19)] <- paste0("N", 0:18)

# reorder
neworder_timevars <- unlist(lapply(0:(n.obs-1), function(t)
  paste0(if (t < n.obs-1) c("L","A","Y","N") else c("L","A","Y"), t)
))
neworder <- c("id", neworder_timevars)
dat.tmle <- dat.tmle[, neworder]

# ----------------
# set up for ltmle
# ----------------

# ------------------
# nodes
# ------------------

Lnodes <- unlist(lapply(0:(n.obs-1), function(t)
  paste0(if (t < n.obs-1) c("L","N") else c("L"), t)
))
Anodes <- paste0("A",0:{n.obs-1})
Ynodes <- paste0("Y",0:{n.obs-1})


SL.lib.g <- c("SL.glm")
SL.lib.Q <- c("SL.glm")


# ------------------
# Qform / gform
# ------------------

gform <- sapply(0:(n.obs - 1), function(i) {
  if (i == 0) {
    paste0("A0 ~ L0 + I(L0^2)")
  } else {
    paste0(
      "A", i,
      " ~ L", i,
      " + I(L", i, "^2) + A", i - 1,
      " + N", i - 1
    )
  }
})



Qform <- unlist(lapply(0:(n.obs - 1), function(i) {
  
  # L model
  out <- c(
    setNames("Q.kplus1 ~ 1", paste0("L", i))
  )
  
  # Y model
  y_formula <-
    if (i == 0) {
      "Q.kplus1 ~ A0 + L0"
      
    } else if (i == 1) {
      "Q.kplus1 ~ A0 + A1 + L0 + L1 + N0 + I(L1*N0)"
      
    } else {
      paste0(
        "Q.kplus1 ~ ",
        paste0("A", (i-2):i, collapse = " + "),
        " + ",
        paste0("L", (i-2):i, collapse = " + "),
        " + ",
        paste0("N", (i-2):(i-1), collapse = " + "),
        " + ",
        paste0(
          sprintf("I(L%d*N%d)", (i-1):i, (i-2):(i-1)),
          collapse = " + "
        )
      )
    }
  
  out <- c(
    out,
    setNames(y_formula, paste0("Y", i))
  )
  
  # N model (skip final time point)
  if (i < n.obs - 1) {
    out <- c(
      out,
      setNames("Q.kplus1 ~ 1", paste0("N", i))
    )
  }
  
  out
}))


# ------------------
# regimes
# ------------------
# array: n x numAnodes x numRegimes
regimes <- array(NA, dim = c(n, n.obs, 2)) 
regimes[,,1] <- rep(1,n.obs)
regimes[,,2] <- rep(0,n.obs)
head(regimes[,,])

# ------------------
# summary measures
# ------------------
# array: num.regimes x num.summary.measures x num.final.Ynodes
summary.measures <- array(NA, dim = c(2,2,n.obs))
# first summary measure = time
summary.measures[1,1,] <- c(1:n.obs)
summary.measures[2,1,] <- c(1:n.obs)
# second summary measure = cumA
summary.measures[1,2,] <- c(1:n.obs)
summary.measures[2,2,] <- rep(0,n.obs)

colnames(summary.measures) <- c("time","cumA")
summary.measures

# ------------------
# msm
# ------------------
#working.msm <- 'Y ~ time + cumA'
working.msm <- 'Y ~ log(time) + time + cumA + log(time):cumA'
final.Ynodes <- Ynodes

# ------------------
# run ltmleMSM
# ------------------

tmle.msm <- ltmleMSM(data = dat.tmle, 
                     Anodes = Anodes,
                     Cnodes = NULL,
                     Lnodes = Lnodes, 
                     Ynodes = Ynodes,
                     survivalOutcome = TRUE, 
                     Qform = Qform,
                     gform = gform,
                     gbounds = c(0.01,1),
                     Yrange = NULL,
                     deterministic.g.function = NULL,
                     SL.library = "glm",
                     SL.cvControl = list(),
                     regimes = regimes,
                     working.msm = working.msm,
                     summary.measures = summary.measures,
                     final.Ynodes = final.Ynodes,
                     stratify = FALSE, # TRUE = stratify on abar when estimating Q and g. FALSE = pool over abar
                     msm.weights = NULL, # options are NULL or "empirical" or user supplied
                     estimate.time = TRUE,
                     gcomp = FALSE,
                     iptw.only = FALSE,
                     deterministic.Q.function = NULL,
                     variance.method = "tmle",
                     observation.weights = NULL,
                     id = NULL
)

# ---------------------------------
# extract coefficients / covariance
# ---------------------------------

#rhs <- ~ time + cumA
rhs <- ~ log(time) + time + cumA + log(time):cumA

beta  <- coef(tmle.msm$msm)
covar <- tmle.msm$variance.estimate / n

# check: design matrix columns must line up with coefficients
nd_check <- data.frame(time = 1:n.obs, cumA = 0)
stopifnot(ncol(model.matrix(rhs, nd_check)) == length(beta),
          all(dim(covar) == length(beta)))

# ---------------------------
# predict survival under a regime
# ---------------------------

predict_surv <- function(cumA_fun) {
  nd <- data.frame(time = 1:n.obs)
  nd$cumA <- cumA_fun(nd$time)
  
  X   <- model.matrix(rhs, nd)
  eta <- as.vector(X %*% beta)
  se_eta <- sqrt(rowSums((X %*% covar) * X))
  
  p <- plogis(eta) 
  
  data.frame(
    time      = nd$time,
    estimate  = 1 - p,
    std.error = p * (1 - p) * se_eta, 
    conf.low  = 1 - plogis(eta + 1.96 * se_eta),
    conf.high = 1 - plogis(eta - 1.96 * se_eta)
  )
}

# ---------------------------
# store results
# ---------------------------

# never treat: cumA = 0 at all times
tmle.msm.never  <- predict_surv(function(t) rep(0, length(t)))

# always treat: cumA = time
tmle.msm.always <- predict_surv(function(t) t)


