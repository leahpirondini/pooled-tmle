# --------------------------------------------------------------------
# --------------------------------------------------------------------
# Standard longitudinal TMLE
# --------------------------------------------------------------------
# --------------------------------------------------------------------

# notes:
# - directly targets mean counterfactual outcome separately at each time-point
# - in this example the goal is to estimate the effect of static treatment strategies: always treat vs never treat
# - uses lmtp package (rather than ltmle) 
# - as the lmtp_survival wrapper allows for estimation of mean counterfactual outcome at every time-point (rather than final time-point only)

library(lmtp)

# ------------------
# general set up
# ------------------

n.obs <- 20

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


# ---------------
# setup for lmtp
# ---------------

# SuperLearner libraries
SL.lib.outcome <- "SL.glm.interaction"
SL.lib.treat <- "SL.glm"

# variables
# "A" nodes (i.e. intervention variables)
Anodes <- paste0("A",0:{n.obs-1})


# "L" nodes (i.e. confounders)
for (i in 2:n.obs){
  assign(paste0("Lnodes.",i-1),
         c(paste0("N",i-2),
           paste0("L",i-1)
         ))
}

Lnodes.baseline <- "L0"
Lnodes.nonbaseline <- mget(paste0("Lnodes.", 1:{n.obs-1}))
Lnodes <- c(list(Lnodes.baseline), Lnodes.nonbaseline)
names(Lnodes)[1] <- "Lnodes.0"

# Y nodes
Ynodes <- paste0("Y",0:{n.obs-1})

# -----------------------------
# -----------------------------
# treatment strategies
# -----------------------------
# -----------------------------

# always and never can use default options "static_binary_on" and "static_binary_off"

# ---------------
# run lmtp
# ---------------

tmle_always <-
  lmtp_survival(
    data = dat.tmle,
    trt = Anodes,
    outcome = Ynodes,
    baseline = NULL,
    time_vary = Lnodes,
    cens = NULL,
    compete = NULL,
    shift = static_binary_on,
    shifted = NULL,
    estimator = "lmtp_tmle",
    k = 3,
    mtp = FALSE,
    id = "id",
    learners_outcome = SL.lib.outcome,
    learners_trt = SL.lib.treat,
    folds = 1,
    weights = NULL,
    control = lmtp_control()
  )


tmle_never <-
  lmtp_survival(
    data = dat.tmle,
    trt = Anodes,
    outcome = Ynodes,
    baseline = NULL,
    time_vary = Lnodes,
    cens = NULL,
    compete = NULL,
    shift = static_binary_off,
    shifted = NULL,
    estimator = "lmtp_tmle",
    k = 3,
    mtp = FALSE,
    id = "id",
    learners_outcome = SL.lib.outcome,
    learners_trt = SL.lib.treat,
    folds = 1,
    weights = NULL,
    control = lmtp_control()
  )


# ----------------------------
# ----------------------------
# extract results
# ----------------------------
# ----------------------------

tmle.lmtp.always <- tidy(tmle_always)
tmle.lmtp.never <- tidy(tmle_never)





