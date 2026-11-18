# -----------------------------------
# -----------------------------------
# ltmleMSM 
# implementation: manual to check
# -----------------------------------
# -----------------------------------

# ----------------------------------------------------------
# manual implementation steps for "pooled TMLE"
# steps correspond to algorithm in section 4.2.2 of
# Peterson et al. 2014
# ----------------------------------------------------------

dat.tmle <- dat.wide.tmle

# ----------------------------------
# ----------------------------------
# Step 1: estimate propensity scores
# ----------------------------------
# ----------------------------------

modA0 <- glm(A0 ~ L0 + I(L0^2), family = "binomial", data = dat.tmle)
modA1 <- glm(A1 ~ L1 + I(L1^2) + A0 + I(L1*A0), family = "binomial", data = dat.tmle)

dat.tmle <- dat.tmle %>%
  mutate(
    predA0 = predict(modA0, type = "response", newdata = dat.tmle),
    predA1 = predict(modA1,
                     type = "response",
                     newdata = dat.tmle),
    g0 = ifelse(A0==1,predA0,1-predA0),
    g1 = ifelse(
      Y0 == 1,
      g0,
      g0 * ifelse(A1==1,predA1,1-predA1)
    )
  )

# ------------------------------------
# ------------------------------------
# Step 2: targeting Y1, time-point k=1
# ------------------------------------
# ------------------------------------

# ---------------
# 2a. initial fit
# ---------------

modY1 <- glm(Y1 ~ L0 + L1 + A0 + A1 + I(L0*A0) + I(L1*A1), family = "binomial", data = dat.tmle)

predY1.always <- predict(modY1, type = "response", newdata = dat.tmle %>% mutate(A0=1,A1=1))
predY1.never <- predict(modY1, type = "response", newdata = dat.tmle %>% mutate(A0=0,A1=0))

# -----------------
# 2b. stack results
# -----------------

always <- data.frame(
  i = 1:n,
  regime = "always",
  Y1 = dat.tmle$Y1,
  Q_Y1_k1 = predY1.always
)

never <- data.frame(
  i = 1:n,
  regime = "never",
  Y1 = dat.tmle$Y1,
  Q_Y1_k1 = predY1.never
)

# ---------------------------------------
# 2c. multidimensional clever covariate
# ---------------------------------------

# for each subject i we create a matrix of covariates
# with dimensions d x no. parameters in MSM = 2 x 4
# i.e. two because our two regimes are always and never
# and 4 because the MSM has 4 parameters

# we need to populate it with h(=1) * derivative wrt beta * I(Abar=dbar) / g1
# msm <- beta0 + beta1*time + beta2*cumA + beta3*time*cumA
# => derivative wrt each beta is: (1,time,cumA,time*cumA)

I_never <- (dat.tmle$A0==0) * (dat.tmle$A1==0)
I_always <- (dat.tmle$A0==1) * (dat.tmle$A1==1)

h <- 1

covariate.matrix.always <-  (h * I_always) / dat.tmle$g1
covariate.matrix.never <- (h * I_never) / dat.tmle$g1

wts <- c(covariate.matrix.always,covariate.matrix.never)

pooled_data_Q_Y1 <- rbind(always,never) %>%
  mutate(intercept = 1,
         time = 2,
         cumA = c(rep(2,n), rep(0,n)))

# deterministic rows where P=1 because event already occurred
pooled_data_Q_Y1 <- pooled_data_Q_Y1 %>%
  mutate(Q_Y1_k1 = ifelse(is.na(Q_Y1_k1),1,Q_Y1_k1))

# ---------------------------------------
# 2d. fit the update model
# ---------------------------------------

Y1.update.mod <- glm(Y1 ~ intercept + time + cumA + time*cumA - 1, 
                      family = "binomial", 
                      offset = qlogis(Q_Y1_k1),
                      weights = wts,
                      data = pooled_data_Q_Y1)
summary(Y1.update.mod)

# ---------------------------------------
# 2e. obtain predictions to get Q*_Y2
# ---------------------------------------

pooled_data_Q_Y1 <- pooled_data_Q_Y1 %>%
  mutate(Q_Y1_k1_update = predict(Y1.update.mod, type = "response", newdata = pooled_data_Q_Y1))

# ------------------------------------
# ------------------------------------
# Step 3: targeting Y1, time-point k=0
# ------------------------------------
# ------------------------------------

# ---------------
# 3a. initial fit
# ---------------

pooled_data_Q_Y1 <- pooled_data_Q_Y1 %>%
  mutate(A0 = rep(dat.tmle$A0,2),
         L0 = rep(dat.tmle$L0,2))

init_fit <- glm(Q_Y1_k1_update ~ L0 + A0 + I(L0*A0), family = "binomial", data = pooled_data_Q_Y1)

newdat_always <- pooled_data_Q_Y1 %>% 
  filter(regime=="always") %>% 
  mutate(A0=1)

newdat_never <- pooled_data_Q_Y1 %>% 
  filter(regime=="never") %>% 
  mutate(A0=0)

Q_Y1_k0_always <- predict(init_fit, type = "response", newdata = newdat_always)
Q_Y1_k0_never <- predict(init_fit, type = "response", newdata = newdat_never)

pooled_data_Q_Y1 <- pooled_data_Q_Y1 %>%
  mutate(Q_Y1_k0 = c(Q_Y1_k0_always,Q_Y1_k0_never)
  )

# ---------------------------------------
# 3b. multidimensional clever covariate
# ---------------------------------------

# NB:
# covariates of MSM relate to time = 2
# but the I() and g relate to time = 1

I_never <- (dat.tmle$A0==0)
I_always <- (dat.tmle$A0==1)

h <- 1

covariate.matrix.always <- (h * I_always) / dat.tmle$g0
covariate.matrix.never <- (h * I_never) / dat.tmle$g0

wts <- c(covariate.matrix.always,covariate.matrix.never)

pooled_data_Q_Y1 <- pooled_data_Q_Y1 %>%
  select(i,regime,Q_Y1_k0,Q_Y1_k1_update,intercept,time,cumA)


# ---------------------------------------
# 3c. fit the update model
# ---------------------------------------

update.mod <- glm(Q_Y1_k1_update ~ intercept + time + cumA + time*cumA - 1, 
                     family = "binomial", 
                     offset = qlogis(Q_Y1_k0),
                     weights = wts,
                     data = pooled_data_Q_Y1)

summary(update.mod)

# ---------------------------------------
# 3d. obtain predictions to get Q*_Y2
# ---------------------------------------

pooled_data_Q_Y1 <- pooled_data_Q_Y1 %>%
  mutate(Q_Y1_k0_update = ifelse(
                  is.na(Q_Y1_k1_update),NA,
                  predict(update.mod, type = "response", newdata = pooled_data_Q_Y1)))

pooled_results <- pooled_data_Q_Y1 %>%
  select(i,regime,Q_Y1_k0_update,intercept,time,cumA) %>%
  rename(Qstar = Q_Y1_k0_update)

# ------------------------------------
# ------------------------------------
# Step 4: targeting Y0, time-point k=0
# ------------------------------------
# ------------------------------------

# ---------------
# 4a. initial fit
# ---------------

modY0 <- glm(Y0 ~ L0 + A0 + I(L0*A0), family = "binomial", data = dat.tmle)

predY0.always <- predict(modY0, type = "response", newdata = dat.tmle %>% mutate(A0=1))
predY0.never <- predict(modY0, type = "response", newdata = dat.tmle %>% mutate(A0=0))

# -----------------
# 4b. stack results
# -----------------

always <- data.frame(
  i = 1:n,
  regime = "always",
  Y0 = dat.tmle$Y0,
  Q_Y0_k0 = predY0.always
)

never <- data.frame(
  i = 1:n,
  regime = "never",
  Y0 = dat.tmle$Y0,
  Q_Y0_k0 = predY0.never
)

# ---------------------------------------
# 4c. multidimensional clever covariate
# ---------------------------------------

h <- 1

I_never <- (dat.tmle$A0==0) 
I_always <- (dat.tmle$A0==1)

covariate.matrix.always <- (h * I_always) / dat.tmle$g0
covariate.matrix.never <- (h * I_never) / dat.tmle$g0

wts <- c(covariate.matrix.always,covariate.matrix.never)

pooled_data_Q_Y0 <- rbind(always,never) %>%
  mutate(intercept = 1,
         time = 1,
         cumA = c(rep(1,n), rep(0,n)))


# ---------------------------------------
# 4d. fit the update model
# ---------------------------------------

Y0.update.mod <- glm(Y0 ~ intercept + time + cumA + time*cumA - 1, 
                     family = "binomial", 
                     offset = qlogis(Q_Y0_k0),
                     weights = wts,
                     data = pooled_data_Q_Y0)

summary(Y0.update.mod)

# ---------------------------------------
# 4e. obtain predictions to get Q*_Y2
# ---------------------------------------

pooled_data_Q_Y0 <- pooled_data_Q_Y0 %>%
  mutate(Q_Y0_k0_update = predict(Y0.update.mod, type = "response", newdata = pooled_data_Q_Y0))

pooled_data_Q_Y0 <- pooled_data_Q_Y0 %>%
  select(i,regime,Q_Y0_k0_update,intercept,time,cumA) %>%
  rename(Qstar = Q_Y0_k0_update)

# ------------------------------------
# ------------------------------------
# Step 5: fitting the MSM
# ------------------------------------
# ------------------------------------

# stack Q* results to get a single vector of length n x regimes x time-points = n x 2 x 2
# add columns time (= 1,2) and cumA (=0,1,2)
# fit MSM with weights h (=1)

pooled_results <- rbind(pooled_results,pooled_data_Q_Y0)

msm <- glm(Qstar ~ time*cumA, family = binomial(), data = pooled_results)
summary(msm)

# comparison of results
summary(msm)
table(tmle.msm$msm[1])

table(predict(msm, type = "response"))
table(tmle.msm$msm[3])
