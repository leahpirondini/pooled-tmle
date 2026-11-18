
library(tidyverse)
library(magrittr)

# -------------------
# data generation
# -------------------

# sample size
n <- 1000

# first generate the time-dependent A and L
L <- matrix(nrow=n,ncol=n.obs)
A <- matrix(nrow=n,ncol=n.obs)
Y <- matrix(nrow=n,ncol=n.obs)


U <- rnorm(n,0,0.1)

L[,1] <- rnorm(n,0+U,Lmod.sd)
A[,1] <- rbinom(n,1,plogis(Amod.beta.0 + Amod.beta.L*L[,1] + Amod.beta.Lsq*(L[,1])^2))

L[,2] <- rnorm(n, Lmod.beta.L*L[,1] + Lmod.beta.Lsq*(L[,1])^2 + Lmod.beta.A*A[,1] + Lmod.beta.0*(1) + Lmod.beta.U*U, Lmod.sd)
A[,2] <- rbinom(n, 1, plogis(Amod.beta.0 + Amod.beta.L*L[,2] + Amod.beta.Lsq*(L[,2])^2 + Amod.beta.A*A[,1] + Amod.beta.AL*A[,1]*L[,2]))

# generate event indicators Y
Y[,1] <- rbinom(n, 1, plogis(hazmod.beta.0 + hazmod.beta.U*U + hazmod.beta.A*A[,1] + hazmod.beta.L*L[,1] + hazmod.beta.AL*A[,1]*L[,1]))

Y[,2] <- rbinom(n, 1, plogis(hazmod.beta.0 + hazmod.beta.U*U + 
                               hazmod.beta.A*A[,1] + hazmod.beta.L*L[,1] + hazmod.beta.AL*A[,1]*L[,1] +
                               hazmod.beta.A*A[,2] + hazmod.beta.L*L[,2] + hazmod.beta.AL*A[,2]*L[,2]
                              ))


# ----------------------
# combine
# ----------------------

colnames(A) <- paste0("A",0:{n.obs-1})
colnames(L) <- paste0("L",0:{n.obs-1})
colnames(Y) <- paste0("Y",0:{n.obs-1})

dat <- data.frame(id=1:n,A,L,U,Y)

dat %$% mean(Y0)
dat %$% mean(Y1)

# ----------------------
# convert to long format
# ----------------------

for (i in 1:(n.obs-1)) {
  prev_var <- if (i == 1) "Y0" else paste0("Y", i - 1)
  curr_var <- paste0("Y", i)
  
  dat <- dat %>%
    mutate(!!curr_var := ifelse(.data[[prev_var]] == 1, 1, .data[[curr_var]]))
}


dat.long <- dat %>%
  pivot_longer(cols=c(paste0("A",0:{n.obs-1}),
                      paste0("L",0:{n.obs-1}),
                      paste0("Y",0:{n.obs-1})),
               names_to="Variable"
  ) %>%
  mutate(time=as.numeric(str_extract(Variable, '[0-9]+'))) %>% 
  separate(Variable, into='Variable', sep='[0-9]+') %>%
  pivot_wider(id_cols=c(id,time,U),
              names_from=Variable,
              values_from=value) %>%
  arrange(id, time)

dat.long <- dat.long %>%
  group_by(id) %>%
  mutate(Y_sum = cumsum(Y)) %>%
  ungroup()
dat.long <- dat.long %>%
  filter(Y_sum < 2) %>% # get rid of rows after event has occurred
  select(-Y_sum) %>%
  rename(event=Y)

# -----------------
# add lag variables
# -----------------

lag_vars <- c("A", "L")
n_lags <- 1

dat.long <- dat.long %>%
  group_by(id)

for (var in lag_vars) {
  for (i in 1:n_lags) {
    dat.long <- dat.long %>%
      mutate(!!paste0(var, "lag", i) := lag(.data[[var]], n = i))
  }
}

dat.long <- dat.long %>%
  ungroup()

for (var in lag_vars) {
  for (i in 1:n_lags) {
    var_name <- paste0(var, "lag", i)
    dat.long <- dat.long %>%
      mutate(!!var_name := ifelse(is.na(.data[[var_name]]), 0, .data[[var_name]]))
  }
}


# -----------------------------
# format the wide data for TMLE
# -----------------------------

new_order <- as.vector(
  rbind(
    paste0("L", 0:{n.obs-1}),
    paste0("A", 0:{n.obs-1}),
    paste0("Y", 0:{n.obs-1})
  )
)
dat.wide.tmle <- dat[, new_order]

for(i in 1:nrow(dat.wide.tmle)) {
  
  for(k in 0:{n.obs-1}) {
    
    if(dat.wide.tmle[i, paste0("Y", k)] == 1) {
      
      if(k < {n.obs-1}) {
        
        for(j in (k+1):{n.obs-1}) {
          
          dat.wide.tmle[i, paste0("L", j)] <- NA
          dat.wide.tmle[i, paste0("A", j)] <- NA
          
        }
      }
      
      break
    }
  }
}

head(dat.wide.tmle)

table(dat.wide.tmle$Y0)
table(dat.wide.tmle$Y1)

# ------------------------
# remove things not needed
# ------------------------

keep_objs <- c("dat", "dat.long", "dat.wide.tmle", "n", "n.obs")
rm(list = setdiff(ls(), keep_objs))

# ----------
# checks
# ----------

table(dat.long$A,useNA="always")
summary(dat.long$L)



