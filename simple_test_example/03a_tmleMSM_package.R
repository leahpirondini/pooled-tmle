# -----------------------------------
# -----------------------------------
# Pooled TMLE
# implementation: ltmle package
# -----------------------------------
# -----------------------------------

library(ltmle)

# ---------------------------
# data
# ---------------------------

dat.ltmle <- dat.wide.tmle

n <- dim(dat.ltmle)[1]

# ------------------
# nodes
# ------------------

Anodes <- c("A0","A1")
Lnodes <- c("L0","L1")
Ynodes <- c("Y0","Y1")

# ------------------
# Qform / gform
# ------------------

Qform = c(L0 = 'Q.kplus1 ~ 1',  
          Y0 = 'Q.kplus1 ~ L0 + A0 + L0*A0',
          L1 = 'Q.kplus1 ~ 1', 
          Y1 = 'Q.kplus1 ~ L0 + L1 + A0 + A1 + L0*A0 + L1*A1')

gform <- c(A0 = 'A0 ~ L0 + L0^2', 
           A1 = 'A1 ~ A0 + L1 + L1^2 + A0*L1')

# ------------------
# regimes
# ------------------
# array: n x numAnodes x numRegimes = n x time.points x ?
regimes <- array(NA, dim = c(n, n.obs, 2)) 
regimes[,,1] <- c(1,1)
regimes[,,2] <- c(0,0)
head(regimes[,,])

# ------------------
# summary measures
# ------------------
# array: num.regimes x num.summary.measures x num.final.Ynodes
summary.measures <- array(NA, dim = c(2,2,2))
# first summary measure = time
summary.measures[1,1,] <- c(1,2)
summary.measures[2,1,] <- c(1,2)
# second summary measure = cumA
summary.measures[1,2,] <- c(1,2)
summary.measures[2,2,] <- c(0,0)

colnames(summary.measures) <- c("time","cumA")
summary.measures

# ------------------
# msm
# ------------------
working.msm <- 'Y ~ time*cumA'
final.Ynodes <- Ynodes

# ------------------
# run ltmleMSM
# ------------------

tmle.msm <- ltmleMSM(data = dat.ltmle, 
                     Anodes = Anodes,
                     Cnodes = NULL,
                     Lnodes = Lnodes, 
                     Ynodes = Ynodes,
                     survivalOutcome = TRUE, 
                     Qform = Qform,
                     gform = gform,
                     gbounds = c(0,1),
                     Yrange = NULL,
                     deterministic.g.function = NULL,
                     SL.library = "glm",
                     SL.cvControl = list(),
                     regimes = regimes,
                     working.msm = working.msm,
                     summary.measures = summary.measures,
                     final.Ynodes = final.Ynodes,
                     stratify = FALSE, # TRUE = stratify on abar when estimating Q and g. FALSE = pool over abar
                     msm.weights = NULL, # options are NULL or "empirical" or user supplied. 
                     estimate.time = TRUE,
                     gcomp = FALSE,
                     iptw.only = FALSE,
                     deterministic.Q.function = NULL,
                     variance.method = "tmle",
                     observation.weights = NULL,
                     id = NULL
)

unlist(tmle.msm$msm[1])

# note: msm.weights
# if "empirical", weight by empirical proportions of rows matching each regime for each final.Ynode, with duplicate regimes given zero weight 
# if NULL, no weights
# this is the h(d,t) component of the multidimensional clever covariate

