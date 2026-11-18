# ------------------
# ------------------
# parameter values
# ------------------
# ------------------

# time points
n.obs <- 2

# ---------------
# model for L
# ---------------
Lmod.beta.0 <- 0.1
Lmod.beta.L <- 0.7
Lmod.beta.Lsq <- 0.1
Lmod.beta.A <- -0.5
Lmod.beta.U <- 1
Lmod.sd <- 1

# ---------------
# model for A
# ---------------
Amod.beta.0 <- -0.5
Amod.beta.L <- 0.6
Amod.beta.Lsq <- 0.05
Amod.beta.A <- 0.7 
Amod.beta.AL <- 0.05
Amod.beta.N <- 3

# ---------------
# model for Y
# ---------------
hazmod.beta.0 <- -0.5
hazmod.beta.L <- 0.3
hazmod.beta.A <- -0.7 
hazmod.beta.AL <- -0.1
hazmod.beta.U <- 0.05
