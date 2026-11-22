# ===== Workspace preparation
rm(list = ls()) # clear all variables
set.seed(1) # fix random seed for debugging and reproduction of results

source('SANEst_helpers.R')

# ===== SAN experiment parameters
theta0 <- c(1,1,1,1,1)  # target (5 x 1)
m <- 1e1                # samples to estimate theta0 on each step
n <- 5e1                # simulation reps on each step
n.true <- 1e7           # simulation reps for true value
steps <- 1e2            # number of data collection steps
macro <- 1e1            # number of macro replications

# ===== Placeholders for results
theta <- array(0, c(steps, 5, macro))   # sample paths of thetahat
Est <- array(0, c(steps, 12, macro))    # 12 different estimators

# ===== Simulation experiment
# ----- True value
if(all.equal(c(1,1,1,1,1), theta0)){
    # numerical integration of known CDF
    Est.True <- integrate(function(x) 1-SAN.CDF(x), 0, Inf,
                          rel.tol = 1e-10, subdivisions = 1000L)$value
} else{
    # MC with a large sample size
    Est.True <- mean(SAN(n.true, theta0)$Y)
}

# ----- Estimates in differet macro reps
for(i in 1:macro){
    print(sprintf("Running the %i-th macro replication.", i))
    Out <- Target.Theta(theta0, m, n, steps)
    theta[, , i] <- Out$theta
    Est[, , i] <- Out$muhat
}

MSE <- apply(Est, c(1,2), function(x) mse(x, Est.True))
colnames(MSE) <- c("Crude", "Brute", "CumAvg", 
                   "PW.BLR", "PW.SLR", "PW.RLR",
                   "OW.BLR", "OW.SLR", "OW.RLR",
                   "SM.BLR", "SM.SLR", "SM.RLR")

save.image(file = "SAN_10K_Macro.RData")