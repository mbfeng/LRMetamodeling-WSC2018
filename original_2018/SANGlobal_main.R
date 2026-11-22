# ===== Workspace preparation
rm(list = ls()) # clear all variables
graphics.off() 
seed <- 1
set.seed(seed) # fix random seed for debugging and reproduction of results

source('SANGlobal_helpers.R')

load("SANGlobal_True0.RData")
colnames(Est)<-c("theta1","theta5","True")

# ===== SAN example
# ===== 2D example, the mean activity time for the first & fifth links vary
theta0 <- c(1,1,1,1,1)  # center theta
# other points to look at in a 2D square
theta1 <- seq(0.5, 1.5, length.out = 20)
theta5 <- seq(0.5, 1.5, length.out = 20)
grid15 <- expand.grid(x = theta1, y = theta5)
n.combo <- nrow(grid15)

n.samp <- 1e2       # simulation reps (sample size) at the center
n.true <- 1e7       # simulation reps for true value
n.macro <- 1e4      # number of macro reps to find out quantities of interest

# # ===== Loop to get "true estimates" (time consuming, run & save)
# Est <- data.frame(theta1 = grid15[,1], theta5 = grid15[,2])
# Est$val <- 0   # accurate estimate at each theta
# theta.cur <- theta0
# for(i in 1:n.combo){
#     print(sprintf("Running the %i-th combination.", i))
#     theta.cur[1] <- grid15[i,1]
#     theta.cur[5] <- grid15[i,2]
# 
#     Est$val[i] <- mean(SAN(n.true, theta.cur)$Y)
# }
#
# save.image("SANGlobal_True.RData")


# ===== Loop to get LR estimates
# ----- In each macro rep, draw a random sample at the center of the design
#       space. Use this random sample to estimate time to completion at other
#       design points. At each design point, we are intersted in:
#       * (1) bias & variance of the estimated value
#       * (2) bias, variance, & skewness of the likelihood ratios
#       There are three different LR estimators (BLR, SLR, RLR)
Est.samp.mean <- matrix(0, nrow = n.combo, ncol = n.macro)   # random sample var
Est.samp.var <- matrix(0, nrow = n.combo, ncol = n.macro)   # random sample var
Est.TTC.mean <- array(0, dim = c(n.combo, 3, n.macro))  # LR metamodel mean
Est.TTC.MSE <- array(0, dim = c(n.combo, 3, n.macro))   # LR metamodel MSE
Est.LR.mean <- matrix(0, nrow = n.combo, ncol = n.macro)
Est.LR.var <- matrix(0, nrow = n.combo, ncol = n.macro)

theta.cur <- theta0
for(j in 1:n.macro){
    print(sprintf("Running the %i-th macro replication", j))
    for(i in 1:n.combo){
        theta.cur[1] <- Est[i, 1]
        theta.cur[5] <- Est[i, 2]
        
        OutSamp <- SAN(n.samp, theta.cur)
        Est.samp.mean[i, j] <- mean(OutSamp$Y)
        Est.samp.var[i, j] <- var(OutSamp$Y)
        
        Out <- SAN(n.samp, theta0)
        
        theta.cur[1] <- grid15[i,1]
        theta.cur[5] <- grid15[i,2]
        
        POut <- Predict.LR(theta.cur, Out$Y, Out$X, Out$L)
        
        Est.LR.mean[i, j] <- mean(POut$LR)
        Est.LR.var[i, j] <- var(POut$LR)
        
        Est.TTC.mean[i, 1, j] <- mean(POut$Y[,1])   # BLR mean
        Est.TTC.mean[i, 2, j] <- mean(POut$Y[,2])   # SLR mean
        Est.TTC.mean[i, 3, j] <- mean(POut$Y[,3])   # RLR mean
        
        Est.TTC.MSE[i, 1, j] <- (Est.TTC.mean[i, 1, j]-Est$True[i])^2   # BLR mean
        Est.TTC.MSE[i, 2, j] <- (Est.TTC.mean[i, 2, j]-Est$True[i])^2   # SLR mean
        Est.TTC.MSE[i, 3, j] <- (Est.TTC.mean[i, 3, j]-Est$True[i])^2   # RLR mean
    }
}

filename <- sprintf("SANGlobal_%iSamp_%iMacro.RData", n.samp, n.macro)
save.image(filename)
