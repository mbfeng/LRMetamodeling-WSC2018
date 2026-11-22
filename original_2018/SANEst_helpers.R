mse <- function(Y, truth){
    # calculate the mean sqaured errors (MSE) between Y & truth
    # input:
    #   Y = scalar, vector, or matrix of estimated values
    #   truth = scalar, vector, or matrix of true values
    # output:
    #   MSE = mean squared error
    
    return(MSE = mean((Y - truth)^2))
}

colMeansVars <- function(data, na.rm = TRUE){
    # compute the column-wise means and variances
    # input
    #   data = 2D numeric matrix
    # output
    #   colmean = column means
    #   colvar  = column variances
    
    if (na.rm) {
        n <- colSums(!is.na(data))
    } else {
        n <- nrow(data)
    }
    
    colmean <- colMeans(data, na.rm = na.rm)
    colvar <- (colMeans(data*data, na.rm = na.rm) - (colmean)^2) * n/(n-1)
    
    list(colmean = colmean, colvar = colvar)
}

SAN.CDF <- function(t){
    # calculate cdf Pr(Y \le t) of SAN with theta = c(1,1,1,1,1)
    
    return((t^2/2 - 3*t - 3)*exp(-2*t) + (-t^2/2 - 3*t + 3)*exp(-t) + 1 - exp(-3*t))
}

SAN.PDF.X <- function(X, theta){
    # Joint pdf prod_i f(X_i|theta_i) of SAN example
    # input
    #   X = n x 5 matrix of generated activity times
    #   theta = vector of exponential RATES for activities 1,2,...,5
    # output
    #   L = n x 1 vector of liklihoods, i.e., prod_i f(X_i | theta_i)
    L <- apply(cbind(dexp(X[, 1], rate = theta[1]), 
                     dexp(X[, 2], rate = theta[2]), 
                     dexp(X[, 3], rate = theta[3]),
                     dexp(X[, 4], rate = theta[4]), 
                     dexp(X[, 5], rate = theta[5])), 1, prod)
    return(L = L)
}

SAN <- function(n, theta){
    # generates reps of the SAN with exponential activity times
    # input
    #   n = number of replications
    #   theta = vector of exponential RATES for activities 1,2,...,5
    # output
    #   Y = vector of times to complete the project
    #   X = n x 5 matrix of generated activity times
    #   L = n x 1 vector of liklihoods, i.e., prod_i(f(X_i | theta))
    
    X <- cbind(rexp(n, rate=theta[1]), 
               rexp(n, rate=theta[2]), 
               rexp(n, rate=theta[3]),
               rexp(n, rate=theta[4]), 
               rexp(n, rate=theta[5]))
    
    A <- cbind(X[,1]+X[,4], X[,1]+X[,3]+X[,5], X[,2]+X[,5])
    
    Y <- as.matrix(apply(A, 1, max))
    
    L <- SAN.PDF.X(X, theta) 
    
    return(list(Y = Y, X = X, L = L))
}

Target.Theta <- function(theta0, m, n, steps){
    # Experiment with thetahat converging to target theta0
    # input
    #   theta0 = 5 x 1 vector of target theta
    #   m = samples to estimate theta0 on each step
    #   n = simulation reps on each step
    #   steps = number of data collection steps
    # output
    #   theta = sequence of thetahats
    #   muhat = steps x 12 matrix 
    #       the 12 estimators are
    #       -- [Crude, Brute, Cumulative Avg] (3 benchmark estimators)
    #       -- [BLR, SLR, RLR] x [PW, OW, SM] (9 LR estimators)
    #
    # NOTE: All estimators (except Brute) share the same outputs 
    
    # auxliary quantities
    nsteps <- n*steps
    msteps <- m*steps
    
    # data for estimating thetahat in each step
    Data <- cbind(rexp(msteps, rate=theta0[1]), rexp(msteps, rate=theta0[2]), 
                  rexp(msteps, rate=theta0[3]), rexp(msteps, rate=theta0[4]), 
                  rexp(msteps, rate=theta0[5]))
    
    # placeholders for results
    theta <- matrix(0, nrow = steps, ncol = 5)
    muhat <- matrix(0, nrow = steps, ncol = 12)
    colnames(muhat) <- c("Crude", "Brute", "CumAvg", 
                         "PW.BLR", "PW.SLR", "PW.RLR",
                         "OW.BLR", "OW.SLR", "OW.RLR",
                         "SM.BLR", "SM.SLR", "SM.RLR")
    
    # storage vectors/matrices for output reusage
    Y.Out <- matrix(0, nrow = nsteps, ncol= 1)
    X.Out <- matrix(0, nrow = nsteps, ncol= 5)
    L.Out <- matrix(0, nrow = nsteps, ncol= 1)  # sampling LR
    L.Mix <- matrix(0, nrow = nsteps, ncol= 1)  # mixture LR
    
    for (i in 1:steps){
        # lower and upper indices for storage vectors/matrices in current step
        li <- (i-1)*n + 1
        ui <- i*n
        
        # estimate thetahat using all data so far
        thetahat <- 1/apply(Data[1:(m*i), ], 2, mean)   
        theta[i, ] <- thetahat  # store current estimate
        
        # SAN simulation in the current step & store the outputs
        Out <- SAN(n, thetahat)
        Y.Out[li:ui] <- Out$Y
        X.Out[li:ui, ] <- Out$X
        L.Out[li:ui] <- Out$L
        L.Mix[li:ui] <- Out$L   # mixture LR requires futher calculations
        
        if(i == 1){
            # first experiment, same estimate for all (save some computations)
            muhat[i, ] = mean(Out$Y)
        }else{
            # subsequent experiments, estimators start to differ from each other
            muhat[i, 1] = mean(Out$Y)                   # curde MC
            muhat[i, 2] = mean(SAN(n*i, thetahat)$Y)    # brute MC
            muhat[i, 3] = mean(Y.Out[1:ui])             # cumulative avg
            
            # PW-LR & OW-LR estimators
            LR.mean <- matrix(0, nrow = 3, ncol = i)
            LR.var <- matrix(0, nrow = 3, ncol = i)
            for(j in 1:i){
                lj <- (j-1)*n + 1
                uj <- j*n
                
                mu.temp <- Predict.LR(thetahat, Y.Out[lj:uj], 
                                      X.Out[lj:uj, ], L.Out[lj:uj])
                
                mvj <- colMeansVars(mu.temp$Y)
                LR.mean[, j] <- mvj$colmean
                LR.var[, j] <- mvj$colvar
            }
            muhat[i, 4] <- mean(LR.mean[1,])    # PW.BLR
            muhat[i, 5] <- mean(LR.mean[2,])    # PW.SLR
            muhat[i, 6] <- mean(LR.mean[3,])    # PW.RLR
            
            LR.varinv <- LR.var^(-1)
            LR.varinvsum <- rowSums(LR.varinv)
            wt.OW <- LR.varinv / LR.varinvsum
            LR.mean.OW <- wt.OW * LR.mean

            muhat[i, 7] <- sum(LR.mean.OW[1,])    # OW.BLR
            muhat[i, 8] <- sum(LR.mean.OW[2,])    # OW.SLR
            muhat[i, 9] <- sum(LR.mean.OW[3,])    # OW.RLR
            
            # compute the likelihoods required for mixture LRs
            L.Mix[1:(li - 1), ] <- L.Mix[1:(li - 1), ]*(i - 1) +
                SAN.PDF.X(X.Out[1:(li - 1), ], thetahat)
            
            for(j in 1:(i - 1)){
                L.Mix[li:ui, ] <- L.Mix[li:ui, ] +
                    SAN.PDF.X(X.Out[li:ui, ], theta[j, ])
            }
            L.Mix[1:ui, ] <- L.Mix[1:ui, ]/i
            
            mu.temp <- Predict.LR(thetahat, Y.Out[1:ui], 
                                  X.Out[1:ui, ], L.Mix[1:ui])$Y
            
            muhat[i, 10] <- mean(mu.temp[, 1])    # SM.BLR
            muhat[i, 11] <- mean(mu.temp[, 2])    # SM.SLR
            muhat[i, 12] <- mean(mu.temp[, 3])    # SM.RLR
        }
    }
    
    return(list(muhat = muhat, theta = theta))
}

Predict.LR <- function(theta0, Y, X, L){
    # LR estimates to predict SAN output at theta0 given SAN output from theta1
    # NOTE: theta1 is not needed, likelihoods are characterized by L
    # input
    #   theta0 = new theta value
    #   theta1 = vector of theta1 values where SAN output was simulated
    #   Y = output at theta1
    #   X = inputs at theta1
    #   L = matrix of likelihoods at theta1
    # output
    #   LR = n x 1 vector of likelihood ratios
    #   Y.BLR = Y*LR, n x 1 vector of Baseline LR estimates
    #   Y.SLR = Y.BLR / (mean LR), n x 1 vector of Self-normalizing LR estimates
    #   Y.RLR = Y.BLR - beta *(mean LR - 1), n x 1 vector of Regression LR estimates
    
    L0 <- SAN.PDF.X(X, theta0)  # likelihood of output X at theta0
    LR <- L0/L                  # likelihood ratios
    
    Y.BLR <- Y*LR               # basline LR
    
    LR.mean <- mean(LR)
    Y.SLR <- Y.BLR / LR.mean    # self-normalizing LR
    
    LR.dif <- LR - LR.mean
    LR.dif2 <- sum(LR.dif^2)
    if(LR.dif2 == 0){
        # All LR are the same, Y.RLR fails, so use Y.BLR
        Y.RLR <- Y.BLR
    } else {
        beta <- sum(LR.dif * Y.BLR)/ sum(LR.dif^2)
        Y.RLR <- Y.BLR - beta*(LR - 1)  # regression LR
    }
    
    return(list(LR = LR, Y = cbind(Y.BLR, Y.SLR, Y.RLR)))
}