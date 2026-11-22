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
    
    return(list(LR = LR, beta = beta, Y = cbind(Y.BLR, Y.SLR, Y.RLR)))
}