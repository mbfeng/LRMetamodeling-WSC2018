################################################################################
# san_model.R                                                                  #
#                                                                              #
# The stochastic activity network (SAN) of Section 4 and the three            #
# basic LR metamodels of Section 3.1:                                          #
#   BLR (baseline), SLR (self-normalized), RLR (regression / control variate). #
#                                                                              #
# Notation (code -> paper):                                                    #
#   theta    -> theta = (theta_1, ..., theta_5), exponential activity RATES    #
#   X        -> activity times X_1, ..., X_5 (one row per replication)         #
#   Y        -> time to completion Y = max(X1+X4, X1+X3+X5, X2+X5)             #
#   density  -> f(X | theta) = prod_k theta_k exp(-theta_k X_k)                #
#   lr       -> likelihood ratio f(X | prediction point) / f(X | design point)  #
################################################################################

simulate_san <- function(n_reps, theta) {
    # Simulate the SAN with exponential activity times
    # Input(s):
    #   -- n_reps: number of replications
    #   -- theta: vector of the 5 exponential rates
    # Output(s): list with
    #   -- X: n_reps x 5 matrix of activity times
    #   -- Y: vector of n_reps times to completion
    #   -- density: vector of n_reps joint densities f(X | theta)
    # NOTE: the five columns are drawn one after another with rexp(n_reps, rate),
    #       exactly as in the 2018 code, so a given seed gives the same numbers.

    X <- cbind(rexp(n_reps, rate = theta[1]),
               rexp(n_reps, rate = theta[2]),
               rexp(n_reps, rate = theta[3]),
               rexp(n_reps, rate = theta[4]),
               rexp(n_reps, rate = theta[5]))
    Y <- completion_time(X)
    density <- san_density(X, theta)
    return(list(X = X, Y = Y, density = density))
}

completion_time <- function(X) {
    # Time to complete the project: the longest of the three paths
    # Input(s):
    #   -- X: n x 5 matrix of activity times
    # Output(s): vector of n completion times
    # formula: Y = max(X1 + X4, X1 + X3 + X5, X2 + X5)
    path_a <- X[, 1] + X[, 4]
    path_b <- X[, 1] + X[, 3] + X[, 5]
    path_c <- X[, 2] + X[, 5]
    return(pmax(path_a, path_b, path_c))
}

san_density <- function(X, theta) {
    # Joint density of the activity times, one value per row of X
    # Input(s):
    #   -- X: n x 5 matrix of activity times
    #   -- theta: vector of the 5 exponential rates
    # Output(s): vector of n densities
    # formula: f(X | theta) = prod_k theta_k exp(-theta_k X_k)
    #                       = prod(theta) * exp(-sum_k theta_k X_k)
    # sum_k theta_k X_k, one value per row. Written out column by column (not
    # as X %*% theta) so each row gets exactly the same rounding whatever the
    # size of X: with an optimized BLAS, X %*% theta can round a row differently
    # in a large matrix than in a small one, and target_macro_run() relies on
    # f(X | theta) / f(X | theta) being exactly 1 for the current step.
    exponent <- theta[1] * X[, 1] + theta[2] * X[, 2] + theta[3] * X[, 3] +
        theta[4] * X[, 4] + theta[5] * X[, 5]
    return(prod(theta) * exp(-exponent))
}

lr_estimates <- function(Y, lr) {
    # The three basic LR metamodel estimates from one sample (Section 3.1)
    # Input(s):
    #   -- Y: vector of n outputs simulated at the design point
    #   -- lr: vector of n likelihood ratios f(X | theta0) / f(X | design point)
    # Output(s): vector c(BLR, SLR, RLR) of estimates of mu(theta0)
    # formula: BLR = mean(Y lr)
    #          SLR = mean(Y lr) / mean(lr)
    #          RLR = mean(Y lr) - beta (mean(lr) - 1), with beta the least-squares
    #                slope of Y lr on lr (lr has known mean 1: control variate)

    y_lr <- Y * lr
    blr <- mean(y_lr)
    lr_mean <- mean(lr)
    slr <- blr / lr_mean
    lr_centered <- lr - lr_mean
    sum_squares <- sum(lr_centered^2)
    if (sum_squares == 0) {
        rlr <- blr # all lr equal (theta0 = design point): no regression possible
    } else {
        beta <- sum(lr_centered * y_lr) / sum_squares
        rlr <- blr - beta * (lr_mean - 1)
    }
    return(c(BLR = blr, SLR = slr, RLR = rlr))
}
