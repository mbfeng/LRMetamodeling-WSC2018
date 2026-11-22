################################################################################
# experiments.R                                                                #
#                                                                              #
# One macro replication of each experiment of Section 4.                       #
#   global_macro_run(): Section 4.1, LR metamodels from one design point       #
#                       (the centre theta* = (1,1,1,1,1)) over a 20 x 20 grid  #
#   target_macro_run(): Section 4.2, LR metamodels from the design points      #
#                       thetahat_1, ..., thetahat_t that converge to theta0     #
# Both draw their random numbers in the same order as the 2018 code.           #
################################################################################

# ===== Section 4.1: global metamodels

global_macro_run <- function(grid, theta_center, n_reps) {
    # One macro replication of the global experiment (Figure 1)
    # Input(s):
    #   -- grid: data frame with columns theta1, theta5 (the 400 grid points)
    #   -- theta_center: design point theta* (rates of all 5 activities)
    #   -- n_reps: replications per sample (100)
    # Output(s): 400 x 4 matrix, one row per grid point, columns
    #   -- SMC: sample mean of n_reps replications simulated at the grid point
    #   -- BLR, SLR, RLR: LR metamodel predictions at the grid point from a
    #                     new sample of n_reps replications at theta_center
    # NOTE: as in the 2018 run, each grid point gets its own sample at the
    #       centre (the paper's MSEs are per grid point, so this does not change
    #       them); the order of the random draws is the 2018 order.

    n_points <- nrow(grid)
    estimates <- matrix(0, nrow = n_points, ncol = 4,
                        dimnames = list(NULL, c("SMC", "BLR", "SLR", "RLR")))
    for (i in 1:n_points) {
        theta <- theta_center           # prediction point theta = (theta1,1,1,1,theta5)
        theta[1] <- grid$theta1[i]
        theta[5] <- grid$theta5[i]

        #-- standard Monte Carlo (SMC) at the grid point itself
        smc_sample <- simulate_san(n_reps, theta)
        estimates[i, "SMC"] <- mean(smc_sample$Y)

        #-- LR metamodels from a sample at the centre
        center_sample <- simulate_san(n_reps, theta_center)
        lr <- san_density(center_sample$X, theta) / center_sample$density
        estimates[i, c("BLR", "SLR", "RLR")] <- lr_estimates(center_sample$Y, lr)
    }
    return(estimates)
}

# ===== Section 4.2: target theta0

target_estimator_names <- c("Crude", "Brute", "CumAvg",
                            "PW.BLR", "PW.SLR", "PW.RLR",
                            "OW.BLR", "OW.SLR", "OW.RLR",
                            "SM.BLR", "SM.SLR", "SM.RLR")

block_means_variances <- function(values, n_reps) {
    # Mean and variance of each block of n_reps consecutive values
    # Input(s):
    #   -- values: vector of length n_reps * number of blocks
    #   -- n_reps: block size
    # Output(s): list with vectors mean and variance (one entry per block)
    # formula: variance = (mean(x^2) - mean(x)^2) n / (n - 1), as in 2018
    blocks <- matrix(values, nrow = n_reps) # column j = block j
    block_mean <- colMeans(blocks)
    block_variance <- (colMeans(blocks^2) - block_mean^2) * n_reps / (n_reps - 1)
    return(list(mean = block_mean, variance = block_variance))
}

target_macro_run <- function(theta0, n_data, n_reps, n_steps) {
    # One macro replication of the target-theta0 experiment (Figure 2)
    # Input(s):
    #   -- theta0: target rates (1,1,1,1,1)
    #   -- n_data: new observations of the activity times per step (m = 10)
    #   -- n_reps: simulation replications per step (n = 50)
    #   -- n_steps: number of steps t = 1, ..., 100
    # Output(s): list with
    #   -- estimates: n_steps x 12 matrix of estimates of mu(theta0); columns
    #        Crude  = mean of the n_reps outputs of step t
    #        Brute  = mean of n_reps * t new outputs at thetahat_t
    #        CumAvg = mean of all outputs so far (not in the paper's figures)
    #        PW.*, OW.*, SM.* = BLR, SLR, RLR combined over the t design points by
    #                 proportional weights (PW: n_s / n, equal here because
    #                 every step has n_reps outputs), optimal weights (OW:
    #                 inverse of the estimated variances) or the mixture LR
    #                 ("M" in the paper; prefix SM, as in the 2018 code)
    #   -- theta: n_steps x 5 matrix of the estimates thetahat_t
    # NOTE: all LR estimators are formed at the CURRENT estimate thetahat_t,
    #       from the outputs of all steps s <= t.

    n_total <- n_reps * n_steps
    #-- input data for estimating theta0 (all drawn at the start, as in 2018)
    input_data <- cbind(rexp(n_data * n_steps, rate = theta0[1]),
                        rexp(n_data * n_steps, rate = theta0[2]),
                        rexp(n_data * n_steps, rate = theta0[3]),
                        rexp(n_data * n_steps, rate = theta0[4]),
                        rexp(n_data * n_steps, rate = theta0[5]))

    estimates <- matrix(0, nrow = n_steps, ncol = 12,
                        dimnames = list(NULL, target_estimator_names))
    theta_hat <- matrix(0, nrow = n_steps, ncol = 5)
    X_all <- matrix(0, nrow = n_total, ncol = 5) # outputs of all steps, stacked
    Y_all <- numeric(n_total)
    density_all <- numeric(n_total)              # f(X | thetahat_s) of its own step s
    mixture_sum <- numeric(n_total)              # sum_{s <= t} f(X | thetahat_s)

    for (t in 1:n_steps) {
        first <- (t - 1) * n_reps + 1            # rows of step t
        last <- t * n_reps
        new_rows <- first:last
        rows_so_far <- 1:last

        #-- estimate theta0 from all data so far: rate = 1 / sample mean
        theta_t <- 1 / colMeans(input_data[1:(n_data * t), , drop = FALSE])
        theta_hat[t, ] <- theta_t

        #-- simulate step t at thetahat_t
        step_sample <- simulate_san(n_reps, theta_t)
        X_all[new_rows, ] <- step_sample$X
        Y_all[new_rows] <- step_sample$Y
        density_all[new_rows] <- step_sample$density

        if (t == 1) {
            estimates[t, ] <- mean(step_sample$Y) # only one design point: all equal
            mixture_sum[new_rows] <- step_sample$density
            next
        }
        estimates[t, "Crude"] <- mean(step_sample$Y)
        brute_sample <- simulate_san(n_reps * t, theta_t)
        estimates[t, "Brute"] <- mean(brute_sample$Y)
        estimates[t, "CumAvg"] <- mean(Y_all[rows_so_far])

        #-- likelihood ratios of all outputs so far, at the current thetahat_t
        density_now <- san_density(X_all[rows_so_far, , drop = FALSE], theta_t)
        lr <- density_now / density_all[rows_so_far]
        y_lr <- Y_all[rows_so_far] * lr

        #-- PW and OW: BLR, SLR, RLR within each step s (a block of n_reps rows),
        #   then combined over the t blocks
        lr_blocks <- block_means_variances(lr, n_reps)
        blr_values <- y_lr                                       # BLR per replication
        slr_values <- y_lr / rep(lr_blocks$mean, each = n_reps)  # SLR per replication
        #   RLR slope per block: beta_s = sum((lr - mean lr) y_lr) / sum((lr - mean lr)^2)
        lr_centered <- lr - rep(lr_blocks$mean, each = n_reps)
        cross <- colSums(matrix(lr_centered * y_lr, nrow = n_reps))
        squares <- colSums(matrix(lr_centered^2, nrow = n_reps))
        beta <- ifelse(squares == 0, 0, cross / squares)  # beta = 0 -> RLR = BLR
        rlr_values <- y_lr - rep(beta, each = n_reps) * (lr - 1)

        basic_values <- list(BLR = blr_values, SLR = slr_values, RLR = rlr_values)
        for (name in c("BLR", "SLR", "RLR")) {
            blocks <- block_means_variances(basic_values[[name]], n_reps)
            estimates[t, paste0("PW.", name)] <- mean(blocks$mean)
            #   OW: weights proportional to 1 / (sample variance of the per-
            #   replication values of step s), as in the 2018 code
            weights <- (1 / blocks$variance) / sum(1 / blocks$variance)
            estimates[t, paste0("OW.", name)] <- sum(weights * blocks$mean)
        }

        #-- SM ("M"): mixture LR f(X | theta_t) / ((1/t) sum_{s <= t} f(X | thetahat_s))
        old_rows <- 1:(first - 1)
        mixture_sum[old_rows] <- mixture_sum[old_rows] + density_now[old_rows]
        for (s in 1:t) { # densities of the new outputs at all design points so far
            mixture_sum[new_rows] <- mixture_sum[new_rows] +
                san_density(step_sample$X, theta_hat[s, ])
        }
        mixture_density <- mixture_sum[rows_so_far] / t
        mixture_lr <- density_now / mixture_density
        estimates[t, c("SM.BLR", "SM.SLR", "SM.RLR")] <-
            lr_estimates(Y_all[rows_so_far], mixture_lr)
    }
    return(list(estimates = estimates, theta = theta_hat))
}
