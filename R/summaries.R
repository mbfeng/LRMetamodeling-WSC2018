################################################################################
# summaries.R                                                                  #
#                                                                              #
# Summaries of squared errors over macro replications, shared by the          #
# simulation scripts and tools/extract_published_runs.R.                       #
################################################################################

squared_error_summary <- function(error2) {
    # MSE, its standard error and the share of the worst 1% of macro runs
    # Input:
    #   -- error2: matrix of squared errors, [point, macro replication]
    # Output: data frame with one row per point and columns
    #   -- mse: mean squared error over the macro replications
    #   -- se_mse: standard error of mse (sd of squared errors / sqrt(K))
    #   -- top1_share: share of mse from the largest 1% of squared errors
    #                  (> 0.5 means the MSE is dominated by a few runs)
    # formula: mse = (1/K) sum_k (muhat_k - mu)^2

    n_macro <- ncol(error2)                        # K
    n_top <- max(1, round(n_macro / 100))          # number of runs in the worst 1%
    mse <- rowMeans(error2)                        # MSE per point
    second_moment <- rowMeans(error2^2)
    variance <- (second_moment - mse^2) * n_macro / (n_macro - 1) # sample variance
    se_mse <- sqrt(pmax(variance, 0) / n_macro)
    top1_share <- numeric(nrow(error2))
    for (i in 1:nrow(error2)) {
        largest <- sort(error2[i, ], decreasing = TRUE)[1:n_top]
        top1_share[i] <- sum(largest) / sum(error2[i, ])
    }
    return(data.frame(mse = mse, se_mse = se_mse, top1_share = top1_share))
}
