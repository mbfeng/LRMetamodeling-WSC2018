################################################################################
# test_against_original_code.R                                                 #
#                                                                              #
# The rewritten functions must give the same numbers as the 2018 code in       #
# original_2018/ on the same random numbers (up to floating-point rounding:    #
# the densities are computed as prod(theta) exp(-X theta) instead of a product #
# of dexp() values).                                                           #
#                                                                              #
# Usage (from the Code folder):  Rscript tests/test_against_original_code.R    #
################################################################################

source("R/san_model.R")
source("R/experiments.R")
original <- new.env() # the 2018 functions, kept apart from the new ones
sys.source("original_2018/SANEst_helpers.R", envir = original)

tolerance <- 1e-10
n_failed <- 0
report <- function(name, new, old) {
    # print one comparison and count failures
    difference <- max(abs(new - old) / pmax(1, abs(old)))
    passed <- difference <= tolerance
    cat(sprintf("%-55s max rel. difference %.1e  %s\n", name, difference,
                ifelse(passed, "ok", "FAILED")))
    if (!passed) {
        n_failed <<- n_failed + 1
    }
}

# ===== The SAN and the basic LR metamodels
theta <- c(0.7, 1, 1, 1, 1.3)
set.seed(11)
new_sample <- simulate_san(1000, theta)
set.seed(11)
old_sample <- original$SAN(1000, theta)
report("simulate_san(): activity times X", new_sample$X, old_sample$X)
report("simulate_san(): completion times Y", new_sample$Y, as.vector(old_sample$Y))
report("simulate_san(): densities", new_sample$density, old_sample$L)

theta_new <- c(1.2, 1, 1, 1, 0.6)
lr <- san_density(new_sample$X, theta_new) / new_sample$density
old_lr <- original$Predict.LR(theta_new, old_sample$Y, old_sample$X, old_sample$L)
report("lr_estimates(): BLR, SLR, RLR", lr_estimates(new_sample$Y, lr),
       colMeans(old_lr$Y))

# ===== Section 4.1: one macro replication on the first 40 grid points
grid <- read.csv("data/global_true_values.csv")[1:40, c("theta1", "theta5")]
theta_center <- c(1, 1, 1, 1, 1)
set.seed(12)
new_global <- global_macro_run(grid, theta_center, 100)
set.seed(12)
old_global <- matrix(0, nrow = 40, ncol = 4)
theta_cur <- theta_center
for (i in 1:40) { # the loop body of original_2018/SANGlobal_main.R
    theta_cur[1] <- grid[i, 1]
    theta_cur[5] <- grid[i, 2]
    out_samp <- original$SAN(100, theta_cur)
    out <- original$SAN(100, theta_center)
    p_out <- original$Predict.LR(theta_cur, out$Y, out$X, out$L)
    old_global[i, ] <- c(mean(out_samp$Y), colMeans(p_out$Y))
}
report("global_macro_run(): SMC, BLR, SLR, RLR (40 points)", new_global, old_global)

# ===== Section 4.2: one macro replication with 30 steps
set.seed(13)
new_target <- target_macro_run(c(1, 1, 1, 1, 1), 10, 50, 30)
set.seed(13)
old_target <- original$Target.Theta(c(1, 1, 1, 1, 1), 10, 50, 30)
report("target_macro_run(): thetahat_t", new_target$theta, old_target$theta)
report("target_macro_run(): 12 estimators x 30 steps", new_target$estimates,
       old_target$muhat)

cat(ifelse(n_failed == 0, "All tests passed\n", sprintf("%d tests FAILED\n", n_failed)))
if (n_failed > 0) {
    quit(status = 1)
}
