################################################################################
# run_target_theta.R                                                           #
#                                                                              #
# Section 4.2, Figure 2: LR metamodels from the design points thetahat_1, ..., #
# thetahat_t that converge to theta0 = (1,1,1,1,1); m = 10 new observations    #
# and n = 50 replications per step, t = 1, ..., 100; 10^4 macro replications.  #
# Writes <out>/target_theta_summary.csv (MSE of the 12 estimators per step).   #
#                                                                              #
# Usage (from the Code folder):                                                #
#   Rscript scripts/run_target_theta.R --out results/rerun                     #
#           [--macro 10000] [--seed 0]                                         #
#   --seed 0 (default) is the seed of the published 2018 run, which it         #
#   reproduces exactly.                                                        #
# The run saves a checkpoint every 500 macro replications and resumes from it. #
################################################################################

source("R/san_model.R")   # simulate_san(), san_density(), lr_estimates()
source("R/experiments.R") # target_macro_run(), target_estimator_names
source("R/summaries.R")   # squared_error_summary()

# ===== Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
get_argument <- function(name, default) {
    # value after "--name" on the command line, or the default
    position <- match(paste0("--", name), args)
    if (is.na(position)) {
        return(default)
    }
    return(args[position + 1])
}
out_folder <- get_argument("out", "results/rerun")
n_macro <- as.integer(get_argument("macro", "10000"))
seed <- as.integer(get_argument("seed", "0"))

# ===== Settings (Section 4.2)
theta0 <- c(1, 1, 1, 1, 1) # target theta0
n_data <- 10               # m: new observations of the activity times per step
n_reps <- 50               # n: simulation replications per step
n_steps <- 100             # t = 1, ..., 100
#-- true value mu(theta0): for theta0 = (1,...,1) the distribution of Y is known,
#   P(Y <= y) = (y^2/2 - 3y - 3) e^{-2y} + (-y^2/2 - 3y + 3) e^{-y} + 1 - e^{-3y},
#   and E[Y] = integral of P(Y > y) = 83/24 (the 2018 run integrated numerically)
true_value <- 83 / 24

# ===== Macro replications (resumable)
set.seed(seed)
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
checkpoint_file <- file.path(out_folder, sprintf("target_checkpoint_seed%d.rds", seed))
estimates <- array(0, dim = c(n_steps, 12, n_macro)) # [step t, estimator, macro]
first_macro <- 1
if (file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    estimates[, , 1:checkpoint$done] <- checkpoint$estimates[, , 1:checkpoint$done]
    first_macro <- checkpoint$done + 1
    assign(".Random.seed", checkpoint$random_seed, envir = globalenv())
    cat("Resuming after macro replication", checkpoint$done, "\n")
}

start <- Sys.time()
if (first_macro <= n_macro) {
    for (k in first_macro:n_macro) {
        run <- target_macro_run(theta0, n_data, n_reps, n_steps)
        estimates[, , k] <- run$estimates
        if (k %% 500 == 0 || k == n_macro) {
            checkpoint <- list(done = k, estimates = estimates[, , 1:k, drop = FALSE],
                               random_seed = get(".Random.seed", envir = globalenv()))
            saveRDS(checkpoint, checkpoint_file)
            cat(sprintf("%d macro replications done (%.1f minutes)\n", k,
                        as.numeric(difftime(Sys.time(), start, units = "mins"))))
        }
    }
}

# ===== Summaries: squared errors [step, macro replication] per estimator
summary <- NULL
for (e in seq_along(target_estimator_names)) {
    error2 <- matrix(estimates[, e, ], nrow = n_steps)
    error2 <- (error2 - true_value)^2
    rows <- cbind(step = 1:n_steps, estimator = target_estimator_names[e],
                  squared_error_summary(error2))
    summary <- rbind(summary, rows)
}
summary$n_macro <- n_macro
summary$true_value <- true_value
summary_file <- file.path(out_folder, "target_theta_summary.csv")
write.csv(summary, summary_file, row.names = FALSE)
cat("Wrote", summary_file, "\n")
