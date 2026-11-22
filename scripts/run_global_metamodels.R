################################################################################
# run_global_metamodels.R                                                      #
#                                                                              #
# Section 4.1, Figure 1: LR metamodels from one design point theta* over the   #
# 20 x 20 grid of (theta1, theta5) in [0.5, 1.5]^2, 10^4 macro replications.   #
# Writes <out>/global_metamodel_summary.csv (MSE of SMC, BLR, SLR, RLR and of  #
# the averages SMC + BLR, SMC + SLR, SMC + RLR at each grid point).            #
#                                                                              #
# Usage (from the Code folder):                                                #
#   Rscript scripts/run_global_metamodels.R --out results/rerun                #
#           [--macro 10000] [--seed 2018]                                      #
#   --seed 2018 (default) starts from the random-number state of the published #
#   2018 run and reproduces it exactly; any other value calls set.seed(value). #
# The run saves a checkpoint every 500 macro replications and resumes from it. #
################################################################################

source("R/san_model.R")   # simulate_san(), san_density(), lr_estimates()
source("R/experiments.R") # global_macro_run()
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
seed <- get_argument("seed", "2018")

# ===== Settings (Section 4.1)
theta_center <- c(1, 1, 1, 1, 1) # design point theta*
n_reps <- 100                    # replications per sample
true_values <- read.csv("data/global_true_values.csv") # mu(theta), 10^7 reps each
grid <- true_values[, c("theta1", "theta5")]           # 20 x 20 grid, theta1 fastest
n_points <- nrow(grid)

# ===== Random numbers
if (seed == "2018") {
    state <- as.integer(readLines("data/global_rng_state_2018.txt"))
    RNGkind("Mersenne-Twister", "Inversion", "Rejection") # kinds stored in state[1]
    assign(".Random.seed", state, envir = globalenv())     # state of the 2018 run
} else {
    set.seed(as.integer(seed))
}

# ===== Macro replications (resumable)
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
checkpoint_file <- file.path(out_folder, sprintf("global_checkpoint_seed%s.rds", seed))
estimates <- array(0, dim = c(n_points, 4, n_macro)) # [grid point, SMC/BLR/SLR/RLR, macro]
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
    for (j in first_macro:n_macro) {
        estimates[, , j] <- global_macro_run(grid, theta_center, n_reps)
        if (j %% 500 == 0 || j == n_macro) {
            checkpoint <- list(done = j, estimates = estimates[, , 1:j, drop = FALSE],
                               random_seed = get(".Random.seed", envir = globalenv()))
            saveRDS(checkpoint, checkpoint_file)
            cat(sprintf("%d macro replications done (%.1f minutes)\n", j,
                        as.numeric(difftime(Sys.time(), start, units = "mins"))))
        }
    }
}

# ===== Summaries: squared errors [grid point, macro replication]
mu <- true_values$true_value
smc <- estimates[, 1, ]
summary <- data.frame(theta1 = grid$theta1, theta5 = grid$theta5, true_value = mu)
add_summary <- function(summary, error2, name) {
    # append mse_<name>, se_mse_<name>, top1_share_<name>
    columns <- squared_error_summary(error2)
    names(columns) <- paste0(names(columns), "_", name)
    return(cbind(summary, columns))
}
summary <- add_summary(summary, (smc - mu)^2, "SMC")
lr_names <- c("BLR", "SLR", "RLR")
for (k in 1:3) {
    lr <- estimates[, 1 + k, ]
    summary <- add_summary(summary, (lr - mu)^2, lr_names[k])
    combined <- 0.5 * smc + 0.5 * lr # average of SMC and LR metamodel (Figure 1f-h)
    summary <- add_summary(summary, (combined - mu)^2, paste0("SMC+", lr_names[k]))
}
summary$n_macro <- n_macro
summary_file <- file.path(out_folder, "global_metamodel_summary.csv")
write.csv(summary, summary_file, row.names = FALSE)
cat("Wrote", summary_file, "\n")
