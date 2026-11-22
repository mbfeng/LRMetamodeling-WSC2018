################################################################################
# extract_published_runs.R                                                     #
#                                                                              #
# One-off provenance script (run once, November 2026). It reduces the two      #
# saved R workspaces of the published 2018 runs to compact summaries in        #
# results/published_2018/, from which reproduce_paper.R redraws the figures.   #
#                                                                              #
#   SANGlobal_100Samp_10000Macro.RData  (300 MB, 2018-04-08) -> Figure 1       #
#   SAN_10K_Macro.RData                 (123 MB, 2018-03-16) -> Figure 2       #
#                                                                              #
# Usage (from the Code folder):                                                #
#   Rscript tools/extract_published_runs.R <folder with the two .RData files>  #
################################################################################

source("R/summaries.R") # squared_error_summary()

args <- commandArgs(trailingOnly = TRUE)
data_folder <- args[1] # folder holding the original workspaces

# ===== Section 4.1: global metamodels (Figure 1)
global <- new.env()
load(file.path(data_folder, "SANGlobal_100Samp_10000Macro.RData"), envir = global)

true_value <- global$Est$True # mu(theta) at the 400 grid points (10^7 reps each)
n_macro <- global$n.macro     # 10,000 macro replications

#-- the grid used for the LR predictions must equal the grid of the true values
stopifnot(all(global$grid15[, 1] == global$Est$theta1))
stopifnot(all(global$grid15[, 2] == global$Est$theta5))

#-- squared errors [grid point i, macro run j] of the seven estimators
smc_error2 <- (global$Est.samp.mean - true_value)^2
estimator_names <- c("SMC", "BLR", "SLR", "RLR", "SMC+BLR", "SMC+SLR", "SMC+RLR")
global_summary <- data.frame(theta1 = global$Est$theta1,
                             theta5 = global$Est$theta5,
                             true_value = true_value)
summary_smc <- squared_error_summary(smc_error2)
global_summary <- cbind(global_summary, setNames(summary_smc,
                        paste0(names(summary_smc), "_SMC")))
for (k in 1:3) {
    lr_estimate <- global$Est.TTC.mean[, k, ]                 # 400 x 10,000
    lr_error2 <- (lr_estimate - true_value)^2
    combo_error2 <- (0.5 * global$Est.samp.mean + 0.5 * lr_estimate - true_value)^2
    summary_lr <- squared_error_summary(lr_error2)
    summary_combo <- squared_error_summary(combo_error2)
    global_summary <- cbind(global_summary,
        setNames(summary_lr, paste0(names(summary_lr), "_", estimator_names[1 + k])),
        setNames(summary_combo, paste0(names(summary_combo), "_", estimator_names[4 + k])))
}

#-- check: same numbers as the original plotting script (SANGlobal_plot.R)
original_ratio_rlr <- apply(global$Est.TTC.MSE, c(1, 2), mean)[, 3] /
    rowMeans((global$Est.samp.mean - true_value)^2)
stopifnot(isTRUE(all.equal(global_summary$mse_RLR / global_summary$mse_SMC,
                           original_ratio_rlr, tolerance = 1e-14)))

global_summary$n_macro <- n_macro
write.csv(global_summary, "results/published_2018/global_metamodel_summary.csv", row.names = FALSE)
rm(global)
gc()

# ===== Section 4.2: target theta0 (Figure 2)
target <- new.env()
load(file.path(data_folder, "SAN_10K_Macro.RData"), envir = target)

target_names <- colnames(target$MSE) # Crude, Brute, CumAvg, PW.BLR, ..., SM.RLR
target_summary <- NULL
for (k in seq_along(target_names)) {
    error2 <- (target$Est[, k, ] - target$Est.True)^2  # 100 steps x 10,000 macro runs
    summary_k <- squared_error_summary(error2)
    summary_k <- cbind(step = 1:target$steps, estimator = target_names[k], summary_k)
    target_summary <- rbind(target_summary, summary_k)
}
target_summary$n_macro <- target$macro
target_summary$true_value <- target$Est.True

#-- check: same numbers as the MSE matrix saved with the run
saved_mse <- as.vector(target$MSE) # column-major: all steps of Crude, then Brute, ...
stopifnot(isTRUE(all.equal(target_summary$mse, saved_mse, tolerance = 1e-14)))

write.csv(target_summary, "results/published_2018/target_theta_summary.csv", row.names = FALSE)

cat("Wrote the two summaries to results/published_2018/\n")

# ===== Inputs of the global experiment (data/)
# The true values mu(theta) (10^7 replications at each grid point) and the
# random-number state the 2018 run started from: SANGlobal_main.R called
# set.seed(1) and then load("SANGlobal_True0.RData"), which replaced the seed
# with the .Random.seed saved in that file.
true_file <- new.env()
load(file.path(data_folder, "SANGlobal_True0.RData"), envir = true_file)
true_values <- data.frame(theta1 = true_file$Est[, 1], theta5 = true_file$Est[, 2],
                          true_value = true_file$Est[, 3])
stopifnot(all(true_values$true_value == global_summary$true_value))
dir.create("data", showWarnings = FALSE)
write.csv(true_values, "data/global_true_values.csv", row.names = FALSE)
writeLines(as.character(get(".Random.seed", envir = true_file)),
           "data/global_rng_state_2018.txt")
cat("Wrote data/global_true_values.csv and data/global_rng_state_2018.txt\n")
cat(sprintf("Settings: global n.macro = %d; target m = %d, n = %d, steps = %d, macro = %d\n",
            n_macro, target$m, target$n, target$steps, target$macro))
