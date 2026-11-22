################################################################################
# reproduce_paper.R                                                            #
#                                                                              #
# Dong, J., Feng, M. B., and Nelson, B. L. (2018). Unbiased Metamodeling via   #
# Likelihood Ratios. Proceedings of the 2018 Winter Simulation Conference,     #
# 1778-1789. https://doi.org/10.1109/WSC.2018.8632506                          #
#                                                                              #
# Entry point. By default it redraws every figure of the paper with numerical  #
# results (Figures 1 and 2) from the summaries in results/ (about a minute),   #
# then checks the summaries against the published runs.                        #
#                                                                              #
# Usage (from the Code folder):  Rscript reproduce_paper.R                     #
################################################################################

# ===== Settings
RERUN_SIMULATIONS <- FALSE # TRUE: rerun both experiments first (~1.5 h; see README)
results_folder <- "results/published_2018" # summaries of the published runs
figures_folder <- "figures"

source("R/plots.R")   # plot_global_metamodels(), plot_target_theta()
source("R/verify.R")  # verify_against_paper()

# ===== Optional rerun of the simulations (new random numbers)
if (RERUN_SIMULATIONS) {
    results_folder <- "results/rerun"
    rscript <- file.path(R.home("bin"), "Rscript")
    system2(rscript, c("scripts/run_global_metamodels.R", "--out", results_folder))
    system2(rscript, c("scripts/run_target_theta.R", "--out", results_folder))
}

# ===== Figures
dir.create(figures_folder, showWarnings = FALSE)
global_summary <- read.csv(file.path(results_folder, "global_metamodel_summary.csv"),
                           check.names = FALSE)
target_summary <- read.csv(file.path(results_folder, "target_theta_summary.csv"))

plot_global_metamodels(global_summary, figures_folder) # Figure 1 (a)-(h)
plot_target_theta(target_summary, figures_folder)      # Figure 2 (a)-(c)
cat("Figures 1 and 2 written to", figures_folder, "\n")

# ===== Verification
verify_against_paper(results_folder)
