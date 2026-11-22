################################################################################
# plots.R                                                                      #
#                                                                              #
# Figures 1 and 2 of the paper, drawn from the summaries in results/.          #
# Figure 1: lattice wireframes over the (theta1, theta5) grid (Section 4.1).   #
# Figure 2: MSE vs. iteration t on log-log axes, ggplot2 (Section 4.2).        #
# Each figure is written as a vector PDF and as a 2018-style JPEG with the     #
# file names used by the paper's LaTeX source (Figures/<name>.jpeg).           #
################################################################################

library(lattice) # wireframe() for the response surfaces (Figure 1)
library(ggplot2) # MSE curves (Figure 2)

# ===== Figure 1: global metamodels

draw_surface <- function(values, summary, title, colors, z_distance, file_stem,
                         out_folder) {
    # Draw one surface of Figure 1 as <file_stem>.pdf and <file_stem>.jpeg
    # Input(s):
    #   -- values: 400 values to plot, in the row order of summary
    #   -- summary: data frame with columns theta1, theta5 (the 20 x 20 grid)
    #   -- title: panel title
    #   -- colors: color palette for the drape (200 colors)
    #   -- z_distance: distance of the z-axis labels (5 or 4, as in 2018)
    #   -- file_stem: file name without extension, e.g. "SANGlobalRLR"
    #   -- out_folder: folder for the figure files
    # Output(s): none (writes two files)

    surface <- data.frame(theta1 = summary$theta1, theta5 = summary$theta5,
                          value = values)
    plot <- wireframe(value ~ theta1 * theta5, data = surface, drape = TRUE,
                      aspect = c(1, 1),
                      screen = list(z = -60, y = 0, x = -60),
                      colorkey = list(labels = list(cex = 2)),
                      scales = list(arrows = FALSE, cex = 2, tick.number = 4,
                                    z = list(arrows = FALSE),
                                    distance = c(1.5, 1.5, z_distance)),
                      col.regions = colors,
                      xlab = list(expression(theta[1]), cex = 2.5),
                      ylab = list(expression(theta[5]), cex = 2.5),
                      zlab = list(""),
                      main = list(title, cex = 2))

    pdf(file.path(out_folder, paste0(file_stem, ".pdf")), width = 500 / 72,
        height = 500 / 72) # same layout as the 500 x 500 pixel JPEG
    print(plot)
    dev.off()
    jpeg(file.path(out_folder, paste0(file_stem, ".jpeg")), width = 500,
         height = 500, units = "px")
    print(plot)
    dev.off()
    return(invisible(NULL))
}

plot_global_metamodels <- function(summary, out_folder) {
    # The eight panels of Figure 1
    # Input(s):
    #   -- summary: results/global_metamodel_summary.csv (one row per grid point)
    #   -- out_folder: folder for the figure files
    # Output(s): none (writes 16 files)
    # formula: relative MSE = MSE(estimator) / Var(SMC), where Var(SMC) is the
    #          MSE of the SMC estimator about the true values (10^7 reps each)

    n_samp <- 100 # replications per SMC estimate (title of panel b)
    var_smc <- summary$mse_SMC

    draw_surface(summary$true_value, summary, "True Response Surface",
                 topo.colors(200), 5, "SANGlobalTrue", out_folder)
    draw_surface(var_smc, summary,
                 sprintf("Var(SMC) of %i independent reps.", n_samp),
                 topo.colors(200), 5, "SANGlobalSampVar", out_folder)

    #-- panels c-e: LR metamodels alone; BLR is truncated at 100 for display
    ratio_blr <- pmin(summary$mse_BLR / var_smc, 100)
    ratio_slr <- summary$mse_SLR / var_smc
    ratio_rlr <- summary$mse_RLR / var_smc
    draw_surface(ratio_blr, summary, "MSE of BLR relative to Var(SMC)",
                 rainbow(200), 5, "SANGlobalBLR", out_folder)
    draw_surface(ratio_slr, summary, "MSE of SLR relative to Var(SMC)",
                 rainbow(200), 5, "SANGlobalSLR", out_folder)
    draw_surface(ratio_rlr, summary, "MSE of RLR relative to Var(SMC)",
                 rainbow(200), 5, "SANGlobalRLR", out_folder)

    #-- panels f-h: average of SMC and LR estimates; SMC + BLR is truncated at
    #   10 for display (as in the 2018 figure)
    ratio_blr_two <- pmin(summary[["mse_SMC+BLR"]] / var_smc, 10)
    ratio_slr_two <- summary[["mse_SMC+SLR"]] / var_smc
    ratio_rlr_two <- summary[["mse_SMC+RLR"]] / var_smc
    draw_surface(ratio_blr_two, summary, "MSE of SMC + BLR relative to Var(SMC)",
                 rainbow(200), 4, "SANGlobalBLRTwo", out_folder)
    draw_surface(ratio_slr_two, summary, "MSE of SMC + SLR relative to Var(SMC)",
                 rainbow(200), 4, "SANGlobalSLRTwo", out_folder)
    draw_surface(ratio_rlr_two, summary, "MSE of SMC + RLR relative to Var(SMC)",
                 rainbow(200), 4, "SANGlobalRLRTwo", out_folder)
    return(invisible(NULL))
}

# ===== Figure 2: target theta0

draw_mse_curves <- function(summary, estimators, labels, y_limits, file_stem,
                            out_folder) {
    # Draw one panel of Figure 2: MSE_t vs. t (log-log) for five estimators
    # Input(s):
    #   -- summary: results/target_theta_summary.csv (one row per step and estimator)
    #   -- estimators: names in the summary, in legend order
    #   -- labels: legend labels, in the same order
    #   -- y_limits: limits of the y-axis, or NULL for automatic limits
    #   -- file_stem: file name without extension, e.g. "SANEst_PW"
    #   -- out_folder: folder for the figure files
    # Output(s): none (writes two files)

    keep <- summary$estimator %in% estimators
    curves <- summary[keep, c("step", "estimator", "mse")]
    curves$estimator <- factor(curves$estimator, levels = estimators,
                               labels = labels) # fixes the legend order

    plot <- ggplot(curves, aes(x = step, y = mse, colour = estimator,
                               linetype = estimator, shape = estimator)) +
        geom_point(size = 5) +
        geom_line(linewidth = 2) +
        scale_x_continuous(name = "Number of Iterations", trans = "log10",
                           breaks = c(1, 10, 100), labels = c("1", "10", "100")) +
        scale_y_continuous(trans = "log10", limits = y_limits,
                           breaks = c(0.01, 0.1), labels = c("0.01", "0.10")) +
        ggtitle("Mean Squared Errors") +
        theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              axis.text = element_text(size = 16),
              axis.title = element_text(size = 16, face = "bold"),
              axis.title.y = element_blank(),
              panel.background = element_rect(fill = "grey95"),
              axis.line = element_line(linewidth = 1, colour = "black"),
              legend.justification = c(-0.1, -0.1),
              legend.position = "inside",
              legend.position.inside = c(0, 0),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              legend.background = element_rect(colour = "black"),
              legend.key.width = unit(2, "cm")) +
        labs(colour = "Estimators", linetype = "Estimators", shape = "Estimators")

    #-- points outside y_limits are dropped, as in 2018; suppress the warning
    suppressWarnings(ggsave(file.path(out_folder, paste0(file_stem, ".pdf")),
                            plot, width = 6, height = 5, units = "in"))
    suppressWarnings(ggsave(file.path(out_folder, paste0(file_stem, ".jpeg")),
                            plot, device = "jpeg", width = 6, height = 5,
                            units = "in"))
    return(invisible(NULL))
}

plot_target_theta <- function(summary, out_folder) {
    # The three panels of Figure 2
    # Input(s):
    #   -- summary: results/target_theta_summary.csv
    #   -- out_folder: folder for the figure files
    # Output(s): none (writes 6 files)

    y_limits <- c(0.002, 0.5) # panels a and c (as in 2018); panel b automatic
    draw_mse_curves(summary, c("Crude", "Brute", "PW.BLR", "PW.SLR", "PW.RLR"),
                    c("Crude", "Brute", "BLR-PW", "SLR-PW", "RLR-PW"),
                    y_limits, "SANEst_PW", out_folder)
    draw_mse_curves(summary, c("Crude", "Brute", "PW.RLR", "OW.RLR", "SM.RLR"),
                    c("Crude", "Brute", "RLR-PW", "RLR-OW", "RLR-M"),
                    NULL, "SANEst_RLR", out_folder)
    draw_mse_curves(summary, c("Crude", "Brute", "SM.BLR", "SM.SLR", "SM.RLR"),
                    c("Crude", "Brute", "BLR-M", "SLR-M", "RLR-M"),
                    y_limits, "SANEst_SM", out_folder)
    return(invisible(NULL))
}
