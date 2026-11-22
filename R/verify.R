################################################################################
# verify.R                                                                     #
#                                                                              #
# verify_against_paper(): checks a set of summaries against the paper.         #
# The paper reports its results only as figures, so two kinds of checks:      #
#   (1) the statements that Sections 4.1 and 4.2 make about the figures       #
#       (ranges, regions where an estimator beats SMC, orderings of curves);  #
#   (2) for a rerun: every plotted value against the published 2018 run,      #
#       within Monte Carlo error (z-scores).                                   #
################################################################################

check_statement <- function(report, description, passed) {
    # Append one checked statement to the report
    # Input(s):
    #   -- report: data frame of checks so far
    #   -- description: the statement, as in the paper
    #   -- passed: TRUE/FALSE
    # Output(s): the report with one more row
    row <- data.frame(check = description, status = ifelse(passed, "ok", "UNEXPECTED"))
    return(rbind(report, row))
}

check_paper_statements <- function(global, target) {
    # Statements of Sections 4.1 and 4.2 about Figures 1 and 2
    # Input(s):
    #   -- global: global_metamodel_summary.csv (Figure 1)
    #   -- target: target_theta_summary.csv (Figure 2)
    # Output(s): data frame (check, status)

    report <- data.frame(check = character(0), status = character(0))
    var_smc <- global$mse_SMC
    ratio_blr <- global$mse_BLR / var_smc       # relative MSEs (Figure 1c-h)
    ratio_slr <- global$mse_SLR / var_smc
    ratio_rlr <- global$mse_RLR / var_smc
    ratio_blr_two <- global[["mse_SMC+BLR"]] / var_smc
    ratio_slr_two <- global[["mse_SMC+SLR"]] / var_smc
    ratio_rlr_two <- global[["mse_SMC+RLR"]] / var_smc

    #-- Section 4.1, Figure 1
    report <- check_statement(report,
        "Fig 1b: Var(SMC) is largest where the response surface is largest",
        which.max(var_smc) == which.max(global$true_value))
    report <- check_statement(report,
        "Fig 1c: BLR relative MSE exceeds 100 somewhere (truncated at 100)",
        max(ratio_blr) > 100)
    report <- check_statement(report,
        "Fig 1c-e: ranges ordered RLR < SLR < BLR",
        max(ratio_rlr) < max(ratio_slr) && max(ratio_slr) < max(ratio_blr))
    report <- check_statement(report,
        "Fig 1e: RLR beats SMC (ratio < 1) in a large area (>= 40% of the grid)",
        mean(ratio_rlr < 1) >= 0.40)
    low_rate <- global$theta1 <= 0.8 | global$theta5 <= 0.8
    report <- check_statement(report,
        "Fig 1c-e: all three worst where a rate is low (theta1 or theta5 <= 0.8)",
        low_rate[which.max(ratio_blr)] && low_rate[which.max(ratio_slr)] &&
            low_rate[which.max(ratio_rlr)])
    center <- which.min(abs(global$theta1 - 1) + abs(global$theta5 - 1))
    report <- check_statement(report,
        "Fig 1e: RLR relative MSE somewhere lower than next to theta* = (1,1)",
        min(ratio_rlr) < ratio_rlr[center])
    report <- check_statement(report,
        "Fig 1f-h: combining with SMC lowers each LR metamodel's relative MSE",
        max(ratio_blr_two) < max(ratio_blr) && max(ratio_slr_two) < max(ratio_slr) &&
            max(ratio_rlr_two) < max(ratio_rlr))
    report <- check_statement(report,
        "Fig 1h: SMC + RLR beats SMC at every grid point",
        all(ratio_rlr_two < 1))
    both_above <- function(rate) global$theta1 > rate & global$theta5 > rate
    report <- check_statement(report,
        "Fig 1g: SMC + SLR beats SMC when both rates are at least 0.75",
        all(ratio_slr_two[both_above(0.75)] < 1))
    report <- check_statement(report,
        "Fig 1f: SMC + BLR beats SMC when both rates > 0.9",
        all(ratio_blr_two[both_above(0.9)] < 1))

    #-- Section 4.2, Figure 2 (iterations t = 2, ..., 100; at t = 1 all are equal)
    mse_at <- function(name) target$mse[target$estimator == name & target$step >= 2]
    report <- check_statement(report,
        "Fig 2: true value mu(theta0) = 83/24 (closed form for theta0 = 1)",
        abs(target$true_value[1] - 83 / 24) < 1e-12)
    report <- check_statement(report,
        "Fig 2a: self-normalization helps (PW: SLR <= BLR at every t)",
        all(mse_at("PW.SLR") <= mse_at("PW.BLR")))
    report <- check_statement(report,
        "Fig 2c: self-normalization helps (M: SLR <= BLR at every t)",
        all(mse_at("SM.SLR") <= mse_at("SM.BLR")))
    report <- check_statement(report,
        "Fig 2b: OW-RLR worse than PW-RLR and M-RLR at t = 100",
        tail(mse_at("OW.RLR"), 1) > tail(mse_at("PW.RLR"), 1) &&
            tail(mse_at("OW.RLR"), 1) > tail(mse_at("SM.RLR"), 1))
    brute <- mse_at("Brute")
    m_based_close <- TRUE
    for (name in c("SM.BLR", "SM.SLR", "SM.RLR")) {
        relative_gap <- mse_at(name) / brute
        m_based_close <- m_based_close && all(relative_gap < 1.10)
    }
    report <- check_statement(report,
        "Fig 2c: M-based estimators almost as accurate as Brute (within 10%, all t)",
        m_based_close)
    report <- check_statement(report,
        "Fig 2: Crude is the least accurate estimator at t = 100",
        tail(mse_at("Crude"), 1) > max(target$mse[target$step == 100 &
                                                  target$estimator != "Crude"]))
    return(report)
}

compare_with_published <- function(new, published, columns, label) {
    # z-scores of new MSEs against the published 2018 run
    # Input(s):
    #   -- new, published: summaries with columns mse*, se_mse*, top1_share*
    #   -- columns: estimator suffixes to compare (e.g. "_BLR"), or "" for long format
    #   -- label: name of each compared value (vector, one per row of new)
    # Output(s): data frame (value, new, published, z, status)
    # formula: z = (x - v) / sqrt(se_x^2 + se_v^2); heavy-tailed if the worst 1% of
    #          macro runs make up more than half of either MSE (SE unreliable)

    rows <- NULL
    for (suffix in columns) {
        x <- new[[paste0("mse", suffix)]]
        v <- published[[paste0("mse", suffix)]]
        se_x <- new[[paste0("se_mse", suffix)]]
        se_v <- published[[paste0("se_mse", suffix)]]
        z <- (x - v) / sqrt(se_x^2 + se_v^2)
        heavy <- new[[paste0("top1_share", suffix)]] > 0.5 |
            published[[paste0("top1_share", suffix)]] > 0.5
        status <- ifelse(abs(z) <= 2.5, "ok", ifelse(heavy, "heavy-tailed", "flagged"))
        rows <- rbind(rows, data.frame(value = paste0(label, suffix), new = x,
                                       published = v, z = z, status = status))
    }
    return(rows)
}

verify_against_paper <- function(results_folder,
                                 published_folder = "results/published_2018") {
    # Run all checks on the summaries in results_folder and print a report
    # Input(s):
    #   -- results_folder: folder with the two summaries to check
    #   -- published_folder: summaries of the published 2018 runs
    # Output(s): (invisibly) list of the two reports; also written to
    #            <results_folder>/verification_report.csv

    global <- read.csv(file.path(results_folder, "global_metamodel_summary.csv"),
                       check.names = FALSE)
    target <- read.csv(file.path(results_folder, "target_theta_summary.csv"))

    statements <- check_paper_statements(global, target)
    cat("\n===== Statements of the paper about Figures 1 and 2\n")
    print(statements, right = FALSE, row.names = FALSE)

    comparison <- NULL
    if (normalizePath(results_folder) != normalizePath(published_folder)) {
        global_2018 <- read.csv(file.path(published_folder, "global_metamodel_summary.csv"),
                                check.names = FALSE)
        target_2018 <- read.csv(file.path(published_folder, "target_theta_summary.csv"))
        stopifnot(all(global$true_value == global_2018$true_value)) # same true values

        grid_label <- sprintf("Fig1 (%.3f,%.3f) ", global$theta1, global$theta5)
        global_cmp <- compare_with_published(global, global_2018,
            c("_SMC", "_BLR", "_SLR", "_RLR", "_SMC+BLR", "_SMC+SLR", "_SMC+RLR"),
            grid_label)
        #-- Figure 2: iterations t >= 2 (at t = 1 all estimators are the crude one)
        keep <- target$step >= 2 & target$estimator != "CumAvg" # CumAvg not plotted
        keep_2018 <- target_2018$step >= 2 & target_2018$estimator != "CumAvg"
        target_new <- target[keep, ]
        target_old <- target_2018[keep_2018, ]
        stopifnot(all(target_new$step == target_old$step))
        target_cmp <- compare_with_published(target_new, target_old, "",
            sprintf("Fig2 %s t=%d", target_new$estimator, target_new$step))
        comparison <- rbind(global_cmp, target_cmp)

        n_checks <- nrow(comparison)
        cat("\n===== Rerun vs. published 2018 run (Monte Carlo error)\n")
        cat(sprintf("%d values compared; expected chance flags (|z| > 2.5): %.1f\n",
                    n_checks, 0.0124 * n_checks))
        print(table(comparison$status))
        flagged <- comparison[comparison$status != "ok", ]
        if (nrow(flagged) > 0) {
            print(head(flagged[order(-abs(flagged$z)), ], 20), row.names = FALSE)
        }
        write.csv(comparison, file.path(results_folder, "comparison_with_2018.csv"),
                  row.names = FALSE)
    }

    write.csv(statements, file.path(results_folder, "verification_report.csv"),
              row.names = FALSE)
    n_unexpected <- sum(statements$status == "UNEXPECTED")
    if (!is.null(comparison)) {
        #-- single flags are expected by chance (1.24% of values with |z| > 2.5);
        #   the rerun disagrees only if there are clearly too many flags (more
        #   than twice the expected number) or a flag with |z| > 4.5
        n_flagged <- sum(comparison$status == "flagged")
        n_expected <- 0.0124 * nrow(comparison)
        largest_z <- max(abs(comparison$z[comparison$status != "heavy-tailed"]))
        cat(sprintf("Flagged: %d (expected by chance: %.1f); largest |z| (not heavy-tailed): %.2f\n",
                    n_flagged, n_expected, largest_z))
        if (n_flagged > 2 * n_expected || largest_z > 4.5) {
            n_unexpected <- n_unexpected + 1
            cat("UNEXPECTED: the rerun does not agree with the 2018 run within Monte Carlo error\n")
        }
    }
    cat(sprintf("\nVerification finished: %d UNEXPECTED\n", n_unexpected))
    return(invisible(list(statements = statements, comparison = comparison)))
}
