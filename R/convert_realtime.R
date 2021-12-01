.convert_realtime <- function(result) {
    # Extract reported LH values
    reported_lh <- result[["expression_data"]][c("ID", "LH")]
    reported_lh["LH"] <- reported_lh[["LH"]] + stats::runif(257, -0.5, 0.5)

    # Extract pseudotimes
    pseudo_lh <- as.data.frame(result[["model_results"]])
    pseudo_lh["ID"] <- rownames(pseudo_lh)
    rownames(pseudo_lh) <- NULL
    pseudo_lh <- pseudo_lh[c("ID", "new_assignments")]
    colnames(pseudo_lh) <- c("ID", "pLH")

    # Arrange samples by pseudotime and pair with reported LHs but ignoring sample IDs
    timing_df <- pseudo_lh[order(pseudo_lh[["pLH"]]), ]
    timing_df["rLH"] <- sort(reported_lh[["LH"]])

    # ----------------------------------------------------------------------
    # Extrapolate data to allow for binning
    bin_radius <- 10 # Half of the smallest bin size used

    # Pseudotimes are already monotonic, so can be extrapolated via linear model
    pseudo_df <- data.frame("rank" = rank(timing_df[["pLH"]]), "pLH" = timing_df[["pLH"]])
    pseudo_lm <- stats::lm(pLH ~ rank, pseudo_df)
    predicted_pseudo <- data.frame("rank" = c(seq(-9, 0), seq(258, 267)))
    predicted_pseudo["pLH"] <- stats::predict(pseudo_lm, newdata = predicted_pseudo)

    # Extrapolate bottom end of rLH
    bottom_rlh <- data.frame("rLH" = timing_df[seq(1, 10), "rLH"])
    bottom_rlh["rank"] <- seq(1, 10)
    bottom_lm <- stats::lm(rLH ~ rank, bottom_rlh)
    predicted_bottom <- data.frame("rank" = seq(-9, 0))
    predicted_bottom["rLH"] <- stats::predict(bottom_lm, newdata = predicted_bottom)

    # Extrapolate top end of rLH
    top_rlh <- data.frame("rLH" = timing_df[seq(248, 257), "rLH"])
    top_rlh["rank"] <- seq(248, 257)
    top_lm <- stats::lm(rLH ~ rank, top_rlh)
    predicted_top <- data.frame("rank" = seq(258, 267))
    predicted_top["rLH"] <- stats::predict(top_lm, newdata = predicted_top)

    # Merge extrapolations
    predicted_rlh <- rbind(predicted_bottom, predicted_top)
    predicted <- merge(predicted_pseudo, predicted_rlh, on = "rank")
    predicted["rank"] <- NULL
    predicted["ID"] <- paste0("extra_", seq(1, 20))

    # Join with existing timing data
    extrapolated_data <- rbind(timing_df, predicted)
    extrapolated_data <- extrapolated_data[order(extrapolated_data[["pLH"]]), ]

    # Bin extrapolated data
    binned_timing <- lapply(seq(1, nrow(timing_df)), function(bin_num) {
        dat <- extrapolated_data[seq(bin_num, bin_num + 20), ]
        sample_id <- dat[11, "ID"]

        # Calculate weighted mean of rLH values within bin
        half_weights <- c(seq(0, 10 - 1) * (1 / 10))
        weights <- c(half_weights, 1, rev(half_weights))
        weighted_mean <- stats::weighted.mean(dat[["rLH"]], weights)

        data.frame("ID" = sample_id, "eLH" = weighted_mean)
    })
    binned_timing <- do.call(rbind, binned_timing)

    pseudo_results <- as.data.frame(result[["model_results"]])
    pseudo_results["ID"] <- rownames(pseudo_results)
    pseudo_results <- pseudo_results[c("ID", "new_assignments")]
    colnames(pseudo_results) <- c("ID", "pLH")

    timing_results <- merge(pseudo_results, binned_timing, on = "ID")

    result[["model_results"]] <- timing_results

    result
}
