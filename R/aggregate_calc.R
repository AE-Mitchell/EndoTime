.aggregate_calc <- function(extrapolated_densities, aggregate_window, panel_genes) {
    density_bins <- lapply(seq(1, nrow(extrapolated_densities) - aggregate_window), function(bin_num) {
        bin_dat <- extrapolated_densities[seq(bin_num, bin_num + aggregate_window - 1), ]
        bin_lh <- stats::median(bin_dat[["medians_LH"]])
        bin_exps <- bin_dat[panel_genes]
        mean_exp <- mean(unlist(bin_exps))
        data.frame("median_lh" = bin_lh, "mean_exp" = mean_exp)
    })
    do.call(rbind, density_bins)
}
