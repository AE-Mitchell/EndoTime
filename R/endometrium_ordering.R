.endometrium_ordering<- function(expression_data, genes, even_tps, iteration = 1, gene_window, aggregate_window, euc, spread = 33,
                                 decreasing_window, new_extrap) {

    samples <- expression_data[, "ID"]
    current_assignments <- expression_data[, "LH_rn"]
    new_assignments <- c()

    if (new_extrap) {
        fake_points <- .extrapolate_data(expression_data, genes, gene_window, iteration)
    } else {
        fake_points <- .extrapolation(dat = expression_data, genes = genes, size = gene_window, dp_per_LH = spread)
    }

    extrapolated_data <- rbind(expression_data[, c("ID", "LH_rn", genes)], fake_points[, c("ID", "LH_rn", genes)])

    new_assignments <- c()

    for (i in samples) {
        a <- .all_genes(unknown_sample = i, genes = genes, data_dCT = extrapolated_data, gene_window = gene_window)

        a <- as.data.frame(a)

        # NB: This process extends expression into the negative values. However, as we're taking the average per bin,
        # it should balance out
        extrapolated_densities <- .extrapolate_aggregate(a, genes, aggregate_window)

        density_bins <- .aggregate_calc(extrapolated_densities, aggregate_window, genes)

        new_assignments <- c(new_assignments, density_bins[which.max(density_bins[["mean_exp"]]), "median_lh"])
    }

    raw_assignments <- new_assignments
    new_assignments[order(new_assignments)] <- even_tps
    euclidean_distance <- stats::dist(rbind(current_assignments, new_assignments))
    print(paste0(gene_window, ", ", aggregate_window))
    print(paste("iteration", iteration, ":", euclidean_distance, sep = " "))

    while (euclidean_distance > euc) {
        expression_data[, "LH_rn"] <- new_assignments
        iteration <- iteration + 1

        if (decreasing_window) {
            gene_window <- ifelse((gene_window - 10) >= 20, gene_window - 10, 20)
            aggregate_window <- ifelse((aggregate_window - 10) >= 20, aggregate_window - 10, 20)
        }

        c <- .endometrium_ordering(expression_data = expression_data, genes = genes, iteration = iteration,
                                   gene_window = gene_window, even_tps = even_tps,
                                   aggregate_window = aggregate_window, euc = euc, spread = spread, decreasing_window = decreasing_window,
                                   new_extrap = new_extrap)
        return(c)
    }
    names(new_assignments) <- samples
    return(cbind(new_assignments, raw_assignments))
}
