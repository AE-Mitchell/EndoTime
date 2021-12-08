.endometrium_ordering<- function(expression_data, genes, even_tps, iteration = 1, gene_window, aggregate_window, euc) {

    samples <- expression_data[, "ID"]
    current_assignments <- expression_data[, "LH_rn"]
    new_assignments <- c()

    fake_points <- .extrapolate_data(expression_data, genes, gene_window)

    extrapolated_data <- rbind(expression_data[, c("ID", "LH_rn", genes)], fake_points[, c("ID", "LH_rn", genes)])

    new_assignments <- c()

    for (unknown_sample in samples) {
        gene_curves <- .all_genes(unknown_sample, genes, extrapolated_data, gene_window)
        gene_curves <- as.data.frame(gene_curves)

        # NB: This process extends expression into the negative values. However, as we're taking the average per bin,
        # it should balance out
        extrapolated_densities <- .extrapolate_aggregate(gene_curves, genes, aggregate_window)
        density_bins <- .aggregate_calc(extrapolated_densities, aggregate_window, genes)

        new_assignments <- c(new_assignments, density_bins[which.max(density_bins[["mean_exp"]]), "median_lh"])
    }

    raw_assignments <- new_assignments
    new_assignments[order(new_assignments)] <- even_tps
    euclidean_distance <- stats::dist(rbind(current_assignments, new_assignments))
    print(paste("iteration", iteration, ":", euclidean_distance, sep = " "))

    while (euclidean_distance > euc) {
        expression_data[, "LH_rn"] <- new_assignments
        iteration <- iteration + 1

        gene_window <- ifelse((gene_window - 10) >= 20, gene_window - 10, 20)
        aggregate_window <- ifelse((aggregate_window - 10) >= 20, aggregate_window - 10, 20)

        return(.endometrium_ordering(expression_data, genes, even_tps, iteration, gene_window, aggregate_window, euc))
    }
    names(new_assignments) <- samples

    return(cbind(new_assignments, raw_assignments))
}
