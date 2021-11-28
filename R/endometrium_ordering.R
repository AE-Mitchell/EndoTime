.endometrium_ordering<- function(expression_data, genes, even_tps, iteration = 1, gene_window, aggregate_window, euc, spread = 33,
                                 decreasing_window) {

    samples <- expression_data[, "ID"]
    current_assignments <- expression_data[, "LH_rn"]
    new_assignments <- c()

    fake_points <- .extrapolation(dat = expression_data, genes = genes, size = gene_window, dp_per_LH = spread)
    extrapolated_data <- rbind(expression_data[, c("ID", "LH_rn", genes)], fake_points[, c("ID", "LH_rn", genes)])

    for (i in samples) {
        a <- .all_genes(unknown_sample = i, genes = genes, data_dCT = extrapolated_data, gene_window = gene_window)
        b <- .all_genes_aggregate(unknown_sample = i, all_genes_result = a, aggregate_window = aggregate_window)
        new_assignments <- c(new_assignments, b)
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
                                   aggregate_window = aggregate_window, euc = euc, spread = spread, decreasing_window = decreasing_window)
        return(c)
    }
    names(new_assignments) <- samples
    return(cbind(new_assignments, raw_assignments))
}
