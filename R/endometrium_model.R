.endometrium_model <- function(expression_data, genes, even_tps, gene_window, aggregate_window, iteration = 1, euc, spread = 33) {

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

    new_assignments_scaled <- .rescaling(new_assignments, min(expression_data[, "LH"]), max(expression_data[, "LH"]))
    euclidean_distance <- stats::dist(rbind(current_assignments, new_assignments_scaled))
    print(paste0("iteration", iteration, ": ", euclidean_distance))

    while (euclidean_distance > euc) {
        expression_data["LH_rn"] <- new_assignments_scaled
        iteration <- iteration + 1
        c <- .endometrium_model(expression_data = expression_data, genes = genes, even_tps = even_tps, gene_window = gene_window, aggregate_window = aggregate_window,
                                spread = spread, iteration = iteration, euc = euc)
        return(c)
    }
    names(new_assignments_scaled) <- samples
    return(new_assignments_scaled)
}
