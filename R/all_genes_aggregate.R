.all_genes_aggregate <- function(unknown_sample, all_genes_result, aggregate_window) {
    average_densities <- c()
    median_LH <- c()

    start_position <- 1
    end_position <- aggregate_window

    while (end_position <= nrow(all_genes_result)){
        current_subset <- all_genes_result[start_position:end_position,]
        median_LH <- c(median_LH, stats::median(current_subset[, 1]))
        average_densities <- c(average_densities, mean(current_subset[, 2:ncol(current_subset)]))
        start_position <- start_position + 1
        end_position <- end_position + 1
    }

    return(median_LH[which.max(average_densities)])
}
