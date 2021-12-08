.scale_sample <- function(sample_data, genes, gene_scaling_info) {
    scaled_values <- c()

    for (i in genes) {
        result <- (sample_data[i] - gene_scaling_info[i, "mins"]) / (gene_scaling_info[i, "maxs"] - gene_scaling_info[i, "mins"])
        scaled_values <- unlist(c(scaled_values, result))
    }
    names(scaled_values) <- genes
    scaled_values
}
