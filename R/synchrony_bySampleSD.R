.synchrony_bySampleSD <- function(all_genes_densities, interval) {
    medians <- all_genes_densities[, "medians_LH"]
    genes <- colnames(all_genes_densities)[2:ncol(all_genes_densities)]

    maxima <- medians[apply(all_genes_densities[, genes], 2, which.max)]
    maxima.sd <- stats::sd(maxima)

    maxima.sd
}
