.synchrony_bySampleSD <- function(all_genes_densities, interval) {
    medians <- all_genes_densities[, "medians_LH"]
    genes <- colnames(all_genes_densities)[2:ncol(all_genes_densities)]
    scores <- apply(all_genes_densities[, genes], 2, function(x) .synchrony_score_sub(x, medians, interval))
    scores <- sort(scores, decreasing = TRUE)


    maxima <- medians[apply(all_genes_densities[, genes], 2, which.max)]
    maxima.sd <- stats::sd(maxima)

    return(maxima.sd)
}
