.all_genes <- function(unknown_sample, genes, data_dCT, gene_window) {

    medians_LH <- .single_gene(unknown_sample, genes[1], data_dCT, gene_window)[, 1]

    densities <- sapply(genes, function(x) .single_gene(unknown_sample, x, data_dCT, gene_window)[, 2])

    result <- cbind(medians_LH, densities)

    return(result)
}
