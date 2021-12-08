.all_genes <- function(unknown_sample, genes, expression_data, gene_window) {

    medians_LH <- .single_gene(unknown_sample, genes[1], expression_data, gene_window)[, 1]

    densities <- sapply(genes, function(x) .single_gene(unknown_sample, x, expression_data, gene_window)[, 2])

    cbind(medians_LH, densities)
}
