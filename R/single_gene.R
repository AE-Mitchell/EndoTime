.single_gene <- function(unknown_sample, gene_name, expression_data, gene_window) {
    unknown_row <- which(expression_data[, "ID"] == unknown_sample)
    unknown_expression <- expression_data[unknown_row, gene_name]
    rest_rows <- which(expression_data[, "ID"] != unknown_sample)
    rest_LH <- expression_data[rest_rows, "LH_rn"][order(expression_data[rest_rows, "LH_rn"])]
    rest_gene <- expression_data[rest_rows, gene_name][order(expression_data[rest_rows, "LH_rn"])]
    medians <- c()
    densities <- c()

    start_position <- 1
    end_position <- gene_window
    while (end_position <= length(rest_LH)) {
        a <- .single_gene_sub(rest_LH, rest_gene, unknown_expression, start_position, end_position)
        medians <- c(medians, a[1])
        densities <- c(densities, a[2])
        start_position <- start_position + 1
        end_position <- end_position + 1
    }

    total <- sum(densities)
    densities <- densities / total

    return(cbind(medians, densities))
}
