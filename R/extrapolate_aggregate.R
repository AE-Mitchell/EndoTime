.extrapolate_aggregate <- function(a, genes, aggregate_window) {
    bin_radius <- aggregate_window / 2

    extrapolated_densities <- lapply(genes, function(gene_name) {
        gene_dat <- a[c("medians_LH", gene_name)]
        colnames(gene_dat) <- c("timing", "expression")
        bottom_gene_dat <- gene_dat[seq(2, (bin_radius + 1)), ]
        top_gene_dat <- gene_dat[seq((nrow(a) - bin_radius), (nrow(a) - 1)), ]

        bottom_extrap <- data.frame(
            "timing" = gene_dat[1, "timing"] + rev(gene_dat[1, "timing"] - bottom_gene_dat[["timing"]]),
            "expression" = gene_dat[1, "expression"] + rev(gene_dat[1, "expression"] - bottom_gene_dat[["expression"]])
        )

        top_extrap <- data.frame(
            "timing" = gene_dat[nrow(a), "timing"] + rev(gene_dat[nrow(a), "timing"] - top_gene_dat[["timing"]]),
            "expression" = gene_dat[nrow(a), "expression"] + rev(gene_dat[nrow(a), "expression"] - top_gene_dat[["expression"]])
        )

        extrapolated_gene_data <- rbind(bottom_extrap, gene_dat, top_extrap)
        colnames(extrapolated_gene_data) <- c("medians_LH", gene_name)
        extrapolated_gene_data
    })
    Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "medians_LH", all.x = TRUE), extrapolated_densities)
}
