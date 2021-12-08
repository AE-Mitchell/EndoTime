.extrapolate_aggregate <- function(gene_curves, genes, aggregate_window) {
    bin_radius <- aggregate_window / 2

    extrapolated_densities <- lapply(genes, function(gene_name) {
        gene_dat <- gene_curves[c("medians_LH", gene_name)]
        colnames(gene_dat) <- c("timing", "expression")
        bottom_gene_dat <- gene_dat[seq(2, (bin_radius + 1)), ]
        top_gene_dat <- gene_dat[seq((nrow(gene_curves) - bin_radius), (nrow(gene_curves) - 1)), ]

        bottom_extrap <- data.frame(
            "timing" = gene_dat[1, "timing"] + rev(gene_dat[1, "timing"] - bottom_gene_dat[["timing"]]),
            "expression" = gene_dat[1, "expression"] + rev(gene_dat[1, "expression"] - bottom_gene_dat[["expression"]])
        )

        top_extrap <- data.frame(
            "timing" = gene_dat[nrow(gene_curves), "timing"] + rev(gene_dat[nrow(gene_curves), "timing"] - top_gene_dat[["timing"]]),
            "expression" = gene_dat[nrow(gene_curves), "expression"] + rev(gene_dat[nrow(gene_curves), "expression"] - top_gene_dat[["expression"]])
        )

        extrapolated_gene_data <- rbind(bottom_extrap, gene_dat, top_extrap)
        colnames(extrapolated_gene_data) <- c("medians_LH", gene_name)
        extrapolated_gene_data
    })
    Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "medians_LH", all.x = TRUE), extrapolated_densities)
}
