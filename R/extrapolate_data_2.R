.extrapolate_data_2 <- function(dat, genes, size, iteration) {
    # Extrapolate times
    time_df <- data.frame("rank" = rank(dat[["LH_rn"]]), "LH_rn" = dat[["LH_rn"]])
    time_lh <- stats::lm(LH_rn ~ rank, time_df)

    predict_time_df <- data.frame("rank" = c(seq(-39, 0), seq(258, 297)))
    predict_time_df["LH_rn"] <- stats::predict(time_lh, newdata = predict_time_df)

    extrapolated_expressions <- sapply(genes, function(gene_name) {
        expressions_df <- dat[c("LH_rn", gene_name)]
        colnames(expressions_df) <- c("timing", "expression")

        expression_lm <- stats::lm(expression ~ poly(timing, 2), expressions_df)

        predict_expressions_df <- data.frame("timing" = predict_time_df[["LH_rn"]])
        predict_expressions_df["expression"] <- stats::predict(expression_lm, newdata = predict_expressions_df)
        colnames(predict_expressions_df) <- c("LH_rn", gene_name)

        predict_expressions_df

    }, USE.NAMES = TRUE, simplify = FALSE)

    extrapolated_expressions <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "LH_rn", all.x = TRUE), extrapolated_expressions)

    plot_expressions <- pivot_longer(extrapolated_expressions, cols = all_of(genes), names_to = "gene", values_to = "expression")
    plot_dat <- pivot_longer(dat, cols = all_of(genes), names_to = "gene", values_to = "expression")

    pdf(paste0("~/Documents/espam/analysis/paper/extrapolate_plots/decreasing_extrap_iteration_", iteration, ".pdf"), width = 3, height = 10)
        print(
            ggplot(plot_dat, aes(LH_rn, expression)) +
                geom_point(color = "orange", size = 0.6) +
                geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
                geom_point(color = "red", size = 0.6, data = plot_expressions) +
                facet_wrap(~ gene, ncol = 1, scales = "free") +
                theme_minimal() +
                theme(
                    panel.grid = element_blank(),
                    axis.line = element_line(color = "black")
                )
        )
    dev.off()

    extrapolated_expressions["ID"] <- paste0("extra_", seq(1, 80))

    extrapolated_expressions
}
