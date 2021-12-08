.extrapolate_data <- function(expression_data, genes, size) {
    width <- size / 2

    # Extrapolate times
    time_df <- data.frame("rank" = rank(expression_data[["LH_rn"]]), "LH_rn" = expression_data[["LH_rn"]])
    time_lh <- stats::lm(LH_rn ~ rank, time_df)

    predict_time_df <- data.frame("rank" = c(seq(1 - width, 0), seq(nrow(expression_data) + 1, nrow(expression_data) + width)))
    predict_time_df["LH_rn"] <- stats::predict(time_lh, newdata = predict_time_df)

    extrapolated_expressions <- sapply(genes, function(gene_name) {
        expressions_df <- expression_data[c("LH_rn", gene_name)]
        colnames(expressions_df) <- c("timing", "expression")

        expression_lm <- stats::lm(expression ~ timing, expressions_df)

        predict_expressions_df <- data.frame("timing" = predict_time_df[["LH_rn"]])
        predict_expressions_df["expression"] <- stats::predict(expression_lm, newdata = predict_expressions_df)
        colnames(predict_expressions_df) <- c("LH_rn", gene_name)

        predict_expressions_df

    }, USE.NAMES = TRUE, simplify = FALSE)

    extrapolated_expressions <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "LH_rn", all.x = TRUE), extrapolated_expressions)

    extrapolated_expressions["ID"] <- paste0("extra_", seq(1, size))

    extrapolated_expressions
}
