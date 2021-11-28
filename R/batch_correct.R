.batch_correct <- function(expression_data, genes) {
    batch_means <- c()

    for (i in as.numeric(levels(factor(expression_data$BATCH)))) {
        subset_df <- expression_data[expression_data[["BATCH"]] == i, ]
        batch_means <- c(batch_means, mean(unlist(subset_df[, genes])))
    }

    correction <- c()
    correction[as.numeric(levels(factor(expression_data$BATCH)))] <- max(batch_means) - batch_means
    corrected <- expression_data[, genes] + correction[expression_data[, "BATCH"]]
    expression_data <- cbind(expression_data[, c("ID", "BATCH", "LH", "LH_rn")], corrected)
    rownames(expression_data) <- expression_data[, "ID"]

    expression_data
}
