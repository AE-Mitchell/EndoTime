#' Title
#'
#' @param expression_data ...
#' @param genes ...
#' @param batch_correct ...
#' @param check_asynchrony ...
#' @param gene_window ...
#' @param aggregate_window ...
#' @param euc ...
#'
#' @return ...
#' @export
endotime_train <- function(expression_data, genes, batch_correct = TRUE, check_asynchrony = TRUE, gene_window = 80, aggregate_window = 80, euc = 2) {
    set.seed(101)

    result <- list()

    # 1) Read in the data, add random noise to LH+ values and invert delta-CT values to reflect expression
    expression_data <- as.data.frame(expression_data)
    expression_data[, "LH_rn"] <- expression_data[, "LH"] + stats::runif(nrow(expression_data), -0.5, 0.5)
    expression_data[, "ID"] <- as.vector(expression_data[, "ID"])
    expression_data <- cbind(expression_data[, c("ID", "BATCH", "LH", "LH_rn")], expression_data[, genes] * -1)

    # 2) Scale genes
    gene_scaling_info <- cbind(apply(expression_data[, genes], 2, min), apply(expression_data[, genes], 2, max))
    colnames(gene_scaling_info) <- c("mins", "maxs")
    result[["gene_scaling"]] <- gene_scaling_info
    expression_data <- cbind(expression_data[, 1:4], t(sapply(1:dim(expression_data)[1], function(x) .scale_sample(expression_data[x, ], genes, gene_scaling_info))))

    if (batch_correct) {
        # 3) Correct for batch effect
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

        result[["batch_correction"]] <- max(batch_means)
    }

    # 3) run the training model
    even_tps <- seq(min(expression_data[, "LH"]), max(expression_data[, "LH"]), length.out = nrow(expression_data))
    training_model_results <- .endometrium_ordering(expression_data = expression_data, even_tps = even_tps,
                                                   genes = genes, gene_window = gene_window, aggregate_window = aggregate_window,
                                                   spread = 11, euc = 2)

    result[["expression_data"]] <- expression_data
    result[["model_results"]] <- training_model_results

    if (check_asynchrony) {
        expression_data["LH_rn"] <- training_model_results[, 1]
        fake_points <- .extrapolation(dat = expression_data, genes = genes, size = gene_window, dp_per_LH = 11)
        data <- rbind(expression_data[, c("ID", "LH_rn", genes)], fake_points[, c("ID", "LH_rn", genes)])
        asynchrony <- sapply(expression_data[, "ID"], function(x) {
            dat <- .all_genes(unknown_sample = x, genes = genes, data_dCT = data, gene_window = gene_window)
            .synchrony_bySampleSD(dat, 20)
        })
        x <- 1:length(asynchrony)
        y <- sort(asynchrony)
        linmod <- stats::lm(y ~ x)
        segmod <- segmented::segmented(linmod, seg.Z = ~x)
        break1 <- round(segmod$psi[[2]])
        threshold <- sort(asynchrony)[break1 - 1]
        result[["asynchronous_scores"]] <- asynchrony
        result[["asynchrony_threshold"]] <- threshold
    }

    return(result)
}
