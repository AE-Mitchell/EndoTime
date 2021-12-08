.single_gene_sub <- function(rest_LH, rest_gene, unknown_expression, start_position, end_position) {
    current_LH <- rest_LH[start_position:end_position]
    median_LH <- stats::median(current_LH)
    current_points <- rest_gene[start_position:end_position]

    weights1 <- (abs(current_LH[1] - current_LH[which(current_LH <= median_LH)])) / (abs(current_LH[1] - median_LH))
    weights2 <- (abs(current_LH[length(current_LH)] - current_LH[which(current_LH > median_LH)])) / (abs(current_LH[length(current_LH)] - median_LH))
    weights <- c(weights1, weights2)
    weights[current_LH == median_LH] <- 1

    l_sd <- radiant.data::weighted.sd(current_points, weights)
    if (is.na(l_sd)) {
        l_sd <- 0
    }

    l_mean <- mean(current_points)

    current_dens <- stats::dnorm(unknown_expression, l_mean, l_sd)

    c(median_LH, current_dens)
}
