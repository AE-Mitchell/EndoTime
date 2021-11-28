.synchrony_score_sub <- function(densities, medians, interval) {
    mode <- which.max(densities)
    lower_limit <- mode - interval
    if (lower_limit < 1) {
        mirrored1 <- abs(lower_limit) + 1
        lower_vicinity <- c(densities[1:mode], densities[1:mirrored1])
    } else {
        lower_vicinity <- densities[lower_limit:mode]
    }
    upper_limit <- mode + interval
    if (upper_limit > length(medians)) {
        mirrored2 <- upper_limit - length(medians) - 1
        upper_vicinity <- c(densities[(mode):length(medians)], densities[(length(medians) - mirrored2):length(medians)])
    } else {
        upper_vicinity <- densities[(mode):upper_limit]
    }

    vicinity <- c(lower_vicinity, upper_vicinity)
    score <- sum(vicinity)
    return(score)
}
