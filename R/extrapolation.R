.extrapolation <- function(dat, genes, size, dp_per_LH = 33) {
    addition <- 1 / dp_per_LH
    upper_boundary_LH <- max(dat[, "LH_rn"]) + (c(1:size) * addition)
    lower_boundary_LH <- min(dat[, "LH_rn"]) - (c(1:size) * addition)
    labs <- c(paste0("a", seq(1:size)), paste0("b", seq(1:size)))

    df <- data.frame(ID = labs, LH_rn = c(upper_boundary_LH, lower_boundary_LH))

    for (i in genes) {
        lm_reg <- stats::lm(dat[, i] ~ dat[, "LH_rn"])
        upper_boundary_dCT <- stats::coef(lm_reg)[2] * upper_boundary_LH + stats::coef(lm_reg)[1]
        lower_boundary_dCT <- stats::coef(lm_reg)[2] * lower_boundary_LH + stats::coef(lm_reg)[1]
        df[, i] <- c(upper_boundary_dCT, lower_boundary_dCT)
    }
    return(df)
}
