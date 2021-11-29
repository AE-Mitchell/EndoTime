.extrapolate_data <- function(endo_data, genes, size, dp_per_LH = 33) {
    # This current version generates the bin size *radius*, or the input size/2, number of extrapolated points
    # at either end of the timescale

    addition <- 1 / dp_per_LH
    upper_boundary_LH <- max(endo_data[, "LH_rn"]) + (c(1:(size/2)) * addition)
    lower_boundary_LH <- min(endo_data[, "LH_rn"]) - (c(1:(size/2)) * addition)
    labs <- c(paste("a", seq(1:(size/2)), sep = ""), paste("b", seq(1:(size/2)), sep = ""))

    model_results <- sapply(genes, function(gene) {
        gene_expressions <- data.frame("LH_rn" = endo_data[["LH_rn"]], "expression" = endo_data[[gene]])

        # Generate models linear and local regression models to describe the data
        lm_model <- stats::lm(gene_expressions[["expression"]] ~ gene_expressions[["LH_rn"]])
        loess_model <- stats::loess(gene_expressions[["expression"]] ~ gene_expressions[["LH_rn"]])

        # Create a data frame from the fitted local regression model, arrange in LH order and identify the expression of points at either end
        loess_curve <- data.frame("LH_rn" = loess_model[["x"]][, 1], "expression" = loess_model[["fitted"]])
        loess_curve <- dplyr::arrange(loess_curve, loess_curve[["LH_rn"]])
        lowest_curve_point <- loess_curve[["expression"]][1]
        highest_curve_point <- loess_curve[["expression"]][nrow(loess_curve)]

        ###
        ### LOWER END
        ###

        # Generate an extrapolated series of points at the lower end, identify the point to connect with the loess curve
        lower_boundary_dCT_slope <- stats::coef(lm_model)[2] * lower_boundary_LH
        lower_boundary_highest_point <- lower_boundary_dCT_slope[1]

        # Generate a linear model using the lowest-LH points, width equal to half the size of the bin and get its slope
        lower_end_data <- loess_curve[seq(1, size/2), ]
        lower_end_slope <- stats::coef(stats::lm(expression ~ LH_rn, data = lower_end_data))[2]

        # If the slope of the overall linear model is positive ...
        if (sign(stats::coef(lm_model)[2]) == 1) {
            # Is the slope of the lower end positive? If so, place the extrapolated points extending the loess curve
            # If not, place the extrapolated points so that they follow the linear model, which will be below the loess curve
            if (sign(lower_end_slope) == 1) {
                lower_loess_slope_difference <- lowest_curve_point - lower_boundary_highest_point
                adjusted_lower_boundary_dCT <- lower_boundary_dCT_slope + lower_loess_slope_difference
            } else {
                adjusted_lower_boundary_dCT <- lower_boundary_dCT_slope + stats::coef(lm_model)[1]
            }
        } else {
            if (sign(lower_end_slope) == 1) {
                adjusted_lower_boundary_dCT <- lower_boundary_dCT_slope + stats::coef(lm_model)[1]
            } else {
                lower_loess_slope_difference <- lowest_curve_point - lower_boundary_highest_point
                adjusted_lower_boundary_dCT <- lower_boundary_dCT_slope + lower_loess_slope_difference
            }
        }

        ###
        ### UPPER END
        ###
        # Generate an extrapolated series of points at the upper end, identify the point to connect with the loess curve
        upper_boundary_dCT_slope <- stats::coef(lm_model)[2] * upper_boundary_LH
        upper_boundary_lowest_point <- upper_boundary_dCT_slope[1]

        # Generate a linear model using the highest-LH points, width equal to half the size of the bin and get its slope
        upper_end_data <- loess_curve[seq(nrow(endo_data) - (size/2) + 1, nrow(endo_data)), ]
        upper_end_slope <- stats::coef(stats::lm(expression ~ LH_rn, data = upper_end_data))[2]

        # If the slope of the overall linear model is positive ...
        if (sign(stats::coef(lm_model)[2]) == 1) {
            # Is the slope of the upper end positive? If so, place the extrapolated points extending the loess curve
            # If not, place the extrapolated points so that they follow the linear model, which will be above the loess curve
            if (sign(upper_end_slope) == 1) {
                upper_loess_slope_difference <- highest_curve_point - upper_boundary_lowest_point
                adjusted_upper_boundary_dCT <- upper_boundary_dCT_slope + upper_loess_slope_difference
            } else {
                adjusted_upper_boundary_dCT <- upper_boundary_dCT_slope + stats::coef(lm_model)[1]
            }
        } else {
            if (sign(upper_end_slope) == 1) {
                adjusted_upper_boundary_dCT <- upper_boundary_dCT_slope + stats::coef(lm_model)[1]
            } else {
                upper_loess_slope_difference <- highest_curve_point - upper_boundary_lowest_point
                adjusted_upper_boundary_dCT <- upper_boundary_dCT_slope + upper_loess_slope_difference
            }
        }

        results <- list("gene_extrapolation" = c(adjusted_upper_boundary_dCT, adjusted_lower_boundary_dCT))
        return(results)
    }, USE.NAMES = TRUE, simplify = FALSE)

    # Construct data frame of gene extrapolations
    df_1 <- data.frame(sapply(model_results, function(gene) {
        return(gene[["gene_extrapolation"]])
    }, USE.NAMES = TRUE))
    df_1["ID"] <- labs
    df_1["LH_rn"] <- c(upper_boundary_LH, lower_boundary_LH)
    df_1 <- df_1[c("ID", "LH_rn", dplyr::all_of(genes))]

    extrapolated_data <- rbind(endo_data[c("ID", "LH_rn", genes)], df_1[c("ID", "LH_rn", genes)])

    return(extrapolated_data)
}
