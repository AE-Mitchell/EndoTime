#' Title
#'
#' @param expression_data ...
#'
#' @return ...
#' @export
.enforce_monotony <- function(expression_data) {
    expression_data["LH_rn"] <- expression_data[["LH"]] + stats::runif(nrow(expression_data), -0.5, 0.5)
    expression_data["rank"] <- rank(expression_data[["LH_rn"]])
    time_lm <- stats::lm(LH_rn ~ rank, data = expression_data)
    expression_data["LH_rn"] <- stats::predict(time_lm, newdata = data.frame("rank" = expression_data[["rank"]]))
    expression_data["rank"] <- NULL
    expression_data
}
