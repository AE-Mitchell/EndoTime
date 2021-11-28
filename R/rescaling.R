.rescaling <- function(values, desired_min, desired_max){
    current_min <- min(values)
    current_max <- max(values)

    result <- (((desired_max - desired_min) * (values - current_min)) / (current_max - current_min)) + desired_min
    return(result)
}
