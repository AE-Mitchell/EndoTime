.rescaling <- function(values, desired_min, desired_max){
    (((desired_max - desired_min) * (values - min(values))) / (max(values) - min(values))) + desired_min
}
