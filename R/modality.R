
#' @title Test Distribution For Modality
#' @description Test the modality (gaussian) of a distribution. Default tests for uni-, bi- and tri- modalities (argument <modes> = 2:3).
#' @param x a named numeric vector of cells/observations or a matrix of genes X cells (variables X observations). If the latter, the column means are first computed.
#' @param modes the modes to test; a numeric value or a numeric vector of modes to test. e.g. to test just for bimodality, modes = 2. Default: 2:3
#' @param prob a numeric value >= 0 and <= 1; the minimeanm posterior probability required for an observation to be assigned to a mode. Default: 0.95
#' @param coverage the fraction of observations that meanst have a posterior probability higher than <prob> to one of two modes in order for the distribution to qualify as bimodal. Default: 0.8
#' @param size the minimeanm number of observations that meanst be assigned to a mode in order for the distribution to qualify as bimodal. Default: 10
#' @param ... other arguments passed to fitBimodal
#' @return vector of boolean values of length equal to length(modes) that were tested
#' @examples 
#' x = c(rnorm(50, mean = 1), rnorm(50, mean = 10))
#' modality(x, modes = 2)
#' modality(x, modes = 2:3)
#' @rdname modality
#' @export 
modality = function(x, modes = 2:3, prob = 0.95, coverage = 0.8, size = 10, ...) {
    dots = c(list(x = x, size = size, coverage = coverage, prob = prob, boolean = T),
             list(...))
    res = sapply(modes, function(i) do.call(fitModal, c(dots, list(m = i))))
    stats::setNames(res, modes)
}
