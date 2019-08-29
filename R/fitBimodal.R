#' @title Fit a Bimodal Gaussian Distribution 
#' @description Fit a bimodal gaussian distribution to a set of observations. 
#' @param x a named numeric vector of cells/observations or a matrix of genes X cells (variables X observations). If the latter, the column means are first computed.
#' @param prob a numeric value >= 0 and <= 1; the minimum posterior probability required for an observation to be assigned to a mode. Default: 0.95
#' @param coverage the fraction of observations that must have a posterior probability higher than <prob> to one of two modes in order for the distribution to qualify as bimodal. Default: 0.8
#' @param size the minimum number of observations that must be assigned to a mode in order for the distribution to qualify as bimodal. Default: 10
#' @param assign if set to TRUE, returns a list of length two containing the vector names that were assigned to each mode. Default: FALSE
#' @param boolean if set to TRUE, returns a boolean value indicating whether the distribution is bimodal. Default: FALSE
#' @param verbose print progress messages. Default: TRUE 
#' @param maxrestarts the maximum number of restarts allowed. See \code{\link[mixtools]{normalmixEM}} for details. Default: 200
#' @param maxit the maximum number of iterations. Default: 10000
#' @param k number of components. Default: 2
#' @return The posterior probabilities of each observation to one of two modes. If boolean = TRUE, return a boolean value indicating whether bimodality was found. If assign = TRUE, return a list of length two with the observations (IDs) in each mode.
#' @examples 
#'  cna = infercna(m = useData(), reference = reference)
#'  # Malignant cells only (remove columns corresponding to reference)
#'  cna = cna[, !colnames(cna) %in% unlist(reference)] 
#'  cnaByChr = splitGenes(cna, by = 'chr')
#'  sapply(cnaByChr, fitBimodal, assign = TRUE)
#'  sapply(cnaByChr, fitBimodal, boolean = TRUE)
#'  sapply(cnaByChr, fitBimodal, boolean = TRUE, coverage = 0.5)
#' @seealso 
#'  \code{\link[mixtools]{normalmixEM}}
#' @rdname fitBimodal 
#' @export 
#' @importFrom stats setNames
#' @importFrom mixtools normalmixEM
#' @importFrom dplyr mutate
fitBimodal = function(x,
                      prob = 0.95,
                      coverage = 0.8,
                      size = 10,
                      assign = FALSE,
                      boolean = FALSE,
                      verbose = TRUE,
                      maxrestarts = 200,
                      maxit = 10000,
                      k = 2) {

    if (!is.null(dim(x))) x = colMeans(x)

    if (length(x) < size * 2) {
        stop('Number of observations is too small for 2 modes >= ' , size)
    }

    result = try(mixtools::normalmixEM(x, maxrestarts = maxrestarts, maxit = maxit)$posterior)
    passed = result >= prob
    failed = is.null(dim(passed)) || ncol(passed) != 2 || nrow(passed) != length(x)
    
    if (failed) {
        if (verbose) message('No bimodal distribution found.')
        return(FALSE)
    }

    if (any(colSums(passed) < size)) {
        if (verbose) message('At least one mode contains < ', size, ' obs.')
        return(FALSE)
    }

    observed = (sum(colSums(passed))) / length(x)
    if (observed < coverage) {
        if (verbose) message('Less than ', coverage * 100, '% of obs. could be assigned to a mode.')
        return(FALSE)
    }

    if (verbose) message('Success!')

    if (assign) {
        a = names(x)[which(passed[,1])]
        b = names(x)[which(passed[,2])]
        return(list(a = a, b = b))
    }

    if (boolean) return(TRUE)

    result
}
