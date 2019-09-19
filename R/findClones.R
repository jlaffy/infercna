
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param m PARAM_DESCRIPTION
#' @param prob PARAM_DESCRIPTION, Default: 0.95
#' @param coverage PARAM_DESCRIPTION, Default: 0.8vto
#' @param size PARAM_DESCRIPTION, Default: 10
#' @param by PARAM_DESCRIPTION, Default: 'chr'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname fetchModes
#' @export 
fetchModes = function(m,
                      prob = 0.95,
                      coverage = 0.8,
                      size = 10,
                      by = 'chr',
                      bySampling = FALSE,
                      nsamp = 2000,
                      minGenes = 50,
                      ...) {

    mats = splitGenes(m, by = by)
    Rows = unlist(sapply(mats, nrow))
    Rows[is.null(Rows)] = 0
    mats = mats[Rows >= minGenes]
    modes = sapply(mats, fitBimodal, prob = prob, coverage = coverage, size = size, assign = T, bySampling = bySampling, nsamp = nsamp, ...)
    modes[!sapply(modes, isFALSE)]
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param modes PARAM_DESCRIPTION
#' @param greaterThan PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{setNames}}
#' @rdname expandToClones
#' @export 
#' @importFrom stats setNames
expandToClones = function(modes,
                          greaterThan = NULL) {

    modes = sapply(modes, function(x) {
                       stats::setNames(unlist(x), rep(names(x), lengths(x)))},
                       simplify = F)

    modes = unlist(modes)
    clones = sapply(split(names(modes), modes), paste0, collapse = '--')
    clones = split(names(clones), clones)
    if (is.null(greaterThan)) return(clones)
    clones[lengths(clones) > greaterThan]
}


#' @title Find Clones 
#' @description Assign cells to genetic subclones from their inferred CNA profiles. You can compute their CNA profiles using infercna::infercna().
#' @param m a matrix of genes X cells (variables X observations) containing CNA values.
#' @param prob a numeric value >= 0 and <= 1; the minimum posterior probability required for an observation to be assigned to a mode. Default: 0.95
#' @param coverage the fraction of observations that must have a posterior probability higher than <prob> to one of two modes in order for the distribution to qualify as bimodal. Default: 0.8
#' @param mode.size the minimum number of observations required to define a mode. Default: 10
#' @param clone.size the minimum number of cells required to define a clone. Default: 3
#' @param by PARAM_DESCRIPTION, Default: 'chr'
#' @return OUTPUT 
#' @rdname findClones
#' @export 
findClones = function(m,
                      prob = 0.95,
                      coverage = 0.8,
                      mode.size = 10,
                      clone.size = 3,
                      by = 'chr',
                      bySampling = FALSE,
                      nsamp = 2000,
                      force.tries = FALSE,
                      verbose = FALSE,
                      ...) {

    if (is.list(as.matrix(m))) {
        if (verbose) {
            modes = sapply(m, fitBimodal, prob = prob, coverage = coverage, size = mode.size, assign = T, bySampling = bySampling, nsamp = nsamp, force.tries = force.tries,...)
        } else {
            modes = suppressWarnings(sapply(m, fitBimodal, prob = prob, coverage = coverage, size = mode.size, assign = T, bySampling = bySampling, nsamp = nsamp, force.tries = force.tries,...))
        }
    } else {
        if (verbose) {
            modes = fetchModes(m, prob = prob, coverage = coverage, size = mode.size, by = by, bySampling = bySampling, nsamp = nsamp,force.tries = force.tries, ...)
        } else {
            modes = suppressWarnings(fetchModes(m, prob = prob, coverage = coverage, size = mode.size, by = by, bySampling = bySampling, nsamp = nsamp, force.tries = force.tries,...))
        }
    }
    
    expandToClones(modes, greaterThan = clone.size)
}
