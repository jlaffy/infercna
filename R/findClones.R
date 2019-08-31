
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param m PARAM_DESCRIPTION
#' @param prob PARAM_DESCRIPTION, Default: 0.95
#' @param coverage PARAM_DESCRIPTION, Default: 0.8
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
                      by = 'chr') {

    mats = splitGenes(m, by = by)
    modes = sapply(mats, fitBimodal, prob = prob, coverage = coverage, size = size, assign = T)
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
#' @description FUNCTION_DESCRIPTION
#' @param m a matrix of genes X cells (variables X observations) containing CNA values.
#' @param prob a numeric value >= 0 and <= 1; the minimum posterior probability required for an observation to be assigned to a mode. Default: 0.95
#' @param coverage the fraction of observations that must have a posterior probability higher than <prob> to one of two modes in order for the distribution to qualify as bimodal. Default: 0.8
#' @param mode.size the minimum number of observations required to define a mode. Default: 10
#' @param clone.size the minimum number of cells required to define a clone. Default: 3
#' @param by PARAM_DESCRIPTION, Default: 'chr'
#' @return 
#' @examples 
#' cna = infercna(useData(), reference = reference)
#' # malignant cells only
#' cna = cna[, !colnames(cna) %in% unlist(reference)]
#' findClones(cna, by = 'chr')
#' findClones(cna, by = 'arm')
#' @rdname findClones
#' @export 
findClones = function(m,
                      prob = 0.95,
                      coverage = 0.8,
                      mode.size = 10,
                      clone.size = 3,
                      by = 'chr') {

    modes = fetchModes(m,
                       prob = prob,
                       coverage = coverage,
                       size = mode.size,
                       by = by)

    expandToClones(modes, greaterThan = clone.size)
}
