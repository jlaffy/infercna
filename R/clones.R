
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param m PARAM_DESCRIPTION, Default: NULL
#' @param L PARAM_DESCRIPTION, Default: NULL
#' @param by PARAM_DESCRIPTION, Default: 'arm'
#' @param prob PARAM_DESCRIPTION, Default: 0.95
#' @param coverage PARAM_DESCRIPTION, Default: 0.8
#' @param size PARAM_DESCRIPTION, Default: 10
#' @param boolean PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname clones
#' @export 
clones = function(m = NULL,
                  L = NULL,
                  by = 'arm',
                  prob = 0.95,
                  coverage = 0.8,
                  size = 10,
                  boolean = FALSE) {

    if (is.null(L)) {
        L = genesplit(m = m, by = by)
        L = L[sapply(L, nrow) >= 5]
    }
    modemat = sapply(L, modes, prob = prob, coverage = coverage, size = size, boolean = boolean, simplify = T)
    browser()
    stopifnot(is.matrix(modemat))
    modemat = modemat[rowSums(modemat) != 0, colSums(modemat) != 0, drop = FALSE]
    gr1 = apply(modemat, 1, function(r) paste0(colnames(modemat)[which(r == 1)], collapse = '_'))
    gr2 = apply(modemat, 1, function(r) paste0(colnames(modemat)[which(r == 2)], collapse = '_'))
    cloneIDs = mapply(function(x, y) paste(x, y, sep = '-'), x = gr1, y = gr2)
    split(names(cloneIDs), cloneIDs)
}
