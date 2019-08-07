#' @title Find Genetic Clones
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
#' @rdname findClones
#' @export 
findClones = function(m, prob = 0.95, rate = 0.8, size = 10) {

    .fetchData = function(m, splitBy) {
        mats = splitGenes(m, by = splitBy)
        modes = sapply(mats, fitBimodal, assign = T)
        modes[!sapply(modes, isFALSE)]}

    armModes = .fetchData(m, splitBy = 'arm')
    chrModes = .fetchData(m, splitBy = 'chr')



    modemat = modemat[rowSums(modemat) != 0, colSums(modemat) != 0, drop = FALSE]
    gr1 = apply(modemat, 1, function(r) paste0(colnames(modemat)[which(r == 1)], collapse = '_'))
    gr2 = apply(modemat, 1, function(r) paste0(colnames(modemat)[which(r == 2)], collapse = '_'))
    cloneIDs = mapply(function(x, y) paste(x, y, sep = '-'), x = gr1, y = gr2)
    split(names(cloneIDs), cloneIDs)
}
