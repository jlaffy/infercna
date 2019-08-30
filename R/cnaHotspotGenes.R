
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cna PARAM_DESCRIPTION
#' @param threshold PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname cnaHotspotGenes
#' @export 
cnaHotspotGenes = function(cna, threshold) {
    stopifnot(is.numeric(threshold))
    msq = rowMeans(cna^2)
    names(msq)[msq >= stats::quantile(msq, threshold)]
}
