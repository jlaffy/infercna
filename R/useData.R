
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION, Default: '125'
#' @param random PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[scrabble]{rowCenter}}
#' @rdname useData
#' @export 
#' @importFrom scrabble rowCenter
useData = function(x = '125', random = FALSE) {

    if (random) x = sample(c('125', 'bt771'), 1)

    if (x == '125') {
        m = mgh125
        cellnm = cells125
    }

    else if (x == "771") {
        m = bt771
        cellnm = cells771
    }

    m = as.matrix(m)[genes, ]
    colnames(m) = cellnm
    scrabble::rowCenter(m)
}
