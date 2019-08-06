#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param m PARAM_DESCRIPTION, Default: NULL
#' @param X PARAM_DESCRIPTION, Default: NULL
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
#' @seealso 
#'  \code{\link[stats]{setNames}}
#'  \code{\link[mixtools]{normalmixEM}}
#'  \code{\link[dplyr]{mutate}}
#' @rdname modes
#' @export 
#' @importFrom stats setNames
#' @importFrom mixtools normalmixEM
#' @importFrom dplyr mutate
modes = function(m = NULL,
                 X = NULL,
                 prob = 0.95,
                 coverage = 0.8,
                 size = 10,
                 boolean = FALSE) {

    if (is.null(X)) X = colMeans(m)
    negativeRes = stats::setNames(rep(0, length(X)), names(X))
    probs = try(mixtools::normalmixEM(X)$posterior)
    if (class(probs) == 'try-error' | is.null(probs)) {
        if (boolean) return(FALSE)
        else return(negativeRes)
    }
    probs = as.data.frame(probs)
    probs = dplyr::mutate(probs,
                          bool.1 = comp.1 >= prob,
                          bool.2 = comp.2 >= prob,
                          bool = bool.1 | bool.2)
    modesexist = mean(probs$bool) > coverage && sum(probs$bool.1) >= size && sum(probs$bool.2) >= size
    if (boolean) return(modesexist)
    if (!modesexist) return(negativeRes)
    probs = dplyr::mutate(probs, modeID = ifelse(comp.1 >= prob, 1, ifelse(comp.2 >= prob, 2, 0)))
    stats::setNames(probs$modeID, names(X))
}

