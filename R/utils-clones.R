
.envmatch = function(genes) {
    match(genes, g.env$gene, nomatch = 0)
}

.envglob = function(genes = NULL) {
    L = as.list(g.env)
    L = L[sapply(L, function(l) is.atomic(l) && length(l) > 1)]
    if (!is.null(genes)) L = sapply(L, `[`, .envmatch(genes), simplify = F)
    L
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param m PARAM_DESCRIPTION, Default: NULL
#' @param genes PARAM_DESCRIPTION, Default: NULL
#' @param by PARAM_DESCRIPTION, Default: 'arm'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname splitGenes
#' @export 
splitGenes = function(m = NULL, genes = NULL, by = 'arm') {
    if (!is.null(m)) genes = rownames(m)
    L = .envglob(genes = genes)
    stopifnot(by %in% names(L))
    splut = split(L$gene, L[[by]])
    if (!is.null(m)) {
        splut = sapply(splut, function(group) m[group, ], simplify = F)}
    splut
}

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

