#' @title Center a matrix column-wise 
#' @description Center a matrix column-wise 
#' @param m a matrix
#' @return a matrix with colMeans == 0
#' @rdname colCenter
#' @export 
colCenter = function(m) {
    scale(m, center = T, scale = F)
}


#' @title Center a matrix row-wise 
#' @description Center a matrix row-wise 
#' @param m a matrix
#' @return a matrix with rowMeans == 0
#' @rdname rowCenter
#' @export 
rowCenter = function(m) {
    t(scale(t(m), center = T, scale = F))
}

#' @title Squish matrix values into range
#' @description Squish matrix values into range
#' @param m matrix to manipulate
#' @param range numeric vector of length two giving desired output range. Default: c(-3, 3)
#' @rdname clip
#' @export 
clip <- function(m, range = c(-3, 3)) {
    m = as.matrix(m)
    m[m < range[[1]]] <- range[[1]]
    m[m > range[[2]]] <- range[[2]]
    m
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param m PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 100
#' @param endrule PARAM_DESCRIPTION, Default: 'mean'
#' @param align PARAM_DESCRIPTION, Default: 'center'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[caTools]{runmean}}
#' @rdname runMean
#' @export 
#' @importFrom caTools runmean
runMean <- function(m,
                    k = 100,
                    endrule = 'mean',
                    align = 'center') {

    m = as.matrix(m)
    mout = caTools::runmean(m,
                            k = k,
                            endrule = endrule,
                            align = align)
    colnames(mout) = colnames(m)
    rownames(mout) = rownames(m)
    as.data.frame(mout)
}

