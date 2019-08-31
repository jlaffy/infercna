#' @title Center a matrix column-wise 
#' @description Center a matrix column-wise 
#' @param m a matrix
#' @param method a character string. one of "mean" or "median" definining how to do the centering. Default: "mean"
#' @param centerBy a numeric vector of length equal to ncol(m) of values to be subtracted from each row.
#' @return a matrix centered by columns such that the statistic defined in method is equal to 0 for each column.
#' @rdname colCenter
#' @export 
colCenter = function(m, method = 'mean', centerBy = NULL) {
    if (!is.null(centerBy)) {
        stopifnot(length(centerBy) == ncol(m))
        center = centerBy
    }
    else if (method == 'median') center = apply(m, 2, stats::median)
    else if (method == 'mean') center = T
    else stop('method not recognised')
    scale(m, center = center, scale = F)
}


#' @title Center a matrix row-wise 
#' @description Center a matrix row-wise 
#' @param m a matrix
#' @param method a character string. one of "mean" or "median" definining how to do the centering. Default: "mean"
#' @param centerBy a numeric vector of length equal to nrow(m) of values to be subtracted from each row.
#' @return a matrix centered by rows such that the statistic defined in method is equal to 0 for each column.
#' @rdname rowCenter
#' @export 
rowCenter = function(m, method = 'mean', centerBy = NULL) {
    if (!is.null(centerBy)) {
        stopifnot(length(centerBy) == nrow(m))
        center = centerBy
    }
    else if (method == 'median') center = apply(m, 1, stats::median)
    else if (method == 'mean') center = T
    else stop('method not recognised')
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

#' @title Rolling Means
#' @description Apply a rolling window mean to a matrix or vector.
#' @param m a numeric vector or matrix. If the latter, each column will be processed separately.
#' @param k width of rolling window. Default: 100
#' @param endrule character string indicating how the values at the beginning and the end of the data should be treated. One of "mean", "trim", "keep", "constant". See caTools::runmean for more details. Default: 'mean'
#' @param align specifies whether result should be centered (default), left-aligned or right-aligned. See caTools::runmean for more details. Default: 'center'
#' @param verbose print progress messages. Default: TRUE
#' @return a numeric vector or matrix of the same size as <m>. Only in case of endrule=trim, the output vectors will be shorter and output matrices will have fewer rows. 
#' @seealso 
#'  \code{\link[caTools]{runmean}}
#' @rdname runMean
#' @export 
#' @importFrom caTools runmean
runMean <- function(m,
                    k = 100,
                    endrule = 'mean',
                    align = 'center',
                    verbose = T) {

    if (!is.null(dim(m))) {
        m = as.matrix(m)
    }
    if (is.null(m) || nrow(m) < k) {
        k = nrow(m)
        if (verbose) message('Setting <k> to nrow(m): ', k)
    }
    mout = caTools::runmean(m,
                            k = k,
                            endrule = endrule,
                            align = align)
    if (!is.null(dim(m))) {
        colnames(mout) = colnames(m)
        rownames(mout) = rownames(m)
        mout = as.data.frame(mout)
    } else {
        names(mout) = names(m)
    }

    mout
}

