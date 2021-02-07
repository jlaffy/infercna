
#' @importFrom scalop colcenter
#' @export
scalop::colcenter

#' @importFrom scalop rowcenter
#' @export
scalop::rowcenter

#' @importFrom scalop logtpm
#' @export
scalop::logtpm

#' @importFrom scalop unlogtpm
#' @export
scalop::unlogtpm

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

    if (is.null(dim(m))) return(FALSE)
    if (nrow(m) == 0) return(FALSE)

    if (nrow(m) < k) {
        k = nrow(m)
    	chr = splitGenes(rownames(m)) %>% .[lengths(.) > 0] %>% names
        if (verbose) message('Adjusting window to the max. number of genes on chromosome ', chr, ' (', k, ')')
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

