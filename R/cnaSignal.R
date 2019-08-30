#' @title Calculate the Means of Absolute CNA Values
#' @description Calculate the Mean of Absolute CNA Values
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @return a numeric vector of CNA signal values or the Mean of Absolute CNA values
#' @rdname cnaSignal
#' @export 
cnaSignal = function(cna, threshold = NULL) {
    cna = as.matrix(cna)
    if (!is.null(threshold)) {
        hotspotGenes = cnaHotspotGenes(cna, threshold = threshold)
        cna = cna[hotspotGenes, ]
    }
    sqmat = cna^2
    colMeans(sqmat)
}
