#' @title Calculate the Means of Squared CNA Values
#' @description Calculate the Mean of Squared CNA Values
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param threshold calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @return a numeric vector of CNA signal values or the Mean of Squared CNA values
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
