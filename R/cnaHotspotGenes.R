
#' @title Find the Genes with the highest CNA signal
#' @description Find Genes in the top nth quantile for CNA Signal values
#' @param cna a matrix of gene rows by cell columns containing CNA values
#' @param threshold calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. 
#' @return gene names in the top nth quantile, where n is specified via <threshold>
#' @rdname cnaHotspotGenes
#' @export 
cnaHotspotGenes = function(cna, threshold) {
    stopifnot(is.numeric(threshold))
    msq = rowMeans(cna^2)
    names(msq)[msq >= quantile(msq, threshold)]
}
