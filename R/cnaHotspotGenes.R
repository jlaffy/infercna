
#' @title Find the Genes with the highest CNA signal
#' @description Find Genes in the top nth quantile for CNA Signal values
#' @param cna a matrix of gene rows by cell columns containing CNA values
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. 
#' @return gene names in the top nth quantile, where n is specified via <gene.quantile>
#' @rdname cnaHotspotGenes
#' @export 
cnaHotspotGenes = function(cna, gene.quantile, cell.quantile = NULL) {
    cna = as.matrix(cna)
    stopifnot(is.numeric(gene.quantile))
    if (!is.null(cell.quantile)) {
        stopifnot(is.numeric(cell.quantile))
        cna = cna[, cnaHotspotCells(cna, cell.quantile = cell.quantile)]
    }
    msq = rowMeans(cna^2)
    names(msq)[msq >= quantile(msq, gene.quantile)]
}
