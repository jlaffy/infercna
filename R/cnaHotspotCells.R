
#' @title Find the Cells with the highest CNA signal
#' @description Find Cells in the top nth quantile for CNA Signal values
#' @param cna a matrix of cell rows by cell columns containing CNA values
#' @param cell.quantile calculate CNA measures including only top / "hotspot" cells according to their squared CNA values across all genes. Value between 0 and 1 denoting the quantile of cells to include. 
#' @return cell names in the top nth quantile, where n is specified via <cell.quantile>
#' @rdname cnaHotspotCells
#' @export 
cnaHotspotCells = function(cna, cell.quantile, gene.quantile = NULL) {
    cna = as.matrix(cna)
    stopifnot(is.numeric(cell.quantile))
    if (!is.null(gene.quantile)) {
        stopifnot(is.numeric(gene.quantile))
        cna = cna[cnaHotspotGenes(cna, gene.quantile = gene.quantile), ]
    }
    msq = colMeans(cna^2)
    names(msq)[msq >= quantile(msq, cell.quantile)]
}
