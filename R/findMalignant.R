#' @title Find Malignant Subset of Cells
#' @description Find Malignant Subset of Cells according to their CNA parameters for correlation to Tumour of origin (using cnaCor()) and for the overall extent of CNA signal (the mean absolute CNA values; using cnaSignal()).
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @return if bimodality was not found, FALSE. Else return list(nonMalignant, malignant) of column IDs.
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname isMalignant
#' @export 
isMalignant = function(cna, samples = NULL, threshold = .9, sep = "-|_") {
    hotspotGenes = cnaHotspotGenes(cna, threshold = threshold)
    cna = cna[hotspotGenes, ]
    ratios = cnaCor(cna, samples = samples)/cnaSignal(cna)
    stopifnot(isTRUE(fitBimodal(ratios, boolean = T)))
    modes = fitBimodal(ratios, assign = T)
    modeRatioMean = sapply(modes, function(mo) mean(ratios[mo]))
    i_nonmal = which(modeRatioMean == min(modeRatioMean))
    nonmal = modes[i_nonmal]
    mal = modes[-1 * i_nonmal]
    stats::setNames(list(nonmal, mal), c("nonMalignant", "malignant"))
}
