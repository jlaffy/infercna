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
#' @rdname findMalignant
#' @export 
findMalignant = function(cna,
                         threshold = NULL,
                         cor.threshold = threshold,
                         signal.threshold = .9,
                         bySample = FALSE,
                         samples = NULL,
                         sep = "-|_",
                         ...) {
    cors = cnaCor(cna,
                  threshold = cor.threshold,
                  bySample = bySample,
                  samples = samples,
                  sep = sep)

    corGroups = fitBimodal(cors, bySampling = TRUE)
    signals = cnaSignal(cna, threshold = signal.threshold)
    sigGroups = fitBimodal(signals, bySampling = TRUE, nsamp = 10000)
    browser()
}
