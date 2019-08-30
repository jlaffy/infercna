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

    ratios = .cnaStatistic(cna,
                           cor.threshold = cor.threshold,
                           signal.threshold = signal.threshold,
                           samples = samples,
                           sep = sep)

    isBimodal = fitBimodal(ratios,
                           boolean = T,
                           verbose = F,
                           ...)

    if (!isBimodal) {
        stop('Two modes not found.')
    }

    modes = fitBimodal(ratios,
                       assign = T,
                       verbose = F,
                       ...)

    modeMeans = sapply(modes, function(mo) mean(ratios[mo]))
    i_nonmal = which(modeMeans == min(modeMeans))
    nonmal = modes[[i_nonmal]]
    mal = modes[[-1 * i_nonmal]]
    modes = stats::setNames(list(nonmal, mal), c("nonmalignant", "malignant"))

    if (bySample) {
    }

    modes
}

.findMalignantBySample = function(cna, modes, samples, sep) {
    groups = splitCellsBySample(modes$malignant, samples = samples, sep = sep)
    cnaBySample = sapply(groups, function(gr) cna[, c(modes$nonmalignant, gr)], simplify = F)
    sapply(cnaBySample, findMalignant, cor.threshold = cor.threshold, signal.threshold = signal.threshold)

}
