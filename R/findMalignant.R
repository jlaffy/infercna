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
                         signal.threshold = threshold,
                         bySample = FALSE,
                         samples = NULL,
                         sep = "-|_",
                         mode.overlap = 0.05,
                         excludeFromAvg = NULL,
                         cor.bySampling = TRUE,
                         signal.bySampling = TRUE,
                         nsamples = 10000,
                         ...) {
    cors = cnaCor(cna,
                  threshold = cor.threshold,
                  bySample = bySample,
                  samples = samples,
                  excludeFromAvg = excludeFromAvg,
                  sep = sep)
    corGroups = fitBimodal(cors, bySampling = cor.bySampling, nsamp = nsamples)
    
    signals = cnaSignal(cna, threshold = signal.threshold)
    sigGroups = fitBimodal(signals, bySampling = signal.bySampling, nsamp = nsamples)
    
    if (scalop::jaccard(sigGroups$a, corGroups$b) > mode.overlap) return(FALSE)
    if (scalop::jaccard(sigGroups$b, corGroups$a) > mode.overlap) return(FALSE)
    sect = scalop::comply(corGroups, sigGroups, FUN = intersect)
    remove = union(unlist(sect[2,1]), unlist(sect[1,2]))
    message('removing ', length(remove), ' cells that were found in both modes')
    corGroups = sapply(corGroups, function(gr) gr[!gr %in% remove], simplify = F)
    sigGroups = sapply(sigGroups, function(gr) gr[!gr %in% remove], simplify = F)
    sect = scalop::comply(corGroups, sigGroups, FUN = intersect)
    list(a = unlist(sect[1,1], use.names = F), b = unlist(sect[2,2], use.names = F))
}
