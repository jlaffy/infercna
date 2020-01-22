#' @title Find Malignant Subset of Cells
#' @description Find Malignant Subset of Cells according to their CNA parameters for correlation to Tumour of origin (using cnaCor()) and for the overall extent of CNA signal (the mean absolute CNA values; using cnaSignal()).
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @return if bimodality was not found, FALSE. Else return list(malignant, nonmalignant, unassigned) of column IDs.
#' @details The function tries to fit cells to one of two Gaussian modes for each of two measures in turn: 1) CNA correlation, or individual cells' pearson correlations to the average CNA profile of their sample of origin; 2) CNA signal, defined as individual cells' mean squared CNA profiles. Cells assigned to the 2nd (higher) mode by both measures are classified as malignant, cells assigned to the 1st mode (lower) by both measures are classified as non-malignant and cells assigned to one of each remain 'unassigned'. Cells missing from the resulting list are those that could not be assigned to any mode by either measure.
#' @rdname findMalignant
#' @export 
findMalignant = function(cna,
                         refCells = NULL,
                         samples = scalop::unique_sample_names(colnames(cna)),
                         gene.quantile = 0.9,
                         cor.threshold = threshold,
                         signal.threshold = threshold,
                         gauss.bootstraps = TRUE,
                         n.bootstraps = 10000,
                         gauss.prob = 0.8,
                         prob.coverage = 0.8,
                         verbose = TRUE,
                         plot = TRUE,
                         ...) {

    if (verbose) message("Calculating cells' CNA correlations...")
    cors = cnaCor(cna,
                  threshold = cor.threshold,
                  samples = samples,
                  refCells = refCells,
                  sep = samples.sep)

    if (verbose) message("Calculating cells' CNA signal...")
    signals = cnaSignal(cna, threshold = signal.threshold)

    if (verbose) {
        message('Fitting CNA correlations to one of two (CNA-low vs. CNA-high) modes...')
    }
    old <- Sys.time()
    invisible(capture.output(corGroups <- suppressMessages(fitBimodal(cors,
                                                                      bySampling = gauss.bootstraps,
                                                                      nsamp = n.bootstraps,
                                                                      prob = gauss.prob,
                                                                      ...))))
    new <- Sys.time()
    if (isFALSE(corGroups)) return(FALSE)

    if (verbose) {
        timerun = scalop::hms_span(start = old, end = new)
        message('(completed in ', timerun, 's)')
    }
    
    if (verbose) {
        message('Fitting CNA signals to one of two (CNA-low vs. CNA-high) modes...')
    }

    old <- Sys.time()
    invisible(capture.output(sigGroups <- suppressMessages(fitBimodal(signals,
                                                           bySampling = gauss.bootstraps,
                                                           nsamp = n.bootstraps,
                                                           prob = gauss.prob,
                                                           ...))))
    new <- Sys.time()
    if (isFALSE(sigGroups)) return(FALSE)

    if (verbose) {
        timerun = scalop::hms_span(start = old, end = new)
        message('(completed in ', timerun, 's)')
    }
    
    if (verbose) {
        message('Identifying cells classified as CNA-high by both parameters...')
    }

    sect = scalop::comply(corGroups, sigGroups, FUN = intersect)
    unassigned = union(unlist(sect[2,1]), unlist(sect[1,2]))
    len.unassig = length(unassigned)
    corGroups = sapply(corGroups, function(gr) gr[!gr %in% unassigned], simplify = F)
    sigGroups = sapply(sigGroups, function(gr) gr[!gr %in% unassigned], simplify = F)
    sect = scalop::comply(corGroups, sigGroups, FUN = intersect)
    result = list(
                  malignant = unlist(sect[2,2], use.names = F),
                  nonmalignant = unlist(sect[1,1], use.names = F),
                  unassigned = unassigned)

    if (length(unassigned) >= 1) {
        warning(round(100 * len.unassig/ncol(cna),2),
                ' % of cells were assigned to opposing modes,',
                '\ni.e.to high CNA-signal and low CNA-correlation or vice versa.',
                '\nThese cells will remain unasssigned.')
    }

    if (plot) {
        cnaScatterPlot(cna,
                       threshold = threshold,
                       signal.threshold,
                       cor.threshold = cor.threshold,
                       samples = samples,
                       refCells = refCells,
                       group = result)
    }

    result

}
