#' @title Find Malignant Subset of Cells
#' @description Find the malignant and non-malignant subsets of cells from a gene-by-cell matrix of CNA values. 
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param refCells a character vector or list of character vectors denoting the cell IDs of the known non-malignant 'Reference' cells (usually those that were used in the call to infercna::infercna to correct the CNA values). This is relevant only if these cells are in the CNA matrix. Default: NULL
#' @param samples a character vector of unique sample names to group the cells by. This is relevant only if multiple tumours are represented in the cna matrix, as it allows tumour-specific average CNA profiles to be computed in the call to infercna::cnaCor. See Details section for more information. Default: scalop::unique_sample_names(colnames(cna))
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: 0.9
#' @param gene.quantile.for.corr as above but for CNA correlations specifically. Default: gene.quantile 
#' @param gene.quantile.for.signal as above but for CNA signal specifically. Default: gene.quantile
#' @param use.bootstraps logical; if TRUE, the function uses a bootstrapping method to subsample values and identify the malignant and non-malignant groups iteratively. This method is more sensitive to differing group sizes, so will be useful if you believe one group to be much smaller than the other. Default: TRUE
#' @param n.bootstraps number of bootstrap replicates. Relevant only if <use.bootstraps> is TRUE. Default: 10000
#' @param gauss.assignment.prob a numeric value >= 0 and <= 1; the minimum posterior probability required for an cell to be assigned to a mode. Default: 0.8
#' @param gauss.assignment.coverage the fraction of cells that must have a posterior probability higher than <prob> to one of two modes in order for the distribution to qualify as bimodal. Default: 0.8
#' @param verbose print progress messages. Default: TRUE
#' @param plot logical; scatter plot of cells' CNA correlations against CNA signal. This uses infercna::cnaScatterPlot. In order for the infercna::findMalignant function to be meaningful, expect two distinct groups of cells in the plot. Default: TRUE
#' @param ... other arguments passed to infercna::fitBimodal. 
#' @return if bimodality was not found, returns FALSE. If bimodality was found, returns a list containing the malignant and non-malignant cell IDs.
#' @details This function attempts to fit gaussian bimodal distributions to each of two parameters that describe the extent of CNAs per cell. The first of these is CNA correlation (infercna::cnaCor), which measures the pearson correlations between individual cells' CNA profiles and the average profile of their tumours of origin. The second measure is CNA signal (infercna::cnaSignal), which computes the cell averages of squared CNA values. The function then assigns cells residing in the 2nd (higher-value) mode by both measures as malignant, cells residing in the 1st (low-value) mode by both measures as non-malignant, and any remaining cells with conflicting assignments to modes as unassigned.
#' @seealso 
#'  \code{\link[scalop]{unique_sample_names}},\code{\link[scalop]{hms_span}},\code{\link[scalop]{comply}}
#' @rdname findMalignant
#' @export 
#' @importFrom scalop unique_sample_names hms_span comply
findMalignant = function(cna,
                         refCells = NULL,
                         samples = scalop::unique_sample_names(colnames(cna), max.nchar = 6),
                         gene.quantile = 0.9,
                         gene.quantile.for.corr = 0.5,
                         gene.quantile.for.signal = gene.quantile,
                         use.bootstraps = TRUE,
                         n.bootstraps = 10000,
                         gauss.assignment.prob = 0.95,
                         gauss.assignment.coverage = 0.8,
                         verbose = TRUE,
                         plot = TRUE,
                         ...) {

    if (verbose) message("Calculating cells' CNA correlations...")
    cors = cnaCor(cna,
                  threshold = gene.quantile.for.corr,
                  samples = samples,
                  refCells = refCells)

    if (verbose) message("Calculating cells' CNA signal...")
    signals = cnaSignal(cna, threshold = gene.quantile.for.signal)

    if (verbose) {
        message('Fitting CNA correlations to one of two (CNA-low vs. CNA-high) modes...')
    }

    old <- Sys.time()
    invisible(capture.output(corGroups <- suppressMessages(fitBimodal(cors,
                                                                      bySampling = use.bootstraps,
                                                                      nsamp = n.bootstraps,
                                                                      prob = gauss.assignment.prob,
                                                                      coverage = gauss.assignment.coverage,
                                                                      assign = TRUE,
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
                                                                      bySampling = use.bootstraps,
                                                                      nsamp = n.bootstraps,
                                                                      prob = gauss.assignment.prob,
                                                                      coverage = gauss.assignment.coverage,
                                                                      assign = TRUE,
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
                       gene.quantile = gene.quantile,
                       gene.quantile.for.signal,
                       gene.quantile.for.corr = gene.quantile.for.corr,
                       samples = samples,
                       refCells = refCells,
                       groups = result)
    }

    result

}
