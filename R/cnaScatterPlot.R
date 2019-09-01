
#' @title Visualise Malignant and Non-Malignant Subsets
#' @description Visualise Malignant and Non-Malignant Subsets of cells. This is achieved by plotting, for each cell, its CNA signal over its CNA correlation. Please see `infercna::cnaSignal` and `infercna::cnaCor` for details.
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param samples a character vector of sample names list of character vectors (cell/column names). Can be provided if the cna matrix contains multiple samples of cells, i.e. multiple tumours, such that the cell-group correlations are calcualted for each group/tumour in turn. Default: FALSE
#' @param threshold PARAM_DESCRIPTION, Default: NULL
#' @param cor.threshold PARAM_DESCRIPTION, Default: threshold
#' @param signal.threshold PARAM_DESCRIPTION, Default: NULL
#' @param group PARAM_DESCRIPTION, Default: NULL
#' @param group.col PARAM_DESCRIPTION, Default: 'magenta'
#' @param hline PARAM_DESCRIPTION, Default: NULL
#' @param vline PARAM_DESCRIPTION, Default: NULL
#' @param bySample PARAM_DESCRIPTION, Default: FALSE
#' @param samples PARAM_DESCRIPTION, Default: NULL
#' @param sep if bySample is TRUE and samples are NULL, split cell IDs by <sep> and take the first substring to be the sample name. Default: '-|_'
#' @param excludeFromAvg a character vector of cell IDs to exclude from the average tumour profile. Default: NULL
#' @param ... other arguments passed to base plot
#' @return a base R plot
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname cnaScatterPlot
#' @export 
cnaScatterPlot = function(cna,
                          cor.method = 'pearson',
                          threshold = NULL,
                          cor.threshold = threshold,
                          signal.threshold = threshold,
                          group = NULL,
                          group.col = 'magenta',
                          cex = 0.5,
                          hline = NULL,
                          vline = NULL,
                          bySample = FALSE,
                          samples = NULL,
                          sep = "-|_",
                          excludeFromAvg = NULL,
                          ...) {

    cors = cnaCor(cna,
                  threshold = cor.threshold,
                  bySample = bySample,
                  samples = samples,
                  sep = sep,
                  excludeFromAvg = excludeFromAvg)

    signals = cnaSignal(cna, threshold = signal.threshold)

    plot(cors, signals, xlab = 'CNA Correlation', ylab = 'CNA Signal', pch = 1, cex = cex, ...)

    if (!is.null(group)) {
        points(cors[group], signals[group], pch = 20, col = group.col, cex = cex)
    }

    if (!is.null(vline)) abline(v = vline, lty = 2)
    if (!is.null(hline)) abline(h = hline, lty = 2)
}
