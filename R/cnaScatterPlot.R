
#' @title Visualise Malignant and Non-Malignant Subsets
#' @description Visualise Malignant and Non-Malignant Subsets of cells. This is achieved by plotting, for each cell, its CNA signal over its CNA correlation. Please see `infercna::cnaSignal` and `infercna::cnaCor` for details.
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param samples a character vector of sample names list of character vectors (cell/column names). Can be provided if the cna matrix contains multiple samples of cells, i.e. multiple tumours, such that the cell-group correlations are calcualted for each group/tumour in turn. Default: FALSE
#' @param threshold  Default: NULL
#' @param cor.threshold PARAM_DESCRIPTION, Default: threshold
#' @param signal.threshold PARAM_DESCRIPTION, Default: NULL
#' @param group PARAM_DESCRIPTION, Default: NULL
#' @param group.col PARAM_DESCRIPTION, Default: 'magenta'
#' @param hline PARAM_DESCRIPTION, Default: NULL
#' @param vline PARAM_DESCRIPTION, Default: NULL
#' @param excludeFromAvg a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to group cells by, iii) TRUE to extract sample names from cell ids and subsequently group. Default: NULL
#' @param ... other arguments passed to scalop::get_sample_names if samples = TRUE.
#' @return a base R plot. If return value is saved to a variable, instead returns data points for cna correlations and cna signal in list form.
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
                          excludeFromAvg = NULL,
                          samples = NULL,
                          ...) {

    cors = cnaCor(cna,
                  threshold = cor.threshold,
                  excludeFromAvg = excludeFromAvg,
                  samples = samples,
                  ...)

    signals = cnaSignal(cna, threshold = signal.threshold)

    invisible(list(cna.cor = cors, cna.signal = signals))

    plot(cors, signals, xlab = 'CNA Correlation', ylab = 'CNA Signal', pch = 1, cex = cex, ...)

    if (!is.null(group)) {
        points(cors[group], signals[group], pch = 20, col = group.col, cex = cex)
    }

    if (!is.null(vline)) abline(v = vline, lty = 2)
    if (!is.null(hline)) abline(h = hline, lty = 2)
}
