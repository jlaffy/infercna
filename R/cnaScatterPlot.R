
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
#' @param refCells a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to group cells by, iii) TRUE to extract sample names from cell ids and subsequently group. Default: NULL
#' @param ... other arguments passed to scalop::unique_sample_names if samples = TRUE.
#' @return a base R plot. If return value is saved to a variable, instead returns data points for cna correlations and cna signal in list form.
#' @rdname cnaScatterPlot
#' @export 
cnaScatterPlot = function(cna,
                          cor.method = 'pearson',
                          threshold = NULL,
                          cor.threshold = threshold,
                          signal.threshold = threshold,
                          group = NULL,
                          group.col = scalop::discrete_colours[1:length(group)],
                          cex = 0.8,
                          pch = 20,
                          hline = NULL,
                          vline = NULL,
                          refCells = NULL,
                          samples = NULL,
                          ...) {

    group.cols = scales::alpha(group.cols, 0.3)
    cors = cnaCor(cna,
                  threshold = cor.threshold,
                  refCells = refCells,
                  samples = samples,
                  ...)

    signals = cnaSignal(cna, threshold = signal.threshold)

    plot(cors,
         signals,
         xlab = 'CNA Correlation',
         ylab = 'CNA Signal',
         pch = 1,
         cex = cex, ...)

    if (!is.null(group)) {
        if (!is.list(group)) group = list(group)
        .Points = function(cors, signals, group, group.col) {
            points(cors[group],
                   signals[group],
                   pch = pch,
                   col = group.col,
                   cex = cex)
        }

        Map(.Points,
            group = group,
            group.col = group.col,
            MoreArgs = list(cors = cors, signals = signals))
    }

    if (!is.null(vline)) abline(v = vline, lty = 2)
    if (!is.null(hline)) abline(h = hline, lty = 2)

    return(invisible(list(cna.cor = cors, cna.signal = signals)))
}
