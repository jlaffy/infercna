
#' @title Visualise Malignant and Non-Malignant Subsets
#' @description Visualise Malignant and Non-Malignant Subsets of cells. This is achieved by plotting, for each cell, its CNA signal over its CNA correlation. Please see `infercna::cnaSignal` and `infercna::cnaCor` for details.
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param samples a character vector of sample names list of character vectors (cell/column names). Can be provided if the cna matrix contains multiple samples of cells, i.e. multiple tumours, such that the cell-groups correlations are calcualted for each groups/tumour in turn. Default: FALSE
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @param gene.quantile.for.corr as above but for CNA correlations specifically. Default: gene.quantile 
#' @param gene.quantile.for.signal as above but for CNA signal specifically. Default: gene.quantile
#' @param groups character vector of cell IDs to colour. Default: NULL
#' @param groups.col colour to colour groups by. Default: 'magenta'
#' @param hline y-intercept line. Default: NULL
#' @param vline x-intercept line. Default: NULL
#' @param refCells a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroupss, provide i) a list of cell id groupss, ii) a character vector of sample names to groups cells by, iii) TRUE to extract sample names from cell ids and subsequently groups. Default: NULL
#' @param ... other arguments passed to scalop::unique_sample_names if samples = TRUE.
#' @return a base R plot. If return value is saved to a variable, instead returns data points for cna correlations and cna signal in list form.
#' @rdname cnaScatterPlot
#' @export 
cnaScatterPlot = function(cna,
                          cor.method = 'pearson',
                          gene.quantile = NULL,
                          gene.quantile.for.corr = gene.quantile,
                          gene.quantile.for.signal = gene.quantile,
                          groups = NULL,
                          groups.col = scalop::discrete_colours[1:length(groups)],
                          cex = 0.8,
                          pch = 20,
                          hline = NULL,
                          vline = NULL,
                          refCells = NULL,
                          samples = NULL,
                          ...) {

    groups.col = scales::alpha(groups.col, 0.3)
    cors = cnaCor(cna,
                  threshold = gene.quantile.for.corr,
                  refCells = refCells,
                  samples = samples,
                  ...)

    signals = cnaSignal(cna, threshold = gene.quantile.for.signal)

    plot(cors,
         signals,
         xlab = 'CNA Correlation',
         ylab = 'CNA Signal',
         pch = 1,
         cex = cex, ...)

    if (!is.null(groups)) {
        if (!is.list(groups)) groups = list(groups)
        .Points = function(cors, signals, groups, groups.col) {
            points(cors[groups],
                   signals[groups],
                   pch = pch,
                   col = groups.col,
                   cex = cex)
        }

        Map(.Points,
            groups = groups,
            groups.col = groups.col,
            MoreArgs = list(cors = cors, signals = signals))
    }

    if (!is.null(vline)) abline(v = vline, lty = 2)
    if (!is.null(hline)) abline(h = hline, lty = 2)

    return(invisible(list(cna.cor = cors, cna.signal = signals)))
}
