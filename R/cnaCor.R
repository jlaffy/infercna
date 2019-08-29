
.cnaCor = function(cna, cor.method = 'pearson') {
    cna = as.matrix(cna)
    genemeans = rowMeans(cna)
    cellcors = stats::cor(genemeans, cna, method = cor.method)
    unlist(as.data.frame(cellcors))
}

#' @title Cell - Tumour CNA Correlations
#' @description Compute the pairwise correlations between individual cells' CNA values and the average CNA values in their tumour of origin. 
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param samples a list of character vectors (cell/column names). Can be provided if the cna matrix contains multiple samples of cells, i.e. multiple tumours, such that the cell-group correlations are calcualted for each group/tumour in turn. Default: FALSE
#' @return a numeric vector or list of numeric vectors
#' @rdname cnaCor
#' @export 
cnaCor = function(cna,
                  cor.method = 'pearson',
                  threshold = NULL,
                  samples = FALSE,
                  sep = "-|_") {

    if (!is.null(threshold)) {
        cna = cna[cnaHotspotGenes(cna, threshold = threshold), ]
    }

    if (isFALSE(samples)) {
        return(.cnaCor(cna, cor.method = cor.method))
    }

    if (isTRUE(samples)) {
        samples = .splitSamples(colnames(cna), sep = sep)
    }

    stopifnot(!is.null(samples))

    Call = quote(.cnaCor(cna[, Sample], cor.method = cor.method))
    
    # calculate the correlation of each cell to its tumour
    # and column-bind the resulting tumour-specific dataframes

    Reduce(cbind.data.frame,
           sapply(samples, function(Sample) eval(Call), simplify = F))
}

