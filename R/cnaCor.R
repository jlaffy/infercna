
.cnaCor = function(cna, cor.method = 'pearson') {
    cna = as.matrix(cna)
    genemeans = rowMeans(cna)
    cellcors = cor(genemeans, cna, method = cor.method)
    unlist(as.data.frame(cellcors))
}

#' @title Cell - Tumour CNA Correlations
#' @description Compute the pairwise correlations between individual cells' CNA values and the average CNA values in their tumour of origin. 
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param groups a list of character vectors (cell/column names). Can be provided if the cna matrix contains multiple groups of cells, i.e. multiple tumours, such that the cell-group correlations are calcualted for each group/tumour in turn. Default: FALSE
#' @return a numeric vector or list of numeric vectors
#' @rdname cnaCor
#' @export 
cnaCor = function(cna, cor.method = 'pearson', groups = FALSE) {
    if (isFALSE(groups)) {
        return(.cnaCor(cna, cor.method = cor.method))
    }
    
    group_apply(mat = cna,
                groups = groups,
                ungroup = TRUE,
                FUN = .cnaCor,
                cor.method = cor.method)
}

