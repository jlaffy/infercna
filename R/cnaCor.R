
.cnaCor = function(cna, cor.method = 'pearson', threshold = NULL, excludeFromAvg = NULL) {
    cna = as.matrix(cna)
    if (!is.null(excludeFromAvg)) {
        genemeans = rowMeans(cna[, !colnames(cna) %in% unlist(excludeFromAvg)])
    } else {
        genemeans = rowMeans(cna)
    }
    cellcors = suppressWarnings(stats::cor(genemeans, cna, method = cor.method))
    unlist(as.data.frame(cellcors))
}

#' @title Cell - Tumour CNA Correlations
#' @description Compute the pairwise correlations between individual cells' CNA values and the average CNA values in their tumour of origin. 
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param samples a character vector of sample names list of character vectors (cell/column names). Can be provided if the cna matrix contains multiple samples of cells, i.e. multiple tumours, such that the cell-group correlations are calcualted for each group/tumour in turn. Default: FALSE
#' @param sep if bySample is TRUE and samples are NULL, split cell IDs by <sep> and take the first substring to be the sample name. Default: '-|_'
#' @return a numeric vector or list of numeric vectors
#' @rdname cnaCor
#' @export 
cnaCor = function(cna,
                  cor.method = 'pearson',
                  threshold = NULL,
                  bySample = F,
                  samples = NULL,
                  sep = "-|_",
                  excludeFromAvg = NULL) {

    if (!is.null(samples)) {
        bySample = TRUE
    }

    if (!is.null(threshold)) {
        cna = cna[cnaHotspotGenes(cna, threshold = threshold), ]
    }

    if (!bySample) {
        return(.cnaCor(cna, cor.method = cor.method))
    }

    cnaBySample = splitCellsBySample(x = cna,
                                     sep = sep,
                                     samples = samples)

    corBySample = sapply(cnaBySample,
                         .cnaCor,
                         cor.method = cor.method,
                         excludeFromAvg = excludeFromAvg,
                         simplify = F)

    corNames = unlist(sapply(corBySample, names, simplify = F))
    cors = stats::setNames(unlist(corBySample), corNames)
    cors[is.na(cors)] = 0

    if (length(cors) == ncol(cna)) {
        maincors = cors
    }

    else if (length(cors) < ncol(cna)) {
        maincors = .cnaCor(cna, cor.method = cor.method)
        maincors[is.na(cors)] = 0
        maincors[names(cors)] = cors
    }

    maincors[colnames(cna)]
}

