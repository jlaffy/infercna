
.cnaCor = function(cna, cor.method = 'pearson', threshold = NULL, refCells = NULL, na.replace = NULL) {
    cna = as.matrix(cna)
    if (!is.null(refCells)) {
        genemeans = rowMeans(cna[, !colnames(cna) %in% unlist(refCells), drop = F])
    } else {
        genemeans = rowMeans(cna)
    }
    cellcors = suppressWarnings(stats::cor(genemeans, cna, method = cor.method))
    cellcors = unlist(as.data.frame(cellcors))
    if (!is.null(na.replace)) cellcors[is.na(cellcors)] <- na.replace
    cellcors
}

#' @title Cell - Tumour CNA Correlations
#' @description Compute the pairwise correlations between individual cells' CNA values and the average CNA values in their tumour of origin. 
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param refCells a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to group cells by, iii) TRUE to extract sample names from cell ids and subsequently group. Default: NULL
#' @param ... other arguments passed to scalop::unique_sample_names if samples = TRUE.
#' @return a numeric vector or list of numeric vectors
#' @rdname cnaCor
#' @export 
cnaCor = function(cna,
                  cor.method = 'pearson',
                  threshold = NULL,
                  refCells = NULL,
                  samples = NULL,
                  ...) {

    if (!is.null(threshold)) cna = cna[cnaHotspotGenes(cna, threshold = threshold), ]

    if (is.null(samples)) {
        cors = .cnaCor(cna,
                       cor.method = cor.method,
                       refCells = refCells,
                       na.replace = 0)
        return(cors)
    }

    if (isTRUE(samples)) {
        samples = scalop::unique_sample_names(colnames(cna), ...)
        message('Samples identified:\n', paste0(samples, collapse = '\n'))
    }

    stopifnot(is.character(samples))

    cellsBySample = scalop::split_by_sample_names(colnames(cna), samples = samples)

    cnaBySample = sapply(cellsBySample, function(cells) cna[, cells], simplify = F)

    corBySample = sapply(cnaBySample,
                         .cnaCor,
                         cor.method = cor.method,
                         na.replace = 0,
                         refCells = refCells,
                         simplify = F)

    corNames = unlist(sapply(corBySample, names, simplify = F))
    cors = stats::setNames(unlist(corBySample), corNames)
    maincors = cors

    if (length(cors) < ncol(cna)) {
        maincors = .cnaCor(cna, cor.method = cor.method, na.replace = 0)
        maincors[names(cors)] = cors
    }

    maincors[colnames(cna)]
}

