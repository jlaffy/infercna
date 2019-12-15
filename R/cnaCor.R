
.cnaCor = function(cna, cor.method = 'pearson', threshold = NULL, excludeFromAvg = NULL, na.replace = NULL) {
    cna = as.matrix(cna)
    if (!is.null(excludeFromAvg)) {
        genemeans = rowMeans(cna[, !colnames(cna) %in% unlist(excludeFromAvg), drop = F])
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
#' @param excludeFromAvg a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to group cells by, iii) TRUE to extract sample names from cell ids and subsequently group. Default: NULL
#' @param ... other arguments passed to scalop::get_sample_names if samples = TRUE.
#' @return a numeric vector or list of numeric vectors
#' @rdname cnaCor
#' @export 
cnaCor = function(cna,
                  cor.method = 'pearson',
                  threshold = NULL,
                  excludeFromAvg = NULL,
                  samples = NULL,
                  ...) {

    if (!is.null(threshold)) cna = cna[cnaHotspotGenes(cna, threshold = threshold), ]

    if (is.null(samples)) {
        cors = .cnaCor(cna,
                       cor.method = cor.method,
                       excludeFromAvg = excludeFromAvg,
                       na.replace = 0)
        return(cors)
    }

    if (isTRUE(samples)) {
        samples = scalop::get_sample_names(colnames(cna), ...)
        message('Samples identified:\n', paste0(samples, collapse = '\n'))
    }

    stopifnot(is.character(samples))

    cellsBySample = scalop::split_by_sample(colnames(cna), samples = samples)

    cnaBySample = sapply(cellsBySample, function(cells) cna[, cells], simplify = F)

    corBySample = sapply(cnaBySample,
                         .cnaCor,
                         cor.method = cor.method,
                         na.replace = 0,
                         excludeFromAvg = excludeFromAvg,
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

