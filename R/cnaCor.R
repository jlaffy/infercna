
.cnaCor = function(cna, cor.method = 'pearson', gene.quantile = NULL, refCells = NULL, na.replace = NULL) {
    cna = as.matrix(cna)
    if (!is.null(refCells)) {
        genemeans = rowMeans(cna[, !colnames(cna) %in% unlist(refCells), drop = F])
        if (length(genemeans) == 0) {
            message('No cells left after filtering out refCells')
            return(FALSE)
        }
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
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @param refCells a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to group cells by, iii) TRUE to extract sample names from cell ids and subsequently group. Default: NULL
#' @param ... other arguments passed to scalop::unique_sample_names if samples = TRUE.
#' @return a numeric vector or list of numeric vectors
#' @rdname cnaCor
#' @export 
cnaCor = function(cna,
                  cor.method = 'pearson',
                  gene.quantile = NULL,
                  refCells = NULL,
                  samples = NULL,
                  ...) {

    if (is.null(refCells)) tmp = cna
    else tmp = cna[, !colnames(cna) %in% unlist(refCells)]

    if (is.null(samples)) {
        if (!is.null(gene.quantile)) {
            cna = cna[cnaHotspotGenes(tmp, gene.quantile = gene.quantile), ]
        }

        return(.cnaCor(cna, cor.method = cor.method, na.replace = 0, refCells = refCells))
    }

    if (isTRUE(samples)) {
        samples = scalop::unique_sample_names(colnames(cna), ...)
        message('Samples identified:\n', paste0(samples, collapse = '\n'))
    }

    if (is.character(samples)) {
        samples = scalop::split_by_sample_names(colnames(cna), samples = samples)
        samples = samples[lengths(samples) != 0]
    }

    stopifnot(is.list(samples))

    if (!is.null(gene.quantile)) {
        tmplist = sapply(samples, function(cells) tmp[, tmp %in% cells, drop = FALSE], simplify = F)
        genelist = sapply(tmplist, cnaHotspotGenes, gene.quantile = gene.quantile)
        rm(tmplist)
    } else {
        genelist = replicate(length(samples), rownames(cna), simplify = F)
    }

    cnalist = Map(function(m, x, y) m[x, y],
                  x = genelist,
                  y = samples,
                  MoreArgs = list(m = m))
    

    cors = sapply(cnalist,
                  .cnaCor,
                  cor.method = cor.method,
                  refCells = refCells,
                  na.replace = 0,
                  simplify = F)

    cors = scalop::Unlist(cors)
    maincors = cors

    if (length(cors) < ncol(cna)) {
        maincors = .cnaCor(cna, cor.method = cor.method, na.replace = 0)
    # not sure what exception this is supposed to catch...
        maincors[names(cors)] = cors
    }

    maincors[colnames(cna)]
}

