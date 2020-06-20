#' @title Calculate the Means of Squared CNA Values
#' @description Calculate the Mean of Squared CNA Values
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @param refCells a character vector of cell ids (e.g. normal reference cell ids) to exclude from calculation of cnaHotspotGenes. Only relevant if gene.quantile is not NULL. Default: NULL
#' @param samples if cnaHotspotGenes should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to group cells by, iii) TRUE to extract sample names from cell ids and subsequently group. Default: NULL
#' @param ... other arguments passed to scalop::unique_sample_names if samples = TRUE.
#' @return a numeric vector of CNA signal values or the Mean of Squared CNA values
#' @rdname cnaSignal
#' @export 
cnaSignal = function(cna, gene.quantile = NULL, refCells = NULL, samples = NULL) {
    cna = as.matrix(cna)

    .cnaSignal = function(cna) {
        sqmat = cna^2
        colMeans(sqmat)
    }

    if (is.null(gene.quantile)) {
        return(.cnaSignal(cna))
    }

    if (is.null(refCells)) tmp = cna
    else tmp = cna[, !colnames(cna) %in% unlist(refCells)]

    if (is.null(samples)) {
        genes = cnaHotspotGenes(tmp, gene.quantile = gene.quantile)
        return(.cnaSignal(cna[genes, ]))
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
                  MoreArgs = list(m = cna))

    scalop::Unlist(sapply(cnalist, .cnaSignal, simplify = F), nested.names = T)
}

