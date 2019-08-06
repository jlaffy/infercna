
geneOrderMatrix <- function(mat = NULL) {
    genes = infercna::geneOrder(genes = rownames(mat))
    mat[genes, , drop = FALSE]
}

geneOrder = function(genes = NULL) {
    gene.order = infercna:::g.env$genes
    if (is.null(genes)) return(gene.order)
#    missing = genes[!genes %in% gene.order]
#    pasted = paste(missing,collapse = '\n')
#    if (length(missing) > 0) {
#        warning('\nRemoving ', length(missing), ' missing genes:\n', pasted)
#    }
    gene.order[gene.order %in% genes]
}

geneFilter = function(genes = NULL,
                      chr = NULL,
                      arm = NULL,
                      prefix = NULL,
                      split = FALSE) {

    genes = geneOrder(genes = genes)

    if (!is.null(prefix)) {
        genes = genes[startsWith(genes, prefix)]
    }
    if (!is.null(chr)) {
        bool = infercna:::g.env$chrm %in% genes & names(infercna:::g.env$chrm) %in% chr
        genes = infercna:::g.env$chrm[bool]
    }
    if (!is.null(arm)) {
        bool = infercna:::g.env$arm %in% genes & names(infercna:::g.env$arm) %in% arm
        genes = infercna:::g.env$arm[bool]
    }
    if (split && length(unique(names(genes))) > 1) {
        genes = split(genes, names(genes))
    }

    genes
}

geneSplit = function(genes = NULL, by = 'arm') {
    by = tolower(by)
    by = stringr::str_replace(by, "omosome", "")
    stopifnot(by %in% c('arm', 'chr'))
    if (by == 'arm') {
        g = g.env$arm
    }
    else if (by == 'chr') {
        g = g.env$chrm
    }
    if (!is.null(genes)) {
        g = g[g %in% genes]
    }
    split(g, names(g))
}



splitMatrixByChr = function(mat) {
    groups = infercna::geneSplit(genes = rownames(mat), by = 'chr')
    sapply(groups, function(g) mat[g, , drop = FALSE], simplify = F)
}

splitMatrixByArm = function(mat) {
    groups = infercna::geneSplit(genes = rownames(mat), by = 'arm')
    sapply(groups, function(g) mat[g, , drop = FALSE], simplify = F)
}

