
.envmatch = function(genes) {
    match(genes, g.env$gene, nomatch = 0)
}

.envglob = function(genes = NULL) {
    L = as.list(g.env)
    L = L[sapply(L, function(l) is.atomic(l) && length(l) > 1)]
    if (!is.null(genes)) L = sapply(L, `[`, .envmatch(genes), simplify = F)
    L
}

#' @title Split Genes By Chromosome (Arm)
#' @description Split Genes By Chromosome (Arm). Input can be one of a matrix to be split into several matrices or a character vector of gene names to be split into several character vectors. 
#' @param m a matrix to be row-split by chromosome (arm). The matrix must have gene names as row names. Default: NULL
#' @param genes a character vector of gene names. Default: NULL
#' @param by a string; one of 'chr' or 'arm', determining how the matrix or character vector should be split. Default: 'chr'
#' @return if a matrix <m> was provided, a list of matrices. If a character vector <genes> was provided, a list of character vectors. Both lists will be of length # of chromosome (arms) in the genome in question. Ie. length 24 for H.sapiens. 
#' @rdname splitGenes
#' @export 
splitGenes = function(m = NULL, genes = NULL, by = 'chr') {
    if (!is.null(m)) genes = rownames(m)
    L = .envglob(genes = genes)
    stopifnot(by %in% names(L))
    splut = split(L$gene, L[[by]])
    if (!is.null(m)) {
        splut = sapply(splut, function(group) m[group, ], simplify = F)}
    splut
}
