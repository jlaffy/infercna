
.envmatch = function(genes) {
    match(genes, Genv$gene, nomatch = 0)
}

.envglob = function(genes = NULL) {
    L = as.list(Genv)
    L = L[sapply(L, function(l) is.atomic(l) && length(l) > 1)]
    if (!is.null(genes)) L = sapply(L, `[`, .envmatch(genes), simplify = F)
    L
}

#' @title Split Genes By Chromosome (Arm)
#' @description Split Genes By Chromosome (Arm). Input can be one of a matrix to be split into several matrices or a character vector of gene names to be split into several character vectors. 
#' @param x a matrix with gene names as rows or a character vector of gene names.
#' @param by a string; one of 'chr' or 'arm', determining how the matrix or character vector should be split. Default: 'chr'
#' @return if a matrix was provided, a list of matrices. If a character vector was provided, a list of character vectors. Both lists will be of length # of chromosome (arms) in the genome in question. Ie. length 24 for H.sapiens. 
#' @examples 
#' m = infercna::useData()
#' a = splitGenes(x = m)
#' b = splitGenes(x = rownames(m))
#' all(sapply(a, nrow) == lengths(b))
#' names(a)
#' names(splitGenes(x = m, by = 'arm'))
#' @rdname splitGenes
#' @export 
splitGenes = function(x, by = 'chr') {
    if (is.null(dim(x))) genes = x
    else genes = rownames(x)
    L = .envglob(genes = genes)
    stopifnot(by %in% names(L))
    splut = split(L$gene, L[[by]])
    if (is.null(dim(x))) return(splut)
    sapply(splut, function(group) x[group, ], simplify = F)
}
