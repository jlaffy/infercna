.geneMatch = function(genes) {
    match(Genv$gene, genes, nomatch = 0)
}

#' @title Order Genes by their Genomic Positions
#' @description Order Genes by their Genomic Positions
#' @param x a matrix with gene names as rows or a character vector of gene names.
#' @return ordered matrix or character vector
#' @examples 
#' m = infercna::useData()
#' orderGenes(m) %>% rownames %>% head
#' orderGenes(rownames(m)) %>% head
#' @rdname orderGenes 
#' @export 
orderGenes = function(x) {
    if (is.null(dim(x))) return(x[.geneMatch(genes = x)])
    x[.geneMatch(rownames(x)), , drop = FALSE]
}

