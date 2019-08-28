
#' @title Filter Genes by their Genome Features
#' @description Filter genes to keep (the default) or to out according to genomic features. For example, filter genes that are either on chromosome 7 or chromosome arm 2p.
#' @param x a character vector of gene names or a matrix with gene row names.
#' @param value a value or a list of values that are the genome features to filter
#' @param attribute the name(s) of the genome feature(s) in question. One of 'chr', 'arm'.
#' @param out boolean value indicating whether the filtered genes should be thrown rather than kept. Default: F
#' @return filtered vector or a matrix with filtered rows, depending on class of <x>.
#' @examples 
#'  m = useData()
#'  filterGenes(m, 7, 'chr')
#'  filterGenes(rownames(m), 7, 'chr')
#'  filterGenes(m, list(7, c('2p', '3q')), c('chr', 'arm'))
#' @rdname filterGenes
#' @export 
filterGenes = function(x, value, attribute, out = F) {
    if (!is.null(dim(x))) genes = rownames(x)
    else genes = x
    bool = genes %in% genesOn(value = value, attribute = attribute)
    if (out) bool = !bool
    if (!is.null(dim(x))) return(x[bool, ])
    x[bool]
}


#' @title Retrieve Genes by their Genome Features
#' @description Retrieve genes by their genome features. e.g. all genes on chromosome 7. 
#' @param value a value or a list of values that are the genome features to filter
#' @param attribute the name(s) of the genome feature(s) in question. One or more of of 'chr', 'arm', 'start', 'end'.
#' @return a character vector of gene names 
#' @examples 
#' genesOn(7, 'chr')
#' @rdname genesOn 
#' @export 
genesOn = function(value, attribute) {

    .genesOn = function(x, value, attribute) {
        attribute = get(envir = Genv, x = attribute)
        names(attribute)[attribute %in% value]
    }

    if (length(attribute) > 1) {
        genelist = Map(.genesOn, value = value, attribute = attribute)
        genes = Reduce(union, genelist)
    } 
    
    else {
        genes = .genesOn(value = value, attribute = attribute)
    }

    genes
}
