
#' @title Filter Genes by their Genome Features
#' @description Filter genes to keep (the default) or to throw according to genomic features. For example, filter genes that are either on chromosome 7 or chromosome arm 2p.
#' @param x a character vector of gene names or a matrix with gene row names.
#' @param value a value or a list of values that are the genome features to filter
#' @param attribute the name(s) of the genome feature(s) in question. One or more of of 'chr', 'arm', 'start', 'end'.
#' @param throw boolean value indicating whether the filtered genes should be thrown away rather than kept. Default: F
#' @return filtered vector or a matrix with filtered rows, depending on class of <x>.
#' @examples 
#'  m = useData()
#'  filterGenes(rownames(m), value = c(7, '2p'), attribute = c('chr', 'arm'))
#'  filterGenes(m, value = c(7, '2p'), attribute = c('chr', 'arm'))
#'  filterGenes(m, value = c(7, '2p', 3q'), attribute = c('chr', 'arm', 'arm')
#'  filterGenes(m, value = list(7, c('2p', '3q')), attribute = c('chr', 'arm'))
#' @rdname filterGenes
#' @export 
filterGenes = function(x, value, attribute, throw = F) {
    
    if (!is.null(dim(x))) {
        mat = x
        x = rownames(mat)
    }

    bool = x %in% genesOn(value = value, attribute = attribute)
    if (throw) bool = !bool
    if (exists(mat, inherits = F)) return(mat[bool, ])
    x[bool]
}



#' @title Retrieve Genes by their Genome Features
#' @description Retrieve genes by their genome features. e.g. all genes on chromosome 7. 
#' @param value a value or a list of values that are the genome features to filter
#' @param attribute the name(s) of the genome feature(s) in question. One or more of of 'chr', 'arm', 'start', 'end'.
#' @param throw boolean value indicating whether the filtered genes should be thrown away rather than kept. Default: F
#' @return a character vector of gene names 
#' @examples 
#' genesOn(7, 'chr')
#' @rdname filterGenes
#' @export 
genesOn = function(value, attribute, x = NULL) {

    .genesOn = function(x, value, attribute) {
        attribute = get(envir = Genv, x = attribute)
        names(attribute)[attribute %in% value]
    }

    if (length(attribute) > 1) {
        genelist = Map(.genesOn, value = value, attribute = attribute)
        genes = Reduce(union, genelist)
    } else {
        genes = .genesOn(value = value, attribute = attribute)
    }

    if (!is.null(x) && !is.null(dim(x))) {
        bool = rownames(x) %in% genes
        return(x[bool, ])
    }

    else if (!is.null(x)) {
        bool = x %in% genes
        return(x[bool])
    }

    genes
}
