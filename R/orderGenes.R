.geneMatch = function(genes) {
    match(g.env$gene, genes, nomatch = 0)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname orderGenes 
#' @export 
orderGenes = function(x = NULL) {
    if (is.null(dim(x))) return(x[.geneMatch(genes = x)])
    x[.geneMatch(rownames(x)), , drop = FALSE]
}

