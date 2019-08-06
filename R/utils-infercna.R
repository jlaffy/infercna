
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param genes PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname genmatch
#' @export 
genmatch = function(genes) {
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
#' @rdname as_genomic
#' @export 
as_genomic = function(x = NULL) {
    if (is.null(dim(x))) return(x[genmatch(genes = x)])
    x[genmatch(rownames(x)), , drop = FALSE]
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param m PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname colCenter
#' @export 
colCenter = function(m) {
    scale(m, center = T, scale = F)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param range PARAM_DESCRIPTION, Default: c(-3, 3)
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname clip
#' @export 
clip <- function(mat, range = c(-3, 3)) {
    mat = as.matrix(mat)
    mat[mat < range[[1]]] <- range[[1]]
    mat[mat > range[[2]]] <- range[[2]]
    mat
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 100
#' @param endrule PARAM_DESCRIPTION, Default: 'mean'
#' @param align PARAM_DESCRIPTION, Default: 'center'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[caTools]{runmean}}
#' @rdname runMean
#' @export 
#' @importFrom caTools runmean
runMean <- function(mat,
                    k = 100,
                    endrule = 'mean',
                    align = 'center') {

    mat = as.matrix(mat)
    mat.out = caTools::runmean(mat,
                               k = k,
                               endrule = endrule,
                               align = align)
    colnames(mat.out) <- colnames(mat)
    rownames(mat.out) <- rownames(mat)
    as.data.frame(mat.out)
}

