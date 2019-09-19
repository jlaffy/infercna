
#' @title Load Example Data
#' @description Load example dataset of human single-cell RNA seq data. Cells and genes are columns and rows, respectively, and all passed QC. Matrix is row-centered.
#' @param rowcenter center the expression matrix row-wise. Default: F
#' @return centered expression matrix of genes X cells 
#' @rdname useData
#' @export 
useData = function(x = NULL, rowcenter = F) {
    if (is.null(x)) x = cbind(as.matrix(bt771), as.matrix(mgh125))
    else x = as.matrix(x)
    x = x[, !duplicated(colnames(x))]
    if (rowcenter) x = rowCenter(x)
    x
}
