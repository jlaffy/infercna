
#' @title Load Example Data
#' @description Load example dataset of human single-cell RNA seq data. Cells and genes are columns and rows, respectively, and all passed QC. Matrix is row-centered.
#' @return centered expression matrix of genes X cells 
#' @rdname useData
#' @export 
useData = function() {
    m = as.matrix(mgh125)
    rowCenter(m)
}
