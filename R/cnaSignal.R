
#' @title Calculate the Means of Absolute CNA Values
#' @description Calculate the Mean of Absolute CNA Values
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @return a numeric vector of CNA signal values or the Mean of Absolute CNA values
#' @rdname cnaSignal
#' @export 
cnaSignal = function(cna) {
    cna = as.matrix(cna)
    sqmat = abs(cna)
    msq = colMeans(sqmat)
    msq
}
