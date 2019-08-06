
.diprange = function(cna, ...) {
    dots = list(...)
    v = sapply(dots, function(ref) rowMeans(cna[, ref, drop = F]), simplify = F)
    list(min = do.call(pmin, v), max = do.call(pmax, v))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param v PARAM_DESCRIPTION
#' @param Min PARAM_DESCRIPTION
#' @param Max PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname dipcenter
#' @export 
dipcenter = function(v, Min, Max) {
    above = v > Min & v > Max
    below = v < Min & v < Max
    normal = !above & !below
    v[above] <- v[above] - Max
    v[below] <- v[below] - Min
    v[normal] <- 0
    v
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cna PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details Correction with diploid cells' <dipcells> CNAs: the boundaries of dipcells CNA values are the boundaries for what should be considered a CNA of 0. Thus, if the boundary is -0.1 and 0.1, then a cell with CNA = -0.5 will be corrected to -0.4 and a cell with CNA value of 1 will be corrected to 0.9.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname dipcorrect
#' @export 
dipcorrect = function(cna, ...) {
    dots = list(...)
    genes = rownames(cna)
    Args = c(list(cna = cna), dots)
    c(Min, Max) %<-% do.call(.diprange, Args)
    n = nrow(cna)
    cna = t(sapply(1:n, function(i) dipcenter(cna[i, ], Min = Min[i], Max = Max[i])))
    rownames(cna) = genes
    cna
}

#' @title Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @description Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @param m a matrix of genes X cells RNA-seq expression data. The matrix should be row-centered.
#' @param dipcells a list of two or more character vectors, each containg cell IDs of normal cell types (one cell type per list element). Since these cell types are presumed diploid cells, their CNA values can be used to correct the remaining CNA values. Note that all cell IDs should correspond to column names in <m>. See infercna::dipcells (diploid cells) for an example. Default: NULL
#' @param window the size of the window to use for the rolling mean. Units are in number of genes. Default: 100
#' @param range values in <m> above and below this range will be set to values in range. Default: c(-3, 3)
#' @return a matrix of genes X cells of inferred CNA values. Note that n = (window - 1)/2 genes will be lost from either extremity of the genome (ie. n genes lost at the start of chromosome 1 and at the end of chromosome Y, if the genome in question is H.sapiens.)
#' @details Correction with diploid cells' <dipcells> CNAs: the boundaries of dipcells CNA values are the boundaries for what should be considered a CNA of 0. Thus, if the boundary is -0.1 and 0.1, then a cell with CNA = -0.5 will be corrected to -0.4 and a cell with CNA value of 1 will be corrected to 0.9.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  dipcells = infercna::dipcells
#'  m = infercna::useData()
#'  cna = infercna::infercna(m = m, dipcells = dipcells)
#'  }
#' }
#' @rdname infercna
#' @export 
infercna = function(m,
                    dipcells = NULL,
                    window = 100,
                    range = c(-3, 3)) {

    # note: m should be row(gene)-centered
    cna = genorder(m)
    cna = clip(cna, range = range)
    cna = runMean(cna, k = window)
    cna = colCenter(cna)
    if (!is.null(dipcells)) {
        Args = c(list(cna = cna), dipcells)
        cna = do.call(dipcorrect, Args)}
    cna
}

