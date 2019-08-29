
.HQgenes = function(m, n) {
    if (n >= 0 & n <= 1) {
        n = n * nrow(m)
    }

    if (nrow(m) < n) {
        n = nrow(m)
    }

    rom = rowMeans(m)
    hq = names(sort(rowMeans(m), decreasing = T)[1:n])
    m[rownames(m) %in% hq, ]
}

#' @title Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @description Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @param m a matrix of genes X cells RNA-seq expression data. The matrix should be row-centered.
#' @param reference a list of two or more character vectors, each containg cell IDs of normal cell types (one cell type per list element). Since these cell types are presumed refloid cells, their CNA values can be used to correct the remaining CNA values. Note that all cell IDs should correspond to column names in <m>. See infercna::reference (refloid cells) for an example. Default: NULL
#' @param window the size of the window to use for the rolling mean. Units are in number of genes. Default: 100
#' @param range values in <m> above and below this range will be set to values in range. See infercna::clip() for more details. Default: c(-3, 3)
#' @param noise a numeric value, which if given, increases the boundaries within which CNA values are considered 0. Increases by <noise> i.e. the bounds become Minimum - noise and Maximum + noise. Default: 0.2
#' @param n a numeric value indicating if only top genes should be used in the CNA calculation and how many. if n is a fraction between 0 and 1, the number of genes included will n * nrow(m). Else the number of genes included will be n. Default: 5000
#' @return a matrix of genes X cells of inferred CNA values. Note that n = (window - 1)/2 genes will be lost from either extremity of the genome (ie. n genes lost at the start of chromosome 1 and at the end of chromosome Y, if the genome in question is H.sapiens.)
#' @details Correction with reference cells' <reference> CNAs: the boundaries of reference CNA values are the boundaries for what should be considered a CNA of 0. Thus, if the boundary is -0.1 and 0.1, then a cell with CNA = -0.5 will be corrected to -0.4, a cell with CNA value of 1 will be corrected to 0.9 and a cell with CNA value of 0.05 will be corrected to 0.
#' @examples 
#'  m = infercna::useData()
#'  cna = infercna::infercna(m = m, reference = infercna::reference)
#' @rdname infercna
#' @export 
infercna = function(m,
                    reference = NULL,
                    window = 100,
                    range = c(-3, 3),
                    n = 5000,
                    noise = 0.2) {

    # note: m should be row(gene)-centered
    if (!is.null(n)) {
        m = .HQgenes(m, n = n)
    }
    cna = orderGenes(m)
    cna = clip(cna, range = range)
    cna = runMean(cna, k = window)
    cna = colCenter(cna)
    if (!is.null(reference)) {
        Args = c(list(cna = cna, noise = noise), reference)
        cna = do.call(toAbsolute, Args)}
    cna
}

