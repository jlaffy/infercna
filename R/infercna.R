
.mostExpressedGenes = function(m, ngenes) {
    if (ngenes >= 0 & ngenes <= 1) {
        ngenes = ngenes * nrow(m)
    } 

    else if (nrow(m) <= ngenes) {
        ngenes = nrow(m)
    }

    else stop('<ngenes> must be either a value no larger than the number of genes in the matrix <m>, or a fraction between 0 and 1.')

    rom = rowMeans(m)
    names(sort(rowMeans(m), decreasing = T)[1:ngenes])
}

#' @title Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @description Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @param m a matrix of genes X cells containing scRNA-seq expression data. The matrix should NOT be row-centered.
#' @param refCells a list of two or more character vectors, each containg cell IDs of normal cell types (one cell type per list element). Since these cell types are presumed normal cells, their CNA values can be used to correct the remaining CNA values. Note that all cell IDs should correspond to column names in <m>. See infercna::refCells (normal cells) for an example. Default: NULL
#' @param window the size of the window to use for the rolling mean. Units are in number of genes. Default: 100
#' @param range values in <m> above and below this range will be set to values in range. See infercna::clip() for more details. Default: c(-3, 3)
#' @param noise a numeric value, which if given, increases the boundaries within which CNA values are considered 0. Increases by <noise> i.e. the bounds become Minimum - noise and Maximum + noise. Default: 0.1
#' @param n a numeric value indicating if only top genes should be used in the CNA calculation and how many. if n is a fraction between 0 and 1, the number of genes included will n * nrow(m). Else the number of genes included will be n. Default: NULL
#' @return a matrix of genes X cells of inferred CNA values. Note that n = (window - 1)/2 genes will be lost from either extremity of the genome (ie. n genes lost at the start of chromosome 1 and at the end of chromosome Y, if the genome in question is H.sapiens.)
#' @details Correction with reference cells' <refCells> CNAs: the boundaries of reference CNA values are the boundaries for what should be considered a CNA of 0. Thus, if the boundary is -0.1 and 0.1, then a cell with CNA = -0.5 will be corrected to -0.4, a cell with CNA value of 1 will be corrected to 0.9 and a cell with CNA value of 0.05 will be corrected to 0.
#' @examples 
#'  m = infercna::useData()
#'  cna = infercna::infercna(m = m, refCells = infercna::refCells)
#' @rdname infercna
#' @export 
infercna = function(m,
                    refCells = NULL,
                    window = 100,
                    range = c(-3, 3),
                    n = NULL,
                    noise = 0.1,
                    center.method = 'median',
                    isLog = FALSE,
                    verbose = TRUE) {

    # note: m should NOT be row(gene)-centered
    if (all(round(range(rowMeans(m)), 3) == 0)) {
        stop('Matrix is row-centered. Please provide non-centered data.')
    }
    if (isLog) m = TPM(m)
    if (!is.null(n)) m = m[.mostExpressedGenes(m, ngenes = n), ]

    m = logTPM(m)
    m = rowCenter(m)
    m = orderGenes(m)
    m = clip(m, range = range)

    m = TPM(m) 
    ms = splitGenes(m, by = 'chr')
    cna = sapply(ms, function(m) try(runMean(m, k = window, verbose = verbose)), simplify = F)
    cna = cna[sapply(cna, class) != 'try-error' | !sapply(class, isFALSE)]
    cna = Reduce(rbind, cna)

    cna = logTPM(cna, dividebyten = T)
    cna = colCenter(cna, method = center.method) # note: using median centering here
    if (!is.null(refCells)) {
        Args = c(list(cna = cna, noise = noise, isLog = TRUE), refCells)
        cna = do.call(refCorrect, Args)}
    cna
}

