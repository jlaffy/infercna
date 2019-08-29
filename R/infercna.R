
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


.refrange = function(cna, ...) {
    dots = list(...)
    v = sapply(dots, function(ref) rowMeans(cna[, ref, drop = F]), simplify = F)
    list(Min = do.call(pmin, v), Max = do.call(pmax, v))
}

.refcenter = function(v, Min, Max, noise = NULL) {
    if (!is.null(noise)) {
        Min = Min - noise
        Max = Max + noise
    }
    above = v > Min & v > Max
    below = v < Min & v < Max
    normal = !above & !below
    v[above] <- v[above] - Max
    v[below] <- v[below] - Min
    v[normal] <- 0
    v
}

#' @title Convert Relative CNA Values To Absolute
#' @description If the identities of normal cells are known, the expected CNA values of these cells should be 0. Thus, CNA values within the range of normal cell CNA values are corrected to 0; CNA values below the range have the minimum substracted; CNA values above the range have the maximum subtracted.
#' @param cna a matrix of rows (genes) by columns (cells) of CNA values. 
#' @param noise a numeric value, which if given, increases the boundaries within which CNA values are considered 0. Increases by <noise> i.e. the bounds become Minimum - noise and Maximum + noise. Default: NULL
#' @param ... normal cell column IDs. Expects at least two character vectors. 
#' @rdname refCorrect
#' @export 
refCorrect = function(cna, noise = NULL, ...) {
    dots = list(...)
    genes = rownames(cna)
    Args = c(list(cna = cna), dots)
    rg = do.call(.refrange, Args)
    n = nrow(cna)
    cna = t(sapply(1:n, function(i) {
                       .refcenter(cna[i, ],
                                  Min = rg$Min[i],
                                  Max = rg$Max[i],
                                  noise = noise)}))
    rownames(cna) = genes
    cna
}

#' @title Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @description Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @param m a matrix of genes X cells RNA-seq expression data. The matrix should be row-centered.
#' @param reference a list of two or more character vectors, each containg cell IDs of normal cell types (one cell type per list element). Since these cell types are presumed refloid cells, their CNA values can be used to correct the remaining CNA values. Note that all cell IDs should correspond to column names in <m>. See infercna::reference (refloid cells) for an example. Default: NULL
#' @param window the size of the window to use for the rolling mean. Units are in number of genes. Default: 100
#' @param range values in <m> above and below this range will be set to values in range. See infercna::clip() for more details. Default: c(-3, 3)
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
                    n = NULL,
                    noise = 0.1) {

    # note: m should be row(gene)-centered
    if (!is.null(nGenes)) {
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

