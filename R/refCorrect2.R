
.refrange2 = function(cna, ...) {
    dots = list(...)
    v = sapply(dots, function(ref) log2(rowMeans(2^(cna[, ref, drop = F])-1)+1), simplify = F)
    res = list(Min = do.call(pmin, v), Max = do.call(pmax, v))
    sapply(res, logtpm, bulk = T, simplify = F)
}

.refcenter = function(v, Min, Max, noise = NULL) {
#    if (!is.null(noise)) {
#        Min = Min - noise
#        Max = Max + noise
#    }
    
    if (is.null(noise)) {
	    noise = 0
    }

    above = v > (Max + noise)
    below = v < (Min - noise)
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
#' @rdname refCorrect2
#' @export 
refCorrect2 = function(cna, noise = NULL, ...) {
    dots = list(...)
    genes = rownames(cna)
    Args = c(list(cna = cna), dots)
    rg = do.call(.refrange2, Args)
    n = nrow(cna)
    cna = t(sapply(1:n, function(i) {
                       .refcenter(cna[i, ],
                                  Min = rg$Min[i],
                                  Max = rg$Max[i],
                                  noise = noise)}))
    rownames(cna) = genes
    cna
}

