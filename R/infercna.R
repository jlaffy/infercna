
diprange = function(cna, ...) {
    dots = list(...)
    v = sapply(dots, function(ref) rowMeans(cna[, ref, drop = F]), simplify = F)
    list(min = do.call(pmin, v), max = do.call(pmax, v))
}

dipcenter = function(v, Min, Max) {
    above = v > Min & v > Max
    below = v < Min & v < Max
    normal = !above & !below
    v[above] <- v[above] - Max
    v[below] <- v[below] - Min
    v[normal] <- 0
    v
}

dipcorrect = function(cna, ...) {
    dots = list(...)
    genes = rownames(cna)
    Args = c(list(cna = cna), dots)
    rg = do.call(diprange, Args)
    n = nrow(cna)
    cna = t(sapply(1:n, function(i) dipcenter(cna[i, ], Min = rg$min[i], Max = rg$max[i])))
    rownames(cna) = genes
    cna
}

infercna2 = function(m,
                     dipcells = NULL,
                     window = 100,
                     range = c(-3, 3)) {

    # note: m should be row(gene)-centered
    cna = infercna::geneOrderMatrix(m)
    cna = infercna::clip(cna, range = range)
    cna = infercna::runMean(cna, k = window)
    cna = infercna::colCenter(cna)
    if (!is.null(dipcells)) {
        Args = c(list(cna = cna), dipcells)
        cna = do.call(infercna::dipcorrect, Args)}
    infercna::geneOrderMatrix(cna)
}

