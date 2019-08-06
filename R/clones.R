
clones = function(m = NULL,
                  L = NULL,
                  by = 'arm',
                  prob = 0.95,
                  coverage = 0.8,
                  size = 10,
                  boolean = FALSE) {

    if (is.null(L)) {
        L = genesplit(m = m, by = by)
        L = L[sapply(L, nrow) >= 5]
    }
    modemat = sapply(L, modes, prob = prob, coverage = coverage, size = size, boolean = boolean, simplify = T)
    browser()
    stopifnot(is.matrix(modemat))
    modemat = modemat[rowSums(modemat) != 0, colSums(modemat) != 0, drop = FALSE]
    gr1 = apply(modemat, 1, function(r) paste0(colnames(modemat)[which(r == 1)], collapse = '_'))
    gr2 = apply(modemat, 1, function(r) paste0(colnames(modemat)[which(r == 2)], collapse = '_'))
    cloneIDs = mapply(function(x, y) paste(x, y, sep = '-'), x = gr1, y = gr2)
    split(names(cloneIDs), cloneIDs)
}
