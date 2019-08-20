
cnasignal = function(cna) {
    cna = as.matrix(cna)
    sqmat = 2^(cna)
    msq = rowMeans(sqmat)
    msq
}

.cnacor = function(cna, cor.method = 'pearson') {
    cna = as.matrix(cna)
    genemeans = rowMeans(cna)
    cellcors = cor(genemeans, cna, method = cor.method)
    unlist(as.data.frame(cellcors))
}


cnacor = function(cna, cor.method, groups = FALSE) {
    if (isFALSE(groups)) {
        return(.cnacor(cna, cor.method = cor.method))
    }
    
    infercna::group_apply(mat = cna,
                          groups = groups,
                          ungroup = TRUE,
                          FUN = .cnacor,
                          cor.method = cor.method)
}


cnastats = function(cna) {
    infercna::cnasignal(cna)
    infercna::cnacor(cna)
}
