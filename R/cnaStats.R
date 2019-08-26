
cnaSignal = function(cna) {
    cna = as.matrix(cna)
    sqmat = 2^(cna)
    msq = rowMeans(sqmat)
    msq
}

.cnaCor = function(cna, cor.method = 'pearson') {
    cna = as.matrix(cna)
    genemeans = rowMeans(cna)
    cellcors = cor(genemeans, cna, method = cor.method)
    unlist(as.data.frame(cellcors))
}


cnaCor = function(cna, cor.method, groups = FALSE) {
    if (isFALSE(groups)) {
        return(.cnaCor(cna, cor.method = cor.method))
    }
    
    infercna::group_apply(mat = cna,
                          groups = groups,
                          ungroup = TRUE,
                          FUN = .cnaCor,
                          cor.method = cor.method)
}


cnaStats = function(cna) {
    sig = infercna::cnaSignal(cna)
    cor = infercna::cnaCor(cna)
}

isMalignant = function(cna) {

}
