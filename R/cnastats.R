
cnasignal = function(cnamat) {
    cnamat = as.matrix(cnamat)
    sqmat = 2^(cnamat)
    msq = rowMeans(sqmat)
    msq
}

.cnacor = function(cnamat, cor.method = 'pearson') {
    cnamat = as.matrix(cnamat)
    genemeans = rowMeans(cnamat)
    cellcors = cor(genemeans, cnamat, method = cor.method)
    unlist(as.data.frame(cellcors))
}


cnacor = function(cnamat, cor.method, groups = FALSE) {
    if (isFALSE(groups)) {
        return(.cnacor(cnamat, cor.method = cor.method))
    }
    
    infercna::group_apply(mat = cnamat,
                          groups = groups,
                          ungroup = TRUE,
                          FUN = infercna::.cnacor,
                          cor.method = cor.method)
}


cnastats = function(cnamat) {
    infercna::cnasignal(cnamat)
    infercna::cnacor(cnamat)
}
