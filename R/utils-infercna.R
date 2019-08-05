
colCenter = function(m) {
    scale(m, center = T, scale = F)
}

clip <- function(mat, range = c(-3, 3)) {
    mat = as.matrix(mat)
    mat[mat < range[[1]]] <- range[[1]]
    mat[mat > range[[2]]] <- range[[2]]
    mat
}

runMean <- function(mat,
                    k = 100,
                    endrule = 'mean',
                    align = 'center') {

    mat = as.matrix(mat)
    mat.out = caTools::runmean(mat,
                               k = k,
                               endrule = endrule,
                               align = align)
    colnames(mat.out) <- colnames(mat)
    rownames(mat.out) <- rownames(mat)
    as.data.frame(mat.out)
}

