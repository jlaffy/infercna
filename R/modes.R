
modes = function(m = NULL,
                 X = NULL,
                 prob = 0.95,
                 coverage = 0.8,
                 size = 10,
                 boolean = FALSE) {
    if (is.null(X)) X = colMeans(m)
    negativeRes = stats::setNames(rep(0, length(X)), names(X))
    probs = try(mixtools::normalmixEM(X)$posterior)
    if (class(probs) == 'try-error' | is.null(probs)) {
        if (boolean) return(FALSE)
        else return(negativeRes)
    }
    probs = as.data.frame(probs)
    probs = dplyr::mutate(probs,
                          bool.1 = comp.1 >= prob,
                          bool.2 = comp.2 >= prob,
                          bool = bool.1 | bool.2)
    modesexist = mean(probs$bool) > coverage && sum(probs$bool.1) >= size && sum(probs$bool.2) >= size
    if (boolean) return(modesexist)
    if (!modesexist) return(negativeRes)
    probs = dplyr::mutate(probs, modeID = ifelse(comp.1 >= prob, 1, ifelse(comp.2 >= prob, 2, 0)))
    stats::setNames(probs$modeID, names(X))
}


clones = function(cnamat,
                  by_chromosome_arm = FALSE,
                  by_chromosome = FALSE,
                  prob = 0.95,
                  coverage = 0.8,
                  size = 10,
                  boolean = FALSE) {

    if (by_chromosome_arm) mats = infercna::splitMatrixByArm(cnamat)
    else if (by_chromosome) mats = infercna::splitMatrixByChr(cnamat)
    mats = mats[sapply(mats, nrow) >= 5]
    modemat = sapply(mats,
                     infercna::modes,
                     prob = prob, coverage = coverage,
                     size = size, boolean = boolean,
                     simplify = T)
    stopifnot(is.matrix(modemat))
    modemat = modemat[rowSums(modemat) != 0, colSums(modemat) != 0, drop = FALSE]
    cloneIDs1 = apply(modemat, 1, function(r) paste0(colnames(modemat)[which(r == 1)],
                                                     collapse = '_'))
    cloneIDs2 = apply(modemat, 1, function(r) paste0(colnames(modemat)[which(r == 2)],
                                                     collapse = '_'))
    cloneIDs = mapply(function(x, y) paste(x, y, sep = '-'),
                      x = cloneIDs1,
                      y = cloneIDs2)
    split(names(cloneIDs), cloneIDs)
}
