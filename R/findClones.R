
fetchModes = function(m,
                      prob = 0.95,
                      coverage = 0.8,
                      size = 10,
                      by = 'chr') {

    mats = splitGenes(m, by = by)
    modes = sapply(mats, fitBimodal, prob = prob, coverage = coverage, size = size, assign = T)
    modes[!sapply(modes, isFALSE)]
}


expandToClones = function(modes,
                          greaterThan = NULL) {

    modes = sapply(modes, function(x) {
                       stats::setNames(unlist(x), rep(names(x), lengths(x)))},
                       simplify = F)

    modes = unlist(modes)
    clones = sapply(split(names(modes), modes), paste0, collapse = '--')
    clones = split(names(clones), clones)
    if (is.null(greaterThan)) return(clones)
    clones[lengths(clones) > greaterThan]
}


findClones = function(m,
                      prob = 0.95,
                      coverage = 0.8,
                      mode.size = 10,
                      clone.size = NULL,
                      by = 'chr') {

    modes = fetchModes(m,
                       prob = prob,
                       coverage = coverage,
                       size = mode.size,
                       by = by)

    expandToClones(modes, greaterThan = clone.size)
}
