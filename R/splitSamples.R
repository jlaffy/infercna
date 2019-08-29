
splitSamples = function(x, sep = "-|_", samples = NULL) {
    if (!is.null(dim(x))) cells = colnames(x)
    else if (!is.null(names(x)) & !is.character(x)) cells = names(x)
    else cells = x

    if (is.null(samples)) {
        splut = stringr::str_split(cells, sep)
        sampleNames = sapply(splut, `[`, 1)
        samples = split(x, samples)
    }

    if (!is.null(dim(x))) {
        return(sapply(samples, function(Sample) x[, Sample], simplify = F))
    }

    else if (!is.null(names(x))) {
        return(sapply(samples, function(Sample) x[Sample], simplify = F))
    }

    samples
}
