
.splitCellsBySample = function(cells, samples = NULL, sep = "-|_", max.nchar = 6) {

    if (is.null(samples)) {
        splut = stringr::str_split(cells, sep)
        samples = sapply(splut, `[`, 1)
        if (!is.null(max.nchar)) {
            samples = sapply(samples, function(sample) {
                                 if (nchar(sample) <= max.nchar) sample
                                 else substr(sample, 1, max.nchar)})
            samples = unique(samples)
        }
    }

    if (length(unique(toupper(samples))) < length(samples)) samples = toupper(samples)

    if (length(samples) < length(cells) & !any(duplicated(samples))) {
        return(sapply(samples, function(Sample) cells[stringr::str_detect(cells, Sample)], simplify = F))
    }

    split(cells, samples)
}

splitCellsBySample = function(x,
                              samples = NULL,
                              cellsBySample = NULL,
                              sep = "-|_",
                              max.nchar = 6) {

    if (!is.null(dim(x))) {
        cells = colnames(x)
    }

    else {
        if (is.null(names(x))) names(x) = x
        cells = names(x)
    }

    if (is.null(cellsBySample)) {
        cellsBySample = .splitCellsBySample(cells = cells, sep = sep, samples = samples, max.nchar = max.nchar)
    }

    if (!is.null(dim(x))) {
        return(sapply(cellsBySample, function(Cells) x[, Cells, drop = F], simplify = F))
    }
    
    out = sapply(cellsBySample, function(Cells) x[Cells], simplify = F)

    if (all(names(out[[1]]) == out[[1]])) {
        out = sapply(out, function(o) stats::setNames(o, NULL), simplify = F)
    }

    out
}
