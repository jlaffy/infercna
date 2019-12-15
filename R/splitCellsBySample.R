
.splitCellsBySample = function(cells,
                               samples = NULL,
                               sep = "-|_",
                               max.nchar = 6,
                               pos = 1) {

    if (is.null(samples)) {
        samples = scalop::get_sample_names(cells, sep = sep, max.nchar = max.nchar, pos = pos)
        return(split(cells, samples))
    }

    else if (length(samples) < length(cells) & !any(duplicated(samples))) {
        return(sapply(samples, function(Sample) cells[stringr::str_detect(cells, Sample)], simplify = F))
    }

    else stop('cells or samples not recognised or not compatible')
}


splitCellsBySample = function(x,
                              samples = NULL,
                              cellsBySample = NULL,
                              sep = "-|_",
                              max.nchar = 6,
                              sample.pos = 1) {

    if (!is.null(dim(x))) cells = colnames(x)

    else {
        if (is.null(names(x))) names(x) = x
        cells = names(x)
    }

    if (is.null(cellsBySample)) {
        cellsBySample = .splitCellsBySample(cells = cells,
                                            sep = sep,
                                            samples = samples,
                                            max.nchar = max.nchar,
                                            pos = sample.pos)
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

