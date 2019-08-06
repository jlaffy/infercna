
genmatch = function(genes) {
    match(g.env$gene, genes, nomatch = 0)
}

as_genomic = function(x = NULL) {
    if (is.null(dim(x))) return(x[genmatch(genes = x)])
    x[genmatch(rownames(x)), , drop = FALSE]
}

genmissing = function(x, name = FALSE) {
    missing = which(genmatch(x) == 0)
    if (name) missing = x[missing]
    missing
}

# filter gene list by include or exclude attributes
genfilter = function(x = NULL, include = NULL, exclude = NULL) {
    xin = sapply(include, attrifilter, x = x, simplify = F) %>%
        Reduce(intersect, .)
    xou = sapply(exclude, attrifilter, x = xin, simplify = F) %>%
        Reduce(union, .)
    setdiff(xin, xou)
}
