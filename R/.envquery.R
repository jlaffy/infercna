
# genome attributes: genes, chromosomes, chromosome arms

# return all attributes 
# is.atomic and length > 1 to exclude genome name & dataframes

envmatch = function(x) {
    match(x, g.env$gene, nomatch = 0)
}

lglob = function(genes = NULL) {
    L = as.list(g.env)
    L = L[sapply(L, function(l) is.atomic(l) && length(l) > 1)]
    if (!is.null(genes)) L = sapply(L, `[`, envmatch(genes), simplify = F)
    L
}

# find the attribute that a value corresponds to
lidentify = function(v, genes = NULL, reveal = TRUE) {
    L = lglob(genes = genes)
    i = which(sapply(L, function(l) all(v %in% l)))
    if (!reveal) return(i)
    stats::setNames(L[[i]], L[['genes']])
}

attrifetch = function(v, positions = FALSE) {
    attri = lidentify(v, reveal = TRUE)
    matching = match(v, attri, nomatch = 0)
    if (positions) return(matching)
}

# filter G
attrifilter = function(v, genes = NULL, preserve = TRUE) {
    L = attrifetch(v, genes = genes)
    if (preserve) L[names(L) %in% v]
}


valmatch = function(v, attribute) {
    attr = attrimatch(v, genes = genes, fetch = TRUE)
    match(v, attr, nomatch = 0)
}

valfilters = function(..., strict = FALSE) {
    dots = list(...)
    if (strict) FUN = intersect
    else FUN = union
    Reduce(FUN, sapply(dots, valfilter, simplify = F))
}

genefilter = function(...) {

}
