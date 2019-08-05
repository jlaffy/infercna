
guess_groups = function(x, verbose = TRUE, sep = "-|_|\\.") {

    groups = stringr::str_split(x, sep)

    if (identical(groups, as.list(x))) {
        stop('No groups found.')
    }

    groups = sapply(groups, `[`, 1)
    grouped = split(x, groups)

    if (verbose) {
        message(length(grouped), ' groups found:')
        print(lengths(grouped))
    }

    grouped
}

split_by_group = function(mat, groups = NULL) {
    if (is.null(groups)) {
        groups = infercna::guess_groups(colnames(mat))
    }
    if (!is.list(groups)) {
        groups = split(groups, names(groups))
    }
    sapply(groups, function(group) mat[, group], simplify = F)
}


group_apply(mat, FUN, groups = NULL, ungroup = F, ...) {

    mats = group_split(mat, groups = groups)
    result = sapply(mats, FUN, ..., simplify = F)

    if (ungroup) {
        Names = rep(names(result), lengths(result))
        result = unlist(result, use.names = F)
        result = stats::setNames(result, Names)
    }

    result
}
