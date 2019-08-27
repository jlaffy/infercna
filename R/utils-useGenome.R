
.fetchGenome = function(name) {
    stopifnot(is.character(name))
    stopifnot(length(name) == 1)
    result = try(utils::getFromNamespace(x = name, ns = 'infercna'))
    if (class(result) == 'try-error') {
        stop("Genome name '", name, "' does not exist in package data.")
    }
   result 
}

.parseGenome = function(data, name) {
    vars = as.list(as.data.frame(data))
    nam = c('symbol', 'chromosome_name', 'arm', 'start_position', 'end_position')
    nam2 = c('gene', 'chr', 'arm', 'start', 'end')
    vars = stats::setNames(vars[nam], nam2)
    vars[2:5] = sapply(vars[2:5], stats::setNames, vars[['gene']], simplify = F)
    c(list(name = name, data = data), vars)
}

.configureGenome = function(name = NULL, data = NULL) {
    if (is.null(data) & is.null(name)) stop('Please provide <name> or/and <data>.')
    else if (is.null(data)) data = .fetchGenome(name)
    else if (is.null(name)) name = 'userDefined'
    vars = .parseGenome(data = data, name = name)
    list2env(vars, envir = Genv)
    Genv
}

