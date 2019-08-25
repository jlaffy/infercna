
Genv = new.env(parent = emptyenv())

.fetchGenome = function(name) {
    stopifnot(is.character(name))
    stopifnot(length(name) == 1)
    if (!name %in% ls(envir = asNamespace('infercna'))) {
        stop("Genome name '", name, "' does not exist in pkg data.")
    }
    data = try(get(name))
    if (class(data) == 'try-error') {
        stop("Genome name '", name, "' does not exist in pkg data.")
    }
    data
}

.prepGenv = function(name = NULL, data = NULL) {

    if (is.null(name) & is.null(data)) {
        stop("No arguments provided.")
    }

    if (is.null(name) & !is.null(data)) {
        name = 'userDefined'
    }
    
    else if (!is.null(name)) {
        data = .fetchGenome(name)
    }
    
    vars = as.list(as.data.frame(data))
    nam = c('symbol', 'chromosome_name', 'arm', 'start_position', 'end_position')
    nam2 = c('gene', 'chr', 'arm', 'start', 'end')
    vars = stats::setNames(vars[nam], nam2)
    vars[2:5] = sapply(vars[2:5], stats::setNames, vars[['gene']], simplify = F)
    vars = c(list(name = name), vars)
    vars
}

.setGenome = function(name = NULL, data = NULL) {
    vars = .prepGenv(name = name, data = data)
    silent = Map(assign, value = vars, x = names(vars), MoreArgs = list(envir = Genv))
}

assign(x = 'data', value = infercna:::hg19, envir = Genv)
vars = .prepGenv(name = 'hg19')
Map(assign, value = vars, x = names(vars), MoreArgs = list(envir = Genv))
