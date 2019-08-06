
availableGenomes = function() {
    avai = paste(names(g.env$datasets), collapse = '\n')
    message('Available genomes:\n', avai)
}

currentGenome = function() {
    message('Genome: ', g.env$name)
}

retrieveGenome = function(x = NULL) {
    if (is.null(x)) x = g.env$name
    else stopifnot(x %in% names(g.env$datasets))
    message('Retrieving: ', x)
    tibble::as_tibble(g.env$data)
}


addGenome = function(genome, name = 'userDefined') {
    stopifnot(all(c('symbol', 'start_position', 'chromosome_name', 'arm') %in% colnames(genome)))
    genome = dplyr::arrange(genome, chromosome_name, start_position)
    g.env$data = genome
    g.env$name = name
    g.env$genes = g.env$data$symbol
    g.env$order = stats::setNames(g.env$data$symbol, g.env$data$start_position)
    g.env$chrm = stats::setNames(g.env$data$symbol, g.env$data$chromosome_name)
    g.env$arm = stats::setNames(g.env$data$symbol, g.env$data$arm)
    message('Genome has been set to ', name)
}

useGenome = function(x) {
    if (!is.character(x) || length(x) != 1 || !x %in% names(g.env$datasets)) {
        stop("Nope sorry, couldn't find genome <", x, ">.")
    }
    g.env$data = g.env$datasets[[x]]
    g.env$name = x
    g.env$genes = g.env$data$symbol
    g.env$order = stats::setNames(g.env$data$symbol, g.env$data$start_position)
    g.env$chrm = stats::setNames(g.env$data$symbol, g.env$data$chromosome_name)
    g.env$arm = stats::setNames(g.env$data$symbol, g.env$data$arm)
    message('Genome has been set to ', x)
}
