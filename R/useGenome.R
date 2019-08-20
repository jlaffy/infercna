
#' @title List of available genomes
#' @description See which genomes are available for use
#' @return character vector
#' @details more genomes will be added in future dev patches.
#' @rdname availableGenomes
#' @export 
availableGenomes = function() {
    avai = paste(names(g.env$datasets), collapse = '\n')
    message('Available genomes:\n', avai)
}

#' @title List the current genome name
#' @description E.g. "hg19" if the genome being used is hg19
#' @return string
#' @details Default genome is "hg19"
#' @rdname currentGenome
#' @export 
currentGenome = function() {
    message('Genome: ', g.env$name)
}

#' @title Retrieve genome data
#' @description Returns a tibble dataframe of the current genome in use. If on function call a character string is supplied that corresponds to a genome in the package, that genome will instead be returned. 
#' @param x a genome name. Currently one of 'hg19' (human), 'hg38' (latest human), 'mm10' (mouse). Default: NULL
#' @return a tibble 
#' @seealso 
#'  \code{\link[tibble]{as_tibble}}
#' @rdname retrieveGenome
#' @export 
#' @importFrom tibble as_tibble
retrieveGenome = function(x = NULL) {
    if (is.null(x)) x = g.env$name
    else stopifnot(x %in% names(g.env$datasets))
    message('Retrieving: ', x)
    tibble::as_tibble(g.env$data)
}

#' @title Add your own Genome
#' @description Add your own Genome for infercna to use. The
#' @param genome PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION, Default: 'userDefined'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @rdname addGenome
#' @export 
#' @importFrom dplyr arrange
#' @importFrom stats setNames
addGenome = function(genome, name = 'userDefined') {
    stopifnot(all(c('symbol', 'start_position', 'chromosome_name', 'arm') %in% colnames(genome)))
    genome = dplyr::arrange(genome, chromosome_name, start_position)
    g.env$data = genome
    g.env$name = name
    g.env$gene = g.env$data$symbol
    g.env$order = stats::setNames(g.env$data$symbol, g.env$data$start_position)
    g.env$chr = stats::setNames(g.env$data$symbol, g.env$data$chromosome_name)
    g.env$arm = stats::setNames(g.env$data$symbol, g.env$data$arm)
    message('Genome has been set to ', name)
}

#' @title Select a Genome for infercna to use
#' @description Select your genome of choice at the start of an analysis. The available genomes in the current implementation are, for H.sapiens, hg38 (latest) and hg19 (preceding) and for mouse, mm10. You can see which genomes are available via availableGenomes(). 
#' @param x character string of genome name. One of 'hg19', 'hg38', 'mm10'.
#' @return genome variables are set internally. No return.
#' @rdname useGenome
#' @export 
#' @importFrom stats setNames
useGenome = function(x) {
    if (!is.character(x) || length(x) != 1 || !x %in% names(g.env$datasets)) {
        stop("Nope sorry, couldn't find genome <", x, ">.")
    }
    g.env$data = g.env$datasets[[x]]
    g.env$name = x
    g.env$gene = g.env$data$symbol
    g.env$order = stats::setNames(g.env$data$symbol, g.env$data$start_position)
    g.env$chr = stats::setNames(g.env$data$symbol, g.env$data$chromosome_name)
    g.env$arm = stats::setNames(g.env$data$symbol, g.env$data$arm)
    message('Genome has been set to ', x)
}
