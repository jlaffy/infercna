
#' @title List of available genomes
#' @description See which genomes are available for use
#' @return character vector
#' @details more genomes will be added in future dev patches.
#' @rdname availableGenomes
#' @export 
availableGenomes = function() {
    avai = paste(c('hg19', 'hg38', 'mm10'), collapse = '\n')
    message('Available genomes:\n', avai)
}

#' @title List the current genome name
#' @description E.g. "hg19" if the genome being used is hg19
#' @return string
#' @details Default genome is "hg19"
#' @rdname currentGenome
#' @export 
currentGenome = function() {
    message('Genome: ', Genv$name)
}

#' @title Retrieve genome data
#' @description Returns a tibble dataframe of the current genome in use. If on function call a character string is supplied that corresponds to a genome in the package, that genome will instead be returned. 
#' @param name a genome name. Currently one of 'hg19' (human), 'hg38' (latest human), 'mm10' (mouse). Default: NULL
#' @return a tibble 
#' @seealso 
#'  \code{\link[tibble]{as_tibble}}
#' @rdname retrieveGenome
#' @export 
#' @importFrom tibble as_tibble
retrieveGenome = function(name = NULL) {

    if (is.null(name)) {
        name = Genv$name
        data = Genv$data
    } else {
        data = .fetchGenome(name)
    }

    message('Retrieving: ', name)
    tibble::as_tibble(data)
    
}

#' @title Add your own Genome
#' @description Add your own Genome for infercna to use. The
#' @param genome genome data provided as a dataframe. The dataframe should contain columns 'symbol', 'start_position', 'end_position', 'chromosome_name', 'arm'. The columns 'chromosome_name' and 'arm' should be factors, with the chromosome/chromosome arms ordered correctly. 
#' @param name a character string; the name of your genome. Default: 'userDefined'
#' @return no return value. The genome variables will be updated internally.
#' @rdname addGenome
#' @export 
#' @importFrom dplyr arrange
#' @importFrom stats setNames
addGenome = function(genome, name = 'userDefined') {
    columns = c('symbol', 'start_position', 'chromosome_name', 'arm')

    if (!all(columns %in% colnames(genome))) {
        stop('Columns must include: ', columns)
    }

    if (is.null(levels(genome$chromosome_name))) {
        stop('Please add levels to the chromosome_name column')
    }

    if (is.null(levels(genome$arm))) {
        stop('Please add levels to the arm column')
    }

    genome = dplyr::arrange(genome, chromosome_name, start_position)
    .configureGenome(data = genome, name = name)
    message('Genome has been set to ', name)
}

#' @title Select a Genome for infercna to use
#' @description Select your genome of choice at the start of an analysis. The available genomes in the current implementation are, for H.sapiens, hg38 (latest) and hg19 (preceding) and for mouse, mm10. You can see which genomes are available via availableGenomes(). 
#' @param name a character string of genome name. One of 'hg19', 'hg38', 'mm10'.
#' @return genome variables are set internally. No return.
#' @rdname useGenome
#' @export 
useGenome = function(name) {
    .configureGenome(name = name)
    message('Genome has been set to ', name)
}
