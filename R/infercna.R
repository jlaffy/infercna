
.mostExpressedGenes = function(m, ngenes) {
    if (ngenes < 0) stop('<ngenes> cannot be negative.')
    if (ngenes > nrow(m)) stop('<ngenes> cannot be larger than nrow(m).')

    if (ngenes >= 0 & ngenes <= 1) ngenes = ngenes * nrow(m)
    else if (ngenes > nrow(m)) ngenes = nrow(m)

    rom = rowMeans(m)
    names(sort(rowMeans(m), decreasing = T)[1:ngenes])
}

#' @title Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @description Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @param m a matrix of genes X cells containing scRNA-seq expression data. The matrix should NOT be row-centered.
#' @param refCells a list of two or more character vectors, each containg cell IDs of normal cell types (one cell type per list element). Since these cell types are presumed normal cells, their CNA values can be used to correct the remaining CNA values. Note that all cell IDs should correspond to column names in <m>. See infercna::refCells (normal cells) for an example. Default: NULL
#' @param window the size of the window to use for the rolling mean. Units are in number of genes. Default: 100
#' @param range values in <m> above and below this range will be set to values in range. See infercna::clip() for more details. Default: c(-3, 3)
#' @param n a numeric value indicating if only top genes should be used in the CNA calculation and how many. if n is a fraction between 0 and 1, the number of genes included will n * nrow(m). Else the number of genes included will be n. Default: NULL
#' @param noise a numeric value, which if given, increases the boundaries within which CNA values are considered 0. Increases by <noise> i.e. the bounds become Minimum - noise and Maximum + noise. Default: 0.1
#' @param center.method method by which to center the cells after calculating CNA values. One of 'median', 'mean', etc.... Default: 'median'
#' @param isLog boolean value indicating whether the input expression matrix <m> is in log2 form. Default: FALSE
#' @param verbose print progress messages. Default: FALSE
#' @return a matrix of genes X cells of inferred CNA values. Note that n = (window - 1)/2 genes will be lost from either extremity of the genome (ie. n genes lost at the start of chromosome 1 and at the end of chromosome Y, if the genome in question is H.sapiens.)
#' @details Correction with reference cells' <refCells> CNAs: the boundaries of reference CNA values are the boundaries for what should be considered a CNA of 0. Thus, if the boundary is -0.1 and 0.1, then a cell with CNA = -0.5 will be corrected to -0.4, a cell with CNA value of 1 will be corrected to 0.9 and a cell with CNA value of 0.05 will be corrected to 0.
#' @examples 
#'  m = infercna::useData()
#'  cna = infercna::infercna(m = m, refCells = infercna::refCells)
#' @rdname infercna
#' @export 
infercna = function(m,
                    refCells = NULL,
                    window = 100,
                    range = c(-3, 3),
                    n = NULL,
                    noise = 0.1,
                    center.method = 'median',
		    window.break = c('chr', 'arm'),
		    exclude.chr = NULL, #c('Y','X')
                    isLog = TRUE,
                    verbose = TRUE,
		    genome = c('hg19', 'hg38')) {


 
    useGenome(match.arg(genome))

    # note: m should NOT be row(gene)-centered
    if (all(round(range(rowMeans(m)), 3) == 0)) {
        stop('Matrix is row-centered. Please provide non-centered data.')
    }


    if (!is.null(n)) {

    	if (isLog) {
    	    if (verbose) {
    	    	message('Converting <m> from log(2) space...')
    	    }
    	    m = unlogtpm(m)
	}

        if (verbose) {
		message('Filtering the expression matrix to include only top ', n, 'genes...')
	}
    	
    	m = m[.mostExpressedGenes(m, ngenes = n), ]

    	if (isLog) {
    	    if (verbose) {
    	    	message('Converting <m> to log(2) space...')
    	    }

    	    m = logtpm(m, bulk = T)
	}
    }

    if (!isLog) {
        if (verbose) {
		message('Converting <m> to log(2) space...')
	}
	    
	m = logtpm(m)
    }


    if (verbose) {
	    message('Performing mean-centering of the genes...') 
    }
    
    m = rowcenter(m, by = 'mean')

    if (verbose) {
	    message('Ordering the genes by their genomic position...')
    }
    
    m = orderGenes(m)
    
    if (verbose) {
	    message('Restricting relative expression values to between ', range[[1]], ' and ', range[[2]], '..')
    }
    
    m = clip(m, range = range)

    if (verbose) {
            message('Converting <m> from log(2) space...')
    }
    
    m = unlogtpm(m, bulk = F) 
    
    if (verbose) {
	    message('Preparing to calculate CNA values on each chromosome in turn...')
    }
    
    ms = splitGenes(m, by = match.arg(window.break))
    ms = ms[!names(ms) %in% exclude.chr]
    
    if (verbose) {
	    message('Calculating rolling means with a window size of ', window, ' genes...')
    }
    
    cna = sapply(ms, function(m) try(runMean(m, k = window, verbose = verbose)), simplify = F)
    cna = cna[sapply(cna, class) != 'try-error' | !sapply(class, isFALSE)]
    cna = Reduce(rbind, cna)

    if (verbose) {
            message('Converting CNA values to log(2) space...')
    }
    
    cna = logtpm(as.matrix(cna), bulk = T)
    
    if (verbose) {
	    message('Performing ', center.method, '-centering of the cells...')
    }
    
    cna = colcenter(cna, by = center.method) # note: using median centering here

    if (!is.null(refCells)) {
        if (verbose) {
		message('Correcting CNA profiles using CNA values from <refCells>...')
	}
        
    	Args = c(list(cna = cna, noise = noise), refCells)
        cna = do.call(refCorrect, Args)
    }

    if (verbose) {
	    message('Done!')
    }

    cna
}
