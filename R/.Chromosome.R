
chrLevels = function() {
    c(paste0('chr', 1:9), paste0('chr1', 0:9), paste0('chr2', 0:2), 'chrX', 'chrY')
}

fixChrLabels = function(chr) {
    if (all(chr %in% chrLevels())) {
        return(chr)
    }

    if (any(chr == 'x')) {
        chr[which(chr == 'x')] = 'X'
    }

    if (any(chr == 'y')) {
        chr[which(chr == 'y')] = 'Y'
    }
    
    sapply(chr, function(n) paste0('chr', n))
}

# THIS EXCLUDES THE MITOCHONDRIAL CHROMOSOME
genomicPositionInfo = function(gene = NULL, chr = NULL) {
    chrdat = readRDS('~/rds/gene.pos.rds') %>%
        dplyr::filter(Chromosome != 'chrM') %>%
        mutate(Chromosome = factor(as.character(Chromosome), levels = chrLevels())) %>%
        arrange(Chromosome)

    if (!is.null(chr)) {
        chr = fixChrLabels(chr)
        chrdat = chrdat %>% dplyr::filter(Chromosome %in% chr)
    }
    if (!is.null(gene)) {
        chrdat = chrdat %>% dplyr::filter(Gene %in% gene)
    }
    chrdat
}

genesInChromosome = function(chr) {
    chr = fixChrLabels(chr)
    chrdat = genomicPositionInfo(chr = chr) %>%
        dplyr::arrange(Chromosome, Start) %>%
        dplyr::mutate(Chromosome = as.character(Chromosome))
    
    if (length(chr) == 1) {
        return(setNames(chrdat$Gene, chrdat$Chromosome))
    }
    split(chrdat$Gene, chrdat$Chromosome)[chr]
}

isChromosome = function(gene = NULL) {
    chrdat = genomicPositionInfo(gene = gene) %>%
        dplyr::mutate(Chromosome = as.character(Chromosome))

    setNames(chrdat$Chromosome, chrdat$Gene)[gene]
}

orderByChromosome = function(genes = NULL, ind = F) {
    chrdat = genomicPositionInfo(gene = genes) %>%
        dplyr::arrange(Chromosome, Start)
    ord = setNames(chrdat$Gene, chrdat$Chromosome)
    if (ind) {
        return(match(genes, ord))
    }
    ord
}

chromosomeBreakPoints = function(...) {
    chromosomeBreaks(...)
}

chromosomeBreaks = function(genes = NULL, halfway = F, stripped = F) {
    chr = isChromosome(gene = genes)
    n = length(chr)
    chrsum = cumsum(table(chr)[chrLevels()])
    Breaks = chrsum/max(chrsum) * n

    if (halfway) {
        b = Breaks
        Breaks = c(b[1]/2, b[2:length(b)] - diff(b)/2)
    }

    if (stripped) {
        names(Breaks) = stringr::str_replace(names(Breaks), "chr", "")
    }

    Breaks
}

groupBy = function(x, Names = NULL) {
    stopifnot(!is.null(Names) | !is.null(names(x)))
    
    if (!is.null(Names)) {
        names(x) <- Names
    }

    split(x, names(x))
}

groupByChromosome = function(x) {
    x = orderByChromosome(x)
    groupBy(x)[chrLevels()]
}

splitByChromosome = function(mat) {
    x = rownames(mat)
    sapply(groupByChromosome(x), function(X) mat[rownames(mat) %in% X, ],
           simplify = F)
}
