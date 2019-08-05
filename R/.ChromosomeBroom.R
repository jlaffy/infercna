
source('~/R/Chromosome.R')

asTumourNames <- function(cellNames, splitChar = "-") {
    stringr::str_split(cellNames, splitChar) %>%
        sapply(., `[`, 1)
}

.chrBroomTidy = function(data = NULL, X = NULL) {
    if (is.null(X)) {
        X = melt(as.matrix(data)) %>%
            setNames(c('Gene', 'Cell', 'CNA')) %>%
            dplyr::mutate(Tumour = asTumourNames(Cell))
    }

    stopifnot(all(c('Gene', 'Cell', 'CNA', 'Tumour') %in% colnames(X)))

    X = X %>%
        dplyr::mutate(Gene = factor(as.character(Gene),
                                    levels = orderByChromosome(Gene)),
                      Chr = factor(isChromosome(Gene),
                                   levels = chrLevels()))

    X %>% dplyr::select(Tumour, Cell, Chr, Gene, CNA)
}

.chrBroomReduce = function(X, ..., full = F) {
    colvars = rlang::enquos(...)
    colnams = sapply(colvars, rlang::quo_name)

    Xout = X %>%
        mutate(sqCNA = CNA^2) %>%
        group_by(!!!colvars) %>%
        summarize(n = n(),
                  Mean = mean(CNA),
                  SD = sd(CNA),
                  MeanSquares = mean(sqCNA))

    if (full) {
        return(full_join(X, Xout, by = colnams))
    }

    Xout
}

chrBroomCell = function(X, highSD = F, nSD = 1) {
    x = X %>%
        select(Tumour, Cell, Gene, CNA)

    y = .chrBroomReduce(X, Tumour, Gene) %>%
        select(Tumour, Gene, Mean, SD)

    z = .chrBroomReduce(X, Tumour, Cell) %>%
        select(Tumour, Cell, MeanSquares)

    DF = reduce(list(x, y, z), full_join)

    if (!is.logical(highSD) & is.numeric(highSD)) {
        DF = DF %>% dplyr::filter(SD >= highSD)
    }

    else if (highSD) {
        DF = DF %>%
            group_by(Tumour) %>%
            summarize(SDcut = mean(SD) + nSD * sd(SD)) %>%
            full_join(DF, .) %>%
            dplyr::filter(SD >= SDcut)
    }

    DF %>%
        group_by(Cell) %>%
        nest() %>%
        mutate(Cor = purrr::map(data, ~cor(.x$Mean, .x$CNA))) %>%
        unnest(Cor, .drop = F) %>%
        unnest(data, .drop = F) %>%
        distinct(Tumour, Cell, Cor, MeanSquares)
}


chrBroomGene = function(X) {
    .chrBroomReduce(X, Tumour, Chr, Gene)
}

BroomTSNE = function(data,
                     mat = NULL,
                     cnamat = NULL,
                     minPts = 5,
                     eps = NULL,
                     DE = F,
                     CNAsigcut = 0.012,
                     CNAcorcut = 0.25,
                     signatures = NULL,
                     SquishSigScores = NULL,
                     SquishCNAcor = NULL,
                     SquishCNAsig = NULL,
                     ngenes = 100,
                     fc.value = 2,
                     p.value = 10^(-3),
                     returning = 'fc',
                     fc.sort = T) {

    DEgroups = NULL
    if (is.null(dim(data))) {
        X = data$Y
    } else {
        X = data
    }

    X = as.data.frame(X)

    if (is.null(eps)) {
        plot(dbscan::hdbscan(X, minPts = minPts))
        return()
    }

    k = fpc::dbscan(X, eps = eps)$cluster
    X = X %>%
        setNames(., c('tSNE1', 'tSNE2')) %>%
        tibble::rownames_to_column('Cell') %>%
        mutate(dbscan = factor(k, levels = sort(unique(k)))) %>%
        arrange(dbscan) %>%
        mutate(Tumour = asTumourNames(Cell))


    if (!is.null(mat) & DE) {
        if (ncol(mat) == nrow(X)) {
            clusters = split(X$Cell, X$dbscan)
            DEcall = quote(DEgenes(k = k,
                                   mat = mat,
                                   fc.value = fc.value,
                                   p.value = p.value,
                                   fc.sort = fc.sort,
                                   returning = returning))

            DEgroups = sapply(clusters, function(k) {
                                  if (length(k) == 1) {
                                      NULL
                                  } else {
                                      eval(DEcall)
                                  }},
                                  simplify = F)
            DEgroups = sapply(DEgroups, function(G) G[1:min(length(G), ngenes)],
                              simplify = F)
        }
    }

    if (!is.null(signatures) & !is.null(mat)) {
        signatures = sapply(signatures, function(s) s[s %in% rownames(mat)], simplify = F)
        scoremat = score(mat, signatures, bin.control = T)
        if (!is.null(SquishSigScores)) {
            scoremat = scales::squish(scoremat, range = SquishSigScores)
        }
        scoremat = scoremat %>%
            as.data.frame %>%
            tibble::rownames_to_column('Cell')
        X = full_join(X, scoremat)
    }

    if (!is.null(cnamat)) {
        cnamat = .chrBroomTidy(cnamat)
        X = chrBroomCell(cnamat, highSD = T) %>%
            select(Cell, Cor, MeanSquares) %>%
            setNames(., c('Cell', 'CNAcor', 'CNAsignal')) %>%
            full_join(X, .)
        if (!is.null(CNAsigcut) & !is.null(CNAcorcut)) {
            X = X %>%
                mutate(CNAhigh = ifelse(CNAcor >= CNAcorcut & CNAsignal >= CNAsigcut, "High", "Low"))
        }
        if (!is.null(SquishCNAcor)) {
            X$CNAcor = scales::squish(X$CNAcor, range = SquishCNAcor)
        }
        if (!is.null(SquishCNAsig)) {
            X$CNAsignal = scales::squish(X$CNAsig, range = SquishCNAsig)
        }
    }

    return(list(X = X,
                eps = eps,
                signatures = signatures,
                DEgroups = DEgroups))
}
