
.cnaStatistic = function(cna,
                         samples = NULL,
                         threshold = NULL,
                         cor.threshold = threshold,
                         signal.threshold = threshold,
                         sep = "-|_") {
    
    s = cnaSignal(cna, threshold = signal.threshold)
    co = cnaCor(cna, samples = samples, threshold = cor.threshold, sep = sep)
    s/co
}
