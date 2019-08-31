
cnaScatterPlot = function(cna,
                          threshold = NULL,
                          cor.threshold = threshold,
                          signal.threshold = .9,
                          group = NULL,
                          group.col = 'magenta',
                          hline = NULL,
                          vline = NULL,
                          bySample = FALSE,
                          samples = NULL,
                          sep = "-|_",
                          excludeFromAvg = NULL) {

    cors = cnaCor(cna,
                  threshold = cor.threshold,
                  bySample = bySample,
                  samples = samples,
                  sep = sep,
                  excludeFromAvg = excludeFromAvg)

    signals = cnaSignal(cna, threshold = signal.threshold)

    plot(cors, signals, xlab = 'CNA Correlation', ylab = 'CNA Signal', pch = 1, cex = 0.5)

    if (!is.null(group)) {
        points(cors[group], signals[group], pch = 20, col = group.col, cex = 0.5)
    }

    if (!is.null(vline)) abline(v = vline, lty = 2)
    if (!is.null(hline)) abline(h = hline, lty = 2)
}
