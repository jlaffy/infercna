
.fitBimodalBySampling = function(x,
                                 sample.size = 100,
                                 tries = 2000,
                                 force.tries = FALSE,
                                 prob = 0.95,
                                 coverage = 0.8,
                                 size = 10,
                                 verbose = T,
                                 assign = T,
                                 boolean = F,
                                 maxrestarts = 100,
                                 maxit = 5000) {
    a = NULL
    b = NULL
    a.flat = 0
    b.flat = 0
    i = 0

    if (!is.null(dim(x))) x = colMeans(x)

    while (i < tries) {

        assign('i', i + 1)
        res = fitBimodal(x[sample(names(x), sample.size)],
                         assign = T,
                         boolean = F,
                         bySampling = F,
                         coverage = coverage,
                         prob = prob,
                         size = size,
                         verbose = verbose,
                         maxit = maxit,
                         maxrestarts = maxrestarts)

        if (isFALSE(res)) {next()}
        if (i > 3 & length(intersect(a, res$a)) == 0) {next()}
        if (i > 3 & length(intersect(b, res$a)) != 0) {next()}
        if (i > 3 & length(intersect(b, res$b)) == 0) {next()}
        if (i > 3 & length(intersect(a, res$b)) != 0) {next()}

        assign('a', c(a, res$a))
        assign('b', c(b, res$b))

        if (!force.tries) {
            assign('a.flat', mean(table(a) >= 5))
            assign('b.flat', mean(table(b) >= 5))

            if (a.flat >= 0.97 & b.flat >= 0.97) {
                message('Found after ', i, ' tries')
                break()
            }
        }
    }

    if (is.null(a) & is.null(b)) return(FALSE)
    taba = table(a)
    tabb = table(b)
    list(a = names(taba)[taba >= 1], b = names(tabb)[tabb > 1])
}
