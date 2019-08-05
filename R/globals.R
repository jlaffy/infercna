
g.env = new.env()
g.env$datasets = list(hg19 = infercna:::hg19,
                      hg38 = infercna:::hg38,
                      mm10 = infercna:::mm10)
g.env$name = 'hg19'
g.env$data = g.env$datasets[['hg19']]
g.env$genes = g.env$data$symbol
g.env$order = stats::setNames(g.env$data$symbol, g.env$data$start_position)
g.env$chrm = stats::setNames(g.env$data$symbol, g.env$data$chromosome_name)
g.env$arm = stats::setNames(g.env$data$symbol, g.env$data$arm)

