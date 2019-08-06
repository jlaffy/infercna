
g.env = new.env()
g.env$datasets = list(hg19 = infercna:::hg19,
                      hg38 = infercna:::hg38,
                      mm10 = infercna:::mm10)
g.env$name = ' hg19'
g.env$data = g.env$datasets[['hg19']]
g.env$gene = g.env$data$symbol
g.env$chr = stats::setNames(g.env$data$chromosome_name, g.env$data$symbol)
g.env$arm = stats::setNames(g.env$data$arm, g.env$data$symbol)
g.env$start = stats::setNames(g.env$data$start_position, g.env$data$symbol)
g.env$end = stats::setNames(g.env$data$end_position, g.env$data$symbol)
