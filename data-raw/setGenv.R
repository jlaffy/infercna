
vars = as.list(as.data.frame(hg19))
nam = c('symbol', 'chromosome_name', 'arm', 'start_position', 'end_position')
nam2 = c('gene', 'chr', 'arm', 'start', 'end')
vars = stats::setNames(vars[nam], nam2)
vars[2:5] = sapply(vars[2:5], stats::setNames, vars[['gene']], simplify = F)
vars = c(list(name = 'hg19', data = data), vars)
Genv = list2env(vars)
rm(vars)
rm(nam)
rm(nam2)

