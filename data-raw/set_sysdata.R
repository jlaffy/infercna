# 
# source('retrieve_genome_info.R')
# hg19 = quickmart('hg19')
# hg38 = quickmart('hg38')
# mm10 = quickmart('mm10')
# marts = list(hg19, hg38, mm10)
# add_attribute = c('hgnc_symbol', 'hgnc_symbol', 'mgi_symbol')
# dats = Map(getGenome, mart = marts, add_attribute = add_attribute)
# 
# dats = sapply(1:3, function(i) colnames(dats[[i]])[7] <<- 'symbol')
# dats[1:2] = sapply(dats[1:2], function(d) d %>% dplyr::mutate(arm = paste0(chromosome_name,
#                                                                            stringr::str_extract(band, "p|q"))),
#                    simplify = F)
# dats[[3]] = dats[[3]] %>% dplyr::mutate(arm = chromosome_name)
# 
# 
# settingchrlevels = function(dat) {
#     x = unique(dat$chromosome_name)
#     xnum = sort(as.numeric(x[x %in% as.character(1:100)]))
#     xchar = sort(x[!x %in% xnum])
#     chrlev = c(xnum, xchar)
#     dat = dplyr::mutate(dat, chromosome_name = factor(as.character(chromosome_name), levels = chrlev))
#     dat = dplyr::arrange(dat, chromosome_name, start_position)
#     armlev = unique(dat$arm)
#     dat = dplyr::mutate(dat, arm = factor(as.character(arm), levels = armlev))
#     dat
# }
# 
# dats = sapply(dats, function(dat) settingchrlevels(dat), simplify = F)
# hg19 = dats[[1]]
# hg38 = dats[[2]]
# mm10 = dats[[3]]
# usethis::use_data(hg19, hg38, mm10, internal = T)
