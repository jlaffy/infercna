

separate_unassigned = function(modes_by_arm) {

}

cnaclones = function() {
    mats_by_arm = split_by_chromosome_arm(mat)
    modes_by_arm = sapply(mats_by_arm, function(m) {
                              try(modes(m = m,
                                        prob = prob,
                                        coverage = coverage,
                                        size = size,
                                        exclude.unassigned = F)),
                          simplify = F)
    modes_by_arm = modes[!which(sapply(modes_by_arm, class) == 'try-error')]
    message('Found modes for chromosome arms:\n', paste(names(modes_by_arm), sep = '\n'))
    clone_ids = Reduce(paste0, sapply(modes_by_arm, names))
    message('Found ', unique(length(clone_ids)))


}
