
new = function(cnamat) {
    arms = infercna::splitMatrixByArm(cnamat)
    arms = arms[scrabble::nrows(arms) >= 5]
    armodes = sapply(armats, modes, exclude.unassigned = FALSE, simplify = FALSE)
    armodes = armodes[lengths(armodes) == 0]
    armodes = armodes[lengths(unlist(sapply(armodes, unique))) >= 2]
    armodes = armodes[lengths(armodes) != 0]
    msg = paste0('Found CNA modes on chromosome arms:\n',
                 paste(names(armodes), sep = '\n')) 
    message(msg)
    browser()
}
