
cnamalignant = function(cnamat, cor.cutoff = 0.4, signal.cutoff = 0.02) {
    ncells = ncol(cnamat)
    cnacor(cnamat) >= cor.cutoff
    cnasignal(cnamat) >= signal.cutoff
    ismalignant = cor.cutoff & signal.cutoff
    stats::setNames(ismalignant, colnames(cnamat))
}
