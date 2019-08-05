
expMatrix = function(x = '125') {

    if (x = '125') {
        m = infercna::mgh125
        cellnm = infercna::cells125
    }

    else if (x == "771") {
        m = infercna::bt771
        cellnm = infercna::cells771
    }

    m = as.matrix(m)[infercna::genes, ]
    colnames(m) = cellnm
    scrabble::rowCenter(m)
}
