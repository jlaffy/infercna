
.axisSpacer <- function(breaks, labels, limits, levels = NULL) {
    if (!is.null(labels) & !is.null(levels)) {
        breaks = levels %in% labels
        labels = levels[breaks]
    }
    if (is.null(breaks)) breaks = seq(limits[[1]], limits[[2]], limits[[2]])
    if (is.null(labels)) labels = breaks
    list(breaks = breaks, labels = labels)
}

.chromosomeBreaks = function(genes = NULL, halfway = F, hide = NULL) {
    n = length(genes)
    chrsum = cumsum(lengths(splitGenes(genes, by = 'chr')))
    Breaks = chrsum/max(chrsum) * n
    # halfway useful for placing x axis chromosome text labels
    # halfway = F for placing x axis chromosome lines
    if (halfway) {
        b = Breaks
        Breaks = c(b[1]/2, b[2:length(b)] - diff(b)/2)
    }
    if (!is.null(hide)) {
        names(Breaks)[which(names(Breaks) %in% hide)] = rep('', length(hide))
    }
    Breaks
}


.cnaPlot = function(...) {
    list2env(list(...), envir = environment())
    geomtile = ggplot2::geom_tile(col = tile.col, size = tile.size)
    if (any(sapply(list(tile.size, tile.col), is.null))) {
        geomtile = ggplot2::geom_tile()
    }

    G = ggplot2::ggplot(dat, ggplot2::aes(x = as.numeric(Gene),
                                          y = as.numeric(Cell),
                                          fill = CNA)) +
        geomtile +
        ggplot2::scale_fill_gradientn(colors = cols,
                                      limits = limits,
                                      oob = scales::squish,
                                      breaks = legend.breaks,
                                      labels = legend.breaks,
                                      name = legend.title,
                                      guide = ggplot2::guide_colorbar(frame.colour = 'black',
                                                                      ticks = F,
                                                                      barheight = grid::unit(legend.height, "cm"),
                                                                      barwidth = grid::unit(legend.width, "cm"),
                                                                      title.position = legend.title.position)) +
        ggplot2::labs(x = x.name,
                      y = y.name,
                      title = title,
                      subtitle = subtitle,
                      caption = caption) +
        hrbrthemes::theme_ipsum(base_size = base.size,
                                axis_title_size = axis.title.size,
                                axis_text_size = axis.text.size) +
        ggplot2::scale_x_continuous(expand = c(0, 0), breaks = xTextBreaks, labels = names(xTextBreaks)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        y.angle +
        x.angle +
        ggplot2::geom_vline(xintercept = xLineBreaks, size = 0.08, linetype = 2) +
        ggplot2::theme(title = ggplot2::element_text(colour = base.col),
                       rect = ggplot2::element_rect(colour = base.col),
                       line = ggplot2::element_line(colour = base.col),
                       panel.border = ggplot2::element_rect(colour = 'black', fill = NA),
                       legend.text = ggplot2::element_text(size = ggplot2::rel(legend.rel)),
                       legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.rel), angle = legend.title.angle),
                       plot.margin = grid::unit(c(1, 0, 0.3, 0.3), "cm"),
                       legend.margin = ggplot2::margin(0,0.3,0,-0.3,'cm'),
                       legend.justification = legend.justification) 
    
    G
}


cnaPlot = function(cna,
                   limits = c(-1, 1),
                   ratio = 0.7,
                   cols = heatCols, 
                   orderCells = F,
                   x.name = 'Chromosome',
                   y.name = 'Cell',
                   angle = NULL,
                   x.angle = NULL,
                   y.angle = 0,
                   axis.rel = 1,
                   base.size = 12,
                   axis.title.size = 12,
                   axis.text.size = 11,
                   base.col = "#073642",
                   title = NULL,
                   subtitle = NULL,
                   caption = NULL,
                   text.size = 12,
                   x.hide = c('13', '18', '21', 'Y'),
                   y.hide = NULL,
                   tile.size = 0.1,
                   tile.col = NULL,
                   legend.position = 'right',
                   legend.height = 2,
                   legend.width = 0.6,
                   legend.rel = 0.9,
                   legend.colour = 'black',
                   legend.breaks = NULL,
                   legend.labels = NULL,
                   legend.title = 'Inferred CNA\n[log2 ratio]',
                   legend.justification = 'top',
                   legend.title.position = 'bottom',
                   legend.title.angle = NULL,
                   legend.title.rel = 0.9) {

    # order genes by genomic position
    cna = orderGenes(cna)
    # order cells?
    if (orderCells) cna = cna[, scrabble::clusta(cna, ord = T, compute.dist = F)]
    # prepare dataframe
    dat = reshape2::melt(as.matrix(cna))
    colnames(dat) = c('Gene', 'Cell', 'CNA')
    # breaks for vertical lines denoting chromosomes
    xLineBreaks = .chromosomeBreaks(levels(dat$Gene))
    # and breaks for their labels
    xTextBreaks = .chromosomeBreaks(levels(dat$Gene), halfway = TRUE, hide = x.hide)
    if (!is.null(angle)) {
        x.angle = angle
        y.angle = angle
    }
    if (isTRUE(x.angle)) x.angle = 45
    if (!is.null(x.angle)) x.angle = ggpubr::rotate_x_text(angle = x.angle)
    if (isTRUE(y.angle)) y.angle = 90
    if (!is.null(y.angle)) y.angle = ggpubr::rotate_y_text(angle = y.angle)
    legend = .axisSpacer(breaks = legend.breaks, labels = legend.labels, limits = limits)
    legend.breaks = legend$breaks
    legend.labels = legend$labels
    # arguments to be passed to .cnaPlot()
    Args = mget(ls()[-1 * which(ls() %in% c('cna', 'x.hide', 'orderCells'))])
    # make ggplot object
    G = do.call(.cnaPlot, Args)
    list(p = G, data = dat)
}
