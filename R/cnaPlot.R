
.getAttributeName = function(val) {
    if (val %in% levels(Genv$chr)) return('chr')
    if (val %in% levels(Genv$arm)) return('arm')
    if (val %in% Genv$symbol) return('symbol')
    else stop('Value ', val, ' not found in genome attributes.')
}

.orderCells = function(cna, by_groups = NULL, include = NULL, euclid.dist = F, ...) {
    cna.orig = cna

    if (!is.null(include)) {
        atts = sapply(include, .getAttributeName)
        Args = split(include, atts)
        cna = filterGenes(x = cna, value = Args, attribute = names(Args))
    }
    
    if (!is.null(by_groups)) {
        ord = colnames(scalop::grouped_reorder(m = cna, groups = cell.groups, interorder = F))
    } else {
        ord = colnames(scalop::hca_order(x = rowcenter(cna)))
    }

    cna.orig[, ord]
}

.axisSpacer = function(breaks, labels, limits, levels = NULL) {
    if (!is.null(labels) & !is.null(levels)) {
        breaks = levels %in% labels
        labels = levels[breaks]
    }
    if (is.null(breaks)) breaks = seq(limits[[1]], limits[[2]], limits[[2]])
    if (is.null(labels)) labels = breaks
    list(breaks = breaks, labels = labels)
}

.cellBreaks = function(groups, halfway = F, hide = NULL) {
    n = length(unlist(groups))
    chrsum = cumsum(lengths(groups))
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
        ggplot2::scale_y_continuous(expand = c(0, 0), breaks = yTextBreaks, labels = names(yTextBreaks)) +
        y.angle +
        x.angle +
        ggplot2::geom_vline(xintercept = xLineBreaks, size = 0.08, linetype = 2) +
        ggplot2::geom_hline(yintercept = yLineBreaks, size = 0.08, linetype = 1) +
        ggplot2::theme(aspect.ratio = ratio,
                       title = ggplot2::element_text(colour = base.col),
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


#' @title Plot a Heatmap of CNA Values
#' @description Plot a heatmap of CNA values.
#' @param cna matrix of CNA values (genes by cells).
#' @param limits colour range. Cells >= upper limit will be the same colour, and likewise for cells <= lower limit. Default: c(-1, 1)
#' @param ratio aspect ratio of the panel. Default: 0.7
#' @param cols character vector of colours to use. Default: heatCols
#' @param orderCells boolean value indicating whether cells should be ordered by hierarchical clustering. Default: F
#' @param x.name x axis label. Default: 'Chromosome'
#' @param y.name y axis label. Default: 'Cell'
#' @param angle angle of axes tick labels. x.angle and y.angle inherit from angle. Default: NULL
#' @param x.angle angle of x axis tick labels. If left, will inherit from angle. Default: NULL
#' @param y.angle angle of y axis tick labels. If left, will inherit from angle. Default: 0
#' @param axis.rel relative size of axes labels. Default: 1
#' @param base.size base text size. Default: 12
#' @param axis.title.size axes titles text size. Default: 12
#' @param axis.text.size axes labels text size. Default: 11
#' @param base.col base text colour. Default: '#073642'
#' @param title title of plot. Default: NULL
#' @param subtitle subtitle of plot. Default: NULL
#' @param caption caption of plot. Default: NULL
#' @param text.size base text size. Default: 12
#' @param x.hide x axis labels to hide from plot. Default: c("13", "18", "21", "Y")
#' @param y.hide y axis labels to hide from plot. Default: NULL
#' @param tile.size size of tile borders. No border if is equal to NULL. Default: 0.1
#' @param tile.col colour of tile borders. No border if is equal to NULL. Default: NULL
#' @param legend.position legend position on plot. Default: 'right'
#' @param legend.height legend height. Default: 2
#' @param legend.width legend width. Default: 0.6
#' @param legend.rel relative size of legend text. Default: 0.9
#' @param legend.colour legend text colour. Default: 'black'
#' @param legend.breaks breaks to include on colour key. Default: NULL
#' @param legend.labels specify labels on colour key. Default: NULL
#' @param legend.title title of legend. Default: 'Inferred CNA `[log2 ratio`]'
#' @param legend.justification legend justification on plot. Default: 'top'
#' @param legend.title.position position of legend title on plot. Default: 'bottom'
#' @param legend.title.angle angle of legend title. Default: NULL
#' @param legend.title.rel relative size of legend title text. Default: 0.9
#' @return a list containing $p, the ggplot plot, and $data, the data (in dataframe format).
#' @seealso 
#'  \code{\link[reshape2]{melt}}
#'  \code{\link[ggpubr]{rotate_axis_text}}
#' @rdname cnaPlot
#' @export 
#' @importFrom reshape2 melt
#' @importFrom ggpubr rotate_x_text rotate_y_text
cnaPlot = function(cna,
                   limits = c(-1, 1),
                   ratio = 0.5,
                   cols = heatCols, 
                   x.name = 'Chromosome',
                   y.name = 'Cell',
                   legend.title = 'Inferred CNA\n[log2 ratio]',
                   x.hide = c('13', '18', '21', 'Y'),
                   orderCells = F,
                   cell.groups = scalop::split_by_sample_names(colnames(m)),
                   order.with = NULL,
                   euclid.dist = F,
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
                   legend.justification = 'top',
                   legend.title.position = 'bottom',
                   legend.title.angle = NULL,
                   legend.title.rel = 0.9) {

    # order genes by genomic position
    cna = orderGenes(cna)
    # order cells?
    if (orderCells) {
        cna = .orderCells(cna,
                          include = order.with,
                          by_groups = cell.groups,
                          euclid.dist = euclid.dist)
    }

    # prepare dataframe
    dat = reshape2::melt(as.matrix(cna))
    colnames(dat) = c('Gene', 'Cell', 'CNA')
    ylineBreaks = ggplot2::waiver()
    yTextBreaks = ggplot2::waiver()

    if (orderCells) {
        cell.groups = split_by_sample_names(levels(dat$Cell))
        # breaks for horizontal lines denoting cell groups
        yLineBreaks = .cellBreaks(cell.groups)
        # and breaks for their labels
        yTextBreaks = .cellBreaks(cell.groups, halfway = TRUE, hide = x.hide)
    }

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
    list(p = G, data = tibble::as_tibble(dat))
}
