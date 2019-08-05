#!/usr/bin/env Rscript

source('~/R/map.R')
source('~/R/Chromosome.R')
source('~/R/solarized.R')

breakPoints = function(x, n, halfway = F) {
    Breaks = cumsum(table(x))
    BreaksNormalisedToScale = c(0, Breaks/max(Breaks) * n)

    if (halfway) {
        b = BreaksNormalisedToScale
        BreaksNormalisedToScale = c(b[1], b[2:length(b)] - diff(b)/2)
    }

    BreaksAtBoundaries = BreaksNormalisedToScale + 0.5
    BreaksAtBoundaries
}

orderByGroups = function(X, Groups = NULL, euclid.dist = T) {
    if (!is.null(Groups)) {
        Group = stack(Groups) %>%
            setNames(., c('Cell', 'Group'))
        X$Group = Group
    }

    dats = split(X, X$Group)
    mats = sapply(dats, function(d) reshape2::acast(d, Gene ~ Cell, value.var = 'CNA'), simplify = F)
    ord = NULL
    for (i in 1:length(mats)) {
        if (ncol(mats[[i]]) == 1) {
            ord = c(ord, colnames(mats[[i]]))
        } else {
            ordi = scrabble::clusta(mat = mats[[i]], compute.dist = euclid.dist, dist.method = 'euclidean', ORD = T)
            ord = c(ord, ordi)
        }
    }

    ord
}

prep.chrHeatMap = function(X = NULL, mat = NULL, Groups = NULL, euclid.dist = T) {

    if (!is.null(mat)) {
        X = reshape2::melt(as.matrix(mat)) %>% stats::setNames(c('Gene','Cell', 'CNA'))
    }

    if ('Group' %in% colnames(X)|is.null(Groups)) {
        cellOrder = orderByGroups(X = X, Groups = Groups)
        X = X %>%
            dplyr::mutate(Cell = factor(as.character(Cell), levels = cellOrder)) %>%
            dplyr::arrange(Cell)
    }

    X
}

chrHeatMap = function(X = NULL,
                      limits = c(-1, 1),
                      lim.find = F,
                      lim.sym = T,
                      x.name = 'Chromosome',
                      y.name = 'Cell',
                      angle = NULL,
                      x.angle = NULL,
                      y.angle = 0,
                      axis.rel = 1,
                      base.size = 12,
                      axis.title.size = 12,
                      axis.text.size = 11,
                      base_col = solarized()['base02'],
                      title = NULL,
                      subtitle = NULL,
                      caption = NULL,
                      text.size = 12,
                      x.hide = c('13', '18', '21', 'Y'),
                      y.hide = NULL,
                      ratio = NULL,
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
    
    cols = readRDS('~/rds/hotmap.rds')

    if ('Group' %in% colnames(X)) {

        ylinebreaks = breakPoints(x = X$Group, n = length(levels(X$Cell))) 
        ylabelbreaks = breakPoints(x = X$Group, n = length(levels(X$Cell)), halfway = T) 

        if (!is.null(y.hide)) {
            names(ylabelbreaks)[which(names(ylabelbreaks) %in% y.hide)] = rep('', length(y.hide))
        }

        ylabels = names(ylabelbreaks)
    }
    
    else {
        X$Cell = as.numeric(X$Cell)
        ylinebreaks = NULL
        ylabelbreaks = waiver()
        ylabels = waiver()
    }

    xlabelbreaks = chromosomeBreakPoints(levels(X$Gene), stripped = T, halfway = T)
    xlinebreaks = chromosomeBreakPoints(levels(X$Gene))

    if (!is.null(x.hide)) {
        x.hide = stringr::str_replace(x.hide, "chr|Chr|CHR", "")
        names(xlabelbreaks)[which(names(xlabelbreaks) %in% x.hide)] = rep('', length(x.hide))
    }

    if (lim.find) {
        v = X %>% pull(CNA) %>% range
        limits = .limits(v = v, symmetric = lim.sym)
    }

    if (!is.null(angle)) {
        x.angle = angle
        y.angle = angle
    }

    if (isTRUE(x.angle)) {
        x.angle = 45
    }

    if (!is.null(x.angle)) {
        x.angle = ggpubr::rotate_x_text(angle = x.angle)
    }

    if (isTRUE(y.angle)) {
        y.angle = 90
    }

    if (!is.null(y.angle)) {
        y.angle = ggpubr::rotate_y_text(angle = y.angle)
    }

    legend = .axis.spacer(breaks = legend.breaks, labels = legend.labels, limits = limits)
    legend.breaks = legend$breaks
    legend.labels = legend$labels
    
    geomtile = geom_tile(col = tile.col, size = tile.size)
    if (any(sapply(list(tile.size, tile.col), is.null))) {
        geomtile = ggplot2::geom_tile()
    }

    G = ggplot(X, aes(x = as.numeric(Gene), y = as.numeric(Cell), fill = CNA)) +
        geomtile +
        ggplot2::scale_fill_gradientn(colors = cols,
                                      limits = limits,
                                      expand = expand,
                                      oob = scales::squish,
                                      breaks = legend.breaks,
                                      labels = legend.breaks,
                                      name = legend.title,
                                      guide = guide_colorbar(frame.colour = 'black',
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
        scale_x_continuous(expand = c(0, 0), breaks = xlabelbreaks, labels = names(xlabelbreaks)) +
        scale_y_continuous(expand = c(0, 0), breaks = ylabelbreaks, labels = ylabels) +
        y.angle +
        x.angle +
        geom_vline(xintercept = xlinebreaks, size = 0.08, linetype = 2) +
        theme(title = element_text(colour = base_col),
              rect = element_rect(colour = base_col),
              line = element_line(colour = base_col),
              panel.border = element_rect(colour = 'black', fill = NA),
              legend.text = element_text(size = rel(legend.rel)),
              legend.title = element_text(size = rel(legend.title.rel), angle = legend.title.angle),
              plot.margin = unit(c(1, 0, 0.3, 0.3), "cm"),
              legend.margin = margin(0,0.3,0,-0.3,'cm'),
              legend.justification = legend.justification) 
    
    if (!is.null(ylinebreaks)) G = G + geom_hline(yintercept = ylinebreaks[-1], size = 0.1)
    return(list(X = X, P = G))
}
