
.ggcna_prep = function(m,
                       reorder = FALSE,
                       reorder.by = NULL,
                       groups = NULL,
                       interorder = TRUE,
                       genome = 'hg38',
                       dist.method = 'euclidean',
                       cluster.method = 'complete') {

    if (!reorder) {
        return(m)
    }

    if (!is.null(reorder.by)) {
        reorder.by = as.character(reorder.by)
        useGenome(genome)
        geneList = c(splitGenes(rownames(m)), splitGenes(rownames(m), by = 'arm'))
        genes = geneList[names(geneList) %in% reorder.by] %>% unlist
        rm(geneList)
    } else {
        genes = rownames(m)
    }

    if (is.null(groups)) {
        cols = scalop::hca_order(m[genes, ],
                         cluster.method = cluster.method,
                         dist.method = dist.method)
    } else {
        cols = grouped_reorder(m[genes,],
                               groups = groups,
                               interorder = interorder,
                               cluster.method = cluster.method,
                               dist.method = dist.method) %>%
            colnames
    }

    m[, cols]
}


#' @title Plot a CNA heatmap
#' @description Uses `ggplot::geom_raster` to generate a plot of copy number aberration (CNA) values. Delineates chromosomes and chromosome arms, and provides options for clustering cells or groups of cells (e.g., cell types, samples, or subclones) prior to plotting.
#' @param m CNA matrix (genes by cells). <m> can be generated using `infercna::infercna`.
#' @param reorder reorder cells using hierarchical clustering, Default: FALSE
#' @param ... other arguments to pass to `scalop::graster`.
#' @param reorder.by reorder cells by genes on select chromosomes or chromosome arms, e.g., c("7", "10", "1p", "19q"). Default: NULL
#' @param groups groups of cells to delineate and label on the plot (e.g., cell types, samples, or subclones). If <reorder> is TRUE, ordering will be performed within groups, Default: NULL
#' @param interorder reorder by hierarchical clustering betweeen groups, Default: TRUE
#' @param genome set genome to use ('hg19' or 'hg38'), Default: 'hg19'
#' @param dist.method distance metric for reordering, Default: 'euclidean'
#' @param cluster.method linkage method for reordering, Default: 'complete'
#' @param hide hide x-axis labels for specific chromosomes (e.g., small chromosomes whose labels would overlap flanking labels), Default: c("21", "Y")
#' @param limits for the colour key; replaces out of bounds values with the nearest limit. Default: c(-1, 1)
#' @param y.angle y-axis labels' angle, Default: 90
#' @param legend.title legend title, Default: NULL
#' @param title plot title, Default: 'Copy-number aberrations'
#' @param axis.text.size x and y axes label size, Default: 12
#' @param legend.height legend bar height, Default: 0.3
#' @param legend.width legend bar width, Default: 0.5
#' @param cols custom colour palette (character vector), Default: NULL
#' @return a `ggplot` object
#' @examples 
#' m = infercna(mgh125, isLog=T, refCells=refCells)
#' malCells = list(Malignant=setdiff(colnames(mgh125), unlist(refCells)))
#' groups = c(malCells, refCells)
#' p = ggcna(m, reorder=T, groups=refCells)
#' @rdname ggcna
#' @export 
ggcna = function(m,
                 reorder = FALSE,
                 reorder.by = NULL,
                 groups = NULL,
                 interorder = TRUE,
                 genome = 'hg19',
                 dist.method = 'euclidean',
                 cluster.method = 'complete',
                 hide = c('21','Y'),
                 limits = c(-1, 1),
                 y.angle=90,
                 legend.title = NULL,
                 title = 'Copy-number aberrations',
		 axis.text.size = 12,
                 legend.height = 0.3,
                 legend.width = 0.5,
                 cols = NULL,
                 ...) {

    if (is.null(cols)) {
        cols = c(
                 "#053061",
                 "#2166AC",
                 "#4393C3",
                 "#92C5DE",
                 "#D1E5F0",
                 "#F7F7F7",
                 "white",
                 "#F7F7F7",
                 "#FDDBC7",
                 "#F4A582",
                 "#D6604D",
                 "#B2182B",
                 "#67001F")
    }
    # geom vlines,breaks,labels for chr and chr arms
    chr.size = 0.1
    arm.size = 0.05
    genes = rownames(m)
    useGenome(genome)
    m = orderGenes(m)
    L = splitGenes(genes)
    xints = cumsum(lengths(L)) + chr.size
    L2 = splitGenes(genes, by = 'arm')
    xints2 = cumsum(lengths(L2)) + arm.size
    xints2 = xints2[str_detect(names(xints2), "p")]
    breaks = (xints-chr.size) - lengths(L)/2
    labels = names(xints)
    labels[names(xints) %in% hide] <- ''

    m = .ggcna_prep(m,
                    reorder = reorder,
                    reorder.by = reorder.by,
                    groups = groups,
                    genome = genome,
                    interorder = interorder,
                    dist.method = dist.method,
                    cluster.method = cluster.method) 

    if (!is.null(groups)) {
        line.size = 0.1
        ord = colnames(m)
        groupNames = flip(Unlist(groups))[ord]
        groupNames = groupNames[!duplicated(groupNames)]
        yints = cumsum(lengths(groups)) + line.size
        ybreaks = (yints-line.size) - lengths(groups)/2
        ylabels = names(yints)
    } else {
        ylabels = ggplot2::waiver()
        ybreaks = ggplot2::waiver()
    }

    d = reshape2::melt(as.matrix(m)) %>% 
	    dplyr::mutate(Var1=as.numeric(factor(as.character(Var1),
						 levels=rownames(m))),
			  Var2=as.numeric(factor(as.character(Var2),
						 levels=colnames(m))))

    G = graster(d,
		x=Var1,
		y=Var2,
		fill=value,
                limits = limits,
                legend.title = legend.title,
                legend.height = legend.height,
                legend.width = legend.width,
                ...) +
        ggplot2::geom_vline(xintercept = xints,
                   col = 'grey20',
                   linetype = 1,
                   size = chr.size) +
        ggplot2::geom_vline(xintercept = xints2,
                   col = 'grey20',
                   linetype = 2,
                   size = arm.size) +
        ggplot2::scale_x_continuous(breaks = breaks,
                           labels = labels,
                           expand = c(0,0)) +
        ggplot2::scale_y_continuous(breaks = ybreaks,
                           labels = ylabels,
                           expand = c(0,0)) +
    scalop::theme_scalop(legend.text.size = 14) +
    ggplot2::theme(plot.margin = margin(.2,.2,.2,.2,'cm'),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
	           axis.text = ggplot2::element_text(size=axis.text.size),
	           axis.ticks.x = ggplot2::element_blank(),
                   legend.key.height = grid::unit(legend.height, 'cm'),
                   legend.key.width = grid::unit(legend.width,'cm')) +
    ggplot2::scale_fill_gradientn(name = legend.title,
                         	  oob = scales::squish,
                         	  colors = cols,
                         	  breaks = c(limits[1],0,limits[2]),
                         	  limits = limits,
                         	  guide = ggplot2::guide_colorbar(frame.colour = 'black',
								  ticks.colour = 'black')) +
    labs(title = title) 


    if (!is.null(groups)) {
        G = G + 
		ggplot2::geom_hline(yintercept = yints, col = 'grey20', size = line.size,linetype=1) +
		ggplot2::theme(axis.text.y = ggplot2::element_text(angle = y.angle, hjust =0.5),
		      axis.ticks.y=ggplot2::element_blank())
    }

    G
}
