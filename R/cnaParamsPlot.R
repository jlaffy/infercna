
#' @title Plot CNA parameters to evaluate cell CNAs 
#' @description Plot cells according to their CNA correlation and CNA signal values, and display this as a contour plot. If a "sample" column is provided, the resulting plot will be faceted by sample, with the sample name written in the plot. A "group" column is currently required (suggested to just make a group column with all values the same if there are no real groups annotated).
#' @param data data.frame with x (CNA correlations), y (CNA signals), group (e.g. presumed Malignant versus Reference Normal) and (optional) sample columns
#' @param x column name corresponding to x values. Either passed as a string or as an unquoted column name. eg 'x' or x is fine. 
#' @param y column name corresponding to y values. Either passed as a string or as an unquoted column name. eg 'y' or y is fine. 
#' @param group column name corresponding to groups. Either passed as a string or as an unquoted column name. eg 'group' or group is fine. 
#' @param points logical; should cell points be plotted as well as the contour plot? Default: TRUE
#' @param point.col relevant only if points = TRUE. specify point colour. Default: 'black'
#' @param point.alpha relevant only if points = TRUE. specify point transparency. Default: 0.075
#' @param point.size relevant only if points = TRUE. specify point size. Default: 0.5
#' @param xlab x-axis title. Default: 'CNA correlation'
#' @param ylab y-axis title. Default: 'CNA signal'
#' @param xlim x-axis limits. Default: NULL
#' @param ylim y-axis limits. Default: NULL
#' @param xintercept PARAM_DESCRIPTION, Default: NULL
#' @param yintercept PARAM_DESCRIPTION, Default: NULL
#' @param bins number of bins to use for contour levels. Fewer bins averages over more of the data but often results with a 'cleaner' looking contour. Default: 7
#' @param breaks breaks. Default: scales::breaks_pretty(n = 4)
#' @param xbreaks x-axis breaks. Default: breaks
#' @param ybreaks y-axis breaks. Default: ggplot2::waiver()
#' @param contour contour lines as well as filled in contour areas? Default: TRUE
#' @param geom I've only used polygon in this specific function but you can go crazy and try other things. Default: 'polygon'
#' @param legend.justification legend justification. Default: 'right'
#' @param legend.position legend position. set to 'none' if not desired. Default: 'top'
#' @param legend.title legend title Default: ggplot2::waiver()
#' @param cols colours of each group. Note that if more than two groups are present, more colours should be provided than the default two. Default: c("red", "blue")
#' @param text.size base text size. Default: 16
#' @param label.hjust only relevant if sample column in data. Adjust position of geom_text with sample name. Default: 1.25
#' @param label.hjust only relevant if sample column in data. Adjust position of geom_text with sample name. Default: 1.25
#' @param label.size only relevant if sample column in data. geom_text label size. Default: 5
#' @param ... other arguments passed to ggplot2::facet_wrap(.~sample, ...) if sample column in data.
#' @return ggplot2 object 
#' @seealso 
#'  \code{\link[scales]{breaks_pretty}}
#'  \code{\link[ggplot2]{waiver}},\code{\link[ggplot2]{tidyeval}},\code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{geom_density_2d}},\code{\link[ggplot2]{geom_label}},\code{\link[ggplot2]{scale_continuous}},\code{\link[ggplot2]{scale_alpha}},\code{\link[ggplot2]{scale_manual}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{c("guide_bins", "guide_colourbar", "guide_coloursteps", "guide_legend", "guides", "guides")}},\code{\link[ggplot2]{guide_legend}},\code{\link[ggplot2]{facet_wrap}}
#' @rdname cnaParamsPlot
#' @export 
#' @importFrom scales breaks_pretty
#' @importFrom ggplot2 waiver ensym ggplot aes geom_point stat_density_2d geom_text scale_x_continuous scale_y_log10 scale_alpha scale_fill_manual scale_colour_manual theme element_blank element_text guides guide_legend facet_wrap
cnaParamsPlot = function(data,
                         x,
                         y,
                         group = NULL,
                         points = TRUE,
                         point.col = 'black',
                         point.alpha = 0.075,
                         point.size = 0.5,
                         xlab = 'CNA correlation',
                         ylab = 'CNA signal',
                         xlim = NULL,
                         ylim = NULL,
                         xintercept = NULL,
                         yintercept = NULL,
                         bins = 30,
                         breaks = scales::breaks_pretty(n = 4),
                         xbreaks = breaks,
                         ybreaks = ggplot2::waiver(),
                         contour = TRUE,
                         geom = 'polygon',
                         legend.justification = 'right',
                         legend.position = 'top',
                         legend.title = ggplot2::waiver(),
                         cols = c('red', 'blue'),
                         text.size = 16,
                         label.hjust = 1.25,
                         label.vjust = 1.25,
                         label.size = 5,
                         ...) {

    x = ggplot2::ensym(x)
    y = ggplot2::ensym(y)
    group = try(ggplot2::ensym(group))

    G = ggplot2::ggplot(data, ggplot2::aes(x = !!x,
                                           y = !!y,
                                           group = !!group)) 

    if (isTRUE(points)) {
        G = G +
            ggplot2::geom_point(ggplot2::aes(col = group),
                                size = point.size,
                                alpha = point.alpha,
                                shape = 16)
    }

    G = G +
        ggplot2::stat_density_2d(ggplot2::aes(alpha = ..level..,
                                              fill = group,
                                              col = group),
                                 geom = geom,
                                 contour = contour,
                                 bins = bins) +
        ggplot2::geom_text(ggplot2::aes(label = sample),
                           hjust = label.hjust,
                           vjust = label.vjust,
                           x = Inf,
                           y = Inf,
                           size = label.size) + 
        ggplot2::scale_x_continuous(expand = c(0,0),
                                    limits = xlim,
                                    name = xlab,
                                    breaks = xbreaks) +
        ggplot2::scale_y_log10(expand = c(0,0),
                               limits = ylim,
                               name = ylab,
                               breaks = ybreaks) +
        ggplot2::scale_alpha(guide = 'none') +
        ggplot2::scale_fill_manual(values = cols, name = legend.title) +
        ggplot2::scale_colour_manual(values = cols, name = legend.title) +
        theme_scalop() +
        ggplot2::theme(strip.text = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(size = 13),
                       legend.text = ggplot2::element_text(size = 13),
                       legend.justification = legend.justification,
                       legend.position = legend.position) +
        ggplot2::guides(fill = ggplot2::guide_legend(title.position = 'top', hjust = .5)) +
        ggplot2::facet_wrap(.~sample, ...)

    G
}
