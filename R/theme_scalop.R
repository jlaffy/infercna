
#' @title scalop ggtheme 
#' @description theme for ggplot2 
#' @param size base text size. Default: 16
#' @param legend.title.size legend title text size. Default: 13
#' @param legend.text.size legend text size. Default: 12
#' @param legend.position legend position. Default: 'right'
#' @param legend.justification legend justification. Default: 'top'
#' @return function to be added to ggplot2 call 
#' @seealso 
#'  \code{\link[ggplot2]{ggtheme}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}}
#'  \code{\link[grid]{unit}}
#' @rdname theme_scalop
#' @export 
#' @importFrom ggplot2 theme_bw theme element_text element_blank
#' @importFrom grid unit
theme_scalop = function(size = 16,
                        legend.title.size = 13,
                        legend.text.size = 12,
                        legend.position = 'right',
                        legend.justification = 'top') {

    ggplot2::theme_bw() +
        ggplot2::theme(text = ggplot2::element_text(size = size),
                       strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_text(size = 16),
                       legend.text = ggplot2::element_text(size = legend.text.size),
                       legend.title = ggplot2::element_text(size = legend.title.size),
                       legend.position = legend.position,
                       legend.justification = legend.justification,
                       panel.grid.minor = ggplot2::element_blank(),
                       axis.ticks.length = grid::unit(-0.25, 'cm'),
                       axis.text.x = ggplot2::element_text(margin = margin(t = 0.35, unit = 'cm')),
                       axis.text.y = ggplot2::element_text(margin = margin(r = 0.35, unit = 'cm')))
}

