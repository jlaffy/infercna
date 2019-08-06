
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cna PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname diprange
#' @export 
diprange = function(cna, ...) {
    dots = list(...)
    v = sapply(dots, function(ref) rowMeans(cna[, ref, drop = F]), simplify = F)
    list(min = do.call(pmin, v), max = do.call(pmax, v))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param v PARAM_DESCRIPTION
#' @param Min PARAM_DESCRIPTION
#' @param Max PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname dipcenter
#' @export 
dipcenter = function(v, Min, Max) {
    above = v > Min & v > Max
    below = v < Min & v < Max
    normal = !above & !below
    v[above] <- v[above] - Max
    v[below] <- v[below] - Min
    v[normal] <- 0
    v
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cna PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname dipcorrect
#' @export 
dipcorrect = function(cna, ...) {
    dots = list(...)
    genes = rownames(cna)
    Args = c(list(cna = cna), dots)
    c(Min, Max) %<-% do.call(diprange, Args)
    n = nrow(cna)
    cna = t(sapply(1:n, function(i) dipcenter(cna[i, ], Min = Min[i], Max = Max[i])))
    rownames(cna) = genes
    cna
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param m PARAM_DESCRIPTION
#' @param dipcells PARAM_DESCRIPTION, Default: NULL
#' @param window PARAM_DESCRIPTION, Default: 100
#' @param range PARAM_DESCRIPTION, Default: c(-3, 3)
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname infercna
#' @export 
infercna = function(m,
                    dipcells = NULL,
                    window = 100,
                    range = c(-3, 3)) {

    # note: m should be row(gene)-centered
    cna = genorder(m)
    cna = clip(cna, range = range)
    cna = runMean(cna, k = window)
    cna = colCenter(cna)
    if (!is.null(dipcells)) {
        Args = c(list(cna = cna), dipcells)
        cna = do.call(dipcorrect, Args)}
    cna
}

