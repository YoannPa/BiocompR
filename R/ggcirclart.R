#' Circlizes ggplot2 objects.
#'
#' @param gg A \code{gg} plot object built with ggplot2.
#' @param facet_fun A \code{character} specifying the facet function to use for
#'                  the circlized plot, if the input gg object has facets. To
#'                  use facet_grid(): facet_fun = "grid"; To use facet_wrap():
#'                  facet_fun = 'wrap'.
#' @param wrap_to   A \code{character} specifying to which axis the facet should
#'                  be assigned, when the input gg object uses facet_wrap() and
#'                  facet_fun = 'grid'. Supported values: wrap_to = 'x' or 'y'.
#' @param wrap_nrow An \code{integer} specifying on how many rows the facets
#'                  should be wrapped when facet_fun = 'wrap'.
#' @param wrap_ncol An \code{integer} specifying on how many columns the facets
#'                  should be wrapped when facet_fun = 'wrap'.
#' @return A circlized \code{gg} object.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Most basic circlizing on a histogram generated with BiocompR::gghist()
#' distrib <- rnorm(100, 5, 1)
#' original_hist <- BiocompR::gghist(x = distrib)
#' original_hist
#' ggcirclart(gg = original_hist)
#' # Circlizing on a faceted histogram generated with BiocompR::gghist() using
#' # facet_wrap()
#' faceted_hist <- BiocompR::gghist(
#'     x = distrib, facet = rep(x = c("A", "B", "C", "D"), each = 25))
#' faceted_hist
#' ggcirclart(gg = faceted_hist, facet_fun = "wrap")
#' # Wrap on only 1 row
#' ggcirclart(gg = faceted_hist, facet_fun = "wrap", wrap_nrow = 1)
#' # Wrap on 4 columns
#' ggcirclart(gg = faceted_hist, facet_fun = "wrap", wrap_ncol = 4)
#' # Circlizing on a faceted histogram generated with BiocompR::gghist() using
#' # facet_grid()
#' ggcirclart(gg = faceted_hist, facet_fun = "grid")
#' # Map wrapping variable to Y axis (instead of X axis, the default)
#' ggcirclart(gg = faceted_hist, facet_fun = "grid", wrap_to = "y")
#' # Customize your circlized ggplot2 object
#' ggcirclart(gg = faceted_hist, facet_fun = "wrap") +
#'     theme(
#'         plot.title = element_text(hjust = 0.5),
#'         axis.text = element_text(size = 11, color = "black"),
#'         axis.title = element_text(size = 13),
#'         panel.grid.major.y = element_line(color = "black", linewidth = 0.2),
#'         panel.border = element_rect(colour = "black", fill = NA),
#'         strip.background = element_rect(color = "black"),
#'         strip.text = element_text(size = 10.5),
#'         panel.spacing = unit(0, "lines"))

ggcirclart <- function(
    gg, facet_fun = "wrap", wrap_to = "x", wrap_nrow = NULL, wrap_ncol = NULL){
    # Check if gg is a ggplot2
    if(!is.ggplot(gg)){
        stop("gg is not a ggplot2 object.")
    }
    # Get all geoms used in the gg object
    ls_geom <- lapply(X = gg$layers, FUN = function(l){
        class(l$geom)[1]
    })
    if("GeomBar" %in% ls_geom){
        # Barplot or Histogram
        cp <- coord_polar(theta = "x", start = 0)
        cp$is_free <- function() TRUE
        circlized_gg <- gg + cp
        # Check which facet function is used
        if(class(circlized_gg$facet)[1] == "FacetGrid"){
            cat("gg object uses 'facet_grid()'. 'wrap_to' parameter ignored.")
            if(facet_fun == "wrap"){
                if(length(circlized_gg$facet$params$cols) >= 1){
                    facet_var <- names(circlized_gg$facet$params$cols)
                } else if(length(circlized_gg$facet$params$rows) >= 1){
                    facet_var <- names(circlized_gg$facet$params$rows)
                } else {
                    stop("Something went wrong. Please contact developer.")
                }
                circlized_gg <- circlized_gg + theme(aspect.ratio = 1) +
                    facet_wrap(
                        eval(parse(text = paste0("~", facet_var))),
                        scales = "free_x", nrow = wrap_nrow, ncol = wrap_ncol) +
                    scale_y_reverse(expand = c(0, 0))
            } else if(facet_fun == "grid"){
                if(length(circlized_gg$facet$params$cols) >= 1){
                    facet_cols <- names(circlized_gg$facet$params$cols)
                } else { facet_cols <- NULL }

                if(length(circlized_gg$facet$params$rows) >= 1){
                    facet_rows <- names(circlized_gg$facet$params$rows)
                } else { facet_rows <- NULL }

                if(is.null(facet_cols) & is.null(facet_rows)){
                    stop("Something went wrong. Please contact developer.")
                } else if(!is.null(facet_cols) & is.null(facet_rows)){
                    circlized_gg <- circlized_gg + theme(aspect.ratio = 1) +
                        facet_grid(
                            cols = vars(eval(parse(text = facet_cols)))) +
                        scale_y_reverse(expand = c(0, 0))
                } else if(is.null(facet_cols) & !is.null(facet_rows)){
                    circlized_gg <- circlized_gg + theme(aspect.ratio = 1) +
                        facet_grid(
                            rows = vars(eval(parse(text = facet_rows)))) +
                        scale_y_reverse(expand = c(0, 0))
                } else if(!is.null(facet_cols) & !is.null(facet_rows)){
                    circlized_gg <- circlized_gg + theme(aspect.ratio = 1) +
                        facet_grid(
                            cols = vars(eval(parse(text = facet_cols))),
                            rows = vars(eval(parse(text = facet_rows)))) +
                        scale_y_reverse(expand = c(0, 0))
                }
            }
        } else if(class(circlized_gg$facet)[1] == "FacetWrap"){
            if(facet_fun == "wrap"){
                if(length(circlized_gg$facet$params$facets) >= 1){
                    facet_var <- names(circlized_gg$facet$params$facets)
                } else {
                    stop("Something went wrong. Please contact developer.")
                }

                circlized_gg <- circlized_gg + theme(aspect.ratio = 1) +
                    facet_wrap(
                        eval(parse(text = paste0("~", facet_var))),
                        scales = "free_x", nrow = wrap_nrow, ncol = wrap_ncol) +
                    scale_y_reverse(expand = c(0, 0))
            } else if(facet_fun == "grid"){
                if(length(circlized_gg$facet$params$facets) >= 1){
                    facet_var <- names(circlized_gg$facet$params$facets)
                } else {
                    stop("Something went wrong. Please contact developer.")
                }

                if(wrap_to == "x"){
                    circlized_gg <- circlized_gg + theme(aspect.ratio = 1) +
                        facet_grid(cols = vars(eval(parse(text = facet_var)))) +
                        scale_y_reverse(expand = c(0, 0))
                } else if(wrap_to == "y"){
                    circlized_gg <- circlized_gg + theme(aspect.ratio = 1) +
                        facet_grid(rows = vars(eval(parse(text = facet_var)))) +
                        scale_y_reverse(expand = c(0, 0))
                } else { stop("Unsupported value for 'wrap_to'.") }
            }
        } else if(class(circlized_gg$facet)[1] != "FacetNull"){
            stop("Something went wrong. Please contact developer.")
        }
    } else { stop("This geom is not supported yet. Please contact developer.") }
    return(circlized_gg)
}
