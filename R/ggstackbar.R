#' Draws stacked barplots from an annotation table
#'
#' @param dt       A \code{data.table} or a \code{data.frame} containing 1
#'                 annotation per column.
#' @param cols     An \code{integer} or a \code{character} vector matching
#'                 columns in 'dt'. If 'cols' is a vector of integers, integers
#'                 must match columns position in the table. If 'cols' is a
#'                 vector of characters, characters must match columns names.
#'                 If 'cols' is NULL (the default), all columns are retrieved
#'                 from 'dt'.
#' @param y_breaks An \code{integer} specifying the number of graduations on the
#'                 Y axis (Default: y_breaks = 11).
#' @param pal      A \code{character} vector matching existing color codes to be
#'                 used as a palette for the stacked barplot.
#' @param lgd_nrow An \code{integer} specifying the number of rows on which the
#'                 legend keys should be displayed (Default: lgd_nrow = 35).
#' @param rm.null  A \code{logical} specifying whether the NULL elements
#'                 generated in the final list by ggstackbar() should be
#'                 automatically removed (rm.null = TRUE), or kept in the final
#'                 result (Default: rm.null = FALSE).
#' @param verbose  A \code{logical} specifying whether the function should
#'                 verbose information while executing
#'                 (Default: verbose = TRUE).
#' @return A \code{gg} objects list, each element being an individual stacked
#'         barplot.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Basic use of ggstackbar
#' res_stackbar <- ggstackbar(dt = ggplot2::mpg)
#' # Simply explore results: to see the stacked barplot of manufacturers, do:
#' res_stackbar$manufacturer
#' # Only generate stacked barplots for columns 'manufacturer', 'model' & 'year'
#' ls_res <- ggstackbar(
#'     dt = ggplot2::mpg, cols = c('manufacturer', 'model', 'year'))
#' ls_res$manufacturer
#' # Increase the number of graduation on the Y axis
#' ls_res <- ggstackbar(
#'     dt = ggplot2::mpg, cols = c('manufacturer', 'model', 'year'),
#'     y_breaks = 21)
#' ls_res$manufacturer
#' # Set custom color palette
#' ls_res <- ggstackbar(
#'     dt = ggplot2::mpg, cols = c('manufacturer', 'model', 'year'),
#'     y_breaks = 21, pal = ggsci::pal_npg(palette = c("nrc"), alpha = 1)(10))
#' ls_res$manufacturer
#' # Set number of lines on which to display legend's keys to 5 lines
#' ls_res <- ggstackbar(
#'     dt = ggplot2::mpg, cols = c('manufacturer', 'model', 'year'),
#'     y_breaks = 21, pal = ggsci::pal_npg(palette = c("nrc"), alpha = 1)(10),
#'     lgd_nrow = 5)
#' ls_res$manufacturer
#' # Arrange the output result in a single plot
#' ls_res <- ggstackbar(
#'     dt = ggplot2::mpg, cols = c('manufacturer', 'model', 'year'),
#'     y_breaks = 21, pal = ggsci::pal_npg(palette = c("nrc"), alpha = 1)(10),
#'     lgd_nrow = 38)
#' egg::ggarrange(plots = ls_res, nrow = 1)

ggstackbar <- function(
    dt, cols = NULL, y_breaks = 11, pal = grDevices::rainbow(n = 10),
    lgd_nrow = 35, rm.null = FALSE, verbose = TRUE){
    if(!data.table::is.data.table(dt)){ dt <- data.table::as.data.table(dt) }
    if(is.null(cols)){ cols <- colnames(dt) } # Get all colnames from data.table
    ls_plt <- lapply(X = cols, FUN = function(i){
        if(verbose){
            cat(i, "\n")
        }
        tab <- rev(sort(table(dt[[i]])))
        if(length(tab) == 0){
            warning(paste0("Skipping empty column '", i, "'."))
            NULL
        } else {
            tab <- c(tab, nrow(dt)-sum(tab))
            names(tab)[which(names(tab) == "")] <- "NA"
            dt_count_var <- data.table::as.data.table(
                tab, keep.rownames = "Categories")
            dt_count_var <- data.table::melt.data.table(
                data = dt_count_var, id.vars = "Categories")
            dt_count_var[, Categories := as.factor(Categories)]
            dt_count_var[, Categories := factor(
                x = Categories, levels = rev(Categories))]
            if(nrow(dt_count_var) <= 70){
                pal_plt <- sample(
                    grDevices::colorRampPalette(pal)(nrow(dt_count_var)))
                if(dt_count_var[Categories == "NA"]$value != 0){
                    pal_plt[1] <- "grey"
                } else {
                    dt_count_var <- dt_count_var[value != 0]
                    dt_count_var[, Categories := droplevels(Categories)]
                }
                stackbar_plt <- ggplot2::ggplot(
                    data = dt_count_var,
                    mapping = ggplot2::aes(
                        x = variable, y = value, fill = Categories)) +
                    ggplot2::geom_bar(position = "stack", stat = "identity") +
                    ggplot2::scale_fill_manual(values = pal_plt) +
                    ggplot2::scale_y_continuous(breaks = seq(
                        from = 0, to = nrow(dt), length.out = y_breaks),
                        expand = c(0, 0)) +
                    ggplot2::scale_x_discrete(expand = c(0, 0)) +
                    ggplot2::labs(
                        x = i, y = "Cumulative amount of data available") +
                    ggplot2::theme(
                        axis.text.x = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_text(
                            size = 11, colour = "black"),
                        axis.title.y = ggplot2::element_text(size = 13),
                        legend.justification = c(0, 1),
                        plot.margin = ggplot2::margin(
                            0.3, 0.1, 0.1, 0.1, unit = "cm")) +
                    ggplot2::guides(
                        fill = ggplot2::guide_legend(nrow = lgd_nrow))
                stackbar_plt
            } else {
                warning(paste0(
                    "Column '", i,
                    "' skipped; too many categories for ggstackbar()."))
                NULL
            }
        }
    })
    names(ls_plt) <- cols
    if(rm.null){ # Remove NULL elements
        ls_plt <- ls_plt[!vapply(
            X = ls_plt, FUN = is.null, FUN.VALUE = logical(length = 1))]
    }
    return(ls_plt)
}
