#' Plots an histogram using ggplot2 from a numeric or character vector.
#'
#' @param x          A \code{numeric} or \code{character} vector to be used for
#'                   plotting the histogram.
#' @param xmin       A \code{numeric} specifying the minimum value of the
#'                   histogram distribution. All values below xmin will be
#'                   aggregated on the first histogram bin
#'                   (Default: xmin = min(x, na.rm = TRUE)).
#' @param xmax       A \code{numeric} specifying the maximum value of the
#'                   histogram distribution. All values above xmax will be
#'                   aggregated on the last histogram bin
#'                   (Default: xmax = max(x, na.rm = TRUE)).
#' @param binwidth   An \code{integer} specifying the width for the histogram
#'                   bins. If the binwidth is not a divisor of the entire
#'                   distribution, the closest divisor is used. By default, the
#'                   bins width is guessed internally.
#' @param nbins      An \code{integer} specifying the total number of bins to be
#'                   displayed on the plot. The resulting plot may not reflect
#'                   the number of bins specified, but should be close.
#'                   Moreover, if the range of the plot display exceed the range
#'                   of the histogram distribution, the number of visible bins
#'                   with a count > 0, will be lower, because the breaks of bins
#'                   are computed on the plot display range.
#' @param ngrad      An \code{integer} specifying the number of graduations to
#'                   display on the X-axis of the plot (Default: ngrad = 10).
#' @param round.grad An \code{integer} specifying the number significant digits
#'                   to be considered when calculating graduations
#'                   (Default: round.grad = 1).
#' @param bin.col    A \code{character} matching a R color code to be used to
#'                   fill the histogram bins (Default: bin.col = "#0570b0").
#' @param show.annot A \code{character} specifying what annotations should be
#'                   displayed on the histogram (Supoorted: show.annot = c(
#'                   "cutoff", "median", "all", "none");
#'                   Default: show.annot = "none").
#' @param label.size A \code{integer} specifying the size of annotations'
#'                   labels (Default: label.size = 4).
#' @param facet      A \code{numeric} or \code{character} vector of the same
#'                   length as x, optionally with n levels to be used to
#'                   generate n histograms in separate panels.
#' @param frow       An \code{integer} to specify the number of rows to display
#'                   panels for the facet_wrap() display.
#' @param fcol       An \code{integer} to specify the number of columns to
#'                   display panels for the facet_wrap() display.
#' @param verbose    A \code{logical} to verbose details about the histogram
#'                   computation. It may be useful to understand how the
#'                   resulting histogram was drawn.
#' @return A \code{gg} plot of an histogram.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create a numerical distribution
#' distrib <- rnorm(100, 5, 1)
#' # Draw the default gghist() histogram
#' gghist(x = distrib)
#' # Increase display range
#' gghist(x = distrib, xmin = 0, xmax = 10)
#' # Decrease distribution range and aggregate extremas
#' gghist(x = distrib, xmin = 4, xmax = 7)
#' # Set binwidth equals 2
#' gghist(x = distrib, xmin = 4, xmax = 7, binwidth = 2)
#' # Set nbins to 10 (drawing 9) bins
#' gghist(x = distrib, xmin = 4, xmax = 7, nbins = 10)
#' # If binwidth and nbins are both specified, nbins is ignored
#' gghist(x = distrib, binwidth = 1, nbins = 10)
#' # When dislay range exceed distribution range, the number of visible bins is
#' # inferior to the value of nbins specified
#' gghist(x = distrib, nbins = 50, xmin = 0, xmax = 10)
#' # Turn verbose on.
#' gghist(x = distrib, nbins = 10, verbose = TRUE)
#' # Display first valley cut-off and median of the distribution.
#' gghist(x = distrib, show.annot = "all")
#' # Change annotations' labels size
#' gghist(x = distrib, show.annot = "all", label.size = 3)
#' # Set another filling color to the histogram
#' gghist(x = distrib, bin.col = "orange")
#' # Facet the data into 4 panels A, B, C & D
#' gghist(x = distrib, facet = rep(x = c("A", "B", "C", "D"), each = 25))
#' # Display all facets onto 1 row
#' gghist(
#'     x = distrib, facet = rep(x = c("A", "B", "C", "D"), each = 25), frow = 1)
#' # Create a character distribution
#' distrib_char <- rep(x = LETTERS, round(runif(26, 1, 99)))
#' # Draw a default gghist() histogram on character data
#' gghist(x = distrib_char)
#' # Set nbins to 10 bins; The closest divisor of the distribution is used
#' # instead
#' gghist(x = distrib_char, nbins = 10)
#' # Set binwidth to 3; The closest divisor of the distribution is used instead
#' gghist(x = distrib_char, binwidth = 3)
#' # Turn verbose on to understand what happened to binwidth internally.
#' gghist(x = distrib_char, binwidth = 3, verbose = TRUE)
#' # Facet the character data into 2 panels 1 & 2
#' gghist(
#'     x = distrib_char, binwidth = 3,
#'     facet = rep(seq(2), length(distrib_char)/2))
#' # Display all facets onto 1 column
#' gghist(
#'     x = distrib_char, binwidth = 3,
#'     facet = rep(seq(2), length(distrib_char)/2), fcol = 1)

gghist <- function(
    x, xmin = min(x, na.rm = TRUE), xmax = max(x, na.rm = TRUE),
    binwidth = NULL, nbins = NULL, ngrads = 10, round.grad = 1,
    bin.col = "#0570b0", show.annot = "none", label.size = 4, facet = NULL,
    frow = NULL, fcol = NULL, verbose = FALSE){
    # Check if is any NAs
    if(any(is.na(x))){ x <- x[!is.na(x)] }
    if(is.numeric(x)){
        if(verbose){ cat("Numeric vector detected.\n") }
        # Set graduations
        xgrads <- round(
            x = seq(xmin, xmax, length.out = ngrads), digits = round.grad)
        # Set X axis labels
        xlabs <- as.character(xgrads)
        if(max(x, na.rm = TRUE) > xmax){
            # Set Maximum
            x[x > xmax] <- xmax
            xlabs[length(xlabs)] <- paste(xlabs[length(xlabs)], "\nor more")
        }
        if(min(x, na.rm = TRUE) < xmin){
            # Set Minimum
            x[x < xmin] <- xmin
            xlabs[1] <- paste(xlabs[1], "\nor less")
        }
        # Plot the histogram
        if(!is.null(facet)){
            dt.distrib <- data.table::data.table(x, facet = facet)
        } else { dt.distrib <- data.table::data.table(x) }
        if(!is.null(binwidth)){
            if(!is.null(nbins)){
                warning("'binwidth' and 'nbins' parameters are clashing. Ignoring 'nbins'.")
            }
            histo <- ggplot2::ggplot() +
                ggplot2::geom_histogram(
                    data = dt.distrib, mapping = ggplot2::aes(x = x),
                    boundary = 0, binwidth = binwidth, fill = bin.col,
                    alpha = 0.7)
        } else {
            if(!is.null(nbins)){
                histo <- ggplot2::ggplot() +
                    ggplot2::geom_histogram(
                        data = dt.distrib, mapping = ggplot2::aes(x = x),
                        boundary = 0, bins = nbins, fill = bin.col, alpha = 0.7)
            } else {
                histo <- ggplot2::ggplot() +
                    ggplot2::geom_histogram(
                        data = dt.distrib, mapping = ggplot2::aes(x = x),
                        boundary = 0, fill = bin.col, alpha = 0.7)
            }
        }
        if(!is.null(facet)){
            histo <- histo +
                ggplot2::facet_wrap(
                    facets = ggplot2::vars(facet), nrow = frow, ncol = fcol)
        }
    } else if(is.character(x)){
        if(verbose){ cat("Character vector detected.\n") }
        if(!is.null(facet)){
            dt.facet <- data.table::data.table(x = x, facet = facet)
            dt.facet[, N := .N, by = c("x", "facet")]
            dt.count <- unique(dt.facet)
            dt.count[, facet := as.factor(facet)]
            # Check if some facets have missing counts and complete with zeros
            ref.vect <- unique(dt.count$x)
            ls.dt <- split.data.frame(x = dt.count, f = dt.count$facet)
            ls.dt <- lapply(X = names(ls.dt), FUN = function(i){
                if(all(ref.vect %in% ls.dt[[i]]$x)){
                    ls.dt[[i]]
                } else {
                    ref.vect[!(ref.vect %in% ls.dt[[i]]$x)]
                    dt.zero <- data.table::data.table(
                        x = ref.vect[!(ref.vect %in% ls.dt[[i]]$x)],
                        facet = i, N = 0)
                    rbind(ls.dt[[i]], dt.zero)
                }
            })
            dt.count <- data.table::rbindlist(ls.dt)
        } else { dt.count <- data.table::data.table(table(x)) }
        if(!is.null(binwidth)){
            if(!binwidth%%1 == 0){
                # Update binwidth to the closest integer
                binwidth <- round(binwidth)
            }
            # Update nbins using binwidth
            if(!is.null(facet)){
                if(!is.null(nbins)){
                    warning(
                        "'binwidth' and 'nbins' parameters are clashing. Ignoring 'nbins'.")
                    nbins <- nrow(
                        dt.count[facet == levels(dt.count$facet)[1]])/binwidth
                } else {
                    if(verbose){ cat("Updating 'nbins' with 'binwidth'.\n") }
                    nbins <- nrow(
                        dt.count[facet == levels(dt.count$facet)[1]])/binwidth
                }
            } else {
                if(!is.null(nbins)){
                    warning(
                        "'binwidth' and 'nbins' parameters are clashing. Ignoring 'nbins'.")
                    nbins <- nrow(dt.count)/binwidth
                } else {
                    if(verbose){ cat("Updating 'nbins' with 'binwidth'.\n") }
                    nbins <- nrow(dt.count)/binwidth }
            }
        }
        if(!is.null(nbins)){
            # Check if nbins is an integer, if not round to closest integer
            if(!nbins%%1 == 0){ nbins <- round(nbins) }
            # nbins must be a divisor of the length of unique(x)
            if(!is.null(facet)){
                if(!(nrow(dt.count[
                    facet == levels(dt.count$facet)[1]])/nbins)%%1 == 0){
                    if(verbose){
                        cat("Cannot plot histogram with current 'nbins' value.",
                            "Using closest divisor as 'nbins'.\n")
                    }
                    ls_divisor <- lapply(
                        X = seq(nrow(dt.count[
                            facet == levels(dt.count$facet)[1]])-1),
                        FUN = function(i){
                            res <- nrow(dt.count[
                                facet == levels(dt.count$facet)[1]])/i
                            if(res%%1 == 0){ res } else { NA }
                        })
                    dt_div <- data.table::data.table(
                        "divisor" = unlist(ls_divisor),
                        "prox" = abs(nbins - unlist(ls_divisor)))
                    nbins <- dt_div[prox == min(prox, na.rm = TRUE), ]$divisor
                }
            } else {
                if(!(nrow(dt.count)/nbins)%%1 == 0){
                    if(verbose){
                        cat("Cannot plot histogram with current 'nbins' value.",
                            "Using closest divisor as 'nbins'.\n")
                    }
                    ls_divisor <- lapply(
                        X = seq(nrow(dt.count)-1), FUN = function(i){
                            res <- nrow(dt.count)/i
                            if(res%%1 == 0){ res } else { NA }
                        })
                    dt_div <- data.table::data.table(
                        "divisor" = unlist(ls_divisor),
                        "prox" = abs(nbins - unlist(ls_divisor)))
                    nbins <- dt_div[prox == min(prox, na.rm = TRUE), ]$divisor
                }
            }
            if(!is.null(facet)){
                step <- nrow(dt.count[facet == levels(dt.count$facet)[1]])/nbins
                dt.count <- dt.count[order(facet, x)]
                dt.count[, index := rep(seq(nbins), each = step), by = facet]
                dt.count[, N := sum(N), by = c("index", "facet")]
                dt.count[, x := paste(
                    x, collapse = " & "), by = c("index", "facet")]
            } else {
                step <- nrow(dt.count)/nbins
                dt.count[, index := rep(seq(nbins), each = step)]
                dt.count[, N := sum(N), by = "index"]
                dt.count[, x := paste(x, collapse = " & "), by = "index"]
            }
            dt.count <- unique(dt.count)
        }
        histo <- ggplot2::ggplot() +
            ggplot2::geom_bar(
                data = dt.count, mapping = ggplot2::aes(x = x, y = N),
                stat = "identity", fill = bin.col, alpha = 0.7,  width = 1)
        if(!is.null(facet)){
            histo <- histo +
                ggplot2::facet_wrap(
                    facets = ggplot2::vars(facet), nrow = frow, ncol = fcol)
        }
    } else { stop("format not supported for 'x'.") }

    #Get recommended cut-off value and median
    if(is.numeric(x)){
        if(show.annot != "none"){
            if(verbose){ cat("show.annot: On. Computing annotations...") }
            if(show.annot == "all"){
                show.median <- TRUE
                show.cutoff <- TRUE
            } else if(show.annot == "median"){
                show.median <- TRUE
                show.cutoff <- FALSE
            } else if(show.annot == "cutoff"){
                show.median <- FALSE
                show.cutoff <- TRUE
            } else { stop("Unsupported value for 'show.annot'.") }
            if(show.cutoff){
                cutoff.val <- quantmod::findValleys(
                    x = ggplot2::layer_data(histo)$count)[1]
            }
            if(show.median){ median.val <- stats::median(x, na.rm = TRUE) }
            if(verbose){ cat("Done.\n") }
        } else if(show.annot == "none"){
            show.median <- FALSE
            show.cutoff <- FALSE
        } else { stop("Unsupported value for 'show.annot'.") }
    } else if(is.character(x)){
        if(show.annot != "none"){
            warning("'show.annot' ignored for character data.")
        }
    }
    #Define limits of the plots for lower minimum and highter maximum
    if(is.numeric(x)){
        if(xmax > max(x, na.rm = TRUE) | xmin < min(x, na.rm = TRUE)){
            histo <- histo +
                ggplot2::scale_x_continuous(
                    expand = c(0, 0), breaks = xgrads, labels = xlabs,
                    limits = c(xmin, xmax), oob = scales::oob_keep)
        } else {
            histo <- histo +
                ggplot2::scale_x_continuous(
                    expand = c(0, 0), breaks = xgrads, labels = xlabs,
                    oob = scales::oob_keep)
        }
    } else if(is.character(x)){
        histo <- histo + ggplot2::scale_x_discrete(expand = c(0, 0))
    }
    #Default ggplot
    if(verbose){ cat("Creating histogram plot...") }
    histo <- histo +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(x = "Values", y = "Frequency") +
        ggplot2::theme(
            plot.margin = ggplot2::margin(0.1, 1, 0.1, 0.1, unit = "cm"),
            axis.title = ggplot2::element_text(size = 13),
            axis.text = ggplot2::element_text(size = 11),
            panel.background = ggplot2::element_rect(fill = "white"),
            panel.grid.major = ggplot2::element_line(colour = "grey"),
            panel.grid.minor = ggplot2::element_line(colour = "grey")
        )
    #Add annotations
    if(is.numeric(x)){
        #If recommended cut-off inferior or equal to median, plot cut-off
        if(show.annot == "all"){
            if(cutoff.val >= median.val){
                if(verbose){
                    cat("Turning 'show.cutoff' off: cut-off value >= median.\n")
                }
                show.cutoff <- FALSE
            }
        }
        if(show.median){
            histo <- histo +
                ggplot2::geom_vline(
                    xintercept = median.val, color = "#313695",
                    linewidth = 0.7) +
                ggrepel::geom_label_repel(
                    data = data.frame(), mapping = ggplot2::aes(
                        x = median.val, y = Inf, fontface = 2, label = paste0(
                            "median = ", round(x = median.val, digits = 2))),
                    vjust = 0.5, color = "#313695", size = label.size)
        }
        if(show.cutoff){
            histo <- histo +
                ggplot2::geom_vline(
                    xintercept = cutoff.val, color = "#d7191c",
                    linewidth = 0.7) +
                ggrepel::geom_label_repel(
                    data = data.frame(), mapping = ggplot2::aes(
                        x = cutoff.val, y = Inf, fontface = 2, label = paste0(
                            "cut-off = ", round(x = cutoff.val, digits = 2))),
                    vjust = 0.5, color = "#d7191c", size = label.size)
        }
    }
    if(verbose){ cat("Done.\n") }
    #Return plot
    return(histo)
}
