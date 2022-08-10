
#' Draws a sunset plot showing the completeness of a dataset.
#'
#' @param mat              A \code{matrix}.
#' @param title            A \code{character} to specify the title of your plot.
#' @param col.pal          A \code{character} vector of length 3 matching R
#'                         colors to used as a palette for the plot.
#'                         Alternatively, col.pal can use the palettes provides
#'                         within the script of the function:
#'                         pal_sunset, pal_westworld, pal_startrek,
#'                         pal_margesimpson, pal_eighties, pal_insta\cr
#'                         (Default: col.pal = pal_sunset).
#' @param horizontal       A \code{logical} to specify whether the plot should
#'                         be drawn vertically or horizontally\cr
#'                         (Default: horizontal = FALSE).
#' @param reverse          A \code{logical} to specify whether the order
#'                         displayed of the bins should reversed\cr
#'                         (Default: reverse = FALSE).
#' @param keep_2nd_ticks   A \code{logical} to specify whether the ticks of the
#'                         secondary axis should be displayed or not\cr
#'                         (Default: keep_2nd_ticks = FALSE).
#' @param n.grad           An \code{integer} specifying the number of graduation
#'                         to use on the Y-Axis\cr(Default: n.grad = 15).
#' @param display.cutoff   A \code{double} specifying the minimum height of
#'                         a bin to label the number of samples on it. Below
#'                         this cutoff the number of samples is displayed in the
#'                         side margin\cr
#'                         (Default: display.cutoff = 0.03).
#' @param display.num.smpl A \code{double} specifying the minimum height of a
#'                         bin to display the number of samples in the side
#'                         margin\cr(Default: display.num.smpl = 0.01).
#'
#' @return A \code{gg} sunset plot.
#' @author Yoann Pageaud.
#' @importFrom data.table `:=` `.I`
#' @export

sunset <- function(
    mat, title = "Number of rows without missing data",
    col.pal = c("red", "gold", "blue4"), horizontal = FALSE, reverse = FALSE,
    n.grad = 15, display.cutoff = 0.03, display.num.smpl = 0.01,
    lgd.pos = "bottom"){
    #Fix BiocCheck() complaining about these objects initialization
    index <- NULL
    . <- NULL
    sample.amount <- NULL
    data.covered <- NULL
    sample.string <- NULL
    hiden <- NULL
    right.y.data.covered <- NULL
    right.y.hiden <- NULL
    right.y.diff.bins <- NULL
    diff.bins <- NULL
    white.line.hide <- NULL
    label.pos <- NULL
    cumulated <- NULL

    if(!is.matrix(mat)){ stop("'mat' must be a matrix.") }
    #Calculate number values not being NAs for each row in the matrix
    pos.cov <- as.integer(rowSums(!is.na(mat)))
    #Get total number of samples
    N <- dim(mat)[2]
    #Get total number of rows
    N.rows <- dim(mat)[1]
    #Tabulate data
    tab.data <- tabulate(pos.cov + 1L, N + 1L)

    #Compute size of bins
    dt.data <- data.table::data.table(
        "sample.amount" = c(0, seq(N)), "data.covered" = tab.data,
        "sample.string" = paste0(
            c("NONE", seq(N)), " (", round((tab.data/sum(tab.data))*100), "%",
            ")"))
    #Add index to dt.data
    dt.data[, index := .I]

    #Reverse option
    if(reverse){
        dt.data[, c("sample.amount", "data.covered", "sample.string") := .(
            rev(sample.amount), rev(data.covered), rev(sample.string)) ]
    }
    #Remove leading & lagging uncovered data
    dt.data <- dt.data[
        min(which(data.covered != 0)):max(which(data.covered != 0)), ]
    #Add label position, percentage, white_line position, difference bins &
    # cumulated sum
    dt.data[, c(
        "label.pos", "diff.bins", "cumulated", "right.y.data.covered",
        "right.y.sample.string", "right.y.label.pos","right.y.diff.bins") := .(
            cumsum(data.covered)-0.5*data.covered,
            data.covered/sum(data.covered),
            cumsum(data.covered), data.covered, sample.string,
            cumsum(data.covered)-0.5*data.covered,
            data.covered/sum(data.covered))]
    #Apply display cut-off on main labels
    dt.data[data.covered < display.cutoff * N.rows, hiden := TRUE]
    #Apply display cut-off on right y labels
    dt.data[
        right.y.data.covered > display.cutoff * N.rows, right.y.hiden := TRUE]
    dt.data[right.y.diff.bins < display.num.smpl, right.y.hiden := TRUE]
    #Set white lines display
    dt.data[diff.bins < display.num.smpl, white.line.hide := TRUE]
    #Set to FALSE all NAs in hiden, right.y.hiden & white.line.hide
    dt.data[is.na(hiden), hiden := FALSE]
    dt.data[is.na(right.y.hiden), right.y.hiden := FALSE]
    dt.data[is.na(white.line.hide), white.line.hide := FALSE]
    #Hide last white line if one
    if(utils::tail(dt.data$white.line.hide, n = 1) == FALSE){
        dt.data[nrow(dt.data), white.line.hide := TRUE]
    }
    #Keep white line of level before the one having a white line
    invisible(lapply(X = dt.data[, .I], FUN = function(i){
        if(!all(is.na(dt.data[i+1,]))){
            if((dt.data[i+1, ]$hiden == FALSE |
                dt.data[i+1, ]$right.y.hiden == FALSE) &
               dt.data[i, ]$white.line.hide != FALSE){
                dt.data[i, white.line.hide := FALSE]
            }
        }
    }))
    #Get positions of strings that will not be displayed
    pos_str <- data.table::data.table(
        dt.data[hiden == TRUE & right.y.hiden == TRUE]$index)
    #Get groups of following samples
    smpl_intervals <- do.call(
        rbind, by(pos_str, cumsum(c(0, diff(pos_str$V1) != 1)), function(g){
            data.table::data.table(start = min(g$V1), end = max(g$V1),
                                   width = diff(range(g$V1)) + 1)}))
    #Get groups average label.pos cumulative diff_bins
    if(nrow(smpl_intervals) != 0){
        invisible(lapply(smpl_intervals[, .I], function(i){
            dt_grp <- dt.data[index >= smpl_intervals[i, ]$start &
                                  index <= smpl_intervals[i, ]$end]
            data_cov_grp <- sum(dt_grp$data.covered)
            per_grp <- sum(dt_grp$diff.bins)
            str_grp <- paste0(
                dt_grp[index == smpl_intervals[i, ]$start]$sample.amount, " to ",
                dt_grp[index == smpl_intervals[i, ]$end]$sample.amount, " (",
                round(per_grp*100), "%)")
            if(dt_grp[1, ]$index == 1){
                pos_lab_grp <- utils::tail(dt_grp, n = 1)$cumulated/2
                amnt_grp <- utils::tail(dt_grp$sample.amount, n = 1)
            } else {
                pos_lab_grp <- (dt.data[dt_grp[1, ]$index - 1,]$cumulated +
                                    utils::tail(dt_grp, n = 1)$cumulated)/2
                amnt_grp <- utils::tail(
                    dt_grp$sample.amount, n = 1) - utils::head(
                        dt_grp$sample.amount, n = 1) + 1
            }
            if(data_cov_grp > display.cutoff * N.rows){
                hiden_grp <- FALSE
                right.y.hiden_grp <- TRUE
            } else {
                hiden_grp <- TRUE
                if(per_grp < display.num.smpl){ right.y.hiden_grp <- TRUE
                } else { right.y.hiden_grp <- FALSE }
            }
            new.dt <- data.table::data.table(
                "sample.amount" = amnt_grp, "data.covered" = data_cov_grp,
                "sample.string" = str_grp,
                "index" = utils::tail(dt_grp$index, n = 1),
                "label.pos" = pos_lab_grp, "diff.bins" = per_grp,
                "cumulated" = utils::tail(dt_grp$cumulated, n = 1),
                "right.y.data.covered" = data_cov_grp,
                "right.y.sample.string" = str_grp,
                "right.y.label.pos" = pos_lab_grp,
                "right.y.diff.bins" = per_grp, "hiden" = hiden_grp,
                "right.y.hiden" = right.y.hiden_grp)
            dt.data[index == new.dt$index, c(
                "sample.string", "label.pos", "right.y.sample.string",
                "right.y.label.pos", "hiden", "right.y.hiden") := .(
                    new.dt$sample.string, new.dt$label.pos,
                    new.dt$sample.string, new.dt$label.pos, new.dt$hiden,
                    new.dt$right.y.hiden)]
        }))
    }
    #"Sunset" Plot of the Amount of CpGs Covered by Number of Samples
    sun.plt <- ggplot2::ggplot() + ggplot2::theme_gray() +
        ggplot2::theme(legend.title.align = 0.5,
                       legend.text = ggplot2::element_text(size = 12),
                       legend.position = lgd.pos) +
        ggplot2::geom_bar(
            data = dt.data, mapping = ggplot2::aes(
                x = 0, y = data.covered, fill = sample.amount),
            stat = "identity") +
        ggplot2::geom_text(
            data = dt.data[hiden == FALSE], mapping = ggplot2::aes(
                x = 0, y = label.pos, label = sample.string), vjust = 0.35,
            hjust = 0.5, color = "white", size = 5) +
        ggplot2::scale_fill_gradient2(
            low = col.pal[1], mid = col.pal[2], high = col.pal[3],
            midpoint = round(N/2), breaks = round(seq(
                min(dt.data$sample.amount, na.rm = TRUE),
                max(dt.data$sample.amount, na.rm = TRUE), length.out = 4))) +
        ggplot2::scale_y_continuous(
            expand = c(0, 0),
            breaks = seq(0, N.rows, length.out = n.grad),
            labels = function(x) format(x, digits = 2, scientific = TRUE),
            sec.axis = ggplot2::sec_axis(
                trans = ~.,
                breaks = dt.data[right.y.hiden == FALSE]$right.y.label.pos,
                labels = dt.data[right.y.hiden == FALSE]$sample.string)) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::geom_hline(
            data = dt.data[white.line.hide == FALSE], mapping = ggplot2::aes(
                yintercept = cumulated), color = "white") +
        ggplot2::guides(
            fill = ggplot2::guide_colorbar(
                ticks.linewidth = 2, ticks.colour = "black"))
    if(horizontal) {
        sun.plt <- sun.plt +
            ggplot2::ggtitle(title) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(
                    title, size = 15, hjust = 0.5),
                axis.title = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
                legend.title = ggplot2::element_text(size = 14, vjust = 0.8),
                plot.margin = ggplot2::margin(0.1, 0, 0, 20, unit = "cm")) +
            ggplot2::labs(fill = "Number of Samples") +
            ggplot2::coord_flip()
    } else {
        sun.plt <- sun.plt +
            ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(size = 12, face = "bold"),
                axis.title = ggplot2::element_text(size = 15),
                legend.title = ggplot2::element_text(size = 14, vjust = 1),
                legend.title.align = 0,
                plot.margin = ggplot2::margin(0.4, 0.2, 0, 0.2, unit = "cm")) +
            ggplot2::ylab(title) +
            ggplot2::labs(fill = "Number\nof samples")
    }
    sun.plt
}
