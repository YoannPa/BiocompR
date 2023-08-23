
#' Computes density by column on a matrix.
#'
#' @param x A \code{matrix}.
#' @param na.rm  A \code{logical}. Should missing values (including NaN) be
#'               omitted from the calculations?
#' @param from   A \code{numeric} to specify the minimum value from which the
#'               kernel estimation starts
#'               (Default: from = min(x, na.rm = TRUE)).
#' @param to     A \code{numeric} to specify the maximum value to which the
#'               kernel estimation ends (Default: to = max(x, na.rm = TRUE)).
#' @param ncores An \code{integer} to specify the number of cores/threads to be
#'               used to parallel-compute densities on matrix columns.
#' @return A \code{list} of data.tables containing the density data, one
#'         data.table per column.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

colDensity <- function(
    x, na.rm = FALSE, from = min(x, na.rm = TRUE), to = max(x, na.rm = TRUE),
    ncores = 1){
    ls.density <- parallel::mclapply(
        X = seq(ncol(x)), mc.cores = ncores, FUN = function(i){
            dens.res <- stats::density(
                x = x[, i], from = from, to = to, na.rm = na.rm)
            dt.dens <- data.table::data.table(
                variable = colnames(x)[i], value = dens.res$x,
                dens = dens.res$y)
            dt.dens
        })
    return(ls.density)
}

#' Plots a density color map from a matrix or a molten data.frame.
#'
#' @param m        Can be a \code{matrix} of numerical values or a molten
#'                 \code{data.frame} (e.g. that you can generate with
#'                 \code{data.table::melt()}) with a minimum of 4 columns:
#'                 \itemize{
#'                  \item{1 column for rownames.}
#'                  \item{1 column for the variable.}
#'                  \item{1 column for values.}
#'                  \item{and at least 1 column specifying the groups to which
#'                  values belong to.}
#'                 }
#'                 The 'variable' and 'groups' columns can be as factors. If so,
#'                 the order of levels will influence the order of columns and
#'                 panels on the resulting plot.
#' @param rn.col   A \code{character} when m is a molten data.frame to specify
#'                 the name of the rownames column.
#' @param var.col  A \code{character} when m is a molten data.frame to specify
#'                 the name of the variable column.
#' @param val.col  A \code{character} when m is a molten data.frame to specify
#'                 the name of the column containing values.
#' @param grp.col  A \code{character} when m is a molten data.frame to specify
#'                 the name of the groups column.
#' @param col.map  A \code{character} vector of colors specifying a color map to
#'                 use for the density heatmap (Default: col.map = biopalette(
#'                 name = "viridis_C_plasma", mute = TRUE). It is advised to use
#'                 'viridis' palettes from \link{biopalette}, but any color
#'                 gradient is supported.
#' @param from     A \code{numeric} to specify the minimum value from which the
#'                 kernel estimation starts
#'                 (Default: from = min(m, na.rm = TRUE)).
#' @param to       A \code{numeric} to specify the maximum value to which the
#'                 kernel estimation ends (Default: to = max(m, na.rm = TRUE)).
#' @param sort.fun A \code{character} matching a function to apply to the matrix
#'                 or the molten dataframe, to sort columns of the density
#'                 heatmap based on the values calculated
#'                 (Default: sort.fun = NULL ; Supported: "base::mean",
#'                 "base::median", "base::sd", "stats::IQR", ...and more).
#' @param facet.by A \code{character} specifying whether the density map on a
#'                 molten data.frame should be faceted on rows
#'                 (Default: facet.by = "rows") or on columns
#'                 (facet.by = "cols").
#' @param ncores   An \code{integer} to specify the number of cores/threads to
#'                 be used to parallel-compute densities on matrix columns
#'                 (Default: ncores = 1).
#' @param verbose  A \code{logical} to specify whether the function process
#'                 should be verbose (verbose = TRUE), or mute
#'                 (Default: verbose = FALSE).
#' @return A \code{gg} plot of a density heatmap.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create a density heatmap on the iris dataset:
#' ggdensity_map(m = as.matrix(iris[, 1:4])) +
#'   theme(axis.text.x = element_text(size = 12))
#' # Sort density heatmap columns using a function (here mean() from base):
#' ggdensity_map(m = as.matrix(iris[, 1:4]), sort.fun = "base::mean") +
#'   theme(axis.text.x = element_text(size = 12))
#' # Facetting the density heatmap on the iris dataset by Species:
#' dt.iris <- as.data.table(iris)
#' molten.iris <- melt.data.table(data = dt.iris, id.vars = "Species")
#' molten.iris[, rn := .I]
#' ggdensity_map(m = molten.iris, grp.col = "Species") +
#'   theme(axis.text.x = element_text(size = 12))
#' # Plot densities on a specific range of values:
#' ggdensity_map(m = molten.iris, grp.col = "Species", from = 2, to = 5) +
#'   theme(axis.text.x = element_text(size = 12))
#' # Change the color map:
#' ggdensity_map(
#'   m = molten.iris, grp.col = "Species",
#'   col.map = biopalette(name = "viridis_D_viridis", mute = TRUE),
#'   from = 2, to = 5) +
#'   theme(axis.text.x = element_text(size = 12))
#' # Facetting the density heatmap by columns:
#' ggdensity_map(m = molten.iris, grp.col = "Species", facet.by = "cols") +
#'   theme(axis.text.x = element_text(size = 12, angle = -90))

ggdensity_map <- function(
    m, rn.col = "rn", var.col = "variable", val.col = "value",
    grp.col = "group", col.map = biopalette(
        name = "viridis_C_plasma", mute = TRUE), from = NULL, to = NULL,
    sort.fun = NULL, facet.by = "rows", ncores = 1, verbose = FALSE){
    #Fix BiocCheck() complaining about these objects initialization
    variable <- NULL
    group <- NULL
    sorting <- NULL
    value <- NULL
    dens <- NULL

    #Check sort.fun
    if(!is.null(sort.fun)){
        sort.fun <- check_fun(
            fun = sort.fun, param.name = "sort.fun", ncores = ncores)
    }
    #Compute density on matrix or molten dataframe
    if(is.matrix(m)){
        if(verbose){cat("Matrix detected\n")}
        data.type <- "matrix"
        if(is.null(from)){ from <- min(m, na.rm = TRUE) }
        if(is.null(to)){ to <- max(m, na.rm = TRUE) }
        if(verbose){cat("Computing densities...")}
        #Compute density
        ls.density <- colDensity(
            x = m, na.rm = TRUE, from = from, to = to, ncores = ncores)
        if(verbose){cat("Done.\n")}
        dt.density <- data.table::rbindlist(l = ls.density)
        #Sort using function if any set
        if(!is.null(sort.fun)){
            sort.fun.vals <- eval(parse(text = paste0(
                "apply(X = m, MARGIN = 2, FUN = function(i){", sort.fun,
                "(i, na.rm = TRUE)})")))
            sorted.names <- names(sort(sort.fun.vals))
            dt.density[, variable := as.factor(variable)]
            dt.density[, variable := factor(variable, levels = sorted.names)]
        }
    } else if(is.data.frame(m)){
        if(verbose){cat("data.frame detected\n")}
        data.type <- "molten"
        #Convert as a data.table if only data.frame
        if(!data.table::is.data.table(m)){
            m <- data.table::as.data.table(as.data.frame(m))}
        #Check if any other column name match "rn", "variable", "value" or "group"
        cols <- c(rn.col, var.col, val.col, grp.col)
        # if(verbose){cat("Discarding unused columns.\n")}
        # m <- m[, ..cols]
        other_cols <- colnames(m)[!colnames(m) %in% cols]
        if(length(other_cols[
            other_cols %in% c("rn", "variable", "value", "group")]) != 0){
            stop(paste(
                "Some additionnal columns are named 'rn', 'variable', 'value'",
                "or 'group'. Please give these columns another name."))
        }
        #Copy data.table before making modifications to it
        m <- data.table::copy(m)
        #Rename column to be used
        data.table::setnames(x = m, old = rn.col, new = "rn")
        data.table::setnames(x = m, old = var.col, new = "variable")
        data.table::setnames(x = m, old = val.col, new = "value")
        data.table::setnames(x = m, old = grp.col, new = "group")
        #Check which column is a factor
        if(is.factor(m$variable)){ variable_lvls <- levels(m$variable) }
        if(is.factor(m$group)){ group_lvls <- levels(m$group) }
        #Set from & to
        if(is.null(from)){ from <- min(m$value, na.rm = TRUE) }
        if(is.null(to)){ to <- max(m$value, na.rm = TRUE) }
        #Get group names
        grp.names <- unique(m$group)
        #Cast each matrix based on group and compute density
        ls.dt.dens <- lapply(X = grp.names, FUN = function(g){
            if(verbose){cat(paste0("\tCasting sub-matrix '", g, "'..."))}
            #Recast smaller matrix
            dt.subset <- data.table::dcast.data.table(
                data = m[group == g, c("rn", "variable", "value")],
                formula = ...~variable, value.var = "value")
            mat.subset <- as.matrix(dt.subset[, -c("rn"), ])
            rm(dt.subset)
            if(verbose){cat("Done.\n")}
            if(verbose){cat("\tComputing densities...")}
            ls.density <- colDensity(
                x = mat.subset, na.rm = TRUE, from = from, to = to,
                ncores = ncores)
            if(verbose){cat("Done.\n")}
            dt.density <- data.table::rbindlist(l = ls.density)
            rm(ls.density)
            dt.density
        })
        names(ls.dt.dens) <- grp.names
        dt.density <- data.table::rbindlist(l = ls.dt.dens, idcol = "group")
        #Remove list
        rm(ls.dt.dens)
        #Set back factors
        if(is.factor(m$variable)){
            dt.density[, variable := as.factor(variable)]
            dt.density[, variable := factor(
                x = variable, levels = variable_lvls)]
        }
        if(is.factor(m$group)){
            dt.density[, group := as.factor(group)]
            dt.density[, group := factor(x = group, levels = group_lvls)]
        }
        #Sort using function if any set
        if(!is.null(sort.fun)){
            eval(parse(text = paste0(
                "m[, sorting := ", sort.fun, "(value), by = variable]")))
            sorted.names <- as.character(unique(m, by = "sorting")[
                order(sorting)]$variable)
            if(!is.factor(dt.density$variable)){
                dt.density[, variable := as.factor(variable)]
            }
            dt.density[, variable := factor(variable, levels = sorted.names)]
        }
    }
    #Remove m
    rm(m)
    #Plot density map
    if(verbose){cat("Plotting density color map...")}
    ggdensmap <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = dt.density, mapping = ggplot2::aes(
            x = variable, y = value, fill = dens)) +
        # ggplot2::scale_fill_viridis_c(option = col.map) +
        ggplot2::scale_fill_gradientn(colors = col.map) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = 11),
            axis.title = ggplot2::element_text(size = 12),
            legend.position = "none",
            strip.text = ggplot2::element_text(size = 12),
            strip.background = ggplot2::element_rect(
                fill = "white", color = "black"),
            panel.spacing = ggplot2::unit(1, "lines"),
            panel.border = ggplot2::element_rect(fill = NA, color = "black"),
            plot.margin = ggplot2::margin(0.3, 0.1, 0.1, 0.1, unit = "cm")) +
        ggplot2::labs(x = var.col, y = val.col)
    if(data.type == "molten"){
        if(facet.by == "rows"){
            ggdensmap <- ggdensmap +
                ggplot2::facet_grid(rows = ggplot2::vars(group))
        } else if(facet.by == "cols"){
            ggdensmap <- ggdensmap +
                ggplot2::facet_grid(cols = ggplot2::vars(group))
        }
    }
    if(verbose){ cat("Done.\n") }
    rm(dt.density)
    #Return plot
    return(ggdensmap)
}
