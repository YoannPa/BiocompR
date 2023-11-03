
#' Checks matching between annotation groups and annotation palettes.
#'
#' @param data       A \code{matrix} or a \code{data.frame} with column names.
#' @param annot.grps A \code{list} of vectors of groups to which variables
#'                   belongs for the annotation sidebars. Vectors' lengths have
#'                   to match the number of variables.
#' @param annot.pal  A \code{vector} or a list of vectors containing colors as
#'                   characters for the annotation sidebars. The length of
#'                   vectors has to match the number of levels of vectors listed
#'                   in 'annot.grps'. If a list is provided, its length must
#'                   match the length of the list provided to 'annot.grps'.
#' @param verbose    A \code{logical} to specify wether the function should be
#'                   run on verbose mode (verbose = TRUE) or not
#'                   (Default: verbose = FALSE).
#' @return An error message if something goes wrong during annotations checks.
#' @author Yoann Pageaud.
#' @keywords internal

check.annotations <- function(data, annot.grps, annot.pal, verbose = FALSE){
    #Groups checking
    if(is.list(annot.grps)){
        if(any(unlist(lapply(annot.grps, length))!= ncol(data))){
            stop("samples are not all assigned to a group.")
        } else{ if(verbose){
            #Print groups values
            invisible(lapply(seq_along(annot.grps), function(i){ cat(paste0(
                names(annot.grps)[i],": ", paste(
                    unique(annot.grps[[i]]), collapse = ", "), ".\n"))
            }))
        }
        }
    } else {
        stop(paste(
            "'annot.grps' must be a named list. e.g. annot.grps =",
            "list('annotation_name' = annotation_vector)", sep = "\n"))
    }
    #Color checking
    if(is.list(annot.pal)){
        if(length(annot.grps) == length(annot.pal)){
            invisible(lapply(seq_along(annot.pal), function(i){
                if(length(annot.pal[[i]]) != length(levels(as.factor(
                    annot.grps[[i]])))){
                    stop(paste0(
                        "The length of annotation '", names(annot.grps)[i],
                        "' levels do not match the length of the ",
                        "corresponding palette."))
                }
            }))
        } else {
            stop(paste(
                "The number of palettes does not match the number of",
                "annotations provided."))
        }
    } else if(!is.list(annot.pal)){ #if a single palette is provided
        invisible(lapply(seq_along(annot.grps), function(i){
            if(length(annot.pal) == 0){ stop("'annot.pal' is empty.") }
            if(length(levels(as.factor(annot.grps[[i]]))) != length(annot.pal)){
                if(isTRUE(all.equal(
                    target = annot.pal, current = grDevices::rainbow(
                        n = ncol(data))))){
                    stop(paste(
                        "A specific palette must be defined in",
                        "'annot.pal' to match the annotation provided."))
                } else {
                    stop(paste0(
                        "The length of annotation '",names(annot.grps)[i],
                        "' levels do not match the length of the ",
                        "corresponding palette."))
                }
            }
        }))
    } else { #If not a list or a vector
        stop(paste(
            "Unknown type for 'annot.pal'. 'annot.pal' should be either",
            "a list or a vector."))
    }
}

#' Checks if a list's attributes has for class 'element_blank'.
#'
#' @param arg A \code{list}.
#' @return A \code{logical}.
#' @author Yoann Pageaud.
#' @keywords internal

is.elt_blank <- function(arg){
    bool<-attributes(arg)$class[1] == "element_blank"
    return(bool)
}

#' Checks if n.breaks and labels match in a ScaleContinuous object.
#'
#' @param ScaleContinuous_obj A \code{ScaleContinuous} object usually generated
#'                            with \link[ggplot2]{scale_fill_gradient},
#'                            \link[ggplot2]{scale_fill_gradient2} or
#'                            \link[ggplot2]{scale_fill_gradientn}.
#' @author Yoann Pageaud.
#' @keywords internal

chk.breaks.labels <- function(ScaleContinuous_obj){
    if(!is.null(ScaleContinuous_obj$n.breaks)){
        if(ScaleContinuous_obj$n.breaks != length(ScaleContinuous_obj$labels)){
            stop("in 'upper_scale_fill_grad' n.breaks is not equal to the",
                 "length of labels.")
        }
    }
}

#' Extracts legend from a ggplot2 object.
#'
#' @param gg2.obj  A \code{gg} object with legends.
#' @return A \code{gg} object only containing the legends of the plot.
#' @author Yoann Pageaud.
#' @keywords internal
#' @references \href{https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs}{Share a legend between two ggplot2 graphs - Mara Averick}

get.lgd <- function(gg2.obj){
    tmp <- R.devices::suppressGraphics(ggplot2::ggplot_gtable(
        ggplot2::ggplot_build(gg2.obj)))
    leg <- which(vapply(X = tmp$grobs, FUN = function(x) x$name,
                        FUN.VALUE = character(length = 1)) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

#' Creates a dendogram in ggplot2.
#'
#' @param df          A \code{data.frame} containing variables and values to be
#'                    used to create the dendrogram.
#' @param orientation A \code{character} specifying the orientation of the
#'                    dendrogram. Possible values are "top" and "left".
#' @param reverse.x   A \code{logical} specifying whether the X-axis should be
#'                    reversed (reverse.x = TRUE like in a correlation plot), or
#'                    kept unchanged (reverse.x = TRUE like in a heatmap), for a
#'                    dendrogram to be displayed on the left of a plot
#'                    (Default: reverse.x = FALSE).
#' @param theme_dend  A ggplot2 \code{theme} to specify any theme parameter you
#'                    wish to custom the dendrogram
#'                    (Default: theme_dend = NULL). For more information about
#'                    how to define a theme, see \link[ggplot2]{theme}.
#' @return A \code{gg} plot of the dendrogram.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create matrix from mtcars dataset
#' mat <- as.matrix(t(scale(mtcars)))
#' # Compute Pearson's correlation between cars
#' corr.res <- psych::corr.test(
#'     mat, use = "pairwise", method = "pearson", adjust = "BH")$r
#' # Normalize the correlation matrix
#' norm.res <- (corr.res-min(corr.res))/(max(corr.res)-min(corr.res))
#' # Do a complete hierarchy clustering on the normalized matrix
#' hierarchy.clust <- fastcluster::hclust(
#'     d = stats::as.dist(1-norm.res), method = "complete")
#' # Convert the cluster object into a dendrogram object
#' ddgr <- stats::as.dendrogram(hierarchy.clust)
#' ddgr_dat <- ggdendro::dendro_data(ddgr)
#' # Draw the dendrogram vertically
#' ggdend(df = ddgr_dat$segments, orientation = "top")
#' # Same horizontally
#' ggdend(df = ddgr_dat$segments, orientation = "left")
#' # Reverse branches order
#' ggdend(df = ddgr_dat$segments, orientation = "left", reverse.x = TRUE)

ggdend <- function(df, orientation, reverse.x = FALSE, theme_dend = NULL){
    #Fix BiocCheck() complaining about these objects initialization
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    #Set default dendrogram theme
    theme_default_dend <- ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "transparent"),
        plot.background = ggplot2::element_rect(fill = "transparent"))
    #Set forced dendrogram theme
    theme_forced_dend <- ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.ticks.length = ggplot2::unit(0, "pt"))

    #Update theme_dend
    if(is.null(theme_dend)){
        theme_dend <- theme_default_dend + theme_forced_dend
    } else {
        theme_dend <- theme_default_dend + theme_dend + theme_forced_dend
    }
    #Create ggplot of dendrogram
    ddplot <- ggplot2::ggplot() + ggplot2::geom_segment(
        data = df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggplot2::expand_limits(x = c(0.5, max(df$x) + 0.5), y = 0) + theme_dend

    if(orientation == "top"){
        if(is.null(theme_dend$plot.margin)){
            dend_margin <- ggplot2::margin(0.1, 0, 0, 0, unit = "cm")
        } else {
            if(as.numeric(theme_dend$plot.margin[1]) != 0.1){
                dend_margin <- ggplot2::margin(
                    0.1, as.numeric(theme_dend$plot.margin[2]),
                    as.numeric(theme_dend$plot.margin[3]),
                    as.numeric(theme_dend$plot.margin[4]), unit = "cm")
            } else { dend_margin <- theme_dend$plot.margin }
        }
        ddplot <- ddplot +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::theme(plot.margin = dend_margin) +
            ggplot2::scale_y_continuous(expand = c(0, 0))
    } else if(orientation == "left"){
        if(is.null(theme_dend$plot.margin)){
            dend_margin <- ggplot2::margin(0, 0, 0, 0.1, unit = "cm")
        } else {
            if(as.numeric(theme_dend$plot.margin[4]) != 0.1){
                dend_margin <- ggplot2::margin(
                    as.numeric(theme_dend$plot.margin[1]),
                    as.numeric(theme_dend$plot.margin[2]),
                    as.numeric(theme_dend$plot.margin[3]), 0.1, unit = "cm")
            } else { dend_margin <- theme_dend$plot.margin }
        }
        ddplot <- ddplot +
            ggplot2::theme(plot.margin = dend_margin) +
            ggplot2::scale_y_reverse(expand = c(0, 0)) +
            ggplot2::coord_flip()
        if(reverse.x){
            ddplot <- ddplot + ggplot2::scale_x_reverse(expand = c(0, 0))
        } else {
            ddplot <- ddplot + ggplot2::scale_x_continuous(expand = c(0, 0)) }
    } else { stop("dendrogram's orientation value not supported by ggdend().") }
    return(ddplot)
}


#' Draws a triangle ggplot for a basic molten triangle matrix.
#'
#' @param melt_tri         A \code{data.frame} melted triangle containing a
#'                         statistical test values.
#' @param grid_col         A \code{character} specifying the color of the grid.
#' @param grid_linewidth   A \code{double} value for the thickness of the grid.
#' @param ggtri_theme      A ggplot2 \code{theme} to specify any theme parameter
#'                         you wish to custom on the triangle plot
#'                         (Default: ggtri_theme = NULL). For more information
#'                         about how to define a theme, see
#'                         \link[ggplot2]{theme}.
#' @param scale_fill_grad  A \code{ScaleContinous} object generated by ggplot2
#'                         functions such as \link[ggplot2]{scale_fill_gradient}
#'                         ,\link[ggplot2]{scale_fill_gradient2} or
#'                         \link[ggplot2]{scale_fill_gradientn} to customize
#'                         heatmap colors and the associated color bar.
#' @param guide_custom_bar A \code{guide} object generated by the ggplot2
#'                         function \link[ggplot2]{guide_colorbar} to custom the
#'                         triangle plot color bar appearance
#'                         (see also 'scale_fill_grad' option).
#' @param x_axis_pos       A \code{character} to specify the position of the X
#'                         axis on the plot (Default: x_axis_pos = "top";
#'                         Supported: c("top", "bottom")).
#' @param y_axis_pos       A \code{character} to specify the position of the Y
#'                         axis on the plot (Default: y_axis_pos = "right";
#'                         Supported: c("right", "left")).
#' @return A \code{gg} object of a basic triangle plot (a 'geom_tile()').
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create matrix from mtcars dataset
#' mat <- as.matrix(t(scale(mtcars)))
#' # Compute Pearson's correlation between cars
#' corr.res <- psych::corr.test(
#'     mat, use = "pairwise", method = "pearson", adjust = "BH")$r
#' # Order samples following a complete hierarchical clustering
#' correlation.order <- corrplot::corrMatOrder(
#'     corr = corr.res, order = "hclust", hclust.method = "complete")
#' corr.res <- corr.res[correlation.order, correlation.order]
#' # Remove duplicated data
#' corr.res[upper.tri(corr.res)] <- NA
#' # Melt Correlation matrix
#' dt.corr <- data.table::as.data.table(x = corr.res, keep.rownames = "Var1")
#' molt.corr <- data.table::melt.data.table(
#'     data = dt.corr, id.vars = "Var1", variable.name = "Var2", na.rm = TRUE,
#'     measure.vars = colnames(dt.corr)[-c(1)])
#' molt.corr[, Var1 := as.factor(x = Var1)]
#' molt.corr[, Var1 := factor(
#'     x = Var1, levels = as.character(unique(molt.corr$Var1)))]
#' # Invert order of samples
#' molt.corr[, Var2 := factor(x = Var2, levels = rev(levels(Var2)))]
#' # Replace identical correlations by NAs
#' molt.corr[Var1 == Var2, value := NA]
#' # Draw the triangle plot
#' ggtriangle(melt_tri = molt.corr)

ggtriangle <- function(
    melt_tri, grid_col = "white", grid_linewidth = 0.3, ggtri_theme = NULL,
    scale_fill_grad = ggplot2::scale_fill_gradientn(
        colors = BiocompR::biopalette(name = "viridis_C_plasma"),
        na.value = "grey"),
    guide_custom_bar = ggplot2::guide_colorbar(
        ticks.linewidth = 0.5, ticks.colour = "black", frame.colour = "black",
        frame.linewidth = 0.5), x_axis_pos = "top", y_axis_pos = "right"){
    # Fix BiocCheck() complaining about these objects initialization
    Var1 <- NULL
    Var2 <- NULL
    value <- NULL
    # Define default theme
    theme_default_ggtri <- ggplot2::theme(
        axis.text.x.top = ggplot2::element_text(
            color = "black", angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y.right = ggplot2::element_text(color = "black"),
        axis.title = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(size = 12),
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        plot.background = ggplot2::element_rect(fill = "transparent"))
    # Update theme
    ggtri_theme <- theme_default_ggtri + ggtri_theme
    # Create the triangle plot
    ggtri_plt <- ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = melt_tri, ggplot2::aes(x = Var1, y = Var2, fill = value),
            color = grid_col, linewidth = grid_linewidth) +
        ggtri_theme +
        ggplot2::scale_x_discrete(position = x_axis_pos, expand = c(0, 0)) +
        ggplot2::scale_y_discrete(position = y_axis_pos, expand = c(0, 0)) +
        scale_fill_grad +
        ggplot2::guides(fill = guide_custom_bar) +
        ggplot2::labs(x = NULL, y = NULL)
    return(ggtri_plt)
}


#' Maps annotations categories to palettes.
#'
#' @param origin.grps A \code{list} of characters or numerics specifying the
#'                    categories for legends, in their order of input.
#' @param groups      A \code{list} of characters or numerics specifying the
#'                    categories for legends, in their order of appearance on
#'                    the annotation sidebar of the plot.
#' @param annot.pal   A \code{list} of character vectors matching valid codes
#'                    for R colors.
#' @return A \code{list} of data.tables where categories of each annotation are
#'         matched with their color palette.
#' @author Yoann Pageaud.
#' @keywords internal

map.cat2pal <- function(origin.grps, groups, annot.pal){
    if(length(groups) == length(annot.pal)){ #if annotations match palettes
        #Map groups to palettes
        ls.dt.grp.pal <- Map(
            data.table::data.table, "Grps" = origin.grps, "Cols" = annot.pal,
            stringsAsFactors = FALSE)
        col_table <- lapply(seq_along(groups), function(i){
            if(length(groups[[i]]) == length(annot.pal[[i]])){
                #if number of groups match number of colors
                ls.dt.grp.pal[[i]][
                    match(groups[[i]], ls.dt.grp.pal[[i]]$Grps), ]
            } else {
                stop(paste0(
                    "The length of annotation '", names(groups)[i],
                    "' levels does not match the length of the corresponding ",
                    "palette."))
            }
        })
    } else {
        stop(paste(
            "The number of annotations does not match the number of",
            "palettes provided."))
    }
    return(col_table)
}


#' Builds color tables for legends.
#'
#' @param annot.grps A \code{list} of characters or numerics specifying the
#'                   categories for legends, in their order of input.
#' @param annot.pal  A \code{list} of character vectors, or a \code{vector} of
#'                   characters matching valid codes for R colors.
#' @return A \code{data.table} list of categories matching their colors for each
#'         legend.
#' @author Yoann Pageaud.
#' @keywords internal

build.col_table <- function(annot.grps, annot.pal){
    # Create list of groups using levels of annotation groups
    origin.grps <- lapply(X = annot.grps, FUN = levels)
    # origin.grps <- lapply(X = annot.grps, FUN = function(i){
    #     if(is.factor(i)){ levels(i) } else { levels(as.factor(i)) }
    # })

    # Update levels following the new order of the annotation
    groups <- origin.grps
    # groups <- lapply(X = lapply(X = annot.grps, FUN = function(i){
    #     factor(x = i, levels =  unique(i))}), FUN = levels)

    # Create list of color tables
    if(is.list(annot.pal)){
        # Map categories to palettes
        col_table <- BiocompR:::map.cat2pal(
            origin.grps = origin.grps, groups = groups, annot.pal = annot.pal)
    } else if(!is.list(annot.pal)){
        col_table <- lapply(X = origin.grps, FUN = function(grp){
            data.table::data.table(
                "Grps" = grp, "Cols" = annot.pal, stringsAsFactors = FALSE)
        })
    }
    return(col_table)
}


#' Gets legend lengths.
#'
#' @param col_table A \code{data.table} list of categories with their matching
#'                  color palettes for legends.
#' @return An \code{integer} list of sizes, one per legend to be displayed.
#' @author Yoann Pageaud.
#' @keywords internal

get.len.legends <- function(col_table){
    lgd_sizes <- unlist(lapply(X = col_table, FUN = nrow))
    lgd_sizes <- lgd_sizes + 1
    return(lgd_sizes)
}


#' Builds legends layout.
#'
#' @param col_table         A \code{data.table} list of categories with their
#'                          matching color palettes for legends. Each data.table
#'                          contains 2 columns:
#'                          \itemize{
#'                           \item{'Grps': for annotation categories.}
#'                           \item{'Cols': for the color palette matching
#'                           categories.}
#'                          }
#' @param height.lgds.space An \code{integer} specifying the height of legend
#'                          space on final plot.
#' @param ncol.override     An \code{integer} to override the layout build,
#'                          specifying the number of columns on which legends'
#'                          keys should be displayed.
#' @return A \code{matrix} to be used as a layout for legend grobs display on
#'         final plot.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Building the simplest layout with 1 legend and 2 categories
#' build_legends_layout(col_table = list(data.table::data.table(
#'   'Grps' = c('cat 1', 'cat 2'), 'Cols' = c('blue', 'red'))),
#'   height.lgds.space = 29)
#' # If more annotations are passed in the list height.lgds.space let the
#' # function stacks legends up to a total size of 29 units on 1 column, before
#' # generating a new column if there is not enough space.

build_legends_layout <- function(
    col_table, height.lgds.space, ncol.override = NULL){
    # Calculate legends length
    lgd_sizes <- BiocompR:::get.len.legends(col_table = col_table)
    #Compute the number of columns necessary for display of large legends
    if(is.null(ncol.override)){
        ncol.by.lgd <- ceiling(lgd_sizes/height.lgds.space)
    } else if(all.equal(ncol.override, as.integer(ncol.override)) == TRUE){
        ncol.by.lgd <- ncol.override
    } else { stop("Value not supported. 'ncol.override' must be an integer.") }

    # Calculate spatial disposition of legends based on their size and the
    # number of columns they occupy to create layout matrix by column
    legend_ids <- seq(length(col_table))
    columns <- lapply(X = legend_ids, FUN = function(i){
        if(i %in% legend_ids){
            lgd_position <- which(cumsum(lgd_sizes) <= height.lgds.space)
            if(length(lgd_position) == 0){
                lgd_position <- which(ceiling(cumsum(
                    lgd_sizes)/ncol.by.lgd[1]) <= height.lgds.space)
            }
            col_lgd <- legend_ids[lgd_position]
            col_mat <- rep(x = rep(x = col_lgd, times = ceiling(lgd_sizes[
                lgd_position]/ncol.by.lgd[1])), times = ncol.by.lgd[1])
            col_mat <- matrix(data = col_mat, ncol = ncol.by.lgd[1])
            #Make NA matrix to complete missing space
            NAmat <- matrix(
                nrow = height.lgds.space - nrow(col_mat), ncol = ncol.by.lgd[1])
            #Rbind legend matrix with NA matrix
            col_mat <- rbind(col_mat, NAmat)
            #Empty vectors from the legend previously used
            base::`<<-` (lgd_sizes, lgd_sizes[-lgd_position])
            base::`<<-` (legend_ids, legend_ids[-lgd_position])
            base::`<<-` (ncol.by.lgd, ncol.by.lgd[-lgd_position])
            #Return column matrix
            col_mat
        }
    })
    #Cbind all columns
    mlayout <- do.call(cbind, columns)
    #Return legends layout
    return(mlayout)
}


#' Draws a ggplot2 of a basic sidebar.
#'
#' @param data      A \code{data.table} with the column names 'Samples','.id'
#'                  and 'Groups' in this order.
#' @param palette   A \code{character} vector containing R colors like a
#'                  palette.
#' @param annot.sep A \code{double} specifying the width of the separation
#'                  spaces between annotations (Default: annot.sep = 0).
#' @param annot.cut A \code{double} specifying the width of cuts separating
#'                  annotation cells (Default: annot.cut = 0).
#' @param lgd.ncol  An \code{integer} specifying the number of columns to be
#'                  used to display a legend (Default: lgd.ncol = 1).
#' @return A \code{gg} object of the basic sidebar (a 'geom_tile()').
#' @author Yoann Pageaud.
#' @importFrom data.table `:=`
#' @export
#' @examples
#' # Prepare the dataset
#' dt_cars <- as.data.table(mtcars, keep.rownames = "Samples")
#' dt_cars <- dt_cars[, c("Samples", "am", "gear", "carb"), ]
#' dt_cars <- melt.data.table(
#'     data = dt_cars, id.vars = "Samples", variable.name = ".id",
#'     value.name = "Groups")
#' dt_cars[, Groups := paste(Groups, .id)]
#' # Create palette (length matches number of unique value in column 'Groups')
#' pal <- rainbow(n = 11)
#' # Draw the most basic annotation sidebar with 1 annotation only ('am')
#' ggsidebar.basic(data = dt_cars[.id == "am"], palette = pal)
#' # Draw an annotation sidebar with multiple annotations
#' ggsidebar.basic(data = dt_cars, palette = pal)
#' # Custom the annotation sidebar using ggplot2
#' ggsidebar.basic(data = dt_cars, palette = pal) +
#'     scale_y_discrete(expand = c(0,0)) +
#'     scale_x_discrete(expand = c(0,0)) +
#'     theme(
#'         axis.title.y = element_blank(),
#'         axis.text.x = element_text(angle = -45, hjust = 0, size = 10))
#' # Separate annotation with a blank space
#' ggsidebar.basic(data = dt_cars, palette = pal, annot.sep = 0.1)
#' # Separate samples within annotations by a blank space
#' ggsidebar.basic(data = dt_cars, palette = pal, annot.cut = 0.3)
#' # Display legend keys on 2 columns instead of 1
#' ggsidebar.basic(data = dt_cars, palette = pal, lgd.ncol = 2)
#' # Specify 1 of the annotation as the facet of the sidebar
#' ggsidebar.basic(data = dt_cars, palette = pal, facet = "am")

ggsidebar.basic <- function(
    data, palette, annot.sep = 0, annot.cut = 0, lgd.ncol = 1, facet = NULL){
    #Fix BiocCheck() complaining about these objects initialization
    .id <- NULL
    facet.annot <- NULL
    Groups <- NULL
    Samples <- NULL
    # Set palette if facet is on
    if(!is.null(facet)){
        if(is.factor(data$Groups)){
            dt_levcol <- data.table::data.table(
                "levels" = levels(data$Groups), "colors" = palette)
            sublev <- data[.id == facet]$Groups
            palette <- dt_levcol[!levels %in% levels(droplevels(sublev))]$colors
        }
    }
    # Create a basic sidebar plot
    basic <- ggplot2::ggplot() +
        ggplot2::theme(
            legend.justification = c(1, 1),
            legend.spacing.y = ggplot2::unit(0.05, 'cm'),
            legend.margin = ggplot2::margin(0, 0.8, 0, 0, unit = "cm"),
            legend.box.margin = ggplot2::margin(3, 0, 0, 0, unit = "cm"),
            axis.text = ggplot2::element_text(size = 12),
            panel.grid = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 0, 0)) +
        ggplot2::scale_fill_manual(values = as.character(palette)) +
        ggplot2::guides(fill = ggplot2::guide_legend(
            ncol = lgd.ncol, byrow = TRUE))
    if(!is.null(facet)){
        # Add faceting
        dt.facet <- data[.id == facet]
        dt.facet[, `:=`(facet.annot, Groups)]
        data <- merge(
            x = data, y = dt.facet[, c("Samples", "facet.annot"), ],
            by = "Samples", all.x = TRUE)
        data <- data[.id != facet]
        # Plot annotation bar with facet
        basic <- basic + ggplot2::geom_tile(data = data, mapping = ggplot2::aes(
            x = Samples, y = .id, fill = Groups, height = 1 - annot.sep,
            width = 1 - annot.cut)) +
            ggplot2::facet_grid(
                . ~ facet.annot, scales = "free", space = "free") +
            ggplot2::theme(
                panel.spacing = ggplot2::unit(0, "lines"),
                strip.background = ggplot2::element_rect(
                    color = "black", linewidth = 0.5))
    } else {
        basic <- basic + ggplot2::geom_tile(data = data, mapping = ggplot2::aes(
            x = Samples, y = .id, fill = Groups, height = 1 - annot.sep,
            width = 1 - annot.cut))
    }
    return(basic)
}


#' Creates a colored side annotation bars in ggplot2.
#'
#' @param sample.names A \code{character} vector of the labels to be used for
#'                     annotation.
#' @param annot.grps   A \code{list} of vectors of groups to which variables
#'                     belongs for the annotation sidebars. vectors lengths
#'                     have to match the number of variables.
#' @param annot.pal    A \code{vector} or a list of vectors containing colors
#'                     as characters for the annotation sidebars. The length of
#'                     vectors has to match the number of levels of vectors
#'                     listed in 'annot.grps'.
#'                     \itemize{
#'                      \item{If annot.pal is a list: its length must match the
#'                            length of the list provided to 'annot.grps'.}
#'                      \item{If annot.pal is a vector: make sure that the
#'                            levels content of annotations listed in
#'                            'annot.grps' is the same, and that no annotation
#'                            contains less or more levels than another one in
#'                            'annot.grps'.}
#'                     }
#' @param annot.pos    A \code{character} specifying the position of the
#'                     annotation sidebar.\cr Possible values are: 'top',
#'                     'left' or 'both'.
#' @param annot.sep    A \code{double} specifying the width of the separation
#'                     spaces between annotations (Default: annot.sep = 0).
#' @param annot.cut    A \code{double} specifying the width of cuts separating
#'                     annotation cells (Default: annot.cut = 0).
#' @param merge.lgd    A \code{logical} to specify whether annotation legends
#'                     should be merged (annot.lgd.merge = TRUE) or remain
#'                     separated (annot.lgd.merge = FALSE)
#'                     (Default: annot.lgd.merge = FALSE).
#' @param right        A \code{logical} to specify that Y-axis should be
#'                     displayed on the right side of the plot
#'                     (Default: right = FALSE).
#' @param lgd.name      A \code{character} to specify a title to the legend of
#'                     the plot, only if 'merge.lgd' = TRUE
#'                     (Default: lgd.name = "Legends").
#' @param lgd.ncol     An \code{integer} specifying the number of columns to be
#'                     used to display a legend (Default: lgd.ncol = 1).
#' @param theme_legend A ggplot2 \code{theme} to specify any theme parameter you
#'                     wish to custom on legends (Default: theme_legend = NULL).
#'                     For more information about how to define a theme, see
#'                     \link[ggplot2]{theme}.
#' @param theme_annot  A ggplot2 \code{theme} to specify any theme parameter you
#'                     wish to custom on the annotation bar
#'                     (Default: theme_annot = NULL). For more information about
#'                     how to define a theme, see \link[ggplot2]{theme}.
#' @param set.x.title  A \code{character}to be used as the title for the X axis.
#' @param set.y.title  A \code{character}to be used as the title for the Y axis.
#' @param dendro.pos   A \code{character} specifying the position of the
#'                     dendrogram (Supported: dendro.pos = c('top', 'left')).
#' @return A \code{list} of length 2:
#'         \itemize{
#'          \item{'sidebar' contains the colored sidebar plot.}
#'          \item{'legends' lists legends of the matching sidebar as \code{gg}
#'                objects.}
#'         }
#' @author Yoann Pageaud.
#' @export ggsidebar.full
#' @export
#' @examples
#' # Draw a default sidebar
#' annotbar <- ggsidebar.full(
#'     sample.names = rownames(mtcars),
#'     annot.grps = list("Carb" = paste(mtcars$carb, "carb")),
#'     annot.pal = rainbow(n = 6))
#' # Draw the sidebar and customize it with ggplot2
#' annotbar$sidebar +
#'     theme(
#'         axis.text.x = element_text(angle = 45, hjust = 0),
#'         plot.margin = margin(0, 3, 0, 0, "cm"))
#' # Draw the sidebar's first legend
#' grid::grid.draw(annotbar$legends[[1]])

ggsidebar.full <- function(
    sample.names, annot.grps, annot.pal, annot.pos = "top", annot.sep = 0,
    annot.cut = 0, merge.lgd = FALSE, right = FALSE, lgd.name = "Legends",
    lgd.ncol = 1, theme_legend = NULL,
    theme_annot = ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12)),
    set.x.title = NULL, set.y.title = NULL, dendro.pos = "top", facet = NULL){
    # Fix BiocCheck() complaining about these objects initialization
    Grps <- NULL
    Cols <- NULL
    .id <- NULL
    # Create list of groups in their original order
    origin.grps <- lapply(X = annot.grps, FUN = function(i){
        if(is.factor(i)){ levels(i) } else { levels(as.factor(i)) }
    })
    # Update levels following the new order of the annotation
    groups <- origin.grps
    # Create list of color tables
    if(is.list(annot.pal)){ # If a list of palettes is provided
        # Map categories to palettes
        col_table <- BiocompR:::map.cat2pal(origin.grps, groups, annot.pal)
    } else if(!is.list(annot.pal)){ #if a single palette is provided
        # Map groups to the same palette
        ls.df.grp.pal <- lapply(X = origin.grps, FUN = function(grp){
            data.table::data.table(
                "Grps" = grp, "Cols" = annot.pal, stringsAsFactors = FALSE)
        })
        col_table <- lapply(seq_along(groups), function(i){
            if(length(groups[[i]]) == length(annot.pal)){
                # If groups match colors
                ls.df.grp.pal[[i]][
                    match(groups[[i]], ls.df.grp.pal[[i]]$Grps), ]
            } else {
                stop(paste0("The length of annotation '", names(groups)[i],
                            "' levels do not match the length of the ",
                            "corresponding palette."))}
        })
    } else { # If not a list or a vector
        stop(paste(
            "Unknown type for 'annot.pal'. 'annot.pal' should be either",
            "a list or a vector."))
    }
    # Create list of annotation data.tables
    dframe.annot <- lapply(annot.grps, function(i){
        data.table::data.table("Samples" = sample.names, "Groups" = i)
    })
    # Order samples following the correlation order provided
    dframe.annot <- lapply(dframe.annot, function(i){
        dt_annot <- i
        sample_order <- dt_annot$Samples
        dt_annot[, Samples := as.factor(Samples)]
        dt_annot[, Samples := factor(x = Samples, levels = sample_order)]
        dt_annot
    })
    # Rbind all annotations
    dframe.annot <- data.table::rbindlist(dframe.annot, idcol = TRUE)
    unique_id <- unique(dframe.annot$.id)
    # Convert .id as factors
    dframe.annot[, .id :=  as.factor(.id)]
    dframe.annot[, .id :=  factor(x = .id, levels = unique_id)]

    if(annot.pos == "top"){ # Change order of levels
        dframe.annot[, .id :=  factor(x = .id, levels = rev(unique_id))]
    }

    # Check color tables
    col_table <- lapply(X = col_table, FUN = function(tbl){
        if(any(duplicated(tbl$Cols))){ # Check palette consistency
            stop(paste(
                "1 color in a palette has been associated to more than 1",
                "group."))
        }
        if(any(duplicated(tbl$Grps))){ # Check annotation consistency
            warning("Duplicated group name provided. Removing duplicated...")
            tbl <- tbl[!duplicated(Grps)]
            tbl
        } else { tbl }
    })
    col_table <- data.table::rbindlist(col_table, idcol = TRUE)
    if(!is.list(annot.pal)){
        if(any(duplicated(col_table$Cols))){ # Check palette consistency
            col_table <- col_table[!duplicated(x = Cols)]
        }
    }
    # Plot color sidebars
    col_sidebar <- BiocompR::ggsidebar.basic(
        data = data.table::copy(dframe.annot), palette = col_table$Cols,
        annot.sep = annot.sep, annot.cut = annot.cut, facet = facet,
        lgd.ncol = lgd.ncol[1])
    # Add legend theme parameters if some
    col_sidebar <- col_sidebar + theme_legend
    # Modify base plot following its position
    if(annot.pos == "top"){
        theme_annot <- theme_annot + ggplot2::theme(
            axis.text.x.top = theme_annot$axis.text.x,
            axis.text.x = ggplot2::element_blank())
        col_sidebar <- col_sidebar + theme_annot +
            ggplot2::scale_x_discrete(expand = c(0, 0), position = "top") +
            ggplot2::xlab(set.x.title)
        if(right){
            col_sidebar <- col_sidebar +
                ggplot2::scale_y_discrete(position = 'right', expand = c(0, 0))
        } else {
            col_sidebar <- col_sidebar +
                ggplot2::scale_y_discrete(expand = c(0, 0))
        }
        if(dendro.pos != "top"){
            theme_annot <- theme_annot +
                ggplot2::theme(axis.title.y = ggplot2::element_blank())
            col_sidebar <- col_sidebar + theme_annot
        } else {
            theme_annot <- theme_annot +
                ggplot2::theme(axis.title = ggplot2::element_blank())
            col_sidebar <- col_sidebar + theme_annot
        }
    } else if(annot.pos == "left"){
        theme_annot <- theme_annot + ggplot2::theme(
            axis.text.x.top = theme_annot$axis.text.x,
            axis.text.x = ggplot2::element_blank())
        col_sidebar <- col_sidebar + ggplot2::coord_flip() + theme_annot +
            ggplot2::scale_x_discrete(expand = c(0, 0)) +
            ggplot2::scale_y_discrete(expand = c(0, 0), position = "right") +
            ggplot2::xlab(set.y.title)
        if(dendro.pos != "left"){
            col_sidebar <- col_sidebar + theme_annot +
                ggplot2::theme(axis.title.x = ggplot2::element_blank())
        } else {
            col_sidebar <- col_sidebar + theme_annot +
                ggplot2::theme(axis.title = ggplot2::element_blank())
        }
    }
    if(merge.lgd){ # Do not split legends
        sidebar.lgd <- list(
            BiocompR:::get.lgd(col_sidebar + ggplot2::labs(fill = lgd.name)))
    } else { # Split legends and return a list of legends
        if(annot.pos == "top"){
            if(!is.null(facet)){
                dframe.annot <- dframe.annot[.id != facet]
                unique_id <- unique_id[unique_id != facet]
                dframe.annot[, Groups := droplevels(Groups)]
                col_table <- col_table[Grps %in% levels(dframe.annot$Groups)]
                col_table[, .id := as.factor(.id)]
                data.table::setattr(
                    x = col_table$.id, name = "levels",
                    value = as.character(seq(length(levels(col_table$.id)))))
                col_table[, .id := as.integer(as.character(.id))]
            }
            dframe.annot[, .id :=  factor(x = .id, levels = unique_id)]
        }
        # Generate separate legends if more than 1 palette available or if only
        # 1 annotation is used
        if((is.list(annot.pal) & length(annot.pal) > 1) | length(
            levels(dframe.annot$.id)) == 1){
            # Get all legends separately
            sidebar.lgd <- lapply(
                X = seq_along(levels(dframe.annot$.id)), FUN = function(i){
                    BiocompR:::get.lgd(BiocompR::ggsidebar.basic(
                        data = dframe.annot[.id == levels(dframe.annot$.id)[i]],
                        palette = col_table[.id == i]$Cols,
                        lgd.ncol = lgd.ncol[i]) + theme_legend + ggplot2::labs(
                            fill = levels(dframe.annot$.id)[i]))
                })
        } else {
            stop(paste(
                "Cannot generate separated legends if only one",
                "annotation palette is given."))
        }
    }
    return(list("sidebar" = col_sidebar +
                    ggplot2::theme(legend.position = "none"),
                "legends" = sidebar.lgd))
}


#' Resizes heights or widths of multiple grobs based on a given grob dimensions.
#'
#' @param ls.grobs   A \code{grob} list. The list can be named if necessary.
#' @param dimensions A \code{character} specifying the type of dimensions to
#'                   resize, either 'heights' or 'widths'.
#' @param start.unit An \code{integer} specifying at which rank of the unit
#'                   object the dimension comparison between grobs should start.
#' @param end.unit   An \code{integer} specifying at which rank of the unit
#'                   object the dimension comparison between grobs should end.
#' @return A \code{grob} list, all resized with their dimensions modified by the
#'         grid::unit.pmax() function.
#' @author Yoann Pageaud.
#' @export
#' @references \href{https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page}{Align two plots on a page - Mara Averick}

resize.grobs <- function(ls.grobs, dimensions, start.unit, end.unit){
    #Get dimension units from the list of grobs to redimension
    ls.dim <- lapply(X = ls.grobs, FUN = function(i){
        if(length(i[[dimensions]]) < end.unit){
            i[[dimensions]][start.unit:length(i[[dimensions]])]
        } else { i[[dimensions]][start.unit:end.unit] }
    })
    #Calculate maximum of all unit objects including the main grob.
    max.dim <- eval(parse(
        text = paste("grid::unit.pmax(",paste(paste(
            rep("ls.dim[[",length(ls.dim)), seq(length(ls.dim)),"]]", sep = ""),
            collapse = ", "), ")", sep = "")))
    #Apply changes to grobs dimensions
    ls.grobs <- lapply(X = ls.grobs, FUN = function(i){
        i[[dimensions]][start.unit:end.unit] <- as.list(max.dim)
        i
    })
    return(ls.grobs)
}


#' Resizes heights or widths of a grob based on the dimensions of another grob.
#'
#' @param grob1      A \code{grob} to be modified.
#' @param grob2      A \code{grob} to be used as reference for dimensions
#'                   modifications.
#' @param dimensions A \code{character} specifying the type of dimensions to
#'                   resize, either 'heights' or 'widths'.
#' @param positions  An \code{integer} vector specifying indexes of the
#'                   dimensions to change.
#' @return A \code{grob} for which the dimensions have been modified.
#' @author Yoann Pageaud.
#' @export

resize.grob.oneway <- function(grob1, grob2, dimensions, positions){
    #Get dimension units from the list of grobs to redimension
    if(max(positions) > length(grob2[[dimensions]]) | min(positions) < 1){
        stop("positions out of range in grob2.")
    }
    if(max(positions) > length(grob1[[dimensions]])){
        warning("some positions out of range in grob1 ignored.")
        positions <- positions[positions <= length(grob1[[dimensions]])]
    }
    grob1[[dimensions]][positions] <- grob2[[dimensions]][positions]
    return(grob1)
}


#' Rasterize a gg plot into a raster grob.
#'
#' @param gg.plot A \code{gg} plot to be rasterized.
#' @param filter  A \code{character} to be used as a filter for ggplot
#'                rasterization. The list of the supported rasterization
#'                filters is available in magick::filter_types()
#'                (Default: raster = "Lanczos"). Warning: Be aware that
#'                rasterization may take several minutes to process the ggplot.
#' @param size    An \code{integer} specifying the size of a squared raster in
#'                pixels. If size = 1080 -> raster = 1080x1080
#'                (Default: size = 1080).
#' @return A \code{grob} of the rasterized ggplot.
#' @author Yoann Pageaud.
#' @export

raster.gg2grob <- function(
    gg.plot, filter = "Lanczos", size = 1080){
    #Catch heatmap in magick::image_graph()
    fig <- magick::image_graph(width = 2160, height = 2160, res = 96)
    print(gg.plot)
    grDevices::dev.off()
    rastered <- magick::image_resize(
        image = fig, geometry = paste0(size, "x", size), filter = filter)
    #Create raster grob
    raster.grob <- grid::rasterGrob(
        rastered, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
        interpolate = TRUE)
    #Return grob annotation
    return(raster.grob)
}

#' This is a special geom intended for use as static annotations derived from
#' ggplot2::annotation_custom() matching a specific panel on a faceted ggplot.
#'
#' @param grob      A \code{grob} to display.
#' @param xmin,xmax x location (in data coordinates) giving horizontal location
#'                  of raster.
#' @param ymin,ymax y location (in data coordinates) giving vertical location
#'                  of raster.
#' @param data      A subset of a \code{data.table} matching the panel where the
#'                  grob annotation should be displayed.
#' @return A ggplot2 \code{LayerInstance} object.
#' @keywords internal

annotation_custom2 <- function(
    grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
    ggplot2::layer(
        data = data, stat = ggplot2::StatIdentity,
        position = ggplot2::PositionIdentity, geom = ggplot2::GeomCustomAnn,
        inherit.aes = TRUE, params = list(
            grob = grob, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
}


#' Update a guide_colorbar with values from a new one, with exceptions.
#'
#' @param old_guide    A \code{guide} of class 'colorbar' you want to update.
#' @param new_guide    A \code{guide} of class 'colorbar' updating the old_guide
#'                     with its values.
#' @param forced_param A named \code{list} containing guide_colorbar parameters
#'                     associated with their values to be conserved.
#'                     Independently from whether the there are differences or
#'                     not between old and new guides, values passed using
#'                     'forced_param' will override old and new guide values.
#' @return An updated \code{guide} of class 'colorbar'.
#' @author Yoann Pageaud.
#' @keywords internal

update_guide_colorbar <- function(
    old_guide = ggplot2::guide_colorbar(), new_guide, forced_param = NULL){
    new_guide <- lapply(X = names(old_guide), function(i){
        #If any difference for the parameter between previous and new guides
        if(isTRUE(all.equal(
            target = new_guide[[i]], current = old_guide[[i]]))){
            if(!is.null(forced_param)){
                if(i %in% names(forced_param)){
                    #If different, but also in forced param, keep forced value
                    forced_param[[i]]
                } else {
                    #if different, but not in forced param, update value
                    old_guide[[i]]
                }
            } else { old_guide[[i]] } #No difference: Keep the reference value
        } else {
            if(!is.null(forced_param)){
                if(i %in% names(forced_param)){
                    #If different, but also in forced param, keep forced value
                    forced_param[[i]]
                } else {
                    #if different, but not in forced param, update value
                    new_guide[[i]]
                }
            } else { new_guide[[i]] }
        }
    })
    names(new_guide) <- names(ggplot2::guide_colorbar())
    class(new_guide) <- c("guide", "colorbar")
    return(new_guide)
}

#' Builds the gg2heatmap layout using the egg package
#'
#' @return A \code{egg} object of the assembled heatmap.
#' @author Yoann Pageaud.
#' @keywords internal

heatmap_layout <- function(
    dd.rows = FALSE, dd.cols = FALSE, show.annot = FALSE, ddgr_seg_row = NULL,
    ddgr_seg_col = NULL, sidebar, htmp, legends, dend.row.size = 0,
    dend.col.size = 0, annot.size = 0, lgd.layout,
    lgd.space.width, override_grob_list = NULL, override_grob_height = NULL){
    # Create empty ggplot
    ggempty <- ggplot2::ggplot(data.frame()) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "cm"),
            panel.background = ggplot2::element_rect(fill = NA))
    # Create Heatmap layout
    if(dd.rows & dd.cols & show.annot){
        main_grob <- R.devices::suppressGraphics(egg::ggarrange(
            ggempty, ggempty, ddgr_seg_row, ddgr_seg_col, sidebar, htmp,
            heights = c(dend.col.size + 2, annot.size, 30),
            widths = c(dend.row.size + 1, 10), byrow = FALSE, draw = FALSE))
    } else if(dd.rows & dd.cols & !show.annot){
        main_grob <- R.devices::suppressGraphics(egg::ggarrange(
            ggempty, ddgr_seg_row, ddgr_seg_col, htmp,
            heights = c(dend.col.size + 2, 30),
            widths = c(dend.row.size + 1, 10), byrow = FALSE, draw = FALSE))
    } else if(!dd.rows & !dd.cols & show.annot){
        main_grob <- R.devices::suppressGraphics(egg::ggarrange(
            sidebar, htmp, heights = c(annot.size, 30), byrow = FALSE,
            draw = FALSE))
    } else if(!dd.rows & !dd.cols & !show.annot){
        main_grob <- R.devices::suppressGraphics(egg::ggarrange(
            htmp, byrow = FALSE, draw = FALSE))
    } else if(dd.rows & !dd.cols & show.annot){
        main_grob <- R.devices::suppressGraphics(egg::ggarrange(
            ggempty, ddgr_seg_row, sidebar, htmp,
            heights = c(annot.size, 30), widths = c(dend.row.size + 1, 10),
            byrow = FALSE, draw = FALSE))
    } else if(dd.rows & !dd.cols & !show.annot){
        main_grob <- R.devices::suppressGraphics(egg::ggarrange(
            ddgr_seg_row, htmp, widths = c(dend.row.size + 1, 10),
            byrow = FALSE, draw = FALSE))
    } else if(!dd.rows & dd.cols & show.annot){
        main_grob <- R.devices::suppressGraphics(egg::ggarrange(
            ddgr_seg_col, sidebar, htmp,
            heights = c(dend.col.size + 2, annot.size, 30), byrow = FALSE,
            draw = FALSE))
    } else if(!dd.rows & dd.cols & !show.annot){
        main_grob <- R.devices::suppressGraphics(egg::ggarrange(
            ddgr_seg_col, sidebar, htmp,
            heights = c(dend.col.size + 2, annot.size, 30), byrow = FALSE,
            draw = FALSE))
    }
    # Add annotation legends or not
    if(show.annot){
        if(anyNA(lgd.layout)){
            lgd.layout[is.na(lgd.layout)] <- max(lgd.layout, na.rm = TRUE) + 1
            # Add an empty grob in the legend to stack them to the top
            legends <- c(legends, list(grid::textGrob("")))
        }
        right.legends <- gridExtra::arrangeGrob(
            grobs = legends, layout_matrix = lgd.layout)
        # Set default legend width space
        grob_list <- c(list(main_grob, right.legends), list(override_grob_list))
        grob_list <- grob_list[!vapply(
            X = grob_list, FUN = is.null, FUN.VALUE = logical(length = 1L))]
        arranged.grob <- gridExtra::arrangeGrob(
            grobs = grob_list, ncol = 2, widths = c(20, 2 + lgd.space.width),
            heights = override_grob_height)
    } else { arranged.grob <- main_grob }
    # # Plot result
    # grid::grid.draw(arranged.grob)
    return(arranged.grob)
}

#' Prepares gg2heatmap() results
#'
#' @return A \code{list} of grobs and gg objects composing the heatmap.
#' @author Yoann Pageaud.
#' @keywords internal

prepare_plot_export <- function(
    return_plots, dd.rows, dd.cols, show.annot, ddgr_seg_row, ddgr_seg_col,
    col_sidebar, htmp, htmp_legend, final.plot){
    if(return_plots == "all"){
        if(dd.rows & dd.cols & show.annot){
            grob_list <- list(
                "heatmap_legends" = htmp_legend,
                "annotation_legends" = col_sidebar$legends)
            ggplot_list <- list(
                "rows_dendrogram" = ddgr_seg_row,
                "cols_dendrogram" = ddgr_seg_col,
                "annotation_bars" = col_sidebar$sidebar, "heatmap" = htmp)
        } else if(dd.rows & dd.cols & !show.annot){
            grob_list <- list("heatmap_legends" = htmp_legend)
            ggplot_list <- list(
                "rows_dendrogram" = ddgr_seg_row,
                "cols_dendrogram" = ddgr_seg_col, "heatmap" = htmp)
        } else if(!dd.rows & !dd.cols & show.annot){
            grob_list <- list(
                "heatmap_legends" = htmp_legend,
                "annotation_legends" = col_sidebar$legends)
            ggplot_list <- list(
                "annotation_bars" = col_sidebar$sidebar, "heatmap" = htmp)
        } else if(!dd.rows & !dd.cols & !show.annot){
            grob_list <- list("heatmap_legends" = htmp_legend)
            ggplot_list <- list("heatmap" = htmp)
        } else if(dd.rows & !dd.cols & show.annot){
            grob_list <- list(
                "heatmap_legends" = htmp_legend,
                "annotation_legends" = col_sidebar$legends)
            ggplot_list <- list(
                "rows_dendrogram" = ddgr_seg_row,
                "annotation_bars" = col_sidebar$sidebar, "heatmap" = htmp)
        } else if(dd.rows & !dd.cols & !show.annot){
            grob_list <- list("heatmap_legends" = htmp_legend)
            ggplot_list <- list(
                "rows_dendrogram" = ddgr_seg_row, "heatmap" = htmp)
        } else if(!dd.rows & dd.cols & show.annot){
            grob_list <- list(
                "heatmap_legends" = htmp_legend,
                "annotation_legends" = col_sidebar$legends)
            ggplot_list <- list(
                "cols_dendrogram" = ddgr_seg_col,
                "annotation_bars" = col_sidebar$sidebar, "heatmap" = htmp)
        } else if(!dd.rows & dd.cols & !show.annot){
            grob_list <- list("heatmap_legends" = htmp_legend)
            ggplot_list <- list(
                "cols_dendrogram" = ddgr_seg_col, "heatmap" = htmp)
        }
        ls.res <- list(
            "main plot" = final.plot, "ggplots" = ggplot_list,
            "legends" = grob_list)
    } else if(return_plots == "main"){ ls.res <- final.plot }
    return(ls.res)
}

#' Proof of concept for making a multi-panel clustered ggplot2 heatmap.
#'
#' @param ls_col_ggdend A \code{list} of ggplot2 dendrograms on matrix colums.
#' @param ls_ggannot    A \code{list} of ggplot2 geom_tile() serving as color
#'                      annotation sidebars.
#' @param ls_row_ggdend A \code{list} of ggplot2 dendrograms on matrix rows.
#' @param ls_ggheatmap  A \code{list} of ggplot2 geom_tile() serving as
#'                      heatmaps.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @keywords internal

ggclust_layout <- function(
    ls_col_ggdend, ls_ggannot, ls_row_ggdend, ls_ggheatmap){
    cut_range <- sqrt(length(ls_htmp))
    ls_htmp <- lapply(
        X = seq(1, length(ls_htmp) - cut_range + 1, by = cut_range),
        FUN = function(i){
            ls_htmp[seq(i, i+cut_range-1)]
        })
    ls_row_htmp <- Map(list, ls_ddgr_row, ls_htmp)
    ls_row_htmp <- lapply(X = ls_row_htmp, FUN = function(i){
        c(list(i[[1]]), i[[2]])
    })
    ls_row_htmp <- do.call(c, ls_row_htmp)
    # Create empty ggplot
    ggempty <- ggplot2::ggplot(data.frame()) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "cm"),
            panel.background = ggplot2::element_rect(fill = NA))
    ls_plots <- c(
        list(ggempty), ls_ddgr_col, list(ggempty), ls_sidebar, ls_row_htmp)
    clustmap <- egg::ggarrange(plots = ls_plots, nrow = 5, ncol = 4)
    return(clustmap)
}
