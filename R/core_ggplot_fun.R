
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


#' Draws a ggplot for a basic upper or lower triangle.
#'
#' @param melt.tri            A \code{data.frame} melted triangle containing
#'                            test values.
#' @param grid.col            A \code{character} specifying the color of the
#'                            grid.
#' @param grid.thickness      A \code{double} value for the thickness of the
#'                            grid.
#' @param lgd.title           An \code{element_text} object to setup the legend
#'                            title text.
#' @param lgd.text            An \code{element_text} object to setup the legend
#'                            labels text.
#' @param lgd.pal             A \code{character} vector of colors to use for the
#'                            palette of the triangle plot.
#' @param min_tri             A \code{numeric} defining the minimum limit for
#'                            the legend of the values in the triangle plot.
#' @param max_tri             A \code{numeric} defining the maximum limit for
#'                            the legend of the values in the triangle plot.
#' @param lgd.breaks          An \code{integer} defining the number of breaks
#'                            wanted in the legend.
#' @param lgd.round           An \code{integer} indicating the number of decimal
#'                            places to be used for the default legends labels.
#' @param lgd.ticks           A \code{logical} to specify wether ticks should be
#'                            diplayed on both axes.
#' @param lgd.nbin            An \code{integer} specifying the number of bins
#'                            for drawing both colorbars. A smoother colorbar
#'                            results from a larger value.
#' @param lgd.height          A \code{double} specifying the height of the
#'                            colorbar.
#' @param lgd.width           A \code{double} specifying the width of the
#'                            colorbar.
#' @param rasteri             A \code{logical} to specify whether or not colors
#'                            in the legend should be rasterized.
#' @param lgd.ticks.linewidth A \code{double} value to specify the thickness of
#'                            legends ticks.
#' @param lgd.frame.col       A \code{character} defining the color of the
#'                            frame of legends.
#' @param lgd.frame.linewidth A \code{double} defining the thickness of the
#'                            frame of both legends.
#' @param diag.col            A \code{character} defining the color of cells
#'                            with of the empty diagonal.
#' @param set.lgd.title       A \code{character} to used as the name of the
#'                            legend.
#' @return A \code{gg} object of a basic triangle plot (a 'geom_tile()').
#' @author Yoann Pageaud.
#' @keywords internal

basic.ggplot.tri <- function(
    melt.tri, grid.col, grid.thickness, lgd.title, lgd.text, lgd.pal, min_tri,
    max_tri, lgd.breaks, lgd.round, lgd.ticks, lgd.nbin, lgd.height, lgd.width,
    rasteri, lgd.ticks.linewidth, lgd.frame.col, lgd.frame.linewidth, diag.col,
    set.lgd.title){
    #Fix BiocCheck() complaining about these objects initialization
    Var1 <- NULL
    Var2 <- NULL
    value <- NULL
    #Create the triangle plot
    ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = melt.tri, ggplot2::aes(x = Var1, y = Var2, fill = value),
            color = grid.col, size=grid.thickness) +
        ggplot2::theme(
            legend.title = lgd.title, legend.text = lgd.text,
            legend.justification = c(1, 0),
            plot.margin = ggplot2::margin(0, 0, 0, 0),
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "transparent"),
            plot.background = ggplot2::element_rect(fill = "transparent")) +
        ggplot2::scale_x_discrete(position = "top",expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::scale_fill_gradientn(
            colours = lgd.pal,
            breaks = seq(min_tri, max_tri,length.out = lgd.breaks),
            labels = round(
                seq(min_tri, max_tri, length.out = lgd.breaks), lgd.round),
            guide = ggplot2::guide_colourbar(
                ticks = lgd.ticks, nbin = lgd.nbin, barheight = lgd.height,
                label = TRUE, barwidth = lgd.width, raster = rasteri,
                ticks.linewidth = lgd.ticks.linewidth,
                frame.colour = lgd.frame.col,
                frame.linewidth = lgd.frame.linewidth), na.value = diag.col,
            limits = c(min_tri, max_tri), name = set.lgd.title)
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
#' @return A \code{matrix} to be used as a layout for legend grobs display on
#'         final plot.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Building the simplest layout with 1 legend and 2 categories
#' build.layout(col_table = list(data.table::data.table(
#'   'Grps' = c('cat 1', 'cat 2'), 'Cols' = c('blue', 'red'))),
#'   height.lgds.space = 29)
#' #If more annotations are passed in the list height.lgds.space let the
#' # function stacks legends up to a total size of 29 units on 1 column, before
#' # generating a new column if there is not enough space.

build.layout <- function(col_table, height.lgds.space){
    # Calculate legends length
    lgd_sizes <- BiocompR:::get.len.legends(col_table = col_table)
    #Compute the number of columns necessary for display of large legends
    ncol.by.lgd <- ceiling(lgd_sizes/height.lgds.space)
    # Calculate spatial disposition of legends based on their size and the
    # number of columns they occupy to create layout matrix by column
    legend_ids <- seq(length(col_table))
    columns <- lapply(X = seq(length(col_table)), FUN = function(i){
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
#' @param data      A \code{data.frame} with the column names 'Samples','.id'
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

basic.sidebar <- function(
    data, palette, annot.sep = 0, annot.cut = 0, lgd.ncol = 1, facet = NULL){
    #Fix BiocCheck() complaining about these objects initialization
    .id <- NULL
    facet.annot <- NULL
    Groups <- NULL
    Samples <- NULL
    #Create a basic sidebar plot
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
        #Add faceting
        dt.facet <- data[.id == facet]
        dt.facet[, `:=`(facet.annot, Groups)]
        data <- merge(
            x = data, y = dt.facet[, c("Samples", "facet.annot"), ],
            by = "Samples", all.x = TRUE)
        #Plot annotation bar with facet
        basic <- basic + ggplot2::geom_tile(data = data, mapping = ggplot2::aes(
            x = Samples, y = .id, fill = Groups, height = 1 - annot.sep,
            width = 1 - annot.cut)) +
            ggplot2::facet_grid(
                . ~ facet.annot, scales = "free", space = "free") +
            ggplot2::theme(
                panel.spacing = ggplot2::unit(0, "lines"),
                strip.background = ggplot2::element_rect(
                    color = "black", size = 0.5))
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
#' @export plot.col.sidebar
#' @export

plot.col.sidebar <- function(
    sample.names, annot.grps, annot.pal, annot.pos = "top", annot.sep = 0,
    annot.cut = 0, merge.lgd = FALSE, right = FALSE, lgd.name = "Legends",
    lgd.ncol = 1, theme_legend = ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_blank()),
    theme_annot = ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12)),
    set.x.title, set.y.title, dendro.pos, facet = NULL){
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
    # groups <- lapply(X = lapply(X = annot.grps, FUN = function(i){
    #     factor(x = i, levels =  unique(i))}), FUN = levels)

    # Create list of color tables
    if(is.list(annot.pal)){ # If a list of palettes is provided
        # Map categories to palettes
        col_table <- BiocompR:::map.cat2pal(origin.grps, groups, annot.pal)
    } else if(!is.list(annot.pal)){ #if a single palette is provided
        # Map groups to the same palette
        ls.df.grp.pal <- lapply(X = origin.grps, FUN = function(grp){
            data.table::data.table(
                "Grps" = grp, "Cols" = annot.pal, stringsAsFactors = FALSE)
            # data.frame(
            #     "Grps" = grp, "Cols" = annot.pal, stringsAsFactors = FALSE)
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
        # data.frame("Samples" = sample.names, "Groups" = i)
    })
    # Order samples following the correlation order provided
    dframe.annot <- lapply(dframe.annot, function(i){
        dt_annot <- i
        sample_order <- dt_annot$Samples
        dt_annot[, Samples := as.factor(Samples)]
        dt_annot[, Samples := factor(x = Samples, levels = sample_order)]
        dt_annot
    })
    # # Order samples following the correlation order provided and categories by
    # # alphabetical order
    # dframe.annot <- lapply(dframe.annot, function(i){
    #     i[["Samples"]] <- factor(i[["Samples"]], levels = i[["Samples"]])
    #     i[["Groups"]] <- factor(i[["Groups"]], levels = unique(i[["Groups"]]))
    #     i
    # })

    # Rbind all annotations
    dframe.annot <- data.table::rbindlist(dframe.annot, idcol = TRUE)
    unique_id <- unique(dframe.annot$.id)
    # Convert .id as factors
    dframe.annot[, .id :=  as.factor(.id)]
    dframe.annot[, .id :=  factor(x = .id, levels = unique_id)]
    # dframe.annot$.id <- factor(
    #     x = dframe.annot$.id, levels = unique(dframe.annot$.id))

    if(annot.pos == "top"){ # Change order of levels
        dframe.annot[, .id :=  factor(x = .id, levels = rev(unique_id))]
        # dframe.annot$.id <- factor(
        #     x = dframe.annot$.id, levels = rev(levels(dframe.annot$.id)))
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
    #Plot color sidebars
    col_sidebar <- BiocompR::basic.sidebar(
        data = data.table::copy(dframe.annot), palette = col_table$Cols,
        annot.sep = annot.sep, annot.cut = annot.cut, facet = facet,
        lgd.ncol = lgd.ncol[1])
    #Add legend theme parameters if some
    col_sidebar <- col_sidebar + theme_legend
    #Modify base plot following its position
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
            dframe.annot[, .id :=  factor(x = .id, levels = unique_id)]
            # dframe.annot$.id <- factor(
            #     dframe.annot$.id, levels = rev(levels(dframe.annot$.id)))
        }
        #Generate separate legends if more than 1 palette available or if only 1
        # annotation is used
        if((is.list(annot.pal) & length(annot.pal) > 1) | length(
            levels(dframe.annot$.id)) == 1){
            #Get all legends separately
            sidebar.lgd <- lapply(
                X = seq_along(levels(dframe.annot$.id)), FUN = function(i){
                    BiocompR:::get.lgd(BiocompR::basic.sidebar(
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
