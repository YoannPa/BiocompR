
#' Fixes a bug related to corrMatOrder when order parameter is set to
#' 'alphabet'.
#'
#' @param cor.order A \code{character} vector of interest.
#' @param str       A \code{character} vector from which to get elements order.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

fix.corrMatOrder.alphabet <- function(cor.order, str){
    unlist(lapply(X = cor.order, FUN = function(i){
        grep(pattern = paste0("^", i, "$"), x = str)
    }))
}

#' Displays 2 triangle matrices fused together in a single plot.
#'
#' @param param1 A \code{type} parameter description.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export

#TODO: Write documentation!
ggfusion.free <- function(
    sample.names, upper.mat, lower.mat,
    order.select, order.method, hclust.method,
    annot.grps = list("Groups" = seq(length(sample.names))),
    annot.pal = grDevices::rainbow(n = length(sample.names)), annot.size = 1,
    lgd.merge = FALSE, lgd.ncol = NULL,
    lgd.space.height = 26, lgd.space.width = 1, dendro.pos = 'none',
    dendro.size = 0, grid_col = "white", grid_linewidth = 0.3,
    upper_theme = NULL, upper_scale_fill_grad = ggplot2::scale_fill_gradientn(
        colors = BiocompR::biopalette(name = "viridis_C_plasma"),
        name = "Upper\ntriangle values"),
    upper_guide_custom_bar = ggplot2::guide_colorbar(
        barheight = 10, barwidth = 0.7, ticks.linewidth = 0.5,
        ticks.colour = "black", frame.linewidth = 0.5, frame.colour = "black"),
    lower_theme = NULL, lower_scale_fill_grad = ggplot2::scale_fill_gradientn(
        colors = BiocompR::biopalette(name = "viridis_D_viridis"),
        name = "Lower\ntriangle values"),
    lower_guide_custom_bar = ggplot2::guide_colorbar(
        barheight = 0.7, barwidth = 10, ticks.linewidth = 0.5,
        ticks.colour = "black", frame.linewidth = 0.5, frame.colour = "black"),
    diagonal_col = "white", plot_title = NULL, verbose = FALSE
){
    # Order Method
    if(order.method %in% c("AOE", "FPC", "hclust", "alphabet", "default")){
        if(verbose){ cat("Apply", order.method, "ordering method.\n") }
        if(order.method == "hclust"){
            # Check hclust.method
            if(hclust.method %in% c(
                'ward.D', 'ward.D2', 'single', 'complete', 'average',
                'mcquitty', 'median', 'centroid')){
                if(verbose){
                    cat("Apply", hclust.method, "clustering method.\n")
                }
            } else { stop("the clustering method is not supported.") }
            # Check dendrogram
            if(dendro.pos != "none"){
                if(!dendro.pos %in% c('top', 'left', 'both')){
                    stop("dendrogram cannot be put here.")
                }
            }
        } else {
            # Check dendrogram
            if(dendro.pos != "none"){
                stop(paste("order != 'hclust'. Dendrogram cannot be generated",
                           "if rows & cols are not ordered following the",
                           "hierarchical clustering."))
            }
        }
        # Order select
        if(order.method %in% c("AOE", "FPC", "hclust", "alphabet")){
            if(!order.select %in% c('upper', 'lower')){
                stop("Wrong value for order.select - You can only select",
                     "values from the 'upper' triangle or from the 'lower'",
                     "triangle to apply order.")
            }
        }
    } else { stop("the order method is not supported.") }
    # Check if length of labels match length of breaks in upper_scale_fill_grad
    # and lower_scale_fill_grad
    BiocompR:::chk.breaks.labels(ScaleContinuous_obj = upper_scale_fill_grad)
    BiocompR:::chk.breaks.labels(ScaleContinuous_obj = lower_scale_fill_grad)
    # Get order of the correlations for the method used
    if(order.select == 'upper'){
        if(order.method %in% c("AOE", "FPC", "hclust")){
            correlation.order <- corrplot::corrMatOrder(
                upper.mat, order = order.method, hclust.method = hclust.method)
        } else if(order.method == "alphabet"){
            # Fix bug of corrMatOrder when alphabet order
            correlation.order <- BiocompR::fix.corrMatOrder.alphabet(
                cor.order = correlation.order, str = colnames(upper.mat))
        }else if(order.method == 'default'){
            correlation.order <- seq(ncol(upper.mat))}
        if(dendro.pos != "none"){
            # Generate Hierarchy Cluster
            hierarchy.clust <- fastcluster::hclust(
                d = stats::as.dist(1-upper.mat), method = hclust.method)
        }
    } else {
        if(order.method %in% c("AOE","FPC","hclust")){
            correlation.order <- corrplot::corrMatOrder(
                lower.mat, order = order.method, hclust.method = hclust.method)
        } else if(order.method == "alphabet"){
            # Fix bug of corrMatOrder when alphabet order
            correlation.order <- BiocompR::fix.corrMatOrder.alphabet(
                cor.order = correlation.order, str = colnames(upper.mat))
        } else if(order.method == 'default'){
            correlation.order <- seq(ncol(lower.mat))}
        if(dendro.pos != "none"){
            # Generate Hierarchy Cluster
            hierarchy.clust <- fastcluster::hclust(
                d = stats::as.dist(1-lower.mat), method = hclust.method)
        }
    }
    # Generate Dendrogram
    if(dendro.pos != "none"){
        dendrogram <- stats::as.dendrogram(hierarchy.clust)
        ddgr_dat <- ggdendro::dendro_data(dendrogram) #Dendrogram data
        ddgr_seg <- BiocompR::ggdend( #Get dendrogram segments
            df = ddgr_dat$segments, orientation = dendro.pos, reverse.x = TRUE)
    }
    # Re-order rows and columns
    if(order.method != 'default'){
        upper.mat <- upper.mat[correlation.order, correlation.order]
        lower.mat <- lower.mat[correlation.order, correlation.order]
    }
    # Replace half matrices by NAs
    upper.mat[upper.tri(upper.mat)] <- NA
    lower.mat[lower.tri(lower.mat)] <- NA
    # Melt Correlation matrix
    dt.upper <- data.table::as.data.table(x = upper.mat, keep.rownames = "Var1")
    upper.melt <- data.table::melt.data.table(
        data = dt.upper, id.vars = "Var1", variable.name = "Var2", na.rm = TRUE,
        measure.vars = colnames(dt.upper)[-c(1)])
    upper.melt[, Var1 := as.factor(x = Var1)]
    upper.melt[, Var1 := factor(
        x = Var1, levels = as.character(unique(upper.melt$Var1)))]
    dt.lower <- data.table::as.data.table(x = lower.mat, keep.rownames = "Var1")
    lower.melt <- data.table::melt.data.table(
        data = dt.lower, id.vars = "Var1", variable.name = "Var2", na.rm = TRUE,
        measure.vars = colnames(dt.lower)[-c(1)])
    lower.melt[, Var1 := as.factor(x = Var1)]
    lower.melt[, Var1 := factor(
        x = Var1, levels = as.character(unique(lower.melt$Var1)))]
    # Invert order of samples
    upper.melt[, Var2 := factor(x = Var2, levels = rev(levels(Var2)))]
    lower.melt[, Var2 := factor(x = Var2, levels = rev(levels(Var2)))]
    # Replace identical correlations by NA
    upper.melt[Var1 == Var2, value := NA]
    lower.melt[Var1 == Var2, value := NA]
    # Define default theme
    upper_default_fused <- ggplot2::theme(
        axis.text.x.top = ggplot2::element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.length.x.top = ggplot2::unit(0, "pt"),
        axis.ticks.length.y.right = ggplot2::unit(0, "pt"),
        legend.title = ggplot2:: element_text(size = 13),
        legend.text = ggplot2::element_text(size = 9),
        legend.justification = c(0, 0.5),
        plot.margin = ggplot2::margin(0, 0, 0, 0))
    # Update theme
    upper_theme <- upper_default_fused + upper_theme
    # Update upper_scale_fill_grad na.value to add diagonal_col
    upper_scale_fill_grad$na.value <- diagonal_col
    # Upper plot
    upper.ggplot <- ggtriangle(
        melt_tri = upper.melt, grid_col = grid_col,
        grid_linewidth = grid_linewidth, ggtri_theme = upper_theme,
        scale_fill_grad = upper_scale_fill_grad,
        guide_custom_bar = upper_guide_custom_bar, y_axis_pos = "left") +
        ggtitle(plot_title)
    # Define default theme
    lower_default_fused <- ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.x.top = ggplot2::element_blank(),
        axis.text.y.right = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.ticks.length = ggplot2::unit(0, "pt"),
        legend.direction = 'horizontal',
        legend.position = "bottom",
        legend.justification = c(0.4, 0.5),
        plot.margin = ggplot2::margin(0, 0, 0, 0))
    # Update theme
    lower_theme <- lower_default_fused + lower_theme
    # Update lower_scale_fill_grad na.value to add diagonal_col
    lower_scale_fill_grad$na.value <- diagonal_col
    # Lower plot
    lower.ggplot <- ggtriangle(
        melt_tri = lower.melt, grid_col = grid_col,
        grid_linewidth = grid_linewidth, ggtri_theme = lower_theme,
        scale_fill_grad = lower_scale_fill_grad,
        guide_custom_bar = lower_guide_custom_bar)
    # Reorder groups and convert as factors
    annot.grps <- lapply(X = annot.grps, FUN = function(i){
        if(!is.factor(i)){ factor(x = i, levels = unique(i)) } else { i }
    })
    annot.grps <- lapply(
        X = annot.grps, FUN = function(i){ i[correlation.order] })
    # Plot Color Sidebar
    col_sidebar <- BiocompR:::plot.col.sidebar(
        sample.names = sample.names[correlation.order], annot.grps = annot.grps,
        annot.pal = annot.pal, annot.pos = "top", theme_annot = theme(
            axis.ticks.x.top = ggplot2::element_line(color = "black"),
            axis.ticks.y.right = element_blank(),
            axis.title.x = element_blank(), axis.text.x.top = element_text(
                size = 12, angle = 90, vjust = 0.5, hjust = 0, color = "black"),
            axis.text.y.right = element_text(size = 12, color = "black"),
            axis.title.y = axis.title.y,
            plot.margin = margin(0, 0, 0.1, 0, unit = "cm")),
        theme_legend = theme(
            legend.title = element_text(size = 13, color = "black"),
            legend.text = element_text(size = 12, color = "black")),
        dendro.pos = 'top', merge.lgd = lgd.merge, right = TRUE)
    # Remove lower legends
    lower.ggplot.nolgd <- lower.ggplot +
        ggplot2::theme(legend.position = "none")
    sidebar.nolgd <- col_sidebar$sidebar
    # Create grob for lower matrix
    lower.grob <- ggplot2::ggplotGrob(lower.ggplot.nolgd)
    # Add lower ggplot as an annotation to the upper ggplot
    main_gg <- upper.ggplot + ggplot2::annotation_custom(lower.grob)
    # Get lower legends
    lower.legend <- BiocompR:::get.lgd(lower.ggplot)
    # Builds color tables for legends.
    col_table <- BiocompR:::build.col_table(
        annot.grps = annot.grps, annot.pal = annot.pal)
    if(lgd.merge){
        # Rbind list color tables because lgd.merge is TRUE
        col_table <- data.table::rbindlist(col_table, idcol = TRUE)
        if(!(is.list(annot.pal)) | length(annot.pal) == 1){
            # Remove duplicated colors
            col_table <- col_table[!duplicated(x = Cols)]
        }
        col_table <- list(col_table[, c("Grps", "Cols"), ])
    }
    # Calculate legend length
    lgd_sizes <- BiocompR:::get.len.legends(col_table = col_table)
    # Calculate legend columns
    if(is.null(lgd.ncol)){
        lgd.ncol <- ceiling(lgd_sizes/lgd.space.height)
    } else if(all.equal(lgd.ncol, as.integer(lgd.ncol)) == TRUE){
        lgd.ncol <- rep(lgd.ncol, times = length(lgd_sizes))
    }
    # Build legends layout
    lgd.layout <- BiocompR::build_legends_layout(
        col_table = col_table, height.lgds.space = lgd.space.height,
        ncol.override = as.numeric(lgd.ncol[1]))
    # Create fused plot
    arranged.grob <- BiocompR:::heatmap_layout(
        dd.rows = FALSE, dd.cols = FALSE, show.annot = TRUE,
        ddgr_seg_row = NULL, ddgr_seg_col = NULL, sidebar = col_sidebar$sidebar,
        htmp = main_gg, legends = col_sidebar$legends, dend.row.size = 0,
        dend.col.size = 0, annot.size = annot.size, lgd.layout = lgd.layout,
        lgd.space.width = lgd.space.width,
        override_grob_list = lower.legend, override_grob_height = c(10, 1))
    # Plot Final Figure
    grid::grid.newpage()
    grid::grid.draw(arranged.grob)
    # Return res
    return(arranged.grob)
}

#' Creates a plot summarizing results from 2 different pairwise comparisons.
#'
#' @param data                  A \code{matrix} or \code{dataframe}.
#' @param ncores                An \code{integer} to specify the number of
#'                              cores/threads to be used to parallel-run tests.
#' @param upper.comp            The comparison for which results will be
#'                              displayed in the upper triangle of the plot as a
#'                              \code{character} matching one of these:
#'                              'pearson','spearman','kendall'.
#' @param upper.value           A \code{character} matching the type of value to
#'                              display in the upper triangle of the plot
#'                              resulting from the comparison.\cr The possible
#'                              value types are:
#'  \itemize{
#'   \item{'r' - the correlation values.}
#'   \item{'n' - the number of cases used for the comparisons.}
#'   \item{'stat' - the values of the comparison statistic.}
#'   \item{'p' - the P-values of the comparison.}
#'   \item{'se' - the standard errors of the comparison.}
#'  }
#' @param lower.comp            The comparison for which results will be
#'                              displayed in the lower triangle of the plot as a
#'                              \code{character} matching one of these:
#'                              'pearson','spearman','kendall'.
#' @param lower.value           See 'upper.value' help.
#' @param na.rm                 A \code{character} to specify how to handle
#'                              missing values when calculating a correlation.
#'                              Possible types are 'pairwise' and 'complete'.
#'                              'pairwise' is the default value and will do
#'                              pairwise deletion of cases. 'complete' will
#'                              select just complete cases.
#' @param order.method          A \code{character} specifying the ordering
#'                              method to apply.
#'                              \cr Possible ordering methods are:
#' \itemize{
#'  \item{'AOE' - angular order of the eigenvectors.}
#'  \item{'FPC' - first principal component order.}
#'  \item{'hclust' - hierarchical clustering order.}
#'  \item{'alphabet' - alphabetical order}
#'  \item{'default'}
#' }
#' @param order.select          A \code{character} specifying comparison results
#'                              to use for the clustering.\cr Possible values
#'                              are 'upper' or 'lower'.
#' @param hclust.method         A \code{character} specifying the method to use
#'                              for the hierarchical clustering if the ordering
#'                              method is 'hclust'.\cr Possible methods are:
#'                              'ward.D','ward.D2', 'single', 'complete',
#'                              'average' (= UPGMA), 'mcquitty' (= WPGMA),
#'                              'median' (= WPGMC) or 'centroid' (= UPGMC).
#' @param p.adjust              A \code{character} specifying what adjustment
#'                              for multiple tests should be used.\cr
#'                              (Default: p.adjust = "BH"; Supported:
#'                              p.adjust = c("holm", "hochberg", "hommel",
#'                              "bonferroni", "BH", "BY", "fdr", "none")).
#' @param annot.grps            A \code{list} of vectors of groups to which
#'                              variables belongs for the annotation sidebars.
#'                              Vectors' lengths have to match the number of
#'                              variables.
#' @param annot.pal             A \code{vector} or a list of vectors containing
#'                              colors as characters for the annotation
#'                              sidebars. The length of vectors has to match the
#'                              number of levels of vectors listed in
#'                              'annot.grps'. If a list is provided, its length
#'                              must match the length of the list provided to
#'                              'annot.grps'.
#' @param annot.pos             A \code{character} specifying the position of
#'                              the annotation sidebar.\cr Possible values are:
#'                              'top', 'left' or 'both'.
#' @param annot.size            An \code{integer} to increase or decrease the
#'                              size of the annotation side bar.
#' @param dendro.pos            A \code{character} specifying the position of
#'                              dendrograms if the selected order is 'hclust'.
#'                              \cr Possible values are: 'top','left','none'.
#' @param dendro.size           An \code{integer} to increase or decrease the
#'                              size of the dendrogram.
#' @param grid_col              A \code{character} specifying the color of the
#'                              grid.
#' @param grid_linewidth        A \code{double} value for the thickness of the
#'                              grid.
#' @param axis.title            An \code{element_text} object to setup axes
#'                              titles.
#' @param axis.title.x          An \code{element_text} object to setup X axis
#'                              title.
#' @param axis.title.y          An \code{element_text} object to setup Y axis
#'                              title.
#' @param axis.text             An \code{element_text} object to setup axes
#'                              text.
#' @param axis.text.x           An \code{element_text} object to setup X axis
#'                              text.
#' @param axis.text.y           An \code{element_text} object to setup Y axis
#'                              text.
#' @param axis.ticks            An \code{element_line} object to setup the ticks
#'                              of the plot.
#' @param set.x.title           A \code{character}to be used as the title for
#'                              the X axis.
#' @param set.y.title           A \code{character}to be used as the title for
#'                              the Y axis.
#' @param set.lgd1.title        A \code{character}to be used as the title of the
#'                              legend.
#' @param set.lgd2.title        A \code{character}to be used as the title of the
#'                              legend.
#' @param diagonal_col              A \code{character} defining the color of cells
#'                              with of the empty diagonal.
#' @param lgd.pal1              A \code{character} vector of colors to use for
#'                              the palette of the upper triangle.
#' @param lgd.pal2              A \code{character} vector of colors to use for
#'                              the palette of the lower triangle.
#' @param lgd.title             An \code{element_text} object to setup both
#'                              legend titles.
#' @param lgd.title1            An \code{element_text} object to setup the
#'                              legend title of the upper triangle.
#' @param lgd.title2            An \code{element_text} object to setup the
#'                              legend title of the lower triangle.
#' @param lgd.text              An \code{element_text} object to setup the text
#'                              of both legends.
#' @param lgd.text1             An \code{element_text} object to setup the text
#'                              of the upper triangle legend.
#' @param lgd.text2             An \code{element_text} object to setup the text
#'                              of the lower triangle legend.
#' @param lgd.breaks            A \code{double} vector defining the graduation
#'                              position along the possible values.
#' @param lgd.breaks1           A \code{double} vector defining the graduation
#'                              position along the possible values of the upper
#'                              triangle.
#' @param lgd.breaks2           A \code{double} vector defining the graduation
#'                              position along the possible values of the lower
#'                              triangle.
#' @param lgd.labels            A \code{double} vector of values to map to
#'                              legends breaks. Will be displayed on the plot.
#' @param lgd.labels1           A \code{double} vector of values to map to the
#'                              upper triangle legend breaks. Will be displayed
#'                              on the plot.
#' @param lgd.labels2           A \code{double} vector of values to map to the
#'                              lower triangle legend breaks. Will be displayed
#'                              on the plot.
#' @param lgd.round             An \code{integer} indicating the number of
#'                              decimal places to be used for the default
#'                              legends labels.
#' @param lgd.round1            An \code{integer} indicating the number of
#'                              decimal places to be used for the default
#'                              upper triangle legend labels.
#' @param lgd.round2            An \code{integer} indicating the number of
#'                              decimal places to be used for the default
#'                              lower triangle legend labels.
#' @param lgd.limits            A \code{double} vector of length 2, giving the
#'                              upper and lower limits for mapping breaks and
#'                              colors to both triangle legends.
#' @param lgd.limits1           A \code{double} vector of length 2, giving the
#'                              upper and lower limits for mapping breaks and
#'                              colors to the upper triangle legend.
#' @param lgd.limits2           A \code{double} vector of length 2, giving the
#'                              upper and lower limits for mapping breaks and
#'                              colors to the lower triangle legend.
#' @param lgd.ticks.linewidth   A \code{double} value to specify the thickness
#'                              of legends ticks.
#' @param lgd.ticks.linewidth1  A \code{double} value to specify the thickness
#'                              of the upper triangle associated legend ticks.
#' @param lgd.ticks.linewidth2  A \code{double} value to specify the thickness
#'                              of the lower triangle associated legend ticks.
#' @param lgd.nbin              An \code{integer} specifying the number of bins
#'                              for drawing both colorbars. A smoother
#'                              colorbar results from a larger value.
#' @param lgd.nbin1             An \code{integer} specifying the number of bins
#'                              for drawing the upper triangle associated
#'                              colorbar. A smoother colorbar results from a
#'                              larger value.
#' @param lgd.nbin2             An \code{integer} specifying the number of bins
#'                              for drawing the lower triangle associated
#'                              colorbar. A smoother colorbar results from a
#'                              larger value.
#' @param lgd.height1           A \code{double} specifying the height of the
#'                              upper triangle associated colorbar.
#' @param lgd.height2           A \code{double} specifying the height of the
#'                              lower triangle associated colorbar.
#' @param lgd.width1            A \code{double} specifying the width of the
#'                              upper triangle associated colorbar.
#' @param lgd.width2            A \code{double} specifying the height of the
#'                              lower triangle associated colorbar.
#' @param lgd.frame.col         A \code{character} defining the color of the
#'                              frame of legends.
#' @param lgd.frame.linewidth   A \code{double} defining the thickness of the
#'                              frame of both legends.
#' @param lgd.frame.linewidth1  A \code{double} defining the thickness of the
#'                              frame of upper triangle associated legends.
#' @param lgd.frame.linewidth2  A \code{double} defining the thickness of the
#'                              frame of lower triangle associated legends.
#' @param raster                A \code{logical} to specify whether or not
#'                              colors in both legends should be rasterized.
#' @param raster1               A \code{logical} to specify whether or not
#'                              colors in the upper triangle associated legend
#'                              should be rasterized.
#' @param raster2               A \code{logical} to specify whether or not
#'                              colors in the lower triangle associated legend
#'                              should be rasterized.
#' @return A \code{gtable} object plotted automatically and a \code{list} of
#'        results from the 2 comparisons as 2 matrices, with the same gtable
#'        object:
#'        \itemize{
#'         \item{'upper.res' - the matrix result of pairwise comparisons
#'         displayed in the upper triangle.}
#'         \item{'lower.res' - the matrix result of pairwise comparisons
#'         displayed in the lower triangle.}
#'         \item{'fused.plot' - the plot.}
#'        }
#' @author Yoann Pageaud.
#' @export
#' @references \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/#B13}{Why, When and How to Adjust Your P Values? Mohieddin Jafari and Naser Ansari-Pour - Cell J. 2019 Winter; 20(4): 604–607.}

fused.plot <- function(
    data, ncores, upper.comp, upper.value, lower.comp, lower.value,
    na.rm = 'pairwise', order.method, order.select, hclust.method = 'complete',
    p.adjust = "BH", annot.grps = list("Groups" = seq(ncol(data))),
    annot.pal = grDevices::rainbow(n = ncol(data)), annot.pos = 'top',
    annot.size = 0, lgd.merge = FALSE,
    annot.split = FALSE,
    dendro.pos = 'none', dendro.size = 0,
    grid_col = "white",grid_linewidth = 0.3,
    axis.title = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 12),
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(color = "black"),
    set.x.title = NULL, set.y.title = NULL,
    set.lgd1.title = NULL, set.lgd2.title = NULL,
    diag.col = "white",
    lgd.pal1 = NULL, lgd.pal2 = NULL, lgd.title = ggplot2::element_blank(),
    lgd.title1 = ggplot2::element_blank(),
    lgd.title2 = ggplot2::element_blank(),
    lgd.text = ggplot2::element_text(size = 12),
    lgd.text1 = ggplot2::element_blank(), lgd.text2 = ggplot2::element_blank(),
    lgd.breaks = NULL, lgd.breaks1 = NULL, lgd.breaks2 = NULL,
    lgd.labels = NULL, lgd.labels1 = NULL, lgd.labels2 = NULL,
    lgd.round = NULL, lgd.round1 = 2, lgd.round2 = 2,
    lgd.limits = NULL, lgd.limits1 = NULL, lgd.limits2 = NULL,
    lgd.ticks = TRUE, lgd.ticks1 = NULL, lgd.ticks2 = NULL,
    lgd.ticks.linewidth = 2, lgd.ticks.linewidth1 = NULL,
    lgd.ticks.linewidth2 = NULL,
    lgd.nbin = NULL, lgd.nbin1 = NULL, lgd.nbin2 = NULL,
    lgd.height1 = 26, lgd.height2 = 1,
    lgd.width1 = 1, lgd.width2 = 30,
    lgd.frame.col = "grey", lgd.frame.linewidth = 1.5,
    lgd.frame.linewidth1 = NULL, lgd.frame.linewidth2 = NULL,
    raster = TRUE, raster1 = NULL, raster2 = NULL, verbose = FALSE){
    # Fix BiocCheck() complaining about these objects initialization
    statistic <- NULL
    # Data format
    if(!(is.matrix(data))){if(is.data.frame(data)){data<-as.matrix(data)
    } else { stop("data is neither a matrix or a dataframe.") } }
    # Upper value type
    if(!(upper.value %in% c('r', 'n', 'stat', 'p', 'se'))){
        stop("value type unknown for upper.value")
    }
    # Lower value type
    if(!(lower.value %in% c('r', 'n', 'stat', 'p', 'se'))){
        stop("value type unknown for lower.value")
    }
    # Comparison Type
    if(na.rm %in% c("pairwise", "complete")){
        if(verbose){ cat("Apply", na.rm, "deletion of missing data.\n") }
    } else { stop("the value of na.rm is not supported.") }

    # P-value adjustment method
    if(p.adjust %in% c(
        "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")){
        if(verbose){ cat("Apply", p.adjust, "multiple testing adjustment.\n") }
    } else { stop("multiple testing adjustment method is not supported.") }

    # Groups and palettes matching
    BiocompR:::check.annotations(
        data = data, annot.grps = annot.grps, annot.pal = annot.pal,
        verbose = verbose)

    # Position annotation
    if(!annot.pos %in% c("top", "left", "both")){
        stop("annotation sidebar cannot be put here.")
    }

    # Checking comparisons
    if(upper.comp %in% c("pearson", "spearman", "kendall")){
        #Overwrite "stat" by "t" for correlation
        if(upper.value == "stat"){ upper.value <- 't' }
        # Compute correlation
        if(verbose){
            cat("Compute pairwise", upper.comp, "correlation test...")
        }
        upper.correlation.res <- psych::corr.test(
            data, use = na.rm, method = upper.comp, adjust = p.adjust)
        upper.mat <- upper.correlation.res[[upper.value]]
        upper.res <- upper.mat
        if(verbose){ cat("Done.\n") }
        if(lower.comp == upper.comp){
            # Overwrite "stat" by "t" for correlation
            if(lower.value == "stat"){ lower.value <- 't' }
            # Assign same matrix
            lower.correlation.res <- upper.correlation.res
            lower.mat <- lower.correlation.res[[lower.value]]
        } else {
            if(lower.comp %in% c( "pearson", "spearman", "kendall")){
                # Overwrite "stat" by "t" for correlation
                if(lower.value == "stat"){ lower.value <- 't' }
                # Compute correlation
                if(verbose){
                    cat("Compute pairwise", lower.comp, "correlation test...")
                }
                lower.correlation.res <- psych::corr.test(
                    data, use = na.rm, method = lower.comp, adjust = p.adjust)
                lower.mat <- lower.correlation.res[[lower.value]]
                if(verbose){ cat("Done.\n") }
            } else if(lower.comp == "KS"){
                if(lower.value %in% c('n', 'stat', 'p')){
                    if(verbose){
                        cat("Compute pairwise", lower.comp, "test...")
                    }
                    # Compute Kolmogorov-Smirnov test
                    ks.res <- BiocompR::pairwise.ks(
                        data = data, statistic = lower.value, ncores = ncores)
                    lower.mat <- ks.res$res.statistic
                    cat("Done.\n")
                } else if(lower.value == 'r'){
                    stop("a KS test does not compute a correlation value.")
                } else if(lower.value == 'se'){
                    stop("a KS test does not compute a standard error.")
                } else { stop("Unknown statistic for 'lower.value'.") }
            } else { stop("'lower.comp' value not supported yet.") }
        }
        lower.res <- lower.mat
    } else if(upper.comp == "KS"){
        if(upper.value %in% c('n', 'stat', 'p')){
            cat("Compute pairwise", upper.comp, "test...")
            # Compute Kolmogorov Smirnov test
            ks.res <- BiocompR::pairwise.ks(
                data = data, statistic = upper.value, ncores = ncores)
            upper.mat <- ks.res$res.statistic
            cat("Done.\n")
        } else if(upper.value == 'r'){
            stop("a KS test does not compute a correlation value.")
        } else if(upper.value == 'se'){
            stop("a KS test does not compute a standard error.")
        } else { stop("Unknown statistic for 'upper.value'.") }
        upper.res <- upper.mat
        if(lower.comp == upper.comp){
            if(lower.value %in% c('n', 'stat', 'p')){
                lower.mat <- BiocompR::get.ks.stat(
                    table_combinations = ks.res$table_combinations,
                    df.ks.tests = ks.res$res.test, statistic)
            } else if(lower.value == 'r'){
                stop("a KS test does not compute a correlation value.")
            } else if(lower.value == 'se'){
                stop("a KS test does not compute a standard error.")
            } else { stop("Unknown statistic for 'upper.value'.") }
        } else {
            if(lower.comp %in% c("pearson", "spearman", "kendall")){
                #Overwrite "stat" by "t" for correlation
                if(lower.value == "stat"){ lower.value <- 't' }
                #Compute correlation
                cat("Compute pairwise", lower.comp, "correlation test...")
                lower.correlation.res <- psych::corr.test(
                    data, use = na.rm, method = lower.comp, adjust = p.adjust)
                lower.mat <- lower.correlation.res[[lower.value]]
                cat("Done.\n")
            } else { stop("'lower.comp' value not supported yet.") }
        }
        lower.res <- lower.mat
    } else { stop("'upper.comp' value not supported yet.") }

    # Create Fused Plot
    fused.res <- BiocompR::fused.view(
        sample.names = colnames(data), upper.mat = upper.mat,
        lower.mat = lower.mat, order.select = order.select,
        order.method = order.method, hclust.method = hclust.method,
        annot.grps = annot.grps, annot.pal = annot.pal, annot.pos = annot.pos,
        annot.size = annot.size,
        lgd.merge = lgd.merge, annot.split = annot.split,
        dendro.pos = dendro.pos, dendro.size = dendro.size,
        grid_col = grid_col, grid_linewidth = grid_linewidth,
        axis.title = axis.title, axis.title.x = axis.title.x,
        axis.title.y = axis.title.y, axis.text = axis.text,
        axis.text.x = axis.text.x, axis.text.y = axis.text.y,
        axis.ticks = axis.ticks, set.x.title = set.x.title,
        set.y.title = set.y.title, set.lgd1.title = set.lgd1.title,
        set.lgd2.title = set.lgd2.title, diagonal_col = diagonal_col,
        lgd.pal1 = lgd.pal1, lgd.pal2 = lgd.pal2,
        lgd.title = lgd.title, lgd.title1 = lgd.title1,
        lgd.title2 = lgd.title2, lgd.text = lgd.text, lgd.text1 = lgd.text1,
        lgd.text2 = lgd.text2, lgd.breaks = lgd.breaks,
        lgd.breaks1 = lgd.breaks1, lgd.breaks2 = lgd.breaks2,
        lgd.labels = lgd.labels, lgd.labels1 = lgd.labels1,
        lgd.labels2 = lgd.labels2, lgd.round = lgd.round,
        lgd.round1 = lgd.round1, lgd.round2 = lgd.round2,
        lgd.limits = lgd.limits, lgd.limits1 = lgd.limits1,
        lgd.limits2 = lgd.limits2, lgd.ticks = lgd.ticks,
        lgd.ticks1 = lgd.ticks1, lgd.ticks2 = lgd.ticks2,
        lgd.ticks.linewidth = lgd.ticks.linewidth,
        lgd.ticks.linewidth1 = lgd.ticks.linewidth1,
        lgd.ticks.linewidth2 = lgd.ticks.linewidth2,
        lgd.nbin = lgd.nbin, lgd.nbin1 = lgd.nbin1, lgd.nbin2 = lgd.nbin2,
        lgd.height1 = lgd.height1, lgd.height2 = lgd.height2,
        lgd.width1 = lgd.width1, lgd.width2 = lgd.width2,
        lgd.frame.col = lgd.frame.col,
        lgd.frame.linewidth = lgd.frame.linewidth,
        lgd.frame.linewidth1 = lgd.frame.linewidth1,
        lgd.frame.linewidth2 = lgd.frame.linewidth2,
        raster = raster, raster1 = raster1, raster2 = raster2)

    return(list("upper.res" = upper.res, "lower.res" = lower.res,
                "fused.plot" = fused.res))
}
