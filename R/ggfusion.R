
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

#' Draws 2 triangle matrices fused together in a single plot.
#'
#' @param param1 A \code{type} parameter description.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export

#TODO: Write documentation!
ggfusion.free <- function(
    sample.names, upper.mat, lower.mat, order.select, order.method,
    hclust.method, dendrograms = FALSE,
    annot.grps = list("Groups" = seq(length(sample.names))),
    annot.pal = grDevices::rainbow(n = length(sample.names)), annot.size = 1,
    lgd.merge = FALSE, lgd.ncol = NULL,
    lgd.space.height = 26, lgd.space.width = 1,
    dend.size = 1, grid_col = "white", grid_linewidth = 0.3,
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
    # Check logicals for dendrograms
    if(length(dendrograms) == 1){
        dd.rows <- dendrograms
        dd.cols <- dendrograms
    } else if(length(dendrograms) == 2){
        dd.rows <- dendrograms[1]
        dd.cols <- dendrograms[2]
    } else { stop("'dendrograms' length > 2. Too many values.") }
    # Check dendrogram sizes
    if(length(dend.size) == 1){
        dend.row.size <- dend.size
        dend.col.size <- dend.size
    } else if(length(dend.size) == 2){
        dend.row.size <- dend.size[1]
        dend.col.size <- dend.size[2]
    } else { stop("'dend.size' length > 2. Too many values.") }
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
        } else {
            # Check dendrogram
            # if(dendro.pos != "none"){
            if(dd.rows | dd.cols){
                stop(paste(
                    "order != 'hclust' and dendrograms = TRUE. Dendrogram",
                     "cannot be generated if rows & cols are not ordered",
                     "following the hierarchical clustering."))
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
        } else if(order.method == 'default'){
            correlation.order <- seq(ncol(upper.mat))}
        # if(dendro.pos != "none"){
        if(dd.rows | dd.cols){
            # Generate Hierarchy Cluster
            hierarchy.clust <- fastcluster::hclust(
                d = stats::as.dist(1-upper.mat), method = hclust.method)
        }
    } else {
        if(order.method %in% c("AOE", "FPC", "hclust")){
            correlation.order <- corrplot::corrMatOrder(
                lower.mat, order = order.method, hclust.method = hclust.method)
        } else if(order.method == "alphabet"){
            # Fix bug of corrMatOrder when alphabet order
            correlation.order <- BiocompR::fix.corrMatOrder.alphabet(
                cor.order = correlation.order, str = colnames(upper.mat))
        } else if(order.method == 'default'){
            correlation.order <- seq(ncol(lower.mat))}
        # if(dendro.pos != "none"){
        if(dd.rows | dd.cols){
            # Generate Hierarchy Cluster
            hierarchy.clust <- fastcluster::hclust(
                d = stats::as.dist(1-lower.mat), method = hclust.method)
        }
    }
    # Generate Dendrogram
    # if(dendro.pos != "none"){
    if(dd.rows | dd.cols){
        ddgr <- stats::as.dendrogram(hierarchy.clust)
        ddgr_dat <- ggdendro::dendro_data(ddgr) # Dendrogram data
        if(dd.cols){
            ddgr_seg_col <- BiocompR::ggdend( # Get dendrogram segments
                df = ddgr_dat$segments, orientation = "top")
        }
        if(dd.rows){
            ddgr_seg_row <- BiocompR::ggdend( # Get dendrogram segments
                df = ddgr_dat$segments, orientation = "left", reverse.x = TRUE)
        }
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
        axis.text.y = ggplot2::element_text(size = 12, color = "black"),
        axis.text.y.right = ggplot2::element_blank(),
        axis.ticks.length.x.top = ggplot2::unit(0, "pt"),
        axis.ticks.length.y.right = ggplot2::unit(0, "pt"),
        legend.title = ggplot2::element_text(size = 13),
        legend.text = ggplot2::element_text(size = 9),
        legend.justification = c(0, 0.5),
        plot.margin = ggplot2::margin(0, 0, 0, 0))
    # Update theme
    upper_theme <- upper_default_fused + upper_theme
    # Update upper_scale_fill_grad na.value to add diagonal_col
    upper_scale_fill_grad$na.value <- diagonal_col
    # Upper plot
    upper.ggplot <- BiocompR::ggtriangle(
        melt_tri = upper.melt, grid_col = grid_col,
        grid_linewidth = grid_linewidth, ggtri_theme = upper_theme,
        scale_fill_grad = upper_scale_fill_grad,
        guide_custom_bar = upper_guide_custom_bar, y_axis_pos = "left") +
        ggplot2::ggtitle(plot_title)
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
    lower.ggplot <- BiocompR::ggtriangle(
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
    # Plot Color Sidebar
    col_sidebar <- BiocompR:::plot.col.sidebar(
        sample.names = sample.names[correlation.order], annot.grps = annot.grps,
        annot.pal = annot.pal, annot.pos = "top", theme_annot = ggplot2::theme(
            axis.ticks.x.top = ggplot2::element_line(color = "black"),
            axis.ticks.y.right = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x.top = ggplot2::element_text(
                size = 12, angle = 90, vjust = 0.5, hjust = 0, color = "black"),
            axis.text.y.right = ggplot2::element_text(
                size = 12, color = "black"),
            plot.margin = ggplot2::margin(0, 0, 0.1, 0, unit = "cm")),
        theme_legend = ggplot2::theme(
            legend.title = ggplot2::element_text(size = 13, color = "black"),
            legend.text = ggplot2::element_text(size = 12, color = "black")),
        dendro.pos = 'top', merge.lgd = lgd.merge, lgd.ncol = lgd.ncol,
        right = TRUE)
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
    # Create fused plot
    arranged.grob <- BiocompR:::heatmap_layout(
        dd.rows = dd.rows, dd.cols = dd.cols, show.annot = TRUE,
        ddgr_seg_row = ddgr_seg_row, ddgr_seg_col = ddgr_seg_col,
        sidebar = col_sidebar$sidebar, htmp = main_gg,
        legends = col_sidebar$legends, dend.row.size = dend.row.size,
        dend.col.size = dend.col.size*3,
        annot.size = annot.size, lgd.layout = lgd.layout,
        lgd.space.width = 1 + lgd.space.width,
        override_grob_list = lower.legend, override_grob_height = c(10, 1))
    # Plot Final Figure
    grid::grid.newpage()
    grid::grid.draw(arranged.grob)
    # Return res
    return(arranged.grob)
}

#' Draws 2 triangle matrices of computed pairwise correlations' results.
#'
#' @param data                  A \code{matrix} or \code{dataframe}.
#' @param ncores                An \code{integer} to specify the number of
#'                              cores/threads to be used to parallel-run tests.
#' @param upper.corr            The comparison for which results will be
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
#'   \item{'p' - the P-values of the comparison (-log10() transformed).}
#'   \item{'se' - the standard errors of the comparison.}
#'  }
#' @param lower.corr            The comparison for which results will be
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
#' @param order.method        A \code{character} specifying the ordering method
#'                            to apply.
#'                            \cr Possible ordering methods are:
#'                            \itemize{
#'                             \item{'AOE' - angular order of the eigenvectors.}
#'                             \item{'FPC' - first principal component order.}
#'                             \item{'hclust' - hierarchical clustering order.}
#'                             \item{'alphabet' - alphabetical order}
#'                             \item{'default'}
#'                            }
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
#' @param annot.size            An \code{integer} to increase or decrease the
#'                              size of the annotation side bar.
#' @param dendro.size           An \code{integer} to increase or decrease the
#'                              size of the dendrogram.
#' @param grid_col              A \code{character} specifying the color of the
#'                              grid.
#' @param grid_linewidth        A \code{double} value for the thickness of the
#'                              grid.
#' @param diagonal_col          A \code{character} defining the color of cells
#'                              with of the empty diagonal.
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
#' @examples
#' # Default fusion correlation plot showing Pearson correlations & P-values
#' mat <- as.matrix(t(scale(mtcars)))
#' res <- ggfusion.corr(data = mat)
#' # Same thing using the Spearman correlation
#' res <- ggfusion.corr(
#'     data = mat, upper.corr = "spearman", lower.corr = "spearman")
#' # Fuse correlations results from Pearson (top) and Spearman (bottom)
#' res <- ggfusion.corr(
#'     data = mat, upper.corr = "pearson", upper.value = "r",
#'     lower.corr = "spearman", lower.value = "r",
#'     lower_scale_fill_grad = ggplot2::scale_fill_gradientn(
#'         colors = BiocompR::biopalette(name = "RColorBrewer_RdBu8"),
#'         limits = c(-1, 1))) # Setting the scale of the Spearman correlations
#' # Set the missing data removal method to complete before computing
#' # correlations
#' res <- ggfusion.corr(data = mat, na.rm = "complete")
#'
#'
#' @references \href{https://www.scholars.northwestern.edu/en/publications/psych-procedures-for-personality-and-psychological-research}{William R Revelle, psych: Procedures for Personality and Psychological Research. Northwestern University, Evanston, Illinois, USA (2017).}
#' @references \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/#B13}{Why, When and How to Adjust Your P Values? Mohieddin Jafari and Naser Ansari-Pour - Cell J. 2019 Winter; 20(4): 604–607.}

#TODO: fix alphabet order error ggfusion.corr(data = mat, order.method = "alphabet")
ggfusion.corr <- function(
    data, upper.corr = "pearson", upper.value = "r", lower.corr = "pearson",
    lower.value = "p", na.rm = 'pairwise', order.method = 'hclust',
    order.select = 'upper', hclust.method = 'complete', p.adjust = "BH",
    dendrograms = FALSE, dend.size = 1,
    annot.grps = list("Groups" = seq(ncol(data))),
    annot.pal = grDevices::rainbow(n = ncol(data)), annot.size = 1,
    upper_theme = NULL, upper_scale_fill_grad = ggplot2::scale_fill_gradientn(
        colors = BiocompR::biopalette(name = "RColorBrewer_RdBu8"),
        name = "Upper\ntriangle values", limits = c(-1, 1)),
    upper_guide_custom_bar = ggplot2::guide_colorbar(
        barheight = 10, barwidth = 0.7, ticks.linewidth = 0.5,
        ticks.colour = "black", frame.linewidth = 0.5, frame.colour = "black"),
    lower_theme = NULL, lower_scale_fill_grad = ggplot2::scale_fill_gradient2(
        low = "darkblue", mid = "white", high = "darkred",
        midpoint = -log10(0.05), name = "Lower\ntriangle values"),
    lower_guide_custom_bar = ggplot2::guide_colorbar(
        barheight = 0.7, barwidth = 10, ticks.linewidth = 0.5,
        ticks.colour = "black", frame.linewidth = 0.5, frame.colour = "black"),
    plot_title = NULL, grid_col = "white", grid_linewidth = 0.3,
    diagonal_col = "white", lgd.merge = FALSE, lgd.ncol = NULL,
    lgd.space.height = 26, lgd.space.width = 1, ncores = 1, verbose = FALSE){
    # Fix BiocCheck() complaining about these objects initialization
    # statistic <- NULL

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
    # Check logicals for dendrograms
    if(length(dendrograms) > 2){
        stop("'dendrograms' length > 2. Too many values.")
    }
    # Groups and palettes matching
    BiocompR:::check.annotations(
        data = data, annot.grps = annot.grps, annot.pal = annot.pal,
        verbose = verbose)

    # # Position annotation
    # if(!annot.pos %in% c("top", "left", "both")){
    #     stop("annotation sidebar cannot be put here.")
    # }

    # Checking comparisons
    if(upper.corr %in% c("pearson", "spearman", "kendall")){
        #Overwrite "stat" by "t" for correlation
        if(upper.value == "stat"){ upper.value <- 't' }
        # Compute correlation
        if(verbose){
            cat("Compute pairwise", upper.corr, "correlation test...")
        }
        upper.correlation.res <- psych::corr.test(
            data, use = na.rm, method = upper.corr, adjust = p.adjust)
        upper.mat <- upper.correlation.res[[upper.value]]
        # if P-value, auto -log10 transform matrix
        if(upper.value == "p"){ upper.mat <- -log10(upper.mat) }
        upper.res <- upper.mat
        if(verbose){ cat("Done.\n") }
        if(lower.corr == upper.corr){
            # Overwrite "stat" by "t" for correlation
            if(lower.value == "stat"){ lower.value <- 't' }
            # Assign same matrix
            lower.correlation.res <- upper.correlation.res
            lower.mat <- lower.correlation.res[[lower.value]]

        } else {
            if(lower.corr %in% c( "pearson", "spearman", "kendall")){
                # Overwrite "stat" by "t" for correlation
                if(lower.value == "stat"){ lower.value <- 't' }
                # Compute correlation
                if(verbose){
                    cat("Compute pairwise", lower.corr, "correlation test...")
                }
                lower.correlation.res <- psych::corr.test(
                    data, use = na.rm, method = lower.corr, adjust = p.adjust)
                lower.mat <- lower.correlation.res[[lower.value]]
                if(verbose){ cat("Done.\n") }
            # } else if(lower.comp == "KS"){
            #     if(lower.value %in% c('n', 'stat', 'p')){
            #         if(verbose){
            #             cat("Compute pairwise", lower.comp, "test...")
            #         }
            #         # Compute Kolmogorov-Smirnov test
            #         ks.res <- BiocompR::pairwise.ks(
            #             data = data, statistic = lower.value, ncores = ncores)
            #         lower.mat <- ks.res$res.statistic
            #         cat("Done.\n")
            #     } else if(lower.value == 'r'){
            #         stop("a KS test does not compute a correlation value.")
            #     } else if(lower.value == 'se'){
            #         stop("a KS test does not compute a standard error.")
            #     } else { stop("Unknown statistic for 'lower.value'.") }
            } else { stop("'lower.corr' value not supported yet.") }
        }
        # if P-value, auto -log10 transform matrix
        if(lower.value == "p"){ lower.mat <- -log10(lower.mat) }
        lower.res <- lower.mat
        # Update colorbar names
        upper_scale_fill_grad$name <- paste(
            BiocompR:::simplecap(x = upper.corr),
            BiocompR:::metric_alias(upper.value), sep = "\n")
        lower_scale_fill_grad$name <- paste(
            BiocompR:::simplecap(x = lower.corr),
            BiocompR:::metric_alias(lower.value), sep = "\n")

    # } else if(upper.comp == "KS"){
    #     if(upper.value %in% c('n', 'stat', 'p')){
    #         cat("Compute pairwise", upper.comp, "test...")
    #         # Compute Kolmogorov Smirnov test
    #         ks.res <- BiocompR::pairwise.ks(
    #             data = data, statistic = upper.value, ncores = ncores)
    #         upper.mat <- ks.res$res.statistic
    #         cat("Done.\n")
    #     } else if(upper.value == 'r'){
    #         stop("a KS test does not compute a correlation value.")
    #     } else if(upper.value == 'se'){
    #         stop("a KS test does not compute a standard error.")
    #     } else { stop("Unknown statistic for 'upper.value'.") }
    #     upper.res <- upper.mat
    #     if(lower.comp == upper.comp){
    #         if(lower.value %in% c('n', 'stat', 'p')){
    #             lower.mat <- BiocompR::get.ks.stat(
    #                 table_combinations = ks.res$table_combinations,
    #                 df.ks.tests = ks.res$res.test, statistic)
    #         } else if(lower.value == 'r'){
    #             stop("a KS test does not compute a correlation value.")
    #         } else if(lower.value == 'se'){
    #             stop("a KS test does not compute a standard error.")
    #         } else { stop("Unknown statistic for 'upper.value'.") }
    #     } else {
    #         if(lower.comp %in% c("pearson", "spearman", "kendall")){
    #             #Overwrite "stat" by "t" for correlation
    #             if(lower.value == "stat"){ lower.value <- 't' }
    #             #Compute correlation
    #             cat("Compute pairwise", lower.comp, "correlation test...")
    #             lower.correlation.res <- psych::corr.test(
    #                 data, use = na.rm, method = lower.comp, adjust = p.adjust)
    #             lower.mat <- lower.correlation.res[[lower.value]]
    #             cat("Done.\n")
    #         } else { stop("'lower.comp' value not supported yet.") }
    #     }
    #     lower.res <- lower.mat
    } else { stop("'upper.corr' value not supported yet.") }

    # Create Fused Plot
    # fused.res <- BiocompR::ggfusion.free(
    fused.res <- ggfusion.free(
        sample.names = colnames(data), upper.mat = upper.mat,
        lower.mat = lower.mat, order.select = order.select,
        order.method = order.method, hclust.method = hclust.method,
        dendrograms = dendrograms,
        annot.grps = annot.grps, annot.pal = annot.pal, annot.size = annot.size,
        lgd.merge = lgd.merge, lgd.ncol = lgd.ncol,
        lgd.space.height = lgd.space.height, lgd.space.width = lgd.space.width,
        dend.size = dend.size, grid_col = grid_col,
        grid_linewidth = grid_linewidth, upper_theme = upper_theme,
        upper_scale_fill_grad = upper_scale_fill_grad,
        upper_guide_custom_bar = upper_guide_custom_bar,
        lower_theme = lower_theme,
        lower_scale_fill_grad = lower_scale_fill_grad,
        lower_guide_custom_bar = lower_guide_custom_bar,
        diagonal_col = diagonal_col, plot_title = plot_title, verbose = verbose)
    # Return plot and matrices
    return(list("upper.res" = upper.res, "lower.res" = lower.res,
                "fused.plot" = fused.res))
}
