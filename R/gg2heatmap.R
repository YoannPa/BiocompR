#' Creates a custom heatmap with dendrograms and annotations.
#'
#' @param m               A \code{matrix} with row and column names set, or a
#'                        molten \code{dataframe}:
#'                        \itemize{
#'                         \item{Column 1 will be use for row names.}
#'                         \item{Column 2 will be use for column names.}
#'                         \item{Column 3 will be use for values.}
#'                         \item{Additional columns coming after can be use for
#'                         categorical splitting heatmap on rows with the
#'                         'split.by.rows' option.}
#'                        }
#' @param na.handle       A \code{character} to specify how missing values
#'                        should be handled.
#'                        \itemize{
#'                         \item{'keep' do not modify the matrix}
#'                         \item{'impute' make use of the groups provided in
#'                               'imputation.grps' to calculate the mean value
#'                               of a group of columns for a given row where at
#'                               least one missing value has been found. The
#'                               missing values are then substituted by the
#'                               calculated mean of the group.}
#'                         \item{'remove' removes all rows containing at least
#'                               1 missing value.}
#'                        }
#'                        (Default: na.handle = 'remove';
#'                        Supported: na.handle = c('keep','impute','remove')).
#' @param dist.method     A \code{character} vector to specify the distance
#'                        methods to be used on the matrix rows and columns.
#'                        \itemize{
#'                         \item{If the vector is of length 1: the given method
#'                               will be applied on rows and columns of the
#'                               matrix.}
#'                         \item{If the vector is of length 2: the first method
#'                               will be applied on matrix rows, and the second
#'                               method will be applied on matrix columns.}
#'                        }
#'                        (Default: dist.method = 'manhattan';
#'                        Supported: dist.method = c('euclidean', 'maximum',
#'                        'manhattan', 'canberra', 'binary', 'minkowski',
#'                        'none')).
#' @param rank.fun        A \code{character} to specify a function to apply on
#'                        matrix rows to order them on the final heatmap.
#'                        Functions should be specified using the syntax
#'                        'package::function'. If rank.fun = NULL the order of
#'                        rows in the matrix will be kept
#'                        (Default: rank.fun = NULL;
#'                        Supported: rank.fun = 'stats::sd').
#' @param top.rows        An \code{integer} to specify the number of top rows
#'                        you want to display on the heatmap
#'                        (Default: top.rows = NULL).
#' @param dendrograms     A \code{logical} vector to specify whether dendrograms
#'                        should be plotted with the heatmap.
#'                        \itemize{
#'                         \item{If dendrograms = TRUE: both dendrograms on row
#'                               and columns of the matrix will be plotted.}
#'                         \item{If dendrograms = FALSE: both dendrograms on
#'                               rows and columns of the matrix will NOT be
#'                               plotted.}
#'                         \item{If the vector is of length 2: the first logical
#'                               will apply for rows, and the second logical
#'                               will apply for columns.}
#'                         \item{If 'rank.fun' and 'top.rows' are not NULL, and
#'                               dendrograms on rows set to TRUE, the ranking of
#'                               rows and the selection of the top ones will be
#'                               made before the computation of dendrograms.}
#'                        }
#' @param dend.size       A \code{numeric} vector defining row and column
#'                        dendrograms size (Default: dend.size = 1).
#'                        \itemize{
#'                         \item{If dend.size is a \code{numeric}: the value is
#'                               used to set both the width of the dendrogram
#'                               built on rows, and the height of the dendrogram
#'                               built on columns.}
#'                         \item{If dend.size is a \code{numeric} vector of
#'                               length 2: the first numeric will apply to the
#'                               dendrogram built on rows, and the second
#'                               numeric will apply to the dendrogram built on
#'                               columns.}
#'                        }
#' @param theme_dend_top  A ggplot2 \code{theme} to specify any theme parameter
#'                        you wish to custom the top dendrogram
#'                        (Default: theme_dend = NULL). For more information
#'                        about how to define a theme, see
#'                        \link[ggplot2]{theme}.
#' @param theme_dend_left A ggplot2 \code{theme} to specify any theme parameter
#'                        you wish to custom the left dendrogram
#'                        (Default: theme_dend = NULL). For more information
#'                        about how to define a theme, see
#'                        \link[ggplot2]{theme}.
#' @param imputation.grps A \code{character} vector defining groups to which
#'                        columns of the matrix belong. These groups are use to
#'                        make group-wise imputation of missing values between
#'                        columns. The vector has to be of the same length than
#'                        the number of columns in the matrix
#'                        (Default: imputation.grps = NULL).
#' @param ncores          An \code{integer} to specify the number of
#'                        cores/threads to be used to parallel-compute distances
#'                        for dendrograms.
#' @param plot.labs       A \code{labels} ggplot2 object to pass labels to be
#'                        displayed on the final heatmap plot. Currently
#'                        'plot.labs' supports 4 different labels:
#'                        \itemize{
#'                         \item{title - A \code{character} to be used as title
#'                         for the plot.}
#'                         \item{x - A \code{character} to specify X-axis title.
#'                         }
#'                         \item{y - A \code{character} to specify Y-axis title.
#'                         }
#'                         \item{legend - A \code{character} specifying the name
#'                         of annotation side bar legends, when all legends are
#'                         merged into a unique one (lgd.merge = TRUE).}
#'                        }
#'                        For more information on how to set labels for
#'                        'plot.labs' see \link[ggplot2]{labs}.
#' @param row.type        A \code{character} to be used in the plot subtitle
#'                        description as a definition of the rows (e.g. 'loci',
#'                        'samples', 'regions', etc.
#'                        Default: row.type = 'rows').
#' @param facet           A \code{character} matching an annotation name in
#'                        'annot.grps' to be used to split heatmap in separate
#'                        panels following the annotation.
#' @param split.by.rows   A \code{character} matching one of the categorical
#'                        column in m if m is a molten data.frame to split
#'                        heatmap on rows. split.by.rows will not work if m is
#'                        matrix.
#' @param border.col      A \code{character} to specify an R color code for the
#'                        border delimitating the heatmap cells
#'                        (Default: border.col = NA).
#' @param border.size     A \code{numeric} to specify the linewidth of cells
#'                        border (Default: border.size = 0.1).
#' @param cell.size       A \code{numeric} vector of length 2 to specify the
#'                        width and height of cells taking values between 0 and
#'                        1.
#'                        \itemize{
#'                         \item{If cell.size is a \code{numeric}: the value is
#'                               used to set both the width and the height of
#'                               cells.}
#'                         \item{If cell.size is a \code{numeric} vector of
#'                               length 2: the first numeric will apply to
#'                               cell height, and the second numeric will apply
#'                               to cell width.}
#'                        }
#' @param raster.filter   A \code{character} to be used as a filter for matrix
#'                        rasterization. The list of the supported rasterization
#'                        filters is available in magick::filter_types(). If no
#'                        value is specified for 'raster.filter' the original
#'                        ggplot2 heatmap is displayed
#'                        (Default: raster.filter = NULL).
#'                        Warning: Be aware that rasterization may take several
#'                        more minutes than the usual process time.
#' @param raster.size     An \code{integer} specifying the size of a squared
#'                        raster in pixels. If raster.size = 1080 -> raster =
#'                        1080x1080 (Default: raster.size = 1080).
#' @param theme_heatmap   A ggplot2 \code{theme} to specify any theme parameter
#'                        you wish to custom on the heatmap part of the plot
#'                        (Default: theme_heatmap = NULL). For more information
#'                        about how to define a theme, see
#'                        \link[ggplot2]{theme}.
#' @param guide_custom_bar A \code{guide} object generated by the ggplot2
#'                         function \link[ggplot2]{guide_colorbar} to custom the
#'                         heatmap color bar appearance
#'                         (see also 'scale_fill_grad' option).
#' @param scale_fill_grad A \code{ScaleContinous} object generated by ggplot2
#'                        functions such as \link[ggplot2]{scale_fill_gradient},
#'                        \link[ggplot2]{scale_fill_gradient2} or
#'                        \link[ggplot2]{scale_fill_gradientn} to customize
#'                        heatmap colors and the associated color bar.
#' @param annot.grps      A \code{list} of vectors of groups to which variables
#'                        belongs for the annotation sidebars. Vectors' lengths
#'                        have to match the number of variables.
#' @param annot.pal       A \code{vector} or a list of vectors containing colors
#'                        as characters for the annotation sidebars.The length
#'                        of vectors has to match the number of levels of
#'                        vectors listed in 'annot.grps'.
#'                        \itemize{
#'                         \item{If annot.pal is a list: its length must match
#'                               the length of the list provided to
#'                               'annot.grps'.}
#'                         \item{If annot.pal is a vector: make sure that the
#'                               levels content of annotations listed in
#'                               'annot.grps' is the same, and that no
#'                               annotation contains less or more levels than
#'                               another one in 'annot.grps'.}
#'                        }
#' @param annot.size      A \code{numeric} defining the width of the annotation
#'                        bars (Default: annot.size = 1).
#' @param annot.sep       A \code{numeric} vector specifying the width of the
#'                        separations between annotations
#'                        (Default: annot.sep = 0):
#'                        \itemize{
#'                         \item{If annot.sep is a \code{numeric}: the value is
#'                               used to set the width of both horizontal and
#'                               vertical separations of annotations.}
#'                         \item{If annot.sep is a \code{numeric} vector of
#'                               length 2: the first numeric will apply to the
#'                               width of horizontal separations, and the second
#'                               numeric will apply to the
#'                               width of vertical separations.}
#'                        }
#' @param theme_annot     A ggplot2 \code{theme} to specify any theme parameter
#'                        you wish to custom on the annotation bar
#'                        (Default: theme_annot = NULL). For more information
#'                        about how to define a theme, see
#'                        \link[ggplot2]{theme}.
#' @param show.annot      A \code{logical} to specify whether annotations should
#'                        be displayed at the top of the heatmap
#'                        (show.annot = TRUE) or not (show.annot = FALSE).
#' @param lgd.merge       A \code{logical} specifying whether the legends of
#'                        multiple annotation bars should be merged
#'                        (lgd.merge = TRUE) or remain separated
#'                        (lgd.merge = FALSE). lgd.merge is especially useful
#'                        when you want to map the same color palette to
#'                        multiple annotations sharing the same values
#'                        (Default: lgd.merge = FALSE).
#' @param theme_legend    A ggplot2 \code{theme} to specify any theme parameter
#'                        you wish to custom on legends
#'                        (Default: theme_legend = NULL). For more information
#'                        about how to define a theme, see
#'                        \link[ggplot2]{theme}.
#' @param lgd.space.width A \code{numeric} specifying the width of the legend
#'                        space (Default: lgd.space.width = 1).
#' @param lgd.space.height An \code{integer} specifying the height of the legend
#'                         space (Default: lgd.space.height = 29).
#' @param y.axis.right    A \code{logical} to specify whether the heatmap Y axis
#'                        should be displayed on the right (y.axis.right = TRUE)
#'                        or not (y.axis.right = FALSE)
#'                        (Default: y.axis.right = FALSE).
#' @param draw            A \code{logical} to specify whether the final heatmap
#'                        should be drawn automatically when gg2heatmap()
#'                        execution ends (draw = TRUE), or if it shouldn't
#'                        (draw = FALSE).
#' @param verbose         A \code{logical} to display information about the
#'                        step-by-step processing of the data if TRUE
#'                        (Default: verbose = FALSE).
#' @return A \code{grob} list containing the final plot, and also each grob
#'         generated separately.
#' @author Yoann Pageaud.
#' @importFrom data.table `:=` `.I`
#' @export
#' @examples
#' # Create the basic gg2heatmap
#' mat <- as.matrix(t(scale(mtcars)))
#' res <- gg2heatmap(m = mat)
#' # Apply Euclidean distance on rows, and Manhattan distance on columns
#' res <- gg2heatmap(m = mat, dist.method = c("euclidean", "manhattan"))
#' # Rank heatmap rows following their standard deviation
#' res <- gg2heatmap(m = mat, dist.method = c("none", "manhattan"),
#'   rank.fun = "stats::sd", dendrograms = c(FALSE, TRUE))
#' # Keep top 5 rows with highest standard variation
#' res <- gg2heatmap(m = mat, dist.method = c("none", "manhattan"),
#'   rank.fun = "stats::sd", dendrograms = c(FALSE, TRUE), top.rows = 5)
#' # Order rows and columns respectively by euclidean and manhattan distance,
#' # and hide both dendrograms
#' res <- gg2heatmap(m = mat, dist.method = c("euclidean", "manhattan"),
#'   dendrograms = FALSE)
#' # Change both dendrograms size
#' res <- gg2heatmap(m = mat, dist.method = c("euclidean", "manhattan"),
#'   dend.size = c(2, 5))
#' # Custom themes of top and left dendrograms
#' res <- gg2heatmap(
#'   m = mat, theme_dend_top = theme(
#'     panel.background = element_rect(
#'       color = "red", linewidth = 1, fill = "lightblue")),
#'   theme_dend_left = theme(panel.background = element_rect(fill = "gold")))
#' # Set heatmap title, and X and Y axis titles
#' res <- gg2heatmap(m = mat, plot.labs = labs(
#'   title = "mtcars example heatmap", x = "Cars", y = "Cars caracteristics"),
#'   y.axis.right = TRUE, # To enable the display of Y-axis on the right
#'   theme_heatmap = theme(
#'     axis.title.y.right = element_text(size = 14, vjust = -10)))
#' # Specify what 'rows' correspond to in the subtitle.
#' res <- gg2heatmap(m = mat, row.type = "caracteristics")
#' # Add 1 annotation for the top color bar
#' res <- gg2heatmap(
#'   m = mat, annot.grps = list("Carb" = mtcars$carb),
#'   annot.pal = grDevices::rainbow(n = 6))
#' # Add more annotations, adjust colorbar width, and legends positions
#' res <- gg2heatmap(
#'   m = mat,
#'   annot.grps = list( #Separate annotations must all have different values
#'     "Carb" = paste(mtcars$carb, "Carb"),
#'     "Gear" = paste(mtcars$gear, "Gear"),
#'     "Am" = paste(mtcars$am, "Am")),
#'   annot.pal = list( #Palettes order must matches annotations order
#'     grDevices::rainbow(n = 6), c("pink", "red", "darkred"),
#'     c("blue", "grey")), annot.size = 3, # Increases top colorbar width
#'   theme_legend = theme(
#'     legend.justification = c(0, 0.8)), # Justifies all legends on the left
#'     lgd.space.height = 22) # Sets vertical space to avoid legends overlaps
#' # Facetting on 1 of the annotations
#' res <- gg2heatmap(
#'   m = mat,
#'   annot.grps = list(
#'     "Carb" = paste(mtcars$carb, "Carb"),
#'     "Gear" = paste(mtcars$gear, "Gear"),
#'     "Am" = paste(mtcars$am, "Am")),
#'   annot.pal = list(
#'     grDevices::rainbow(n = 6), c("pink", "red", "darkred"),
#'     c("blue", "grey")),
#'     annot.size = 4, # Increases top colorbar width for facets strips
#'   theme_legend = theme(
#'     legend.justification = c(0, 0.8)), # Justifies all legends on the left
#'     lgd.space.height = 22,
#'   facet = "Carb", # Facetting heatmap using annotation 'Carb'
#'   dendrograms = FALSE, # Hide dendrogram since the facetting is On.
#'   dist.method = "none") # Disable clustering to let facetting work.
#' # Splitting on rows following a group annotation
#' molten.cars <- melt.data.table(
#'   data = as.data.table(mat, keep.rownames = TRUE), id.vars = "rn",
#'   variable.name = "cars")
#' molten.cars[
#'   rn %in% c('mpg', 'disp', 'hp', 'drat', 'qsec'), groups := "kinetic"]
#' molten.cars[
#'   rn %in% c('cyl', 'wt', 'vs', 'am', 'gear', 'carb'), groups := "mechanic"]
#' res <- gg2heatmap(
#'   m = molten.cars, split.by.rows = "groups", # Splitting on rows by 'groups'
#'   dist.method = c('none', 'euclidean'), # Disable clustering on rows
#'   dendrograms = c(FALSE, TRUE), # Disable dendrogram display on rows
#'   lgd.space.width = 5, # Increase legend space
#'   theme_heatmap = theme(
#'     axis.text.y.left = element_text(size = 12), # Show Y labels on the left
#'     axis.ticks.y.left = element_line(linewidth = 0.5)), # Show Y ticks on left
#'   theme_legend = theme(legend.justification = c(0.5, 0.5))) # Center legend
#' # Custom heatmap cells borders color and width
#' res <- gg2heatmap(m = mat, border.col = "black", border.size = 0.5)
#' # Custom heatmap cells height and width
#' res <- gg2heatmap(m = mat, cell.size = c(0.9, 0.5))
#' # Plot a rastered heatmap using the Lanczos filter
#' res <- gg2heatmap(m = mat, raster.filter = "Lanczos")
#' # Plot a rastered heatmap using the Lanczos filter changing the raster size
#' res <- gg2heatmap(m = mat, raster.filter = "Lanczos", raster.size = 240)
#' # Plot a rastered faceted heatmap using the Lanczos filter
#' res <- gg2heatmap(
#'   m = mat,
#'   annot.grps = list(
#'     "Carb" = paste(mtcars$carb, "Carb"),
#'     "Gear" = paste(mtcars$gear, "Gear"),
#'     "Am" = paste(mtcars$am, "Am")),
#'   annot.pal = list(
#'     grDevices::rainbow(n = 6), c("pink", "red", "darkred"),
#'     c("blue", "grey")),
#'     annot.size = 4, # Increases top colorbar width for facets strips
#'   theme_legend = theme(
#'     legend.justification = c(0, 0.8)), # Justifies all legends on the left
#'   lgd.space.height = 22,
#'   facet = "Carb", # Facetting heatmap using annotation 'Carb'
#'   dendrograms = FALSE, # Hide dendrogram since the facetting is On
#'   dist.method = "none", # Disable clustering to let facetting work
#'   raster.filter = "Lanczos") # Rasterization using the Lanczos filter
#' # Custom the theme of the heatmap
#' res <- gg2heatmap(
#'   m = mat,
#'   theme_heatmap = theme(
#'     panel.border = element_rect(
#'       color = "black", linewidth = 1, fill = "transparent"), # Add a border
#'     axis.text.y.right = element_text(size = 12), # Add right Y axis labels
#'     axis.ticks.y.right = element_line(color = "black"), # Add Y axis ticks
#'     axis.text.x = element_text(color = "red")), # Color X axis labels in red
#'   y.axis.right = TRUE) # Enable display of Y axis on the right side
#' # Custom the colorbar shape
#' res <- gg2heatmap(m = mat, guide_custom_bar = guide_colorbar(
#'   title = "My custom colorbar", # Set colorbar title
#'   barwidth = 10, # Set colorbar width
#'   barheight = 0.5, # Set colorbar height
#'   ticks.linewidth = 4, # Set ticks width
#'   ticks.colour = "red", # Set ticks color
#'   title.vjust = 1, # Set vertical justification of colorbar title
#'   raster.filter = TRUE, # Rasterize colorbar
#'   nbin = 3, # Set the number of color bins in the colorbar
#'   frame.colour = "green", # Set colorbar frame color
#'   frame.linewidth = 2)) # Set colorbar frame linewidth
#' # Custom colorbar caracteristics using 'scale_fill_gradient'-like functions
#' res <- gg2heatmap(m = mat, scale_fill_grad = scale_fill_gradientn(
#'   colors = c("green", "black", "red"), # Set gradient colors
#'   na.value = "grey", # Set missing value color
#'   n.breaks = 10, # Set number of breaks
#'   labels = c("<=-1", seq(-0.8, 0.8, by = 0.2), ">= 1"), # Map custom labels
#'   limits = c(-1, 1), # Set limits for the colorbar
#'   oob = scales::squish)) # Squish values out of the bound of limits
#' # Custom vertical and horizontal separations width of annotations
#' res <- gg2heatmap(m = mat, annot.grps = list(
#'   "Carb" = paste(mtcars$carb, "Carb"), "Gear" = paste(mtcars$gear, "Gear"),
#'   "Am" = paste(mtcars$am, "Am")),
#'   annot.pal = list(
#'     grDevices::rainbow(n = 6), c("pink", "red", "darkred"),
#'     c("blue", "grey")), annot.size = 3,
#'   theme_legend = theme(legend.justification = c(0, 0.8)),
#'   lgd.space.height = 22,
#'   annot.sep = c(0.3, 0.1)) # Horizontal space of 0.3, and vertical of 0.1
#' # Custom annotation theme
#' res <- gg2heatmap(
#'   m = mat, annot.grps = list("Carb" = mtcars$carb),
#'   annot.pal = grDevices::rainbow(n = 6), theme_annot = theme(
#'     panel.border = element_rect(
#'       color = "black", # Set annotation border color to black
#'       linewidth = 1, # Set annotation border width to 1
#'       fill = "transparent"))) # Make the annotation visible through the rect
#' # Hide the annotation and annotation legends
#' res <- gg2heatmap(m = mat, show.annot = FALSE)
#' # Merge annotations legends into a single legend
#' res <- gg2heatmap(
#'   m = mat, annot.grps = list(
#'     "Carb" = paste(mtcars$carb, "Carb"), "Gear" = paste(mtcars$gear, "Gear"),
#'     "Am" = paste(mtcars$am, "Am")),
#'   annot.pal = list(
#'     grDevices::rainbow(n = 6), c("pink", "red", "darkred"),
#'     c("blue", "grey")), annot.size = 3,
#'   lgd.merge = TRUE) # Merge annotations legends into a single legend
#' # Custom legend theme
#' res <- gg2heatmap(m = mat, theme_legend = theme(
#'   legend.background = element_rect(
#'     fill = "lightblue", # Add lightblue background
#'     color = "red"), # Add red border
#'   legend.key.size = unit(0.5, "lines"), # Decrease legend keys overall size
#'   legend.key.width = unit(1.4, "lines"), # Increase legend keys width
#'   legend.text = element_text(color = "red"))) # Color legend labels in red
#' # Increase legend space width on the left of the heatmap
#' res <- gg2heatmap(m = mat, lgd.space.width = 3)
#' # Change legend space height
#' res <- gg2heatmap(m = mat, lgd.space.height = 15)
#' # Disable automatic drawing of the heatmap (to only retrieve grobs data)
#' res <- gg2heatmap(m = mat, draw = FALSE)
#' # Print step-by-step execution of gg2heatmap (useful for debugging)
#' res <- gg2heatmap(m = mat, verbose = TRUE)

#TODO: Add option to hide subtitle information
#TODO: Add option to disable the return of grobs after execution
#TODO: Rewrite the assembling of plots using egg package functions and handling
# all unsolved remaining cases.
#TODO: Make it possible to facet without showing annotation
gg2heatmap <- function(
    m, na.handle = 'remove', dist.method = 'manhattan', rank.fun = NULL,
    top.rows = NULL, dendrograms = TRUE, dend.size = 1, theme_dend_top = NULL,
    theme_dend_left = NULL, imputation.grps = NULL, ncores = 1,
    plot.labs = NULL, row.type = 'rows', facet = NULL, split.by.rows = NULL,
    border.col = NA, border.size = 0.1, cell.size = 1, raster.filter = NULL,
    raster.size = 1080, theme_heatmap = NULL,
    guide_custom_bar = ggplot2::guide_colorbar(
        title = "Values", barwidth = 15, ticks.linewidth = 1,
        title.vjust = 0.86),
    scale_fill_grad = ggplot2::scale_fill_gradientn(
        colors = c("steelblue", "gray95", "darkorange"), na.value = "black"),
    annot.grps = list("Groups" = seq(ncol(m))),
    annot.pal = grDevices::rainbow(n = ncol(m)), annot.size = 1, annot.sep = 0,
    theme_annot = NULL, show.annot = TRUE, lgd.merge = FALSE,
    theme_legend = NULL, lgd.space.width = 1,
    lgd.space.height = 26, y.axis.right = FALSE, draw = TRUE, verbose = FALSE){
    #Fix BiocCheck() complaining about these objects initialization
    rows <- NULL
    value <- NULL
    variable <- NULL
    Cols <- NULL
    facet.annot <- NULL
    rn <- NULL
    raster.grob <- NULL
    #Check m is a matrix
    if(is.matrix(m)){ m.type <- "matrix"
    } else if(is.data.frame(m)){
        #Check if data.table, if not convert into data.table
        if(!data.table::is.data.table(m)){ m <- data.table::as.data.table(m) }
        #Get number columns of data.frame
        n.col <- ncol(m)
        #Cast matrix from molten data.frame
        colnames(m)[seq(3)] <- c("rows", "cols", "values")
        m <- data.table::dcast(
            data = m, formula = ... ~ cols, value.var = "values")
        #Add index
        m[, I := .I]
        #Store additional columns of splitting
        dt.cat <- m[, c(seq(1, n.col-2), ncol(m)), with = FALSE]
        dt.cat[, I := as.factor(I)]
        #Keep matrix with row names
        m <- m[, -c(seq(1, n.col-2)), with = FALSE]
        m <- data.table:::as.matrix.data.table(x = m, rownames = "I")
        m.type <- "molten"
    } else { stop("m must be a matrix or a molten data.frame.") }
    #Check if na.handle method  supported
    na.method <- c('keep', 'impute', 'remove')
    if(!na.handle %in% na.method){ stop("na.handle method not supported.") }
    #Check dimensions of parameters
    if(length(dist.method) == 1){
        method.rows <- dist.method
        method.cols <- dist.method
    } else if(length(dist.method) == 2){
        method.rows <- dist.method[1]
        method.cols <- dist.method[2]
    } else { stop("'dist.method' length > 2. Too many values.") }
    #Check distance methods
    methods.ls <- c('euclidean', 'maximum', 'manhattan', 'canberra',
                    'binary', 'minkowski', 'none')
    if(!method.rows %in% methods.ls){
        stop("Unknown method for distance computation on rows.")
    }
    if(!method.cols %in% methods.ls){
        stop("Unknown method for distance computation on columns.")
    }
    #Check is method.rows is not none and split.by.rows is not NULL
    if(method.rows != 'none' & !is.null(split.by.rows)){
        stop(paste("Cannot compute distances on rows if heatmap is also",
                   "splitted on rows."))
    }
    #Check if top.rows is an integer
    if(!is.null(top.rows)){
        if(is.numeric(top.rows)){
            top.rows <- as.integer(top.rows)
            if(top.rows <1){
                stop("'top.rows' must be a non-zero positive integer.") }
        } else { stop("'top.rows' must be a non-zero positive integer.") }
    }
    #Check logicals for dendrograms
    if(length(dendrograms) == 1){
        dd.rows <- dendrograms
        dd.cols <- dendrograms
    } else if(length(dendrograms) == 2){
        dd.rows <- dendrograms[1]
        dd.cols <- dendrograms[2]
    } else { stop("'dendrograms' length > 2. Too many values.") }

    #Check dendrogram sizes
    if(length(dend.size) == 1){
        dend.row.size <- dend.size
        dend.col.size <- dend.size
    } else if(length(dend.size) == 2){
        dend.row.size <- dend.size[1]
        dend.col.size <- dend.size[2]
    } else { stop("'dend.size' length > 2. Too many values.") }

    #Check plot.labs provided in input
    if(is.null(plot.labs)){ plot.labs <- ggplot2::labs(
        title = "", x = "Samples", y = "Values", legend = 'Legends')
    } else {
        if(is.null(plot.labs$title)){ plot.labs$title <- "" }
        if(is.null(plot.labs$x)){ plot.labs$x <- "Samples" }
        if(is.null(plot.labs$y)){ plot.labs$y <- "Values" }
        if(is.null(plot.labs$legend)){ plot.labs$legend <- "Legends" }
    }

    #Check annotation separation widths
    if(show.annot){
        if(length(annot.sep) == 1){
            annot.cut <- annot.sep
            annot.sep <- annot.sep
        } else if(length(annot.sep) == 2){
            annot.cut <- annot.sep[2]
            annot.sep <- annot.sep[1]
        } else { stop("'annot.sep' length > 2. Too many values.") }
    }

    #Check cell size
    if(length(cell.size) == 1){
        cell.height <- cell.size
        cell.width <- cell.size
    } else if(length(cell.size) == 2){
        cell.height <- cell.size[1]
        cell.width <- cell.size[2]
    } else { stop("'cell.size' length > 2. Too many values.") }

    #Check if y.axis.right = TRUE when axis.text.y.right or axis.title.y.right
    # or axis.ticks.y.right are not ggplot2::element_blank()
    theme_to_check <- c(
        "axis.text.y.right", "axis.title.y.right", "axis.ticks.y.right")
    thm_val <- theme_to_check[theme_to_check %in% names(theme_heatmap)]
    if(!identical(thm_val, character(0))){
        if(any(unlist(lapply(X = thm_val, FUN = function(t){
            !BiocompR:::is.elt_blank(theme_heatmap[[t]])}))) & !y.axis.right){
            warning(paste("'y.axis.right' has to be set to TRUE in order to",
                          "display the Y-axis on the right side of the",
                          "heatmap."))
        }
    }

    #Check annotations groups and palettes matching
    BiocompR:::check.annotations(
        data = m, annot.grps = annot.grps, annot.pal = annot.pal,
        verbose = verbose)

    #Handle NAs
    if(verbose){ cat("Managing missing values...") }
    m <- BiocompR::manage.na(
        data = m, method = na.handle, groups = imputation.grps, ncores = ncores)
    if(verbose){ cat("Done.\n") }

    #Apply ranking function if any function defined
    if(verbose){ cat("Ranking data by rows...") }
    if(!is.null(rank.fun)){
        # Check the structure of the rank.fun string
        rank.fun <- check_fun(fun = rank.fun, param.name = "rank.fun")
        m <- eval(parse(text = paste0(
            "m[order(apply(m, 1, ", rank.fun,
            ", na.rm = TRUE), decreasing = TRUE), , drop = FALSE]")))
    }
    if(verbose){ cat("Done.\n") }

    #Subset top rows if any value defined
    if(!is.null(top.rows)){ m <- utils::head(x = m, n = top.rows) }

    #Remove NAs if some for dendrogram matrix
    if(method.rows != 'none' | method.cols != 'none'){
        if(verbose){ cat("Clustering data...") }
        m_rname <- rownames(m)[apply(X = m, MARGIN = 1, FUN = anyNA) == FALSE]
        dend_mat <- m[stats::complete.cases(m), ]
        if(is.vector(dend_mat)){
            dend_mat <- matrix(
                dend_mat, nrow = 1, dimnames = list(m_rname, names(dend_mat)))
        }
        if(nrow(dend_mat) != nrow(m)){
            warning(paste(
                "Distance method selected need complete data.",
                nrow(m) - nrow(dend_mat), "incomplete rows removed out of",
                nrow(m), "rows selected."))
        }
        #Check how many rows dend_mat has
        if(nrow(dend_mat) == 0){
            stop(paste("Cannot compute distances on rows. All rows are missing",
                       "values."))
        }
    }
    #Set theme_dend_top plot.margin based on them_heatmap plot.margin
    if(!is.null(theme_heatmap)){
        if(!is.null(theme_heatmap$plot.margin)){
            if(is.null(theme_dend_top)){
                theme_dend_top <- ggplot2::theme(
                    plot.margin = theme_heatmap$plot.margin)
            } else {
                theme_dend_top <- theme_dend_top + ggplot2::theme(
                    plot.margin = theme_heatmap$plot.margin)
            }
        }
    }
    #Compute rows distances & create rows dendrogram
    if(method.rows != 'none'){
        row_dist <- tryCatch(
            parallelDist::parDist(
                dend_mat, method = method.rows, threads = ncores),
            error = function(cond){ if(grepl(
                pattern = "impossible d'allouer un vecteur de taille", x = cond)
                | grepl(pattern = "cannot allocate vector of size", x = cond)){
                cond$message <- paste(
                    "Cannot compute distances on rows.",
                    "Too many rows containing too many values.")
                stop(cond) } else { stop(cond) }
            }, warning = function(cond){ warning(cond$message) }, finally = {})
        row_hclust <- fastcluster::hclust(row_dist)
        rm(row_dist)
        rowclust <- stats::as.dendrogram(row_hclust)
        row.order <- stats::order.dendrogram(rowclust)
        if(dd.rows){
            #Get dendrogram segments and order matrix rows
            ddgr_seg_row <- BiocompR::ggdend(
                df = ggdendro::dendro_data(rowclust)$segments,
                orientation = "left", reverse.x = TRUE,
                theme_dend = theme_dend_left)
        }
    } else if(dd.rows & method.rows == 'none'){
        stop(paste(
            "Cannot plot dendrogram on rows with dist.method = 'none' for",
            "rows. To avoid this error message, set 'dendrograms' to FALSE,",
            "or choose another method for dist.method."))
    }
    #Compute columns distances & create columns dendrogram
    if(method.cols != 'none'){
        ddgr <- stats::as.dendrogram(fastcluster::hclust(parallelDist::parDist(
            t(dend_mat), method = method.cols, threads = ncores)))
        column.order <- stats::order.dendrogram(ddgr)
        if(dd.cols){
            #Get dendrogram data
            ddgr_dat <- ggdendro::dendro_data(ddgr)
            #Get dendrogram segments and order matrix columns
            ddgr_seg_col <- BiocompR::ggdend(
                df = ddgr_dat$segments, orientation = "top",
                theme_dend = theme_dend_top)
        }
    } else if(dd.cols & method.cols == 'none'){
        stop("Cannot plot dendrogram on columns with method.cols = 'none'.")
    }
    #Reorder rows and columns matrix following dendrograms orders
    if(method.rows != 'none' & method.cols != 'none'){
        # Keep rows selected for the method applied on rows
        m <- m[row_hclust$labels, ]
        # All dendrograms on and all methods specified
        dframe <- m[row.order, column.order, drop = FALSE]
    } else if(method.rows != 'none' & is.null(rank.fun) &
              method.cols == 'none'){
        # Keep rows selected for the method applied on rows
        m <- m[row_hclust$labels, ]
        # row.dendrogram on, col.dendrogram off, method.row specified
        dframe <- m[row.order, , drop = FALSE]
    } else if(method.rows == 'none' & method.cols != 'none'){
        # row.dendrogram off, col.dendrogram on, method.col specified
        dframe <- m[, column.order, drop = FALSE]
    } else { #Leave matrix unchanged
        dframe <- m
    }
    rm(m)
    if(method.rows != 'none' | method.cols != 'none'){
        if(verbose){ cat("Done.\n") }
    }

    #Melt Matrix into a data.table
    if(verbose){ cat("Melting matrix...") }
    dt.frame <- data.table::as.data.table(x = dframe, keep.rownames = TRUE)
    dt.frame[, rn := factor(x = rn, levels = rev(rn))]
    melted_mat <- data.table::melt.data.table(
        data = dt.frame, id.vars = "rn",
        measure.vars = colnames(dt.frame)[-c(1)])
    if(m.type == "molten"){
        #Merge original row names and additional categories stored
        melted_mat <- data.table::merge.data.table(
            x = melted_mat, y = dt.cat, by.x = "rn", by.y = "I", all.x = TRUE)
        #replace rn values by rows values and remove rows column
        melted_mat[, rn := rows]
        melted_mat <- melted_mat[, -c("rows"), ]
        #Apply split.by.rows option
        if(!is.null(split.by.rows)){
            #Keep melted_mat with the category specified in split.by.rows
            melted_mat <- melted_mat[, c(
                "rn", "variable", "value", split.by.rows), with = FALSE]
            #Rename splitting categorical column
            data.table::setnames(
                x = melted_mat, old = split.by.rows, new = "split.by.rows")
        }
    }
    if(!is.null(facet)){
        if(method.cols == "none"){
            dt.facet <- data.table::data.table(
                "variable" = colnames(dt.frame)[-c(1)],
                "facet.annot" = annot.grps[[facet]])
            melted_mat <- merge(x = melted_mat, y = dt.facet, by = "variable",
                                all.x = TRUE, sort = FALSE)
            # If no factor defined, define some
            if(!is.factor(melted_mat$facet.annot)){
                melted_mat[, facet.annot := as.factor(facet.annot)]
                melted_mat[, facet.annot := factor(
                    x = facet.annot, levels = sort(levels(facet.annot)))]
            }
            melted_mat[, variable := as.factor(variable)]
            melted_mat[, variable := factor(x = variable, levels = unique(
                melted_mat, by = "variable")[order(facet.annot)]$variable)]
            #Get new column order based on the facetting
            column.order <- match(
                levels(melted_mat$variable), colnames(dt.frame)[-c(1)])
        } else {
            stop(paste("Cannot facet the heatmap if clustering method is",
                       "already applied on columns."))
        }
    }
    rm(dt.frame)
    if(verbose){ cat("Done.\n") }

    if(verbose){ cat("Configure heatmap...") }
    #Set theme_forced_htmp
    theme_forced_htmp <- ggplot2::theme(legend.position = "bottom")
    #Set theme_default_htmp
    theme_default_htmp <- ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 12, color = 'black'),
        axis.text.x = ggplot2::element_text(
            size = 11, angle = -45, hjust = 0, vjust = 0.5, face = 'bold'),
        axis.ticks.x = ggplot2::element_line(color = 'black'),
        axis.title.y.right = ggplot2::element_blank(),
        axis.text.y.right = ggplot2::element_blank(),
        axis.ticks.y.right = ggplot2::element_blank(),
        axis.title.y.left = ggplot2::element_blank(),
        axis.text.y.left = ggplot2::element_blank(),
        axis.ticks.y.left = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "transparent"),
        plot.background = ggplot2::element_rect(fill = "transparent"),
        plot.margin = ggplot2::margin(0, 0, 0, 0),
        legend.text = ggplot2::element_text(size = 11),
        legend.title = ggplot2::element_text(size = 12),
        legend.justification = c(0.4, 0.5),
        strip.text.y = ggplot2::element_text(size = 12),
        strip.background.y = ggplot2::element_blank()
    )
    #Update theme_heatmap
    if(is.null(theme_heatmap)){
        theme_heatmap <- theme_default_htmp + theme_forced_htmp
    } else {
        theme_heatmap <- theme_default_htmp + theme_heatmap + theme_forced_htmp
    }
    #Update guide_colorbar and force order parameter
    guide_custom_bar <- BiocompR:::update_guide_colorbar(
        new_guide = guide_custom_bar, forced_param = list("order" = 1))
    #Set heatmap source parameters
    htmp.source <- ggplot2::ggplot() +
        scale_fill_grad +
        ggplot2::scale_color_manual(values = NA) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::guides(fill = guide_custom_bar) +
        ggplot2::guides(color = ggplot2::guide_legend(
            "NA", override.aes = list(fill = scale_fill_grad$na.value),
            title.vjust = 0.5, order = 2)) +
        ggplot2::labs(x = plot.labs$x, y = plot.labs$y)

    #If facetting is on
    if(!is.null(facet)){
        if(show.annot){
            htmp <- htmp.source +
                ggplot2::facet_grid(
                    . ~ facet.annot, scales = "free", space = "free") +
                ggplot2::theme(panel.spacing = ggplot2::unit(0, "lines"),
                               strip.background = ggplot2::element_blank(),
                               strip.text = ggplot2::element_blank())
        } else {
            htmp <- htmp.source +
                ggplot2::facet_grid(
                    . ~ facet.annot, scales = "free", space = "free") +
                ggplot2::theme(
                    panel.spacing = ggplot2::unit(0, "lines"),
                    strip.background = ggplot2::element_rect(
                        color = "black", linewidth = 0.5))
        }
    } else {
        #If m is a molten data.frame and split.by.rows is on
        if(m.type == "molten" & !is.null(split.by.rows)){
            if(!y.axis.right){
                htmp <- htmp.source + ggplot2::facet_grid(
                    split.by.rows ~ ., scales = "free", space = "free")
            } else {
                stop(paste("Cannot split heatmap on rows when 'y.axis.right'",
                           "TRUE. Please contact the developper."))
            }
        } else { htmp <- htmp.source }
    }

    #Draw heatmap with geom_tile
    if(nrow(melted_mat[is.na(value)]) != 0){
        #Display legend of missing values if any
        htmp <- htmp + ggplot2::geom_tile(data = melted_mat, ggplot2::aes(
            x = variable, y = rn, fill = value, color = " "),
            color = border.col, linewidth = border.size, width = cell.width,
            height = cell.height)
    } else {
        htmp <- htmp + ggplot2::geom_tile(
            data = melted_mat, ggplot2::aes(x = variable, y = rn, fill = value),
            color = border.col, linewidth = border.size,
            width = cell.width, height = cell.height)
    }

    if(dd.rows){
        theme_heatmap <- theme_heatmap +
            ggplot2::theme(axis.title.y.left = ggplot2::element_blank(),
                           axis.ticks.y.left = ggplot2::element_blank(),
                           axis.ticks.length.y.left = ggplot2::unit(0, "pt"),
                           axis.text.y.left = ggplot2::element_blank())
    }

    if(y.axis.right){
        htmp <- htmp +
            ggplot2::scale_y_discrete(position = 'right', expand = c(0, 0))
    } else { htmp <- htmp + ggplot2::scale_y_discrete(expand = c(0, 0)) }

    if(verbose){ cat("Done.\n") }

    if(show.annot){
        if(verbose){ cat("Configure annotations...") }
        #Reorder groups and convert as factors
        annot.grps <- lapply(X = annot.grps, FUN = function(i){
            factor(x = i, levels = unique(i))})
        #If the distance method used on columns is not "none" reorder columns
        if(method.cols != "none" | !is.null(facet)){
            annot.grps <- lapply(
                X = annot.grps, FUN = function(i){ i[column.order] })
        }
        #Get ordered sample names
        if(method.cols != "none"){ sample.names <- colnames(dframe)
        } else { sample.names <- colnames(dframe) }
        #Set theme_forced_annot
        theme_forced_annot <- ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.ticks.length.x.top = ggplot2::unit(0, "pt")
        )
        #Set theme_default_annot
        theme_default_annot <- ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 12, color = "black"),
            axis.ticks.y = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 2, 0),
            panel.background = ggplot2::element_rect(fill = "white")
        )
        #Update theme_heatmap
        if(is.null(theme_annot)){
            theme_annot <- theme_default_annot + theme_forced_annot
        } else{
            theme_annot <- theme_default_annot + theme_annot +
                theme_forced_annot
        }
        #Set number of columns to display annotations legends
        # Builds color tables for legends.
        col_table <- BiocompR:::build.col_table(
            annot.grps = annot.grps, annot.pal = annot.pal)
        if(lgd.merge){
            #Rbind list color tables because lgd.merge is TRUE
            col_table <- data.table::rbindlist(col_table, idcol = TRUE)
            if(!(is.list(annot.pal)) | length(annot.pal) == 1){
                #Remove duplicated colors
                col_table <- col_table[!duplicated(x = Cols)]
            }
            col_table <- list(col_table[, c("Grps", "Cols"), ])
        }
        #Build legends layout
        lgd.layout <- BiocompR::build.layout(
            col_table = col_table, height.lgds.space = lgd.space.height)
        #Calculate legend length
        lgd_sizes <- BiocompR:::get.len.legends(col_table = col_table)
        #Calculate legend columns
        lgd.ncol <- ceiling(lgd_sizes/lgd.space.height)

        #Set theme_default_legend
        theme_default_legend <- ggplot2::theme(
            legend.title = ggplot2::element_text(size = 12),
            legend.text = ggplot2::element_text(size = 11)
        )
        #Update theme_legend
        if(is.null(theme_legend)){
            theme_legend <- theme_default_legend
        } else { theme_legend <- theme_default_legend + theme_legend }
        #Create Color Sidebar
        col_sidebar <- BiocompR::plot.col.sidebar(
            sample.names = sample.names, annot.grps = annot.grps,
            annot.pal = annot.pal, annot.pos = 'top', annot.sep = annot.sep,
            annot.cut = annot.cut, merge.lgd = lgd.merge, right = TRUE,
            lgd.name = plot.labs$legend, lgd.ncol = lgd.ncol,
            theme_legend = theme_legend, theme_annot = theme_annot,
            set.x.title = NULL, set.y.title = NULL, dendro.pos = 'top',
            facet = facet)
        rm(sample.names)
        if(verbose){ cat("Done.\n") }
    }

    #Extract Legend
    if(verbose){ cat("Extracting legends...") }
    #Create a subset plot of the original heatmap
    xtrm.melted_mat <- melted_mat[
        value %in% c(min(value, na.rm = TRUE), max(value, na.rm = TRUE))]
    if(nrow(melted_mat[is.na(value)]) != 0){
        subplot.htmp <- htmp.source +
            ggplot2::geom_tile(
                data = xtrm.melted_mat, ggplot2::aes(
                    x = variable, y = rn, fill = value, color = " "))
    } else {
        subplot.htmp <- htmp.source +
            ggplot2::geom_tile(
                data = xtrm.melted_mat,
                mapping = ggplot2::aes(x = variable, y = rn, fill = value))
    }
    htmp_legend <- BiocompR:::get.lgd(gg2.obj = subplot.htmp + theme_heatmap)
    rm(subplot.htmp)

    # #Create ggplot2 heatmap
    # htmp <- htmp + theme_heatmap + ggplot2::theme(legend.position = "none")
    # #Convert ggplots into grobs
    # if(is.null(facet)){
    #     if(dd.cols){
    #         ddgr_seg_col <- R.devices::suppressGraphics(
    #             ggplot2::ggplotGrob(ddgr_seg_col))
    #     }
    #     if(dd.rows){
    #         ddgr_seg_row <- R.devices::suppressGraphics(
    #             ggplot2::ggplotGrob(ddgr_seg_row))
    #     }
    #     if(show.annot){
    #         sidebar_legend <- col_sidebar$legends
    #         col_sidebar_grob <- R.devices::suppressGraphics(
    #             ggplot2::ggplotGrob(col_sidebar$sidebar))
    #         rm(col_sidebar)
    #     }
    # } else { # Use the egg package
    #     if(!dd.rows & !dd.cols & show.annot){
    #         facet_size <- 0
    #         #Create main grob
    #         main_grob <- R.devices::suppressGraphics(egg::ggarrange(
    #             col_sidebar$sidebar, htmp, ncol = 1, nrow = 2,
    #             heights = c(annot.size + facet_size, 30),
    #             draw = FALSE))
    #         #Get sidebar legends
    #         sidebar_legend <- col_sidebar$legends
    #     } else {
    #         if(dd.cols){
    #             stop("Cannot facet and apply clustering on the data.")
    #         } else if(dd.rows){
    #             #TODO: Handle this case!
    #             stop(paste("Cannot plot heatmap with row dendrogram and",
    #                        "annotation yet. Please contact the developper."))
    #         }
    #     }
    #     if(show.annot){ rm(col_sidebar) }
    # }
    if(verbose){ cat("Done.\n") }

    #Heatmap rasterization
    if(!is.null(raster.filter)){
        if(verbose){ cat("Rasterizing...\n") }
        if(raster.filter %in% magick::filter_types()){
            #Create "empty" theme
            theme_empty <- ggplot2::theme(
                plot.margin = ggplot2::margin(0, 0, 0, 0),
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "transparent"),
                plot.background = ggplot2::element_rect(fill = "transparent"),
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.title.y.right = ggplot2::element_blank(),
                axis.ticks.y.right = ggplot2::element_blank(),
                axis.text.y.right = ggplot2::element_blank(),
                axis.title.y.left = ggplot2::element_blank(),
                axis.ticks.y.left = ggplot2::element_blank(),
                axis.text.y.left = ggplot2::element_blank(),
                axis.ticks.length = ggplot2::unit(0, "pt"),
                legend.position = "none")

            #If facet is used
            if(!is.null(facet)){
                if(verbose){ cat("Facet rasterization:\n") }
                ls.rasters <- lapply(
                    X = levels(melted_mat$facet.annot), FUN = function(i){
                        if(verbose){ cat(paste("\t", i, "\n")) }
                        #Create sub DT
                        sub.melted <- melted_mat[facet.annot == i]
                        #Create sub heatmap and remove all customization
                        sub.htmp <- ggplot2::ggplot() +
                            ggplot2::geom_tile(data = sub.melted, ggplot2::aes(
                                x = variable, y = rn, fill = value),
                                color = border.col, linewidth = border.size,
                                width = cell.width, height = cell.height) +
                            scale_fill_grad +
                            ggplot2::scale_color_manual(values = NA) +
                            ggplot2::scale_x_discrete(expand = c(0, 0)) +
                            ggplot2::scale_y_discrete(expand = c(0, 0)) +
                            ggplot2::guides(fill = guide_custom_bar) +
                            ggplot2::guides(color = ggplot2::guide_legend(
                                "NA", override.aes = list(
                                    fill = scale_fill_grad$na.value),
                                title.vjust = 0.5, order = 2)) +
                            ggplot2::labs(x = plot.labs$x, y = plot.labs$y) +
                            theme_empty +
                            ggplot2::theme(legend.position = "none")
                        rm(sub.melted)
                        #Rasterize ggplot into a grob
                        raster.grob <- BiocompR::raster.gg2grob(
                            gg.plot = sub.htmp, filter = raster.filter,
                            size = raster.size)
                        # raster.grob <- BiocompR::raster.ggplot.to.grob(
                        #     gg.plot = sub.htmp, filter = raster)
                        rm(sub.htmp)
                        #Make grob annotation
                        raster.annot <- BiocompR:::annotation_custom2(
                            grob = raster.grob, xmin = -Inf, xmax = Inf,
                            ymin = -Inf, ymax = Inf,
                            data = melted_mat[facet.annot == i])
                        rm(raster.grob)
                        return(raster.annot)
                    })
                #Fit the list of raster grobs into a ggplot
                htmp <- ggplot2::ggplot(
                    data = melted_mat, ggplot2::aes(
                        x = variable, y = rn, fill = value)) +
                    ggplot2::geom_blank() + theme_heatmap +
                    ggplot2::theme(legend.position = "none") +
                    ggplot2::labs(x = plot.labs$x, y = plot.labs$y) +
                    ggplot2::facet_grid(. ~ facet.annot, scales = "free",
                                        space = "free") +
                    ggplot2::theme(
                        panel.spacing = ggplot2::unit(0, "lines"),
                        strip.background = ggplot2::element_blank(),
                        strip.text = ggplot2::element_blank()) +
                    ggplot2::scale_x_discrete(
                        expand = ggplot2::expansion(add = 0.5)) +
                    ls.rasters

            } else {
                #Remove all customization
                htmp <- htmp + theme_empty
                #Catch heatmap in magick::image_graph()
                raster.grob <- BiocompR::raster.gg2grob(
                    gg.plot = htmp, filter = raster.filter, size = raster.size)
                # raster.grob <- BiocompR::raster.ggplot.to.grob(
                #     gg.plot = htmp, filter = raster)
                #Make grob annotation
                raster.annot <- ggplot2::annotation_custom(
                    raster.grob, -Inf, Inf, -Inf, Inf)
                #Fit the raster grob into a ggplot
                htmp <- ggplot2::ggplot(
                    data = melted_mat, ggplot2::aes(
                        x = variable, y = rn, fill = value)) +
                    ggplot2::geom_blank() + raster.annot + theme_heatmap +
                    ggplot2::theme(legend.position = "none") +
                    ggplot2::labs(x = plot.labs$x, y = plot.labs$y) +
                    ggplot2::scale_x_discrete(
                        expand = ggplot2::expansion(add = 0.5))
            }
        } else { stop("Rasterization filter not supported.") }
        if(verbose){ cat("Done.\n") }
    }

    #Create ggplot2 heatmap
    htmp <- htmp + theme_heatmap + ggplot2::theme(legend.position = "none")
    #Convert ggplots into grobs
    if(is.null(facet)){
        if(dd.cols){
            ddgr_seg_col <- R.devices::suppressGraphics(
                ggplot2::ggplotGrob(ddgr_seg_col))
        }
        if(dd.rows){
            ddgr_seg_row <- R.devices::suppressGraphics(
                ggplot2::ggplotGrob(ddgr_seg_row))
        }
        if(show.annot){
            sidebar_legend <- col_sidebar$legends
            col_sidebar_grob <- R.devices::suppressGraphics(
                ggplot2::ggplotGrob(col_sidebar$sidebar))
            rm(col_sidebar)
        }
    } else { # Use the egg package
        if(!dd.rows & !dd.cols & show.annot){
            facet_size <- 0
            #Create main grob
            main_grob <- R.devices::suppressGraphics(egg::ggarrange(
                col_sidebar$sidebar, htmp, ncol = 1, nrow = 2,
                heights = c(annot.size + facet_size, 30),
                draw = FALSE))
            #Get sidebar legends
            sidebar_legend <- col_sidebar$legends
            # Convert sidebar to a grob
            col_sidebar_grob <- R.devices::suppressGraphics(
                ggplot2::ggplotGrob(col_sidebar$sidebar))
        } else {
            if(dd.cols){
                stop("Cannot facet and apply clustering on the data.")
            } else if(dd.rows){
                #TODO: Handle this case!
                stop(paste("Cannot plot heatmap with row dendrogram and",
                           "annotation yet. Please contact the developper."))
            }
        }
        if(show.annot){ rm(col_sidebar) }
    }

    #Remove melted_mat
    rm(melted_mat)
    if(verbose){ cat("Converting ggplot into grid object...") }
    htmp <- R.devices::suppressGraphics(ggplot2::ggplotGrob(x = htmp))
    if(verbose){ cat("Done.\n") }
    #Resize grobs widths
    if(verbose){ cat("Redimensioning grobs...") }
    if(dd.cols & dd.rows){
        if(show.annot){
            if(is.null(facet)){
                ls.w.grobs <- list(
                    'dd_col' = ddgr_seg_col, 'sidebar' = col_sidebar_grob,
                    'htmp' = htmp)
            }
        } else { ls.w.grobs <- list('dd_col' = ddgr_seg_col, 'htmp' = htmp) }
        upd.grob_w <- BiocompR::resize.grobs(
            ls.grobs = ls.w.grobs, dimensions = "widths", start.unit = 4,
            end.unit = 7)
        rm(ddgr_seg_col)
    } else if(dd.cols & !dd.rows){
        if(show.annot){
            if(is.null(facet)){
                ls.w.grobs <- list(
                    'dd_col' = ddgr_seg_col, 'sidebar' = col_sidebar_grob,
                    'htmp' = htmp)
            }
        } else { ls.w.grobs <- list('dd_col' = ddgr_seg_col, 'htmp' = htmp) }
        rm(ddgr_seg_col)

        if(is.null(facet)){
            upd.grob_w <- BiocompR::resize.grobs(
                ls.grobs = ls.w.grobs, dimensions = "widths", start.unit = 3,
                end.unit = 7)
        }
        # } else {
        #     #TODO: Handle this case!
        #     stop(paste("Cannot plot heatmap with column dendrogram with facet",
        #                "and without annotation yet. Please contact the",
        #                "developper."))
        # }
    } else {
        if(show.annot){
            if(is.null(facet)){
                ls.w.grobs <- list('sidebar' = col_sidebar_grob, 'htmp' = htmp)
            } else {
                htmp_grob <- htmp
            }
        } else { ls.w.grobs <- list('htmp' = htmp) }

        if(is.null(facet)){
            upd.grob_w <- BiocompR::resize.grobs(
                ls.grobs = ls.w.grobs, dimensions = "widths", start.unit = 3,
                end.unit = 7)
        }
        # } else {
        #     upd.grob_w <- BiocompR::resize.grobs(
        #         ls.grobs = ls.w.grobs, dimensions = "widths", start.unit = 1,
        #         end.unit = max(unlist(lapply(X = ls.w.grobs, FUN = function(i){
        #             length(i[["widths"]]) }))))
        # }
    }

    rm(htmp)
    if(show.annot){ if(!(!dd.rows & show.annot)){ rm(col_sidebar_grob) } }
    #Resize grobs heights
    if(dd.rows){
        if(is.null(facet)){
            ls.h.grobs <- list(
                'dd_row' = ddgr_seg_row, 'htmp' = upd.grob_w$htmp)
            upd.grob_h <- BiocompR::resize.grobs(
                ls.grobs = ls.h.grobs, dimensions = 'heights', start.unit = 7,
                end.unit = 9)
            rm(ddgr_seg_row)
        } else {
            #TODO: Handle this case!
            stop(paste("Cannot plot heatmap with row dendrogram without",
                       "annotation yet. Please contact the developper."))
        }
    } else {
        if(!(!is.null(facet) & !dd.rows & show.annot)){
            #TODO: Cannot find upd.grob_w
            upd.grob_h <- list("htmp" = upd.grob_w$htmp)
        }
    }
    if(verbose){ cat("Done.\n") }

    #Create the Right Panel for legends
    if(show.annot){
        if(verbose){ cat("Setting legends layout...") }
        #Update layout with empty grob
        if(anyNA(lgd.layout)){
            lgd.layout[is.na(lgd.layout)] <- max(lgd.layout, na.rm = TRUE) + 1
            #Add an empty grob in the legend to stack them to the top
            sidebar_legend <- c(sidebar_legend, list(grid::textGrob("")))
        }
        right.legends <- gridExtra::arrangeGrob(
            grobs = sidebar_legend, layout_matrix = lgd.layout)
        rm(sidebar_legend)

        if(verbose){ cat("Done.\n") }
    }

    #Combine Dendrogram with Color Sidebar and Heatmap
    if(verbose){ cat("Creating final plot...") }
    if(dd.rows & dd.cols){
        if(show.annot){
            #Create main grob
            main_grob <- gridExtra::arrangeGrob(
                grobs = list(grid::textGrob(""), upd.grob_w$dd_col,
                             grid::textGrob(""), upd.grob_w$sidebar,
                             upd.grob_h$dd_row, upd.grob_h$htmp),
                ncol = 2, nrow = 3, heights = c(
                    dend.col.size + 2, annot.size, 30),
                widths = c(dend.row.size + 1, 10))
            #Set default legend width space
            def.lgd.width <- 2
        } else {
            #Create main grob
            main_grob <- gridExtra::arrangeGrob(
                grobs = list(grid::textGrob(""), upd.grob_w$dd_col,
                             upd.grob_h$dd_row, upd.grob_h$htmp),
                ncol = 2, nrow = 2, heights = c(dend.col.size + 2, 30),
                widths = c(dend.row.size + 1, 10))
        }
    } else if(!dd.rows & !dd.cols){
        if(show.annot){
            #Create main grob
            if(is.null(facet)){
                facet_size <- 0
                main_grob <- gridExtra::arrangeGrob(grobs = list(
                    upd.grob_w$sidebar, upd.grob_h$htmp), ncol = 1, nrow = 2,
                    heights = c(annot.size + facet_size, 30), widths = 10)
            }
            # if(is.null(facet)){ facet_size <- 0 } else { facet_size <- 1 }
            # main_grob <- gridExtra::arrangeGrob(grobs = list(
            #     upd.grob_w$sidebar, upd.grob_h$htmp), ncol = 1, nrow = 2,
            #     heights = c(annot.size + facet_size, 30), widths = 10)
            #Set default legend width space
            def.lgd.width <- 1
        } else {
            #Create main grob
            main_grob <- gridExtra::arrangeGrob(
                grobs = list(upd.grob_h$htmp), ncol = 1, nrow = 1, heights = 30,
                widths = 10)
        }
    } else if(dd.rows & !dd.cols){
        if(show.annot){
            #Create main grob
            main_grob <- gridExtra::arrangeGrob(grobs = list(
                grid::textGrob(""), upd.grob_w$sidebar, upd.grob_h$dd_row,
                upd.grob_h$htmp), ncol = 2, nrow = 2,
                heights = c(annot.size, 30), widths = c(dend.row.size + 1, 10))
            #Set default legend width space
            def.lgd.width <- 2
        } else {
            #Create main grob
            main_grob <- gridExtra::arrangeGrob(grobs = list(
                upd.grob_h$dd_row, upd.grob_h$htmp), ncol = 2, nrow = 1,
                heights = 30, widths = c(dend.row.size + 1, 10))
        }
    } else if(!dd.rows & dd.cols){
        if(show.annot){
            if(is.null(facet)){
                #Create main grob
                main_grob <- gridExtra::arrangeGrob(grobs = list(
                    upd.grob_w$dd_col, upd.grob_w$sidebar, upd.grob_h$htmp),
                    ncol = 1, nrow = 3, heights = c(
                        dend.col.size + 2, annot.size, 30), widths = 10)
            }
            #Set default legend width space
            def.lgd.width <- 1
        } else {
            #Create main grob
            main_grob <- gridExtra::arrangeGrob(grobs = list(
                upd.grob_w$dd_col, upd.grob_h$htmp), ncol = 1, nrow = 2,
                heights = c(dend.col.size + 2, 30), widths = 10)
        }
    }
    #Final plot
    if(show.annot){
        arranged.grob <- gridExtra::arrangeGrob(
            grobs = list(main_grob, right.legends), ncol = 2,
            widths = c(20, def.lgd.width + lgd.space.width))
    } else { arranged.grob <- main_grob }

    final.plot <- gridExtra::arrangeGrob(gridExtra::arrangeGrob(
        top = grid::textGrob(
            plot.labs$title, gp = grid::gpar(fontsize = 15, font = 1)),
        grobs = list(grid::textGrob(paste0(
            "Columns ordered by ", method.cols, " distance; Rows ordered by ",
            method.rows, " distance; ", nrow(dframe), " ", row.type, "."),
            gp = grid::gpar(fontsize = 12, fontface = 3L)), arranged.grob,
            htmp_legend), nrow = 3, heights = c(3, 50, 6)))
    rm(main_grob)
    #Plot final heatmap
    if(draw){
        grid::grid.newpage()
        grid::grid.draw(final.plot)
    } else { if(verbose){
        message(paste("To draw heatmap use grid::grid.draw() on",
                      "your_object$result.grob"))
    }}
    #Prepare results
    if(show.annot){
        if(!dd.rows & !is.null(facet)){
            ls.res <- list(
                "result.grob" = final.plot, "heatmap.grob" = htmp_grob,
                "heatmap.lgd.grob" = htmp_legend,
                "sidebar.grob" = col_sidebar_grob,
                "sidebar.lgds.grob" = right.legends)
        } else {
            ls.res <- list(
                "result.grob" = final.plot, "heatmap.grob" = upd.grob_h$htmp,
                "heatmap.lgd.grob" = htmp_legend,
                "sidebar.grob" = upd.grob_w$sidebar,
                "sidebar.lgds.grob" = right.legends)
        }
        rm(right.legends)
    } else {
        ls.res <- list(
            "result.grob" = final.plot, "heatmap.grob" = upd.grob_h$htmp,
            "heatmap.lgd.grob" = htmp_legend)
    }
    rm(final.plot, htmp_legend)
    if(!(!dd.rows & !is.null(facet))){
        if(dd.cols){ ls.res[["cols.dendrogram.grob"]] <- upd.grob_w$dd_col }
        if(dd.rows){ ls.res[["rows.dendrogram.grob"]] <- upd.grob_h$dd_row }
        rm(upd.grob_w, upd.grob_h)
    }
    if(verbose){ cat("Done.\n") }
    #Return a list of grobs with final plot and separate grobs.
    return(ls.res)
}
