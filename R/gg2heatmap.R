#' Creates a custom heatmap with dendrograms and annotations.
#'
#' @param m               A \code{matrix}.
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
#'                        matrix rows to order them on the final heatmap. If
#'                        rank.fun = NULL the order of rows in the matrix will
#'                        be kept (Default: rank.fun = NULL;
#'                        Supported: rank.fun = 'sd').
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
#' @param raster          A \code{character} to be used as a filter for matrix
#'                        rasterization. The list of the supported rasterization
#'                        filters is available in magick::filter_types(). If no
#'                        value is specified for 'raster' the original ggplot2
#'                        heatmap is displayed (Default: raster = NULL).
#'                        Warning: Be aware that rasterization may take several
#'                        more minutes than the usual process time.
#' @param facet           A \code{character} matching an annotation name in
#'                        'annot.grps' to be used to split heatmap in separate
#'                        panels following the annotation.
#' @param plot.title      A \code{character} to be used as title for the plot.
#' @param row.type        A \code{character} to be used in the plot subtitle
#'                        description as a definition of the rows (e.g. 'loci',
#'                        'samples', 'regions', etc.
#'                        Default: row.type = 'rows').
#' @param imputation.grps A \code{character} vector defining groups to which
#'                        columns of the matrix belong. These groups are use to
#'                        make group-wise imputation of missing values between
#'                        columns. The vector has to be of the same length than
#'                        the number of columns in the matrix
#'                        (Default: imputation.grps = NULL).
#' @param ncores          An \code{integer} to specify the number of
#'                        cores/threads to be used to parallel-compute distances
#'                        for dendrograms.
#' @param heatmap.pal     A \code{character} vector to be used as a color
#'                        palette for the heatmap.
#' @param na.col          A \code{character} to be used as a color for missing
#'                        values on the heatmap (Default: na.col = "black").
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
#' @param verbose         A \code{logical} to display information about the
#'                        step-by-step processing of the data if TRUE
#'                        (Default: verbose = FALSE).
#' @param lgd.scale.name  A \code{character} to be used as legend title for the
#'                        heatmap scale (Default: lgd.scale.name = 'values').
#' @param lgd.bars.name   A \code{character} specifying the name of annotation
#'                        side bar legends, when legends are merged
#'                        (lgd.merge = TRUE)
#'                        (Default: lgd.bars.name = 'Legends').
#' @param lgd.title       An \code{element_text} object to setup legend titles
#'                        (Default: lgd.title = ggplot2::element_blank()).
#' @param lgd.text        An \code{element_text} object to setup legend labels
#'                        (Default: lgd.text = ggplot2::element_blank()).
#' @param lgd.merge       A \code{logical} specifying whether the legends of
#'                        multiple annotation bars should be merged
#'                        (lgd.merge = TRUE) or remain separated
#'                        (lgd.merge = FALSE). lgd.merge is especially useful
#'                        when you want to map the same color palette to
#'                        multiple annotations sharing the same values
#'                        (Default: lgd.merge = FALSE).
#' @param lgd.space.width A \code{numeric} specifying the width of the legend
#'                        space (Default: lgd.space.width = 1).
#' @param y.lab           A \code{character} to specify Y-axis title
#'                        (Default: y.lab = 'Values').
#' @param x.lab           A \code{character} to specify X-axis title
#'                        (Default: x.lab = 'Samples').
#' @param theme_heatmap   A ggplot2 \code{theme} to specify any theme parameter
#'                        you wish to custom on the heatmap part of the plot
#'                        (Default: theme_heatmap = NULL).
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
#'                               cell width, and the second
#'                               numeric will apply to cell height.}
#'                        }
#' @param y.axis.right    A \code{logical} to specify whether the heatmap Y axis
#'                        should be displayed on the right (y.axis.right = TRUE)
#'                        or not (y.axis.right = FALSE)
#'                        (Default: y.axis.right = FALSE).
#' @return A \code{grob} list containing the final plot, and also each grob
#'         generated separately.
#' @author Yoann Pageaud.
#' @export

#TODO: Try to improve display for the annotation when a lot of columns used (try geom_raster()).
#TODO: Create a guide argument using guide().
#TODO: Add option to not display any annotation bar.
#TODO: Remove auto-display of heatmap once it is made
#TODO: Add option to support molten data.frame
#TODO: Add option to split heatmap by rows in n smaller heatmaps based on a categorical variable from the molten data.frame
#TODO: Create an annotation.theme argument using the ggplot2::theme() function.
#TODO: Fix the issue with annotation legends overlapping.
#TODO: merge y.lab and x.lab using the function ggplot2::labs()
gg2heatmap <- function(
  m, na.handle = 'remove', dist.method = 'manhattan', rank.fun = NULL,
  top.rows = NULL, dendrograms = TRUE, dend.size = 1, raster = NULL,
  facet = NULL, plot.title = "", row.type = 'rows', imputation.grps = NULL,
  ncores = 1, heatmap.pal = c("steelblue", "gray95", "darkorange"),
  na.col = "black", annot.grps = list("Groups" = seq(ncol(m))),
  annot.pal = grDevices::rainbow(n = ncol(m)), annot.size = 1, annot.sep = 0,
  lgd.scale.name = 'Values', lgd.bars.name = 'Legends',
  lgd.title = ggplot2::element_text(size = 12),
  lgd.text = ggplot2::element_text(size = 11), lgd.merge = FALSE,
  lgd.space.width = 1, y.lab = "Values", x.lab = "Samples",
  theme_heatmap = NULL, border.col = NA, border.size = 0.1, cell.size = 1,
  y.axis.right = FALSE,
  verbose = FALSE){

  #Check m is a matrix
  if(!is.matrix(m)){ stop("m must be a matrix.") }
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
  #Check if rank.fun is a supported function
  rank.method <- c("sd")
  if(!is.null(rank.fun)){
    if(!rank.fun %in% rank.method){ stop("ranking function not supported.") }
  }
  #Check if top.rows is an integer
  if(!is.null(top.rows)){
    if(is.numeric(top.rows)){
      top.rows <- as.integer(top.rows)
      if(top.rows <1){ stop("'top.rows' must be a non-zero positive integer.") }
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

  #Check annotation separation widths
  if(length(annot.sep) == 1){
    annot.cut <- annot.sep
    annot.sep <- annot.sep
  } else if(length(annot.sep) == 2){
    annot.cut <- annot.sep[2]
    annot.sep <- annot.sep[1]
  } else { stop("'annot.sep' length > 2. Too many values.") }

  #Check cell size
  if(length(cell.size) == 1){
    cell.width <- cell.size
    cell.height <- cell.size
  } else if(length(cell.size) == 2){
    cell.width <- cell.size[1]
    cell.height <- cell.size[2]
  } else { stop("'cell.size' length > 2. Too many values.") }


  #Check if y.axis.right = TRUE when axis.text.y.right or axis.title.y.right or
  # axis.ticks.y.right are not ggplot2::element_blank()
  theme_to_check <- c(
    "axis.text.y.right", "axis.title.y.right", "axis.ticks.y.right")
  thm_val <- theme_to_check[theme_to_check %in% names(theme_heatmap)]
  if(!identical(thm_val, character(0))){
    if(any(unlist(lapply(X = thm_val, FUN = function(t){
      !BiocompR::is.elt_blank(theme_heatmap[[t]])}))) & !y.axis.right){
      warning(
        paste("'y.axis.right' has to be set to TRUE in order to display the",
              "Y-axis on the right side of the heatmap."))
    }
  }
  #Check annotations groups and palettes matching
  BiocompR::check.annotations(
    data = m, annot.grps = annot.grps, annot.pal = annot.pal, verbose = verbose)

  #Handle NAs
  if(verbose){ cat("Managing missing values...") }
  m <- BiocompR::manage.na(
    data = m, method = na.handle, groups = imputation.grps, ncores = ncores)
  if(verbose){ cat("Done.\n") }
  #Apply ranking function if any function defined
  if(verbose){ cat("Ranking data by rows...") }
  if(!is.null(rank.fun)){
    #TODO: improve this part to support more function with eval() & parse()
    m <- m[order(apply(m, 1, sd, na.rm = T), decreasing = TRUE), , drop = FALSE]
  }
  if(verbose){ cat("Done.\n") }

  #Subset top rows if any value defined
  if(!is.null(top.rows)){ m <- utils::head(x = m, n = top.rows) }

  #Remove NAs if some for dendrogram matrix
  if(method.rows != 'none' | method.cols != 'none'){
    if(verbose){ cat("Clustering data...") }
    print(m)
    m_rname <- rownames(m)[apply(X = m, MARGIN = 1, FUN = anyNA) == FALSE]
    dend_mat <- m[stats::complete.cases(m), ]
    if(is.vector(dend_mat)){
      dend_mat <- matrix(
        dend_mat, nrow = 1, dimnames = list(m_rname, names(dend_mat)))
    }
    if(nrow(dend_mat) != nrow(m)){
      warning(paste("Distance method selected need complete data.",
                    nrow(m) - nrow(dend_mat), "incomplete rows removed out of",
                    nrow(m), "rows selected."))
    }
    #Check how many rows dend_mat has
    if(nrow(dend_mat) == 0){
      stop("Cannot compute distances on rows. All rows are missing values.")
    }
  }
  #Compute rows distances & create rows dendrogram
  if(method.rows != 'none'){
    row_dist <- tryCatch(
      parallelDist::parDist(dend_mat, method = method.rows, threads = ncores),
      error = function(cond){
        if(
          grepl(pattern = "impossible d'allouer un vecteur de taille", x = cond)
          | grepl(pattern = "cannot allocate vector of size", x = cond)){
          cond$message <- paste("Cannot compute distances on rows.",
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
        df = ggdendro::dendro_data(rowclust)$segments, orientation = "left",
        reverse.x = TRUE)
    }
  } else if(dd.rows & method.rows == 'none'){
    stop(paste(
      "Cannot plot dendrogram on rows with dist.method = 'none' for rows.",
      "To avoid this error message, set 'dendrograms' to FALSE,",
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
        df = ddgr_dat$segments, orientation = "top")
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
  } else if(method.rows != 'none' & is.null(rank.fun) & method.cols == 'none'){
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

  # dt.frame[, rn := factor(x = rn, levels = rownames(dframe))]
  melted_mat <- data.table::melt.data.table(
    data = dt.frame, id.vars = "rn", measure.vars = colnames(dt.frame)[-c(1)])
  if(!is.null(facet)){
    dt.facet <- data.table::data.table("variable" = colnames(dt.frame)[-c(1)],
                                       "facet.annot" = annot.grps[[facet]])
    melted_mat <- merge(
      x = melted_mat, y = dt.facet, by = "variable", all.x = TRUE, sort = FALSE)
  }
  rm(dt.frame)
  if(verbose){ cat("Done.\n") }

  if(verbose){ cat("Configure heatmap...") }

  #Set theme_forced_htmp
  theme_forced_htmp <- ggplot2::theme(
    plot.margin = ggplot2::margin(0, 0, 0, 0),
    legend.position = "bottom",
  )
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
    legend.text = ggplot2::element_text(size = 11),
    legend.title = ggplot2::element_text(size = 12),
    legend.justification = c(0.4, 0.5)
  )
  #Update theme_heatmap
  if(is.null(theme_heatmap)){
    theme_heatmap <- theme_default_htmp + theme_forced_htmp
  } else{
    theme_heatmap <- theme_default_htmp + theme_heatmap + theme_forced_htmp
  }
  # theme_heatmap <- ggplot2::theme(
  #   plot.margin = ggplot2::margin(0, 0, 0, 0),
  #   panel.grid = ggplot2::element_blank(),
  #   panel.background = ggplot2::element_rect(fill = "transparent"),
  #   plot.background = ggplot2::element_rect(fill = "transparent"),
  #   legend.text = ggplot2::element_text(size = 11),
  #   legend.title = ggplot2::element_text(size = 12), legend.position = "bottom",
  #   legend.justification = c(0.4, 0.5), axis.title.x = axis.title.x,
  #   axis.text.x = axis.text.x, axis.ticks.x = axis.ticks.x,
  #   axis.title.y.right = axis.title.y.right,
  #   axis.ticks.y.right = axis.ticks.y.right,
  #   axis.text.y.right = axis.text.y.right)

  #Set heatmap source parameters
  htmp.source <- ggplot2::ggplot() +
    ggplot2::scale_fill_gradientn(colours = heatmap.pal, na.value = na.col) +
    ggplot2::scale_color_manual(values = NA) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      ticks = TRUE, label = TRUE, barwidth = 15, ticks.linewidth = 1,
      title.vjust = 0.86, order = 1)) +
    ggplot2::guides(color = ggplot2::guide_legend(
      "NA", override.aes = list(fill = na.col), title.vjust = 0.5, order = 2)) +
    ggplot2::labs(x = x.lab, y = y.lab, fill = lgd.scale.name)

  #If facetting is on
  if(!is.null(facet)){
    htmp <- htmp.source +
      ggplot2::facet_grid(. ~ facet.annot, scales = "free", space = "free") +
      ggplot2::theme(panel.spacing = ggplot2::unit(0, "lines"),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_blank())
  } else { htmp <- htmp.source }

  #Draw heatmap with geom_tile
  if(nrow(melted_mat[is.na(value)]) != 0){
    #Display legend of missing values if any
    htmp <- htmp + ggplot2::geom_tile(data = melted_mat, ggplot2::aes(
      x = variable, y = rn, fill = value, color = " "), color = border.col,
      size = border.size, width = cell.width, height = cell.height)
  } else {
    htmp <- htmp + ggplot2::geom_tile(
      data = melted_mat, ggplot2::aes(x = variable, y = rn, fill = value),
      color = border.col, size = border.size,
      width = cell.width, height = cell.height)
  }

  if(dd.rows){
    theme_heatmap <- theme_heatmap +
      ggplot2::theme(axis.title.y.left = ggplot2::element_blank(),
                     axis.ticks.y.left = ggplot2::element_blank(),
                     axis.ticks.length.y.left = ggplot2::unit(0, "pt"),
                     axis.text.y.left = ggplot2::element_blank())
    # ,
    #                  plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"))
  }
  # } else {
  #   theme_heatmap <- theme_heatmap +
  #     ggplot2::theme(axis.title.y.left = axis.title.y.left,
  #                    axis.ticks.y.left = axis.ticks.y.left,
  #                    axis.text.y.left = axis.text.y.left,
  #                    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"))
  # }
  if(y.axis.right){
    htmp <- htmp +
      ggplot2::scale_y_discrete(position = 'right', expand = c(0, 0))
  } else { htmp <- htmp + ggplot2::scale_y_discrete(expand = c(0, 0)) }

  if(verbose){ cat("Done.\n") }

  if(verbose){ cat("Configure annotations...") }
  #Reoder groups and convert as factors
  annot.grps <- lapply(X = annot.grps, FUN = function(i){
    factor(x = i, levels =  unique(i))})
  #If the distance method used on columns is not "none" reorder columns
  if(method.cols != "none"){
    annot.grps <- lapply(X = annot.grps, FUN = function(i){ i[column.order] })
  }

  #Set number of columns to display annotations legends
  if(lgd.merge){
    origin.grps <- lapply(X = annot.grps, FUN = function(i){
      if(is.factor(i)){ levels(i) } else { levels(as.factor(i)) }
    })
    if(is.list(annot.pal)){
      if(length(origin.grps) == length(annot.pal)){
        ls.df.grp.pal <- Map(
          data.frame, "Grps" = origin.grps, "Cols" = annot.pal,
          stringsAsFactors = FALSE)
      } else {
        stop("The number of annotations does not match the number of palettes provided.")
      }
    } else if(is.vector(annot.pal)){
      ls.df.grp.pal <- lapply(X = origin.grps, FUN = function(grp){
        data.frame("Grps" = grp, "Cols" = annot.pal, stringsAsFactors = FALSE)
      })
    }
    #Rbind list color tables
    col_table <- data.table::rbindlist(ls.df.grp.pal, idcol = TRUE)
    if(is.vector(annot.pal) | length(annot.pal) == 1){
      #Remove duplicated colors
      col_table <- col_table[!duplicated(x = Cols)]
    }
    #Calculate legend length
    lgdsizes <- nrow(col_table) + 1
  } else {
    #Calculate legend length
    lgdsizes <- lapply(X = annot.grps, FUN = function(i){ length(unique(i)) })
    lgdsizes <- sum(unlist(lgdsizes)) + length(lgdsizes)
  }
  #Calculate legend columns
  lgd.ncol <- ceiling(lgdsizes/30)

  #Get ordered sample names
  if(method.cols != "none"){ sample.names <- colnames(dframe)
  } else { sample.names <- colnames(dframe) }

  #Create Color Sidebar
  col_sidebar <- BiocompR:::plot.col.sidebar(
    sample.names = sample.names, annot.grps = annot.grps,
    annot.pal = annot.pal, annot.pos = 'top', annot.sep = annot.sep,
    annot.cut = annot.cut, cor.order = seq_along(colnames(dframe)),
    merge.lgd = lgd.merge, right = TRUE, lgd.name = lgd.bars.name,
    lgd.title = lgd.title, lgd.text = lgd.text, lgd.ncol = lgd.ncol,
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 12, color = "black"),
    axis.ticks.y = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    set.x.title = NULL, set.y.title = NULL, dendro.pos = 'top', facet = facet)
  rm(sample.names)
  if(verbose){ cat("Done.\n") }

  #Extract Legend
  if(verbose){ cat("Extracting legends...") }
  #Create a subset plot of the original heatmap
  xtrm.melted_mat <- melted_mat[
    value %in% c(min(value, na.rm = TRUE), max(value, na.rm = TRUE))]
  if(nrow(melted_mat[is.na(value)]) != 0){
    subplot.htmp <- htmp.source +
      ggplot2::geom_tile(data = xtrm.melted_mat,
                         ggplot2::aes(x = variable, y = rn, fill = value, color = " "))
  } else {
    subplot.htmp <- htmp.source +
      ggplot2::geom_tile(
        data = xtrm.melted_mat,
        mapping = ggplot2::aes(x = variable, y = rn, fill = value))
  }
  htmp_legend <- BiocompR::get.lgd(gg2.obj = subplot.htmp + theme_heatmap)
  rm(subplot.htmp)

  sidebar_legend <- col_sidebar$legends
  #Convert ggplots into grobs
  if(dd.cols){ ddgr_seg_col <- ggplot2::ggplotGrob(ddgr_seg_col) }
  if(dd.rows){ ddgr_seg_row <- ggplot2::ggplotGrob(ddgr_seg_row) }
  col_sidebar_grob <- ggplot2::ggplotGrob(col_sidebar$sidebar)
  rm(col_sidebar)
  htmp <- htmp + theme_heatmap + ggplot2::theme(legend.position = "none")
  if(verbose){ cat("Done.\n") }

  #Heatmap rasterization
  if(!is.null(raster)){
    if(verbose){ cat("Rasterizing...\n") }
    if(raster %in% magick::filter_types()){
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
        axis.ticks.length = ggplot2::unit(0, "pt"))

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
                x = variable, y = rn, fill = value, color = " ")) +
              ggplot2::scale_fill_gradientn(
                colours = heatmap.pal, na.value = na.col) +
              ggplot2::scale_color_manual(values = NA) +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_discrete(expand = c(0, 0)) +
              ggplot2::guides(fill = ggplot2::guide_colorbar(
                ticks = TRUE, label = TRUE, barwidth = 15, ticks.linewidth = 1,
                title.vjust = 0.86, order = 1)) +
              ggplot2::guides(color = ggplot2::guide_legend(
                "NA", override.aes = list(fill = na.col), title.vjust = 0.5,
                order = 2)) +
              ggplot2::labs(x = x.lab, y = y.lab, fill = lgd.scale.name) +
              theme_empty + ggplot2::theme(legend.position = "none")
            rm(sub.melted)
            #Rasterize ggplot into a grob
            raster.grob <- BiocompR::raster.ggplot.to.grob(
              gg.plot = sub.htmp, filter = raster)
            rm(sub.htmp)
            #Make grob annotation
            raster.annot <- BiocompR::annotation_custom2(
              grob = raster.grob, xmin = -Inf, xmax = Inf, ymin = -Inf,
              ymax = Inf, data = melted_mat[facet.annot == i])
            rm(raster.grob)
            return(raster.annot)
          })
        #Fit the list of raster grobs into a ggplot
        htmp <- ggplot2::ggplot(
          data = melted_mat, ggplot2::aes(x = variable, y = rn, fill = value)) +
          ggplot2::geom_blank() + theme_heatmap +
          ggplot2::theme(legend.position = "none") +
          ggplot2::labs(x = x.lab, y = y.lab) +
          ggplot2::facet_grid(. ~ facet.annot, scales = "free",
                              space = "free") +
          ggplot2::theme(
            panel.spacing = ggplot2::unit(0, "lines"),
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_blank()) +
          ls.rasters
      } else {
        #Remove all customization
        htmp <- htmp + theme_empty
        #Catch heatmap in magick::image_graph()
        raster.annot <- BiocompR::raster.ggplot.to.grob(
          gg.plot = htmp, filter = raster)
        #Make grob annotation
        raster.annot <- ggplot2::annotation_custom(
          raster.grob, -Inf, Inf, -Inf, Inf)
        #Fit the raster grob into a ggplot
        htmp <- ggplot2::ggplot(
          data = melted_mat, ggplot2::aes(x = variable, y = rn, fill = value)) +
          ggplot2::geom_blank() + raster.annot + theme_heatmap +
          ggplot2::theme(legend.position = "none") +
          ggplot2::labs(x = x.lab, y = y.lab)
      }
    } else { stop("Rasterization filter not supported.") }
    if(verbose){ cat("Done.\n") }
  }
  #Remove melted_mat
  rm(melted_mat)
  if(verbose){ cat("Converting ggplot into grid object...") }
  htmp <- ggplot2::ggplotGrob(x = htmp)
  if(verbose){ cat("Done.\n") }
  #Resize grobs widths
  if(verbose){ cat("Redimensioning grobs...") }
  if(dd.cols & dd.rows){
    ls.w.grobs <- list(
      'dd_col' = ddgr_seg_col, 'sidebar' = col_sidebar_grob, 'htmp' = htmp)
    upd.grob_w <- BiocompR::resize.grobs(
      ls.grobs = ls.w.grobs, dimensions = "widths", start.unit = 4,
      end.unit = 7)
    rm(ddgr_seg_col)
  } else if(dd.cols & !dd.rows){
    ls.w.grobs <- list(
      'dd_col' = ddgr_seg_col, 'sidebar' = col_sidebar_grob, 'htmp' = htmp)
    upd.grob_w <- BiocompR::resize.grobs(
      ls.grobs = ls.w.grobs, dimensions = "widths", start.unit = 3,
      end.unit = 7)
    rm(ddgr_seg_col)
  } else {
    ls.w.grobs <- list('sidebar' = col_sidebar_grob, 'htmp' = htmp)
    if(is.null(facet)){
      upd.grob_w <- BiocompR::resize.grobs(
        ls.grobs = ls.w.grobs, dimensions = "widths", start.unit = 3,
        end.unit = 7)
    } else {
      upd.grob_w <- BiocompR::resize.grobs(
        ls.grobs = ls.w.grobs, dimensions = "widths", start.unit = 1,
        end.unit = max(unlist(lapply(X = ls.w.grobs, FUN = function(i){
          length(i[["widths"]]) }))))
    }
  }
  rm(htmp, col_sidebar_grob)
  #Resize grobs heights
  if(dd.rows){
    ls.h.grobs <- list('dd_row' = ddgr_seg_row, 'htmp' = upd.grob_w$htmp)
    upd.grob_h <- BiocompR::resize.grobs(
      ls.grobs = ls.h.grobs, dimensions = 'heights', start.unit = 7,
      end.unit = 9)
    rm(ddgr_seg_row)
  } else { upd.grob_h <- list("htmp" = upd.grob_w$htmp) }
  if(verbose){ cat("Done.\n") }

  #Create the Right Panel for legends
  if(verbose){ cat("Stacking legends...") }
  right.legends <- BiocompR:::stack.grobs.legends(
    grobs.list = sidebar_legend, annot.grps = annot.grps,
    height.lgds.space = 29)
  rm(sidebar_legend)

  if(verbose){ cat("Done.\n") }

  #Combine Dendrogram with Color Sidebar and Heatmap
  if(verbose){ cat("Creating final plot...") }
  if(dd.rows & dd.cols){
    #Create main grob
    main_grob <- gridExtra::arrangeGrob(
      grobs = list(grid::textGrob(""), upd.grob_w$dd_col,
                   grid::textGrob(""), upd.grob_w$sidebar,
                   upd.grob_h$dd_row, upd.grob_h$htmp),
      ncol = 2, nrow = 3, heights = c(dend.col.size + 2, annot.size, 30),
      widths = c(dend.row.size + 1, 10))
    #Set default legend width space
    def.lgd.width <- 2
  } else if(!dd.rows & !dd.cols){
    #Create main grob
    main_grob <- gridExtra::arrangeGrob(grobs = list(
      upd.grob_w$sidebar, upd.grob_h$htmp), ncol = 1, nrow = 2,
      heights = c(annot.size, 30), widths = 10)
    #Set default legend width space
    def.lgd.width <- 1
  } else if(dd.rows & !dd.cols){
    #Create main grob
    main_grob <- gridExtra::arrangeGrob(grobs = list(
      grid::textGrob(""), upd.grob_w$sidebar, upd.grob_h$dd_row,
      upd.grob_h$htmp), ncol = 2, nrow = 2, heights = c(annot.size, 30),
      widths = c(dend.row.size + 1, 10))
    #Set default legend width space
    def.lgd.width <- 2
  } else if(!dd.rows & dd.cols){
    #Create main grob
    main_grob <- gridExtra::arrangeGrob(grobs = list(
      upd.grob_w$dd_col, upd.grob_w$sidebar, upd.grob_h$htmp), ncol = 1,
      nrow = 3, heights = c(dend.col.size + 2, annot.size, 30), widths = 10)
    #Set default legend width space
    def.lgd.width <- 1
  }
  #rm(upd.grob_w, upd.grob_h)
  #Final plot
  final.plot <- gridExtra::grid.arrange(gridExtra::arrangeGrob(
    top = grid::textGrob(
      plot.title, gp = grid::gpar(fontsize = 15, font = 1)),
    grobs = list(grid::textGrob(paste0(
      "Columns ordered by ", method.cols, " distance; Rows ordered by ",
      method.rows, " distance; ", nrow(dframe), " ", row.type, "."),
      gp = grid::gpar(fontsize = 12, fontface = 3L)),
      gridExtra::arrangeGrob(grobs = list(main_grob, right.legends), ncol = 2,
                             widths = c(20, def.lgd.width + lgd.space.width)),
      htmp_legend), nrow = 3, heights = c(3, 50, 6)))
  rm(main_grob)
  #Prepare results
  ls.res <- list(
    "result.grob" = final.plot, "heatmap.grob" = upd.grob_h$htmp,
    "heatmap.lgd.grob" = htmp_legend, "sidebar.grob" = upd.grob_w$sidebar,
    "sidebar.lgds.grob" = right.legends)
  rm(final.plot, htmp_legend, right.legends)
  if(dd.cols){ ls.res[["cols.dendrogram.grob"]] <- upd.grob_w$dd_col }
  if(dd.rows){ ls.res[["rows.dendrogram.grob"]] <- upd.grob_h$dd_row }
  rm(upd.grob_w, upd.grob_h)
  if(verbose){ cat("Done.\n") }
  #Return a list of grobs with final plot and separate grobs.
  return(ls.res)
}
