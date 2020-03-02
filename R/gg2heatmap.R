
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
#' @param dend.col.size   A \code{numeric} defining the height of the dendrogram
#'                        made on columns (Default: dend.col.size = 1).
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
#' @param annot.lgd.space A \code{numeric} defining the size of the space
#'                        separating the different annotation bar legends
#'                        (Default: annot.lgd.space = 0).
#' @param annot.split     A \code{logical} to specify whether a space should be
#'                        added between each annotation (annot.split = TRUE) or
#'                        not (annot.split = FALSE)
#'                        (Default: annot.split = FALSE).
#' @param lgd.scale.name  A \code{character} to be used as legend title for the
#'                        heatmap scale (Default: lgd.scale.name = 'values').
#' @param lgd.bars.name   A \code{character} specifying the name of annotation
#'                        side bar legends, when legends are merged
#'                        (lgd.merge = TRUE)
#'                        (Default: lgd.bars.name = 'Legends').
#' @param lgd.pos.x       A \code{numeric} vector or unit object specifying the
#'                        x-location of the legends (Default: lgd.pos.x = 0.5).
#' @param lgd.pos.y       A \code{numeric} vector or unit object specifying the
#'                        y-location of the legends (Default: lgd.pos.y = 0.5).
#' @param lgd.merge       A \code{logical} specifying whether the legends of
#'                        multiple annotation bars should be merged
#'                        (lgd.merge = TRUE) or remain separated
#'                        (lgd.merge = FALSE). lgd.merge is especially useful
#'                        when you want to map the same color palette to
#'                        multiple annotations sharing the same values
#'                        (Default: lgd.merge = FALSE).
#' @param lgd.space.width A \code{numeric} specifying the width of the legend
#'                        space (Default: lgd.space.width = 1).
#' @param axis.title.x    An \code{element_text} object to setup X axis title
#'                        (Default: axis.title.x = element_text(size = 12,
#'                        color = 'black')).
#' @param axis.text.x     An \code{element_text} object to setup X axis text
#'                        (Default: axis.text.x = element_text(size = 11,
#'                        angle = -45, hjust = 0, vjust = 0.5, face = 'bold')).
#' @param axis.ticks.x    An \code{element_line} object to setup X axis ticks
#'                        (Default: axis.ticks.x = element_line(color='black')).
#' @param y.axis.right    A \code{logical} to specify whether the heatmap Y axis
#'                        should be displayed on the right (y.axis.right = TRUE)
#'                        or not (y.axis.right = FALSE)
#'                        (Default: y.axis.right = FALSE).
#' @param axis.title.y.right An \code{element_text} object to setup right Y axis
#'                           title
#'                           (Default: axis.title.y.right = element_blank()).
#' @param axis.text.y.right  An \code{element_text} object to setup right Y axis
#'                           text
#'                           (Default: axis.text.y.right = element_blank()).
#' @param axis.ticks.y.right An \code{element_line} object to setup right Y axis
#'                           ticks (Default: axis.ticks.x = element_blank()).
#' @return A \code{grob} of a heatmap plotted with ggplot2.
#' @author Yoann Pageaud.
#' @export

#TODO: Fix how legends are stacked and how spacing is calculated.
#TODO: Support rank.fun alone
#TODO: Support top.rows alone
#TODO: Merge dend.col.size & dend.row.size into dend.size being a tuple
#TODO: Create a theme argument using the theme() function
#TODO: Add some data in return
gg2heatmap<-function(m, na.handle = 'remove', dist.method = 'manhattan',
                     rank.fun = NULL, top.rows = NULL, dendrograms = TRUE,
                     dend.col.size = 1,
                     plot.title = "DMRs in 129 and BL6J WT Vs. 129 Mut DMNT1 HSPCs",
                     row.type = 'rows',
                     imputation.grps = NULL, ncores = 1,
                     heatmap.pal = c("steelblue", "gray95", "darkorange"),
                     annot.grps = list("Groups" = seq(ncol(m))),
                     annot.pal = rainbow(n = ncol(m)), annot.size = 1,
                     annot.lgd.space = 0, annot.split = FALSE,
                     lgd.scale.name = 'values', lgd.bars.name = 'Legends',
                     lgd.pos.x = 0.5, lgd.pos.y = 0.5, lgd.merge = FALSE,
                     lgd.space.width = 1,
                     axis.title.x = element_text(
                       size = 12, color = 'black'), axis.text.x = element_text(
                         size = 11, angle = -45, hjust = 0, vjust = 0.5,
                         face = 'bold'),
                     axis.ticks.x = element_line(color = 'black'),
                     y.axis.right = FALSE,
                     axis.title.y.right = element_text(size = 12),
                     axis.text.y.right = element_text(size = 12),
                     axis.ticks.y.right = element_line(color = 'black')){
  #Check m is a matrix
  if(!is.matrix(m)){ stop("m must be a matrix.") }
  #Check if na.handle method  supported
  na.method <- c('keep','impute','remove')
  if(!na.handle %in% na.method){ stop("na.handle method not supported.") }
  #Check dimensions of parameters
  if(length(dist.method) == 1){
    method.rows<-dist.method
    method.cols<-dist.method
  } else if(length(dist.method) == 2){
    method.rows<-dist.method[1]
    method.cols<-dist.method[2]
  } else { stop("'dist.method' length > 2. Too many values.") }
  #Check distance methods
  methods.ls <- c(
    'euclidean','maximum','manhattan','canberra','binary','minkowski','none')
  if(!method.rows %in% methods.ls){
    stop("Unknown method for distance computation on rows.")
  }
  if(!method.cols %in% methods.ls){
    stop("Unknown method for distance computation on columns.")
  }
  #Check if rank.fun is a supported function
  rank.method<-c("sd")
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
    dd.rows<-dendrograms
    dd.cols<-dendrograms
  } else if(length(dendrograms) == 2){
    dd.rows<-dendrograms[1]
    dd.cols<-dendrograms[2]
  } else { stop("'dendrograms' length > 2. Too many values.") }



  #Check annotations groups and palettes matching
  check.annotations(data = m, annot.grps = annot.grps, annot.pal = annot.pal)

  #Handle NAs
  m<-manage.na(data = m, method = na.handle, groups = imputation.grps,
               ncores = ncores)

  #Apply ranking function if any function defined
  if(!is.null(rank.fun)){
    #TODO: improve this part to support more function with eval() & parse()
    m <- m[order(apply(m, 1, sd, na.rm = T), decreasing = TRUE), , drop = FALSE]
  }
  #Subset top rows if any value defined
  if(!is.null(top.rows)){
    m <- head(x = m, n = top.rows)
  }

  #Remove NAs if some for dendrogram matrix
  dend_mat<- m[complete.cases(m),]

  #Create Dendrograms
  if(dd.rows & method.rows != 'none'){
    #Create Rows Dendrogram
    row_dist <- parallelDist::parDist(
      dend_mat, method = method.rows, threads = ncores)
    row_hclust<-hclust(row_dist)
    rm(row_dist)
    rowclust<-as.dendrogram(row_hclust)
    #Get dendrogram segments and order matrix rows
    ddgr_seg_row <- ggdend(df = ggdendro::dendro_data(rowclust)$segments,
                           orientation = "left", plot.type = 'heatmap')
    row.order<-order.dendrogram(rowclust)
  } else if(dd.rows & method.rows == 'none'){
    stop("Cannot plot rows dendrogram with method.rows = 'none'.")
  }

  if(dd.cols & method.cols != 'none'){
    #Create Columns Dendrogram
    ddgr <- as.dendrogram(hclust(parallelDist::parDist(
      t(dend_mat), method = method.cols, threads = ncores)))
    #Get dendrogram data
    ddgr_dat<-ggdendro::dendro_data(ddgr)
    #Get dendrogram segments and order matrix columns
    ddgr_seg_col <- ggdend(
      df = ddgr_dat$segments, orientation = "top", plot.type = 'heatmap')
    column.order<-order.dendrogram(ddgr)
  } else if(dd.cols & method.cols == 'none'){
    stop("Cannot plot columns dendrogram with method.cols = 'none'.")
  }

  #Reorder rows and columns matrix following dendrograms orders
  if(dd.rows & method.rows != 'none' & dd.cols & method.cols != 'none'){
    dframe <- m[row.order, column.order, drop = FALSE]
  } else if(dd.rows & method.rows != 'none' & is.null(rank.fun) & !dd.cols){
    dframe <- m[row.order, , drop = FALSE]
  } else if(!dd.rows & dd.cols & method.cols != 'none'){
    dframe <- m[, column.order, drop = FALSE]
  }

  melted_mat <- melt(dframe) #Melt Matrix into a dataframe
  #TODO: rm this
  # colnames(melted_mat)[3]<-"Methylation"

  #Plot Heatmap
  htmp <-ggplot() +
    geom_tile(data = melted_mat, aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradientn(colours = heatmap.pal) +
    scale_x_discrete(expand = c(0,0)) +
    theme(legend.justification = 'left', plot.margin = margin(0, 0, 0, 0),
          legend.position = c(0.5,0.5), panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"))
  if(dendrograms[1]){
    htmp <- htmp +
      theme(axis.title.y.left = element_blank(),
            axis.ticks.y.left = element_blank(),
            axis.ticks.length.y.left = unit(0, "pt"),
            axis.text.y.left = element_blank(),
            plot.margin=unit(c(0,0,0,0),"cm"))
  }
  htmp <- htmp +
    theme(legend.text=element_text(size= 11),
          legend.title = element_text(size = 12),
          axis.title.x = axis.title.x, axis.text.x = axis.text.x,
          axis.ticks.x = axis.ticks.x, axis.title.y.right = axis.title.y.right,
          axis.ticks.y.right = axis.ticks.y.right,
          axis.text.y.right = axis.text.y.right) +
    xlab("Samples") +
    labs(fill = lgd.scale.name)
  if(y.axis.right){
    htmp <- htmp + scale_y_discrete(position = 'right', expand = c(0,0))
  } else {
    htmp <- htmp + scale_y_discrete(expand = c(0,0))
  }

  #Reoder groups and convert as factors
  annot.grps <- lapply(X = annot.grps, FUN = function(i){
    factor(x = i, levels =  unique(i))})
  annot.grps <- lapply(X = annot.grps, FUN = function(i){i[column.order]})

  #Create Color Sidebar
  col_sidebar<-plot.col.sidebar(
    sample.names = colnames(m[, column.order]), annot.grps = annot.grps,
    annot.pal = annot.pal, annot.pos = 'top',
    cor.order = seq_along(colnames(dframe)), split.annot = annot.split,
    merge.lgd = lgd.merge, right = TRUE, lgd.lab = lgd.bars.name,
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    set.x.title = NULL, set.y.title = NULL, dendro.pos = 'top')

  #Extract Legend
  htmp_legend <- get.lgd(gg2.obj = htmp)
  sidebar_legend <- col_sidebar$legends
  #Convert ggplots into grobs
  ddgr_seg_col <- ggplotGrob(ddgr_seg_col)
  col_sidebar_grob <- ggplotGrob(col_sidebar$sidebar)
  ddgr_seg_row <- ggplotGrob(ddgr_seg_row)
  htmp <- ggplotGrob(htmp + theme(legend.position="none"))
  #Resize grobs
  upd.grob_w<-resize.grobs(ls.grobs = list(
    'dd_col' = ddgr_seg_col, 'sidebar' = col_sidebar_grob, 'htmp' = htmp),
    dimensions = "widths", start.unit = 4, end.unit = 7)
  upd.grob_h<- resize.grobs(ls.grobs = list(
    'dd_row' = ddgr_seg_row, 'htmp' = upd.grob_w$htmp), dimensions = 'heights',
    start.unit = 7, end.unit = 9)

  #Combine Dendrogram with Color Sidebar and Heatmap
  if((length(dendrograms) == 1 & dendrograms) |
     (length(dendrograms) == 2 & all(dendrograms == TRUE))){
    #Create main grob
    main_grob <- gridExtra::arrangeGrob(
      grobs = list(grid::textGrob(""), upd.grob_w$dd_col,
                   grid::textGrob(""), upd.grob_w$sidebar,
                   upd.grob_h$dd_row, upd.grob_h$htmp),
      ncol = 2, nrow = 3, heights = c(dend.col.size+2, annot.size, 30),
      widths = c(2, 10))
    #Create the Right Panel for legends
    sidebar_legend.grob <- gridExtra::arrangeGrob(
      grobs = sidebar_legend, nrow = 4, heights = c(4,1+annot.lgd.space,4,4))
    right.legends <- gridExtra::arrangeGrob(
      htmp_legend, grid::textGrob(""), sidebar_legend.grob,
      layout_matrix = cbind(c(1,1,1,2), c(2,2,2,2), c(3,3,3,3), c(2,2,2,2)),
      vp = grid::viewport(x= lgd.pos.x-0.6, y = lgd.pos.y))
    #Final plot
    final.plot<-gridExtra::grid.arrange(gridExtra::arrangeGrob(
      top = grid::textGrob(plot.title, gp = grid::gpar(fontsize = 15, font=1)),
      grobs = list(grid::textGrob(paste0(
        "Columns ordered by ", method.cols, " distance; Rows ordered by ",
        method.rows, " distance; ", nrow(m), " ", row.type, "."),
        gp = grid::gpar(fontsize = 12, fontface = 3L)),
        gridExtra::arrangeGrob(grobs = list(main_grob, right.legends), ncol = 2,
                               widths = c(20, 2 + lgd.space.width))),
      nrow = 2, heights = c(3,50)))
  }
  return(final.plot)
}
