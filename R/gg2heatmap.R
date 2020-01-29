
#' WARNING: Unstable! Creates a custom heatmap with dendrograms and annotations.
#'
#' @param m A \code{matrix}.
#' @param na.handle A \code{character} to specify how missing values should be
#'                  handled (Default: na.handle = 'remove';
#'                  Supported: na.handle = c('keep','impute','remove')).
#' @param dist.method A \code{character} vector to specify the distance methods
#'                    to be used on the matrix rows and columns.
#'                    \itemize{
#'                     \item{If the vector is of length 1: the given method will
#'                           be applied on rows and columns of the matrix.}
#'                     \item{If the vector is of length 2: the first method will
#'                           be applied on matrix rows, and the second method
#'                           will be applied on matrix columns.}
#'                    }
#'                    (Default: dist.method = 'manhattan';
#'                    Supported: dist.method = c('euclidean', 'maximum',
#'                    'manhattan', 'canberra', 'binary', 'minkowski', 'none')).
#' @param rank.fun    A \code{character} to specify a function to apply on
#'                    matrix rows to order them on the final heatmap.
#'                    \itemize{
#'                     \item{If rank.fun = NULL the order of rows in the matrix
#'                           will be kept.}
#'                     \item{If a function is given, 'dist.method' must be set
#'                     to 'none' for the method to apply on rows.}
#'                    }
#'                    (Default: rank.fun = NULL; Supported: rank.fun = 'sd').
#' @param top.rows    An \code{integer} to specify the number of top rows you
#'                    want to display on the heatmap. If an integer is given,
#'                    'dist.method' must be set to 'none' for the method to
#'                    apply on rows (Default: top.rows = NULL).
#' @param dendrograms A \code{logical} vector to specify whether dendrograms
#'                    should be plotted with the heatmap.
#'                    \itemize{
#'                     \item{If dendrograms = TRUE: both dendrograms on row
#'                           and columns of the matrix will be plotted.}
#'                     \item{If dendrograms = FALSE: both dendrograms on rows
#'                           and columns of the matrix will NOT be plotted.}
#'                     \item{If the vector is of length 2: the first logical
#'                           will apply for rows, and the second logical will
#'                           apply for columns.}
#'                    }
#' @param heatmap.pal A \code{character} vector to be used as a color palette
#'                    for the heatmap.
#' @param annot.grps  A \code{list} of vectors of groups to which variables
#'                    belongs for the annotation sidebars. 'Vectors' lengths
#'                    have to match the number of variables.
#' @param annot.pal   A \code{list} of vectors to be used as color palettes for
#'                    the annotation sidebars. The length of vectors has to
#'                    match the number of levels of vectors listed in
#'                    'annot.grps'. If a list is provided, its length must match
#'                    the length of the list provided to 'annot.grps'.
#' @return A \code{grob} oh a heatmap plotted with ggplot2.
#' @author Yoann Pageaud.
#' @export
#' @examples

#TODO: Support na.handle
#TODO: Support tuple for dist.method
#TODO: Support rank.fun
#TODO: Support top.rows
#TODO: Check dendrograms
#TODO: Support annot.grps
#TODO: Support annot.pal
#TODO: Add some data in return
gg2heatmap<-function(m, na.handle = 'remove', dist.method = 'manhattan',
                     rank.fun = NULL, top.rows = NULL, dendrograms = TRUE,
                     heatmap.pal = c("steelblue", "gray95", "darkorange"),
                     annot.grps = list("Groups" = seq(ncol(m))),
                     annot.pal = rainbow(n = ncol(m))){
  ##_Handle NAs_################################################################
  m<-manage.na(data = m, method = na.handle, groups = sample_tbl$Cell.types)
  ##_Create Dendrogram_#########################################################
  #Remove NAs if some for dendrogram matrix
  dend_mat<- m[complete.cases(m),]

  if(length(dendrograms) == 1){
    if(dendrograms){
      if(length(dist.method) == 1 & dist.method != 'none'){
        #Create Rows Dendrogram
        row_dist<-parDist(dend_mat, method = dist.method, threads = 7)
        row_hclust<-hclust(row_dist)
        rm(row_dist)
        rowclust<-as.dendrogram(row_hclust)
        #Create Columns Dendrogram
        ddgr <- as.dendrogram(hclust(
          parDist(t(dend_mat), method = dist.method, threads = 7)))
      } else if(length(dist.method) == 2 & !('none' %in% dist.method)){
        #Create Rows Dendrogram
        row_dist<-parDist(dend_mat, method = dist.method[1], threads = 7)
        row_hclust<-hclust(row_dist)
        rm(row_dist)
        rowclust<-as.dendrogram(row_hclust)
        #Create Columns Dendrogram
        ddgr <- as.dendrogram(hclust(
          parDist(t(dend_mat), method = dist.method[2], threads = 7)))
      } else { stop("'dist.method' length > 2. Too many values.") }

      #Get dendrogram data
      ddgr_dat<-dendro_data(ddgr)
      #Get dendrogram segments and order matrix rows and/or columns
      ddgr_seg_row <- ggdend(df = dendro_data(rowclust)$segments,
                             orientation = "left", plot.type = 'heatmap')
      ddgr_seg_col <- ggdend(
        df = ddgr_dat$segments, orientation = "top", plot.type = 'heatmap')

      row.order<-order.dendrogram(rowclust)
      column.order<-order.dendrogram(ddgr)

      dframe<-m[row.order, column.order, drop = FALSE]
    }
  } else if(length(dendrograms) == 2){
    if(dendrograms[1]){
      if(length(dist.method) == 1 & dist.method != 'none'){
        #Create Rows Dendrogram
        row_dist<-parDist(dend_mat, method = dist.method, threads = 7)
      } else if(length(dist.method) == 2 & !('none' %in% dist.method)){
        #Create Rows Dendrogram
        row_dist<-parDist(dend_mat, method = dist.method[1], threads = 7)
      } else { stop("'dist.method' length > 2. Too many values.") }

      row_hclust<-hclust(row_dist)
      rm(row_dist)
      rowclust<-as.dendrogram(row_hclust)
      #Get dendrogram segments and order matrix rows and/or columns
      ddgr_seg_row <- ggdend(df = dendro_data(rowclust)$segments,
                             orientation = "left", plot.type = 'heatmap')
      row.order<-order.dendrogram(rowclust)
    }
    if(dendrograms[2]){
      if(length(dist.method) == 1 & dist.method != 'none'){
        #Create Columns Dendrogram
        ddgr <- as.dendrogram(hclust(
          parDist(t(dend_mat), method = dist.method, threads = 7)))
      } else if(length(dist.method) == 2 & !('none' %in% dist.method)){
        #Create Columns Dendrogram
        ddgr <- as.dendrogram(hclust(
          parDist(t(dend_mat), method = dist.method[2], threads = 7)))
      } else { stop("'dist.method' length > 2. Too many values.") }

      #Get dendrogram data
      ddgr_dat<-dendro_data(ddgr)
      #Get dendrogram segments and order matrix rows and/or columns
      ddgr_seg_col <- ggdend(
        df = ddgr_dat$segments, orientation = "top", plot.type = 'heatmap')
      column.order<-order.dendrogram(ddgr)
    }

    if(all(dendrograms == TRUE)){
      dframe<-m[row.order, column.order, drop = FALSE]
    } else if(dendrograms[1]){ dframe<-m[row.order,, drop = FALSE]
    } else if(dendrograms[2]){ dframe<-m[, column.order, drop = FALSE] }
  } else {
    stop("'dendrograms' length > 2. Too many values.")
  }
  melted_mat <- melt(dframe) #Melt Matrix into a dataframe
  colnames(melted_mat)[3]<-"Methylation"

  #Plot Heatmap
  htmp <-ggplot() +
    geom_tile(data = melted_mat, aes(x = Var2, y = Var1, fill=Methylation)) +
    scale_fill_gradientn(colours = palette.heatmap) +
    scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
    theme(legend.justification = 'left', plot.margin = margin(0, 0, 0, 0),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"))
  if(dendrograms[1]){
    htmp <- htmp +
      theme(axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.length.y = unit(0, "pt"),
            axis.text.y = element_blank(),
            plot.margin=unit(c(0,0,0,0),"cm"))
  }
  htmp <- htmp +
    theme(legend.text=element_text(size= 11),
          legend.title = element_text(size = 12),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 11, angle = -45, hjust = 0,
                                     vjust = 0.5, face = "bold")) +
    xlab("Samples")

  ##_Create Color Side Bar_#####################################################
  Samples<-colnames(dframe)
  #Create Genotypes column
  # Genotypes<-sample_tbl$Short.genotypes[column.order]
  Genotypes<-sample_tbl$Geno.short[column.order]
  # Genotypes<-factor(Genotypes,levels = color_table$Genotypes)
  Genotypes<-factor(Genotypes,levels = color_tbl$Genotypes)
  Genotype_table<-data.frame("Samples" = factor(Samples,levels = Samples),
                             "Genotypes" = Genotypes,row.names = NULL)

  col_sidebar<-plot.col.sidebar(
    sample.names = Genotype_table$Samples,
    annot.grps = list("Genotypes" = Genotypes),
    annot.pal = list(color_tbl$Colors), annot.pos = 'top',
    cor.order = seq_along(colnames(dframe)), split.annot = FALSE,merge.lgd=TRUE,
    right = TRUE, lgd.lab = "Genotypes", axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    set.x.title = NULL, set.y.title = NULL, dendro.pos = 'top')

  ##_Extract Legend_############################################################
  htmp_legend <- get.lgd(gg2.obj = htmp)
  sidebar_legend <- col_sidebar$legends
  ##convert ggplots into grobs
  ddgr_seg_col <- ggplotGrob(ddgr_seg_col)
  col_sidebar_grob <- ggplotGrob(col_sidebar$sidebar)
  ddgr_seg_row <- ggplotGrob(ddgr_seg_row)
  htmp <- ggplotGrob(htmp + theme(legend.position="none"))
  ## Resize grobs
  upd.grob_w<-resize.grobs(ls.grobs = list(
    'dd_col' = ddgr_seg_col, 'sidebar' = col_sidebar_grob, 'htmp' = htmp),
    dimensions = "widths", start.unit = 4, end.unit = 6)
  upd.grob_h<- resize.grobs(ls.grobs = list(
    'dd_row' = ddgr_seg_row, 'htmp' = upd.grob_w$htmp), dimensions = 'heights',
    start.unit = 7, end.unit = 9)

  ##_Combine Dendrogram with Color Sidebar and Heatmap_#########################
  if((length(dendrograms) == 1 & dendrograms) |
     (length(dendrograms) == 2 & all(dendrograms == TRUE))){
    #Create main grob
    main_grob <- arrangeGrob(grobs = list(
      textGrob(""), upd.grob_w$dd_col, textGrob(""), upd.grob_w$sidebar,
      upd.grob_h$dd_row, upd.grob_h$htmp),ncol = 2,nrow = 3,
      heights = c(3,1,30), widths = c(2,10))
    #Create the Right Panel for legends
    sidebar_legend.grob <- arrangeGrob(grobs = sidebar_legend, ncol = 1)
    right.legends <- arrangeGrob(
      textGrob(""), sidebar_legend.grob, htmp_legend, textGrob(""),textGrob(""),
      ncol = 1, vp = viewport(x = -0.5))
    #Final plot
    grid.arrange(arrangeGrob(
      top = textGrob("DMRs in 129 and BL6J WT Vs. 129 Mut DMNT1 HSPCs",
                     gp = gpar(fontsize = 15, font = 1)),
      grobs = list(
        textGrob(paste("(Columns ordered by", dist.method,
                       "distance; Rows ordered by", dist.method, "distance;",
                       nrow(m), "loci.)"), gp = gpar(fontsize=12, fontface=3L)),
        arrangeGrob(grobs = list(main_grob, right.legends), ncol = 2,
                    widths = c(20,1))),
      nrow = 2, heights = c(3,50)))
  }
}
