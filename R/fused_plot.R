
#' Checks if min & max are given by legend:
#' \itemize{
#'  \item{If given: updates min & max of the plot.}
#'  \item{If none: calculates min & max from melted triangle.}
#' }
#'
#' @param lgd.limits A \code{vector} of length 2, containing limit values.
#'                   Can be NULL\cr(Default: lgd.limits=NULL).
#' @param melt.tri   A \code{data.frame} melted triangle containing test values.
#' @param tri.type   A \code{character} specifying the type of triangle provided
#'                   \cr(Supported: tri.type = c("upper","lower"))
#' @return A \code{list} of length 2 containing the updated min and max values.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

min.max.update<-function(lgd.limits=NULL, melt.tri, tri.type){
  if(is.null(lgd.limits)){
    min_tri<-min(melt.tri$value,na.rm = T)
    max_tri<-max(melt.tri$value,na.rm = T)
    if(min_tri == max_tri){
      stop(paste0("min and max values are equals. Please set limits for the ",
                  tri.type," triangle legend."))
    }
  } else {
    if(length(lgd.limits) != 2){
      if(tri.type == "upper"){lgd.name = "lgd.limits1"
      } else if(tri.type == "lower"){lgd.name = "lgd.limits2"}
      stop(paste0(lgd.name," should be a vector of length 2, containing the ",
                  "upper limit and the lower limit."))
    } else {
      min_tri<-lgd.limits[1]
      max_tri<-lgd.limits[2]
    }
  }
  return(list("min_tri"=min_tri,"max_tri"=max_tri))
}

#' Fixes a bug related to corrMatOrder when order parameter is set to 'alphabet'.
#'
#' @param cor.order A \code{character} vector of interest.
#' @param str       A \code{character} vector from which to get elements order.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

fix.corrMatOrder.alphabet<-function(cor.order,str){
  unlist(lapply(cor.order, function(i){
    grep(pattern=paste0("^",i,"$"),x=str)
  }))
}

#' Displays 2 matrices of results as a fused plot.
#'
#' @param param1 A \code{type} parameter description.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export

#TODO: Write documentation!
fused.view<-function(
  sample.names, upper.mat, lower.mat,
  order.select, order.method, hclust.method,
  correlation.order,
  annot.grps = list("Groups"=seq(ncol(data))),annot.pal=rainbow(n = ncol(data)),
  annot.pos = 'top', annot.size = 0, annot.text = NULL, annot.lgd.merge = FALSE,
  annot.split = TRUE,
  dendro.pos = 'none', dendro.size = 0,
  grid.col = "grey", grid.thickness = 0.5,
  axis.title = element_blank(), axis.title.x = element_blank(),
  axis.title.y = element_blank(), axis.text = element_text(size = 12),
  axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
  axis.text.y = element_blank(),
  axis.ticks = element_line(color = "black"),
  set.x.title = NULL, set.y.title = NULL,
  set.lgd1.title = NULL, set.lgd2.title = NULL,
  diag.col = "white",
  lgd.pal1 = NULL, lgd.pal2 = NULL,
  lgd.title = element_blank(), lgd.title1 = element_blank(),
  lgd.title2 = element_blank(),
  lgd.text = element_text(size = 12), lgd.text1 = element_blank(),
  lgd.text2 = element_blank(),
  lgd.breaks = NULL, lgd.breaks1 = NULL, lgd.breaks2 = NULL,
  lgd.labels = NULL, lgd.labels1 = NULL, lgd.labels2 = NULL,
  lgd.round = NULL, lgd.round1 = 2, lgd.round2 = 2,
  lgd.limits = NULL, lgd.limits1 = NULL, lgd.limits2 = NULL,
  lgd.ticks = TRUE, lgd.ticks1 = NULL, lgd.ticks2 = NULL,
  lgd.ticks.linewidth = 2, lgd.ticks.linewidth1 = NULL,
  lgd.ticks.linewidth2 = NULL,
  lgd.nbin = NULL, lgd.nbin1 = NULL, lgd.nbin2 = NULL,
  lgd.height1 = 30, lgd.height2 = 1, lgd.width1 = 1, lgd.width2 = 30,
  lgd.frame.col = "grey",lgd.frame.linewidth = 1.5,lgd.frame.linewidth1 = NULL,
  lgd.frame.linewidth2 = NULL,
  raster = FALSE, raster1 = NULL, raster2 = NULL,
  add.ggplot.arg = NULL
){

  #Get order of the correlations for the method used
  if(order.select == 'upper'){
    correlation.order<-corrMatOrder(upper.mat, order = order.method,
                                    hclust.method = hclust.method)
    if(order.method == "alphabet"){ #Fix bug of corrMatOrder when alphabet order
      correlation.order<-fix.corrMatOrder.alphabet(cor.order=correlation.order,
                                                   str =colnames(upper.mat))}
    if(dendro.pos != "none"){
      #Generate Hierarchy Cluster
      hierarchy.clust<-hclust(d = as.dist(1-upper.mat), method = hclust.method)
    }
  } else {
    correlation.order<-corrMatOrder(lower.mat, order = order.method,
                                    hclust.method = hclust.method)
    if(order.method == "alphabet"){ #Fix bug of corrMatOrder when alphabet order
      correlation.order<-fix.corrMatOrder.alphabet(cor.order=correlation.order,
                                                   str =colnames(upper.mat))}
    if(dendro.pos != "none"){
      #Generate Hierarchy Cluster
      hierarchy.clust<-hclust(d = as.dist(1-lower.mat), method = hclust.method)
    }
  }
  #Generate Dendrogram
  if(dendro.pos != "none"){
    dendrogram<-as.dendrogram(hierarchy.clust)
    ddgr_dat<-dendro_data(dendrogram) #Dendrogram data
    ddgr_seg <- ggdend( #Get dendrogram segments
      df = ddgr_dat$segments, orientation = dendro.pos, plot.type = 'corrplot')
  }
  #Re-order rows and columns
  upper.mat<-upper.mat[correlation.order,correlation.order]
  lower.mat<-lower.mat[correlation.order,correlation.order]
  #Replace half matrices by NAs
  upper.mat[upper.tri(upper.mat)]<-NA
  lower.mat[lower.tri(lower.mat)]<-NA
  #Melt Correlation matrix
  upper.melt<-melt(upper.mat,na.rm=T)
  lower.melt<-melt(lower.mat,na.rm=T)
  #Invert order of samples
  upper.melt$Var2<-factor(upper.melt$Var2,levels = rev(levels(upper.melt$Var2)))
  lower.melt$Var2<-factor(lower.melt$Var2,levels = rev(levels(lower.melt$Var2)))
  #Replace identical correlations by NA
  upper.melt[upper.melt$Var1 == upper.melt$Var2,]$value<-NA
  lower.melt[lower.melt$Var1 == lower.melt$Var2,]$value<-NA
  #Overwrite min and max of upper and lower melt matrices if limits have been
  # given for legends.
  if(is.null(lgd.limits)){
    upper.minmax<-min.max.update(
      lgd.limits = lgd.limits1, melt.tri = upper.melt, tri.type = "upper")
    min_upper<-upper.minmax$min_tri
    max_upper<-upper.minmax$max_tri

    lower.minmax<-min.max.update(
      lgd.limits = lgd.limits2, melt.tri = lower.melt, tri.type = "lower")
    min_lower<-lower.minmax$min_tri
    max_lower<-lower.minmax$max_tri
  } else {
    if(length(lgd.limits) != 2){
      stop("lgd.limits should be a vector of length 2, containing the upper
           limit and the lower limit.")
    } else {
      min_upper<-lgd.limits[1]
      max_upper<-lgd.limits[2]
      min_lower<-lgd.limits[1]
      max_lower<-lgd.limits[2]
      if(min_upper == max_upper){
        stop("min and max values are equals. Please set limits for the upper
             triangle legend.")
      }
      if(min_lower == max_lower){
        stop("min and max values are equals. Please set limits for the lower
             triangle legend.")
      }
    }
  }
  #Upper plot
  upper.ggplot<-basic.ggplot.tri(
    melt.tri = upper.melt, grid.col = grid.col, grid.thickness = grid.thickness,
    lgd.title = lgd.title1, lgd.text = lgd.text1, lgd.pal = lgd.pal1,
    min_tri = min_upper, max_tri = max_upper, lgd.breaks = lgd.breaks1,
    lgd.round = lgd.round1, lgd.ticks = lgd.ticks1, lgd.nbin = lgd.nbin1,
    lgd.height = lgd.height1, lgd.width = lgd.width1, rasteri = raster1,
    lgd.ticks.linewidth = lgd.ticks.linewidth1, lgd.frame.col = lgd.frame.col,
    lgd.frame.linewidth = lgd.frame.linewidth1, diag.col = diag.col,
    set.lgd.title = set.lgd1.title) + theme(axis.text = axis.text) +
    xlab(label = set.x.title) + ylab(label = set.y.title)
  if(annot.pos == "top"){
    upper.ggplot<-upper.ggplot + theme(
      axis.text.x = element_blank(),axis.text.y = axis.text.y,
      axis.ticks.x = element_blank(),axis.ticks.y = axis.ticks,
      axis.title.x = element_blank(),axis.title.y = axis.title.y)
  } else if(annot.pos == "left"){
    upper.ggplot<-upper.ggplot + theme(
      axis.text.x.top = axis.text.x, axis.text.y = element_blank(),
      axis.ticks.x = axis.ticks, axis.ticks.y = element_blank(),
      axis.title.x = axis.title.x, axis.title.y = element_blank())
  } else {stop("'annot.pos' value not supported.")}
  if(!(is.null(add.ggplot.arg))){upper.ggplot <- upper.ggplot + add.ggplot.arg}
  #Lower plot
  lower.ggplot<-basic.ggplot.tri(
    melt.tri = lower.melt, grid.col = grid.col, grid.thickness = grid.thickness,
    lgd.title = lgd.title2, lgd.text = lgd.text2, lgd.pal = lgd.pal2,
    min_tri = min_lower, max_tri = max_lower, lgd.breaks = lgd.breaks2,
    lgd.round = lgd.round2, lgd.ticks = lgd.ticks2, lgd.nbin = lgd.nbin2,
    lgd.height = lgd.height2, lgd.width = lgd.width2, rasteri = raster2,
    lgd.ticks.linewidth = lgd.ticks.linewidth2, lgd.frame.col = lgd.frame.col,
    lgd.frame.linewidth = lgd.frame.linewidth2, diag.col = diag.col,
    set.lgd.title = set.lgd2.title) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"),
          legend.direction = 'horizontal')
  #Plot Color Sidebar
  col_sidebar<-plot.col.sidebar(
    sample.names = sample.names ,annot.grps = annot.grps,annot.pal = annot.pal,
    annot.pos = annot.pos, cor.order = correlation.order,
    axis.ticks.x = axis.ticks, axis.ticks.y = axis.ticks,
    axis.title.x = axis.title.x, axis.title.y = axis.title.y,
    set.x.title = set.x.title, set.y.title = set.y.title,
    dendro.pos = dendro.pos, merge.lgd = annot.lgd.merge,
    split.annot = annot.split)
  if(annot.pos == "top"){
    if(is.null(annot.text)){
      annot.text <- element_text(size = 12, face = 'bold', hjust = 1, vjust=0.5)
    }
    col_sidebar$sidebar <- col_sidebar$sidebar +
      theme(axis.text.x.top = axis.text.x, axis.text.y = annot.text)
  } else if(annot.pos == "left"){
    if(is.null(annot.text)){
      annot.text <- element_text(size=12,angle=90,face='bold',hjust=0,vjust=0.5)
    }
    col_sidebar$sidebar <- col_sidebar$sidebar +
      theme(axis.text.x.top = annot.text, axis.text.y = axis.text.y)
  }
  #Remove all legends
  upper.ggplot.nolgd<-upper.ggplot + theme(legend.position = "none")
  lower.ggplot.nolgd<-lower.ggplot + theme(legend.position = "none")
  sidebar.nolgd<-col_sidebar$sidebar
  #Create grob for lower matrix
  lower.grob<-ggplotGrob(lower.ggplot.nolgd)
  #Add lower ggplot as an annotation to the upper ggplot
  main_grob<-upper.ggplot.nolgd+annotation_custom(lower.grob)
  #Convert ggplots into grobs
  main_grob<-ggplotGrob(main_grob)
  sidebar_grob<-ggplotGrob(sidebar.nolgd)
  if(dendro.pos!='none'){ dendro_grob<-ggplotGrob(ddgr_seg) }
  #Get legends
  upper.legend <- get.lgd(upper.ggplot)
  lower.legend <- get.lgd(lower.ggplot)
  sidebar.legend <- arrangeGrob(grobs = col_sidebar$legends,ncol = 1,
                                vp = viewport(height = 0.3))
  #Assemble grobs
  if(annot.pos == "top"){
    if(dendro.pos == "top"){
      #Get common width of a list of objects widths
      upd_grobs<-resize.grobs(ls.grobs = list(
        "main_grob" = main_grob, "sidebar_grob" = sidebar_grob,
        "dendro_grob" = dendro_grob), dimensions = "widths", start.unit = 2,
        end.unit = 5)

      main_grob<-arrangeGrob(upd_grobs$dendro_grob, upd_grobs$sidebar_grob,
                             upd_grobs$main_grob, ncol=1,
                             heights = c(10+dendro.size,9+annot.size,40))
    } else {
      #Get common width of a list of objects widths
      upd_grobs<-resize.grobs(ls.grobs = list(
        "main_grob" = main_grob, "sidebar_grob" = sidebar_grob),
        dimensions = "widths", start.unit = 2, end.unit = 5)

      main_grob<-arrangeGrob(upd_grobs$sidebar_grob, upd_grobs$main_grob,ncol=1,
                             heights = c(9+annot.size,40))
    }
    #Create the Right Panel for legends
    right.legends<-arrangeGrob(upper.legend,sidebar.legend, nrow = 1)
  } else {
    #Annotation on the left
    if(dendro.pos == "left"){
      #Get common width of a list of objects widths
      upd_grobs<-resize.grobs(ls.grobs = list(
        "main_grob" = main_grob, "sidebar_grob" = sidebar_grob,
        "dendro_grob" = dendro_grob), dimensions = "heights", start.unit = 3,
        end.unit = 8)

      #Make main grob
      main_grob<-arrangeGrob(
        upd_grobs$dendro_grob, upd_grobs$sidebar_grob,upd_grobs$main_grob,
        nrow=1, widths = c(10+dendro.size,9+annot.size,40))
    } else {
      #Get common width of a list of objects widths
      upd_grobs<-resize.grobs(ls.grobs = list(
        "main_grob" = main_grob, "sidebar_grob" = sidebar_grob),
        dimensions = "heights", start.unit = 3, end.unit = 8)
      #Make main grob
      main_grob<-arrangeGrob(upd_grobs$sidebar_grob, upd_grobs$main_grob,ncol=1,
                             widths = c(9+annot.size,40))
    }
    right.legends<-arrangeGrob(upper.legend, sidebar.legend, ncol = 2)
  }
  #Plot Final Figure
  fused.res<-grid.arrange(
    arrangeGrob(grobs = list(main_grob,right.legends,lower.legend),
                ncol = 2,nrow = 2, heights = c(40,3), widths = c(20,7)))
  return(fused.res)
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
#'                              method to apply.\cr Possible ordering methods are:
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
#'                              median' (= WPGMC) or 'centroid' (= UPGMC).
#' @param p.adjust              A \code{character} specifying what adjustment
#'                              for multiple tests should be used.\cr Possible
#'                              values are: "holm", "hochberg", "hommel",
#'                              "bonferroni", "BH", "BY", "fdr", "none".
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
#' @param grid.col              A \code{character} specifying the color of the
#'                              grid.
#' @param grid.thickness        A \code{double} value for the thickness of the
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
#' @param diag.col              A \code{character} defining the color of cells
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
#' @param add.ggplot.args       Ultimate ggplot2 additionnal components that you
#'                              want to pass to the plot (Use only if you are an
#'                              ggplot2 advanced user).
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

#TODO: Add parameter documentation for annot.text.
fused.plot<-function(data,ncores,
                     upper.comp,upper.value,lower.comp,lower.value,
                     na.rm = 'pairwise', order.method, order.select,
                     hclust.method = 'complete', p.adjust,
                     annot.grps = list("Groups" = seq(ncol(data))),
                     annot.pal = rainbow(n = ncol(data)),
                     annot.pos = 'top', annot.size = 0,
                     annot.text = NULL,
                     annot.lgd.merge = FALSE, annot.split = TRUE,
                     dendro.pos = 'none', dendro.size = 0,
                     grid.col = "grey",grid.thickness = 0.5,
                     axis.title = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text = element_text(size = 12),
                     axis.text.x =
                       element_text(angle = 90, hjust = 0, vjust = 0.5),
                     axis.text.y = element_blank(),
                     axis.ticks = element_line(color = "black"),
                     set.x.title = NULL, set.y.title = NULL,
                     set.lgd1.title = NULL,set.lgd2.title = NULL,
                     diag.col = "white",
                     lgd.pal1 = NULL, lgd.pal2 = NULL,
                     lgd.title = element_blank(),
                     lgd.title1 = element_blank(),
                     lgd.title2 = element_blank(),
                     lgd.text = element_text(size = 12),
                     lgd.text1 = element_blank(),lgd.text2 = element_blank(),
                     lgd.breaks = NULL,lgd.breaks1 = NULL,lgd.breaks2 = NULL,
                     lgd.labels = NULL,lgd.labels1 = NULL,lgd.labels2 = NULL,
                     lgd.round = NULL,lgd.round1 = 2,lgd.round2 = 2,
                     lgd.limits = NULL,lgd.limits1 = NULL,lgd.limits2 = NULL,
                     lgd.ticks = TRUE,lgd.ticks1 = NULL,lgd.ticks2 = NULL,
                     lgd.ticks.linewidth = 2,lgd.ticks.linewidth1 = NULL,
                     lgd.ticks.linewidth2 = NULL,
                     lgd.nbin = NULL,lgd.nbin1 = NULL, lgd.nbin2 = NULL,
                     lgd.height1 = 26,lgd.height2 = 1,
                     lgd.width1 = 1,lgd.width2 = 30,
                     lgd.frame.col = "grey", lgd.frame.linewidth = 1.5,
                     lgd.frame.linewidth1 = NULL,lgd.frame.linewidth2 = NULL,
                     raster = FALSE, raster1 = NULL, raster2 = NULL,
                     add.ggplot.arg = NULL){
  #Data format
  if(!(is.matrix(data))){if(is.data.frame(data)){data<-as.matrix(data)
  } else { stop("data is neither a matrix or a dataframe.") } }
  #Upper value type
  if(!(upper.value %in% c('r','n','stat','p','se'))){
    stop("value type unknown for upper.value")
  }
  #Lower value type
  if(!(lower.value %in% c('r','n','stat','p','se'))){
    stop("value type unknown for lower.value")
  }
  #Comparison Type
  if(na.rm %in% c("pairwise","complete")){
    cat(paste("Apply",na.rm,"deletion of missing data.\n"))
  } else { stop("the value of na.rm is not supported.") }
  #Order Method
  if(order.method %in% c("AOE","FPC","hclust","alphabet","default")){
    cat(paste("Apply",order.method,"ordering method.\n"))

    if(order.method == "hclust"){
      #check hclust.method
      if(hclust.method %in% c('ward.D','ward.D2', 'single', 'complete',
                              'average','mcquitty', 'median', 'centroid')){
        cat(paste("Apply",hclust.method,"clustering method.\n"))
      } else { stop("the clustering method is not supported.") }
      #check dendrogram
      if(dendro.pos != "none"){
        if(!(dendro.pos %in% c('top','left','both'))){
          stop("dendrogram cannot be put here.")
        }
      }
    } else {
      #check dendrogram
      if(dendro.pos != "none"){
        stop("order != 'hclust'. Dendrogram cannot be generated if rows & cols
         are not ordered following the hierarchical clustering.")
      }
    }
    #order select
    if(order.method %in% c("AOE","FPC","hclust","alphabet")){
      if(!(order.select %in% c('upper','lower'))){
        stop("Wrong value for order.select - You can only select values from the
         'upper' triangle or from the 'lower' triangle to apply order.")
      }
    }
  } else { stop("the order method is not supported.") }
  #P-value adjustment method
  if(p.adjust %in% c("holm","hochberg","hommel","bonferroni","BH","BY","fdr",
                     "none")){
    cat(paste("Apply",p.adjust,"multiple testing adjustment.\n"))
  } else { stop("multiple testing adjustment method is not supported.") }

  #Groups and palettes matching
  check.annotations(data = data, annot.grps = annot.grps, annot.pal = annot.pal)

  #Position annotation
  if(!(annot.pos %in% c("top","left","both"))){
    stop("annotation sidebar cannot be put here.")
  }
  #Overwrite the plot titles if they take the same parameters
  if(!(is.elt_blank(axis.title))){
    axis.title.x <- axis.title ; axis.title.y <- axis.title
  }
  #Label matching Breaks
  if(!(is.null(lgd.labels))){
    if(length(lgd.labels) != length(seq(lgd.breaks))){
      stop("the number of labels does not match the number of breaks.")
    }
  }
  if(!(is.null(lgd.labels1))){
    if(length(lgd.labels1) != length(lgd.breaks1)){
      stop("the number of labels does not match the number of breaks for upper
         triangle legend.")
    }
  }
  if(!(is.null(lgd.labels2))){
    if(length(lgd.labels2) != length(lgd.breaks2)){
      stop("the number of labels does not match the number of breaks for lower
         triangle legend.")
    }
  }
  #Overwrite the legends rounds if they take the same parameters
  if(!(is.null(lgd.round))){ lgd.round1<-lgd.round;lgd.round2<-lgd.round }
  #Overwrite the legends titles if they take the same parameters
  if(!(is.elt_blank(lgd.title))){ lgd.title1<-lgd.title;lgd.title2<-lgd.title }
  #Overwrite the legends titles if they take the same parameters
  if(!(is.elt_blank(lgd.text))){ lgd.text1<-lgd.text;lgd.text2<-lgd.text }
  #Overwrite the legends breaks if they take the same parameters
  if(!(is.null(lgd.breaks))){lgd.breaks1<-lgd.breaks;lgd.breaks2<-lgd.breaks}
  #Overwrite the legends bins if they take the same value
  if(!(is.null(lgd.nbin))){ lgd.nbin1<-lgd.nbin;lgd.nbin2<-lgd.nbin }
  #Overwrite the legends raster arg if they take the same value
  if(!(is.null(raster))){ raster1 <- raster ; raster2 <- raster }
  #Overwrite the legends ticks arg if they take the same value
  if(!(is.null(lgd.ticks))){ lgd.ticks1<-lgd.ticks ; lgd.ticks2<-lgd.ticks }
  #Overwrite the legends ticks linewidth arg if they take the same value
  if(!(is.null(lgd.ticks.linewidth))){
    lgd.ticks.linewidth1 <- lgd.ticks.linewidth
    lgd.ticks.linewidth2 <- lgd.ticks.linewidth
  }
  #Overwrite the legends frame linewidth arg if they take the same value
  if(!(is.null(lgd.frame.linewidth))){
    lgd.frame.linewidth1 <- lgd.frame.linewidth
    lgd.frame.linewidth2 <- lgd.frame.linewidth
  }
  #Checking comparisons
  if(upper.comp %in% c("pearson","spearman","kendall")){
    #Overwrite "stat" by "t" for correlation
    if(upper.value == "stat"){ upper.value <- 't' }
    #Compute correlation
    cat(paste("Compute pairwise",upper.comp,"correlation test..."))
    upper.correlation.res<-corr.test(data,use = na.rm, method = upper.comp,
                                     adjust = p.adjust)
    upper.mat<-upper.correlation.res[[upper.value]]
    upper.res<-upper.mat
    cat("Done.\n")
    if(lower.comp == upper.comp){
      #Overwrite "stat" by "t" for correlation
      if(lower.value == "stat"){ lower.value <- 't' }
      #Assign same matrix
      lower.correlation.res <- upper.correlation.res
      lower.mat<-lower.correlation.res[[lower.value]]
    } else {
      if(lower.comp %in% c("pearson","spearman","kendall")){
        #Overwrite "stat" by "t" for correlation
        if(lower.value == "stat"){ lower.value <- 't' }
        #Compute correlation
        cat(paste("Compute pairwise",lower.comp,"correlation test..."))
        lower.correlation.res<-corr.test(data,use = na.rm, method = lower.comp,
                                         adjust = p.adjust)
        lower.mat<-lower.correlation.res[[lower.value]]
        cat("Done.\n")
      } else if(lower.comp == "KS"){
        if(lower.value %in% c('n','stat','p')){
          cat(paste("Compute pairwise",lower.comp,"test..."))
          #Compute Kolmogorov-Smirnov test
          ks.res<-pairwise.ks(data=data,statistic=lower.value,ncores=ncores)
          lower.mat<-ks.res$res.statistic
          cat("Done.\n")
        } else if(lower.value == 'r'){
          stop("a KS test does not compute a correlation value.")
        } else if(lower.value == 'se'){
          stop("a KS test does not compute a standard error.")
        } else { stop("Unknown statistic for 'lower.value'.") }
      } else { stop("'lower.comp' value not supported yet.") }
    }
    lower.res<-lower.mat
  } else if(upper.comp == "KS"){
    if(upper.value %in% c('n','stat','p')){
      cat(paste("Compute pairwise",upper.comp,"test..."))
      #Compute Kolmogorov Smirnov test
      ks.res<-pairwise.ks(data = data, statistic = upper.value,ncores=ncores)
      upper.mat<-ks.res$res.statistic
      cat("Done.\n")
    } else if(upper.value == 'r'){
      stop("a KS test does not compute a correlation value.")
    } else if(upper.value == 'se'){
      stop("a KS test does not compute a standard error.")
    } else { stop("Unknown statistic for 'upper.value'.") }
    upper.res<-upper.mat
    if(lower.comp == upper.comp){
      if(lower.value %in% c('n','stat','p')){
        lower.mat<-get.ks.stat(table_combinations = ks.res$table_combinations,
                               df.ks.tests = ks.res$res.test, statistic)
      } else if(lower.value == 'r'){
        stop("a KS test does not compute a correlation value.")
      } else if(lower.value == 'se'){
        stop("a KS test does not compute a standard error.")
      } else { stop("Unknown statistic for 'upper.value'.") }
    } else {
      if(lower.comp %in% c("pearson","spearman","kendall")){
        #Overwrite "stat" by "t" for correlation
        if(lower.value == "stat"){ lower.value <- 't' }
        #Compute correlation
        cat(paste("Compute pairwise",lower.comp,"correlation test..."))
        lower.correlation.res<-corr.test(data,use = na.rm, method = lower.comp,
                                         adjust = p.adjust)
        lower.mat<-lower.correlation.res[[lower.value]]
        cat("Done.\n")
      } else { stop("'lower.comp' value not supported yet.") }
    }
    lower.res<-lower.mat
  } else { stop("'upper.comp' value not supported yet.") }

  #Create Fused Plot
  fused.res<-fused.view(
    sample.names = colnames(data), upper.mat = upper.mat, lower.mat = lower.mat,
    order.select = order.select, order.method = order.method,
    hclust.method = hclust.method, correlation.order = correlation.order,
    annot.grps = annot.grps, annot.pal = annot.pal, annot.pos = annot.pos,
    annot.size = annot.size, annot.text = annot.text, annot.lgd.merge = annot.lgd.merge,
    annot.split = annot.split, dendro.pos = dendro.pos, dendro.size=dendro.size,
    grid.col = grid.col, grid.thickness = grid.thickness,
    axis.title = axis.title, axis.title.x = axis.title.x,
    axis.title.y = axis.title.y, axis.text = axis.text, axis.text.x=axis.text.x,
    axis.text.y = axis.text.y,axis.ticks = axis.ticks,set.x.title = set.x.title,
    set.y.title = set.y.title, set.lgd1.title = set.lgd1.title,
    set.lgd2.title = set.lgd2.title,
    diag.col = diag.col,
    lgd.pal1 = lgd.pal1, lgd.pal2 = lgd.pal2,
    lgd.title = lgd.title, lgd.title1 = lgd.title1, lgd.title2 = lgd.title2,
    lgd.text = lgd.text, lgd.text1 = lgd.text1, lgd.text2 = lgd.text2,
    lgd.breaks = lgd.breaks,lgd.breaks1 = lgd.breaks1,lgd.breaks2 = lgd.breaks2,
    lgd.labels = lgd.labels,lgd.labels1 = lgd.labels1,lgd.labels2 = lgd.labels2,
    lgd.round = lgd.round, lgd.round1 = lgd.round1, lgd.round2 = lgd.round2,
    lgd.limits = lgd.limits,lgd.limits1 = lgd.limits1,lgd.limits2 = lgd.limits2,
    lgd.ticks = lgd.ticks, lgd.ticks1 = lgd.ticks1, lgd.ticks2 = lgd.ticks2,
    lgd.ticks.linewidth = lgd.ticks.linewidth,
    lgd.ticks.linewidth1 = lgd.ticks.linewidth1,
    lgd.ticks.linewidth2 = lgd.ticks.linewidth2,
    lgd.nbin = lgd.nbin, lgd.nbin1 = lgd.nbin1, lgd.nbin2 = lgd.nbin2,
    lgd.height1 = lgd.height1, lgd.height2 = lgd.height2,
    lgd.width1 = lgd.width1, lgd.width2 = lgd.width2,
    lgd.frame.col = lgd.frame.col, lgd.frame.linewidth = lgd.frame.linewidth,
    lgd.frame.linewidth1 = lgd.frame.linewidth1,
    lgd.frame.linewidth2 = lgd.frame.linewidth2,
    raster = raster, raster1 = raster1, raster2 = raster2,
    add.ggplot.arg = add.ggplot.arg)

  return(
    list("upper.res"=upper.res, "lower.res"=lower.res, "fused.plot"=fused.res))
}
