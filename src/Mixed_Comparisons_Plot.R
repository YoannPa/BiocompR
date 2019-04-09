
##IMPORTS
Imports = c('gridExtra','Hmisc','psych','corrplot','reshape2','ggsci',
            'ggdendro','gtools','grid','fastcluster')
lapply(Imports, library, character.only = T)

##FUNCTIONS

# is.in ########################################################################

is.in<-function(arg,vector.val){
  bool<-arg %in% vector.val
  return(bool)
}

# is.elt_blank #################################################################

is.elt_blank<-function(arg){
  bool<-attributes(arg)$class[1] == "element_blank"
  return(bool)
}

# is.none ######################################################################

is.none<-function(arg){
  bool<-arg == 'none'
  return(bool)
}

# get.lgd ######################################################################

#' Extract legend from a ggplot object.
#' 
#' @param gg2.obj  A \code{gg} object with legends.
#' @return a \code{gg} object only containing the legends of the plot.
#' @author Yoann Pageaud.

get.lgd<-function(gg2.obj){
tmp <- ggplot_gtable(ggplot_build(gg2.obj))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
legend
}

# ggdend #######################################################################

#' Create a dendogram in ggplot
#' 
#' @param df          A \code{data.frame} containing variables and value to be
#'                    used to create the dendrogram.
#' @param orientation A \code{character} specifying the orientation of the
#'                    dendrogram. Possible values are "top" and "left".
#' @return A \code{gg} object of the dendrogram that can be plotted.
#' @author Yoann Pageaud.

ggdend <- function(df,orientation) {
  ddplot<- ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    expand_limits(x = c(0.5,max(df$x)+0.5), y = 0)
  if(orientation == "top"){
    ddplot <- ddplot +
      scale_x_continuous(expand = c(0,0)) +
      theme(plot.margin=unit(c(0.1,0.1,0.1,-0.1),"cm")) +
      scale_y_continuous(expand = c(0, 0))
  } else {
    ddplot <- ddplot +
      theme(plot.margin=unit(c(0,0.1,0,0.1),"cm")) +
      scale_y_reverse(expand = c(0,0)) +
      scale_x_reverse(expand = c(0,0)) +
      coord_flip()
  }
  ddplot
}

# mix.comp.plot ################################################################

#' Create a plot summarizing results from 2 different pairwise comparisons.
#' 
#' @param data                  A \code{matrix} or \code{dataframe}.
#' @param upper.comp            The comparison for which results will be
#'                              displayed in the upper triangle of the plot as a
#'                              \code{character} matching one of these:
#'                              'pearson','spearman','kendall'.
#' @param upper.value           A \code{character} matching the type of value to
#'                              display in the upper triangle of the plot
#'                              resulting from the comparison. The possible
#'                              value types are: 'r' - the correlation values,
#'                              'n' - the number of cases used for the
#'                              comparisons, 'stat' - the values of the
#'                              comparison statistic, 'p' - the P-values of the
#'                              comparison, 'se' - the standard errors of the
#'                              comparison.
#' @param lower.comp            The comparison for which results will be
#'                              displayed in the lower triangle of the plot as a
#'                              \code{character} matching one of these:
#'                              'pearson','spearman','kendall'.
#' @param lower.value           A \code{character} matching the type of value to
#'                              display resulting from the comparison. The
#'                              possible value types are: 'r' - the correlation
#'                              value, 'n' - the number of cases used for the
#'                              comparison, 'stat' - the value of the comparison
#'                              statistic, 'p' - the P-value of the comparison,
#'                              'se' - the standard error of the comparison. 
#' @param na.rm                 A \code{character} to specify how to handle
#'                              missing values when calculating a correlation.
#'                              Possible types are 'pairwise' and 'complete'.
#'                              'pairwise' is the default value and will do
#'                              pairwise deletion of cases. 'complete' will
#'                              select just complete cases.  
#' @param order.method          A \code{character} specifying the ordering
#'                              method to apply. Possible ordering methods are:
#'                              'AOE' - angular order of the eigenvectors,
#'                              'FPC' - first principal component order,
#'                              'hclust' - hierarchical clustering order,
#'                              'alphabet' - alphabetical order, 'default'.
#' @param order.select          A \code{character} specifying comparison results
#'                              to use for the clustering. Possible values are
#'                              'upper' or 'lower'.
#' @param hclust.method         A \code{character} specifying the method to use
#'                              for the hierarchical clustering if the ordering
#'                              method is 'hclust'. Possible methods are:
#'                              'ward.D','ward.D2', 'single', 'complete',
#'                              'average' (= UPGMA), 'mcquitty' (= WPGMA),
#'                              median' (= WPGMC) or 'centroid' (= UPGMC).
#' @param p.adjust              A \code{character} specifying what adjustment
#'                              for multiple tests should be used. Possible
#'                              values are: "holm", "hochberg", "hommel",
#'                              "bonferroni", "BH", "BY", "fdr", "none".
#' @param annot.grps            A \code{character} vector of groups to which
#'                              variables belongs for the annotation sidebar.
#'                              The length of this vector has to match the
#'                              number of variables.
#' @param annot.pal             A \code{character} vector of colors for the
#'                              annotation sidebar. The length of this vector
#'                              has to match the number of different groups
#'                              existing.
#' @param annot.pos             A \code{character} specifying the position of
#'                              the annotation sidebar. Possible values are:
#'                              'top', 'left' or 'both'.
#' @param annot.size            An \code{integer} to increase or decrease the
#'                              size of the annotation side bar.
#' @param dendro.pos            A \code{character} specifying the position of
#'                              dendrograms if the selected order is 'hclust'.
#'                              Possible values are: 'top','left','none'.
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
#' @param lgd.nbin             An \code{integer} specifying the number of bins
#'                              for drawing both colourbars. A smoother
#'                              colourbar results from a larger value.
#' @param lgd.nbin1            An \code{integer} specifying the number of bins
#'                              for drawing the upper triangle associated
#'                              colourbar. A smoother colourbar results from a
#'                              larger value.
#' @param lgd.nbin2            An \code{integer} specifying the number of bins
#'                              for drawing the lower triangle associated
#'                              colourbar. A smoother colourbar results from a
#'                              larger value.
#' @param lgd.height1           A \code{double} specifying the height of the
#'                              upper triangle associated colourbar.
#' @param lgd.height2           A \code{double} specifying the height of the
#'                              lower triangle associated colourbar.
#' @param lgd.width1            A \code{double} specifying the width of the
#'                              upper triangle associated colourbar.
#' @param lgd.width2            A \code{double} specifying the height of the
#'                              lower triangle associated colourbar.
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
#' @param add.ggplot.args       Ultimate Ggplot2 additionnal argument that you
#'                              want to pass to the plot (Use only if you are an
#'                              ggplot2 advanced user).
#' @return a \code{gtable} object automatically plot and a \code{list} of the
#'         results of the 2 comparisons made as 2 matrices: 'upper.res' - the
#'         matrix result of pairwise comparisons displayed in the upper
#'         triangle, 'lower.res' - the matrix result of pairwise comparisons
#'         displayed in the lower triangle.
#' 
#' @author Yoann Pageaud.

mix.comp.plot<-function(data,
                        upper.comp,upper.value,lower.comp,lower.value,
                        na.rm = 'pairwise', order.method, order.select,
                        hclust.method = 'complete', p.adjust,
                        annot.grps = seq(ncol(data)),
                        annot.pal = rainbow(n = ncol(data)),
                        annot.pos = 'top', annot.size = 0,
                        dendro.pos = 'none',dendro.size = 0,
                        grid.col = "grey",grid.thickness = 0.5,
                        axis.title = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text = element_text(size = 12),
                        axis.text.x = element_text(angle = 90,hjust = 0),
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
                        lgd.breaks = NULL ,lgd.breaks1 = NULL,lgd.breaks2 = NULL,
                        lgd.labels = NULL,lgd.labels1 = NULL,lgd.labels2 = NULL,
                        lgd.round = NULL,lgd.round1 = 2,lgd.round2 = 2,
                        lgd.limits = NULL,lgd.limits1 = NULL,lgd.limits2 = NULL,
                        lgd.ticks = T,lgd.ticks1 = NULL,lgd.ticks2 = NULL,
                        lgd.ticks.linewidth = 2,lgd.ticks.linewidth1 = NULL,
                        lgd.ticks.linewidth2 = NULL,
                        lgd.nbin = NULL,lgd.nbin1 = NULL, lgd.nbin2 = NULL,
                        lgd.height1 = 33,lgd.height2 = 1,
                        lgd.width1 = 1,lgd.width2 = 33,
                        lgd.frame.col = "grey", lgd.frame.linewidth = 1.5,
                        lgd.frame.linewidth1 = NULL,lgd.frame.linewidth2 = NULL,
                        raster = F, raster1 = NULL, raster2 = NULL,
                        add.ggplot.arg = NULL){
  ##Check arguments
  #data format
  if(is.matrix(data) == F){
    if(is.data.frame(data)){
      data<-as.matrix(data)
      } else {
        stop("data is neither a matrix or a dataframe.")
      }
  }
  #upper value type
  if(!(is.in(upper.value,c('r','n','stat','p','se')))){
    stop("value type unknown for upper.value")
  } 
  #lower value type
  if(!(is.in(lower.value,c('r','n','stat','p','se')))){
    stop("value type unknown for lower.value")
  }
  #comparison type
  if(is.in(na.rm,c("pairwise","complete"))){
    cat(paste("Apply",na.rm,"deletion of missing data.\n"))
  } else {
    stop("the value of na.rm is not supported.")
  }
  #order method
  if(is.in(order.method,c("AOE","FPC","hclust","alphabet","default"))){
    cat(paste("Apply",order.method,"ordering method.\n"))
    
    if(order.method == "hclust"){
      #hclust.method
      if(is.in(hclust.method,c('ward.D','ward.D2', 'single', 'complete',
                               'average','mcquitty', 'median', 'centroid'))){
        cat(paste("Apply",hclust.method,"clustering method.\n"))
      } else {
        stop("the clustering method is not supported.")
      }
      #check dendrogram
      if(!(is.none(dendro.pos))){
        if(!(is.in(dendro.pos,c('top','left','both')))){
          stop("dendrogram cannot be put here.")
        }
      }
    } else {
      #check dendrogram
      if(!(is.none(dendro.pos))){
        stop("order != 'hclust'. Dendrogram cannot be generated if rows & cols
         are not ordered following the hierarchical clustering.")
      }
    }
    #order select 
    if(is.in(order.method,c("AOE","FPC","hclust","alphabet"))){
      if(!(is.in(order.select,c('upper','lower')))){
        stop("Wrong value for order.select - You can only select values from the
         'upper' triangle or from the 'lower' triangle to apply order.")
      }
    }
  } else { stop("the order method is not supported.") }
  
  #P-value adjustment method
  if(is.in(p.adjust,c("holm","hochberg","hommel","bonferroni","BH","BY","fdr",
                      "none"))){
    cat(paste("Apply",p.adjust,"multiple testing adjustment.\n"))
  } else { stop("multiple testing adjustment method is not supported.") }
  
  #Groups checking
  if(length(annot.grps) != ncol(data)){
    stop("samples are not all assigned to a group.")
  } else{
    cat(paste("Groups: ",paste(unique(annot.grps), collapse = ", "),".\n",
              sep = ""))
  }
  #Color checking
  if(length(unique(annot.pal)) != length(unique(annot.grps))) {
    stop("the number of colors does not match the number of groups.")    
  }
  #Position annotation
  if(!(is.in(annot.pos,c("top","left","both")))){
    stop("annotation sidebar cannot be put here.")
  }
  #Overwrite the plot titles if they take the same parameters
  if(!(is.elt_blank(axis.title))){
    axis.title.x <- axis.title ; axis.title.y <- axis.title
  }
  #Label matching Breaks
  if(is.null(lgd.labels) == F){
    if(length(lgd.labels) != length(seq(lgd.breaks))){
      stop("the number of labels does not match the number of breaks.")
    }
  }
  if(is.null(lgd.labels1) == F){
    if(length(lgd.labels1) != length(lgd.breaks1)){
      stop("the number of labels does not match the number of breaks for upper
         triangle legend.")
    }
  }
  if(is.null(lgd.labels2) == F){
    if(length(lgd.labels2) != length(lgd.breaks2)){
      stop("the number of labels does not match the number of breaks for lower
         triangle legend.")
    }
  }
  #Overwrite the legends rounds if they take the same parameters
  if(is.null(lgd.round) == F){ lgd.round1<-lgd.round;lgd.round2<-lgd.round }
  
  #Overwrite the legends titles if they take the same parameters
  if(!(is.elt_blank(lgd.title))){ lgd.title1<-lgd.title;lgd.title2<-lgd.title }
     
  #Overwrite the legends titles if they take the same parameters
  if(!(is.elt_blank(lgd.text))){ lgd.text1<-lgd.text;lgd.text2<-lgd.text }
  
  #Overwrite the legends breaks if they take the same parameters
  if(is.null(lgd.breaks) == F){lgd.breaks1<-lgd.breaks;lgd.breaks2<-lgd.breaks}
  
  #Overwrite the legends bins if they take the same value
  if(is.null(lgd.nbin) == F){ lgd.nbin1<-lgd.nbin;lgd.nbin2<-lgd.nbin }
  #Overwrite the legends raster arg if they take the same value
  if(is.null(raster) == F){ raster1 <- raster ; raster2 <- raster }
  
  #Overwrite the legends ticks arg if they take the same value
  if(is.null(lgd.ticks) == F){ lgd.ticks1<-lgd.ticks ; lgd.ticks2<-lgd.ticks }
  
  #Overwrite the legends ticks linewidth arg if they take the same value
  if(is.null(lgd.ticks.linewidth) == F){
    lgd.ticks.linewidth1 <- lgd.ticks.linewidth
    lgd.ticks.linewidth2 <- lgd.ticks.linewidth
  }
  #Overwrite the legends frame linewidth arg if they take the same value
  if(is.null(lgd.frame.linewidth) == F){
    lgd.frame.linewidth1 <- lgd.frame.linewidth
    lgd.frame.linewidth2 <- lgd.frame.linewidth
  }
  
  #Checking comparison
  if(is.in(upper.comp,c("pearson","spearman","kendall"))){
    #Overwrite "stat" by "t" for correlation
    if(upper.value == "stat"){ upper.value <- 't' }
    #Compute correlation
    cat(paste("Compute pariwise",capitalize(upper.comp),
              "correlation test...\n"))
    upper.correlation.res<-corr.test(data,use = na.rm, method = upper.comp,
                                     adjust = p.adjust)
    upper.mat<-upper.correlation.res[[upper.value]]
    upper.res<-upper.mat
    
    if(lower.comp == upper.comp){
      #Overwrite "stat" by "t" for correlation
      if(lower.value == "stat"){ lower.value <- 't' }
      #Assign same matrix
      lower.correlation.res = upper.correlation.res
    } else {
      if(is.in(lower.comp,c("pearson","spearman","kendall"))){
        #Overwrite "stat" by "t" for correlation
        if(lower.value == "stat"){ lower.value <- 't' }
        #Compute correlation
        cat(paste("Compute pariwise",capitalize(lower.comp),
                  "correlation test...\n"))
        lower.correlation.res<-corr.test(data,use = na.rm, method = lower.comp,
                                         adjust = p.adjust)
      }
    }
    lower.mat<-lower.correlation.res[[lower.value]]
    lower.res<-lower.mat
  } else { stop("comparison not supported yet.") }

  #Get order of the correlations for the method used
  if(order.select == 'upper'){
    correlation.order<-corrMatOrder(upper.mat, order = order.method,
                                    hclust.method = hclust.method)
    if(!(is.none(dendro.pos))){
    #Generate Hierarchy Cluster
    hierarchy.clust<-hclust(d = as.dist(1-upper.mat), method = hclust.method)
    }
  } else {
    correlation.order<-corrMatOrder(lower.mat, order = order.method,
                                    hclust.method = hclust.method)
    if(!(is.none(dendro.pos))){
    #Generate Hierarchy Cluster
    hierarchy.clust<-hclust(d = as.dist(1-lower.mat), method = hclust.method)
    }
  }
  
  if(!(is.none(dendro.pos))){
    #Generate Dendrogram
    dendrogram<-as.dendrogram(hierarchy.clust)
    ddgr_dat<-dendro_data(dendrogram) #Dendrogram data    
    ddgr_seg <- ggdend(ddgr_dat$segments,dendro.pos) #Get dendrogram segments
  }
  upper.mat<-upper.mat[correlation.order,correlation.order]
  lower.mat<-lower.mat[correlation.order,correlation.order]
  
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
    if(is.null(lgd.limits1)){
      min_upper<-min(upper.melt$value,na.rm = T)
      max_upper<-max(upper.melt$value,na.rm = T)
      if(min_upper == max_upper){
        stop("min and max values are equals. Please set limits for the upper
             triangle legend.")
      }
    } else {
      if(length(lgd.limits1) != 2){
        stop("lgd.limits1 should be a vector of length 2, containing the upper
           limit and the lower limit.")
      } else {
        min_upper<-lgd.limits1[1]
        max_upper<-lgd.limits1[2]  
      }
    }
    
    if(is.null(lgd.limits2)){
      min_lower<-min(lower.melt$value,na.rm = T)
      max_lower<-max(lower.melt$value,na.rm = T)
      if(min_lower == max_lower){
        stop("min and max values are equals. Please set limits for the lower
             triangle legend.")
      }
    } else {
      if(length(lgd.limits2) != 2){
        stop("lgd.limits2 should be a vector of length 2, containing the upper
           limit and the lower limit.")
      } else {
        min_lower<-lgd.limits2[1]
        max_lower<-lgd.limits2[2]  
      }
    }
    
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
  
  ##Upper plot
  upper.ggplot<-ggplot() + 
    geom_tile(data = upper.melt, aes(x=Var1, y=Var2, fill = value),
              color = grid.col, size=grid.thickness) +
    theme(axis.text = axis.text,
          panel.grid = element_blank(),
          legend.title = lgd.title1,
          legend.text = lgd.text1,
          legend.justification = c(1, 0.5),
          legend.position = c(1,0.37), 
          plot.margin = margin(0, 0, 0, 0)) +
    scale_x_discrete(position = "top",expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0))+
    scale_fill_gradientn(colours = lgd.pal1,
                         breaks = seq(min_upper,max_upper,
                                      length.out = lgd.breaks1),
                         labels = paste(" ",
                                        round(seq(min_upper,max_upper,
                                                  length.out = lgd.breaks1),
                                              lgd.round1)),
                         guide=guide_colourbar(ticks=lgd.ticks1,
                                               nbin = lgd.nbin1,
                                               barheight=lgd.height1,
                                               label=T, barwidth=lgd.width1,
                                               raster = raster1,
                                               ticks.linewidth
                                               = lgd.ticks.linewidth1,
                                               frame.colour = lgd.frame.col,
                                               frame.linewidth
                                               = lgd.frame.linewidth1),
                         na.value = diag.col,
                         limits=c(min_upper,max_upper),
                         name = set.lgd1.title) +
    xlab(label = set.x.title) +
    ylab(label = set.y.title)
  if(annot.pos == "top"){
    upper.ggplot<-upper.ggplot + theme(axis.text.x = element_blank(),
                                       axis.text.y = axis.text.y,
                                       axis.ticks.x = element_blank(),
                                       axis.ticks.y = axis.ticks,
                                       axis.title.x = element_blank(),
                                       axis.title.y = axis.title.y)
  } else {
    upper.ggplot<-upper.ggplot + theme(axis.text.x = axis.text.x,
                                       axis.text.y = element_blank(),
                                       axis.ticks.x = axis.ticks,
                                       axis.ticks.y = element_blank(),
                                       axis.title.x = axis.title.x,
                                       axis.title.y = element_blank())
  }
  
  if(is.null(add.ggplot.arg) == F){
    upper.ggplot <- upper.ggplot + add.ggplot.arg
  }
  
  ##Lower plot
  lower.ggplot<-ggplot() + 
    geom_tile(data = lower.melt, aes(x=Var1, y=Var2, fill = value),
              color = grid.col, size=grid.thickness) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          legend.title = lgd.title2,
          legend.text = lgd.text2,
          legend.position = 'bottom',
          legend.justification = c(1, 0),
          plot.margin = margin(-0.08, 0, 0, -0.1, "cm")) +
    scale_x_discrete(position = "top",expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0))+
    scale_fill_gradientn(colours = lgd.pal2,
                         breaks = seq(min_lower,max_lower,
                                      length.out = lgd.breaks2),
                         labels = round(seq(min_lower,max_lower,
                                            length.out = lgd.breaks2),
                                        lgd.round2),
                         guide=guide_colourbar(ticks=lgd.ticks2,
                                               nbin = lgd.nbin2,
                                               barheight=lgd.height2, label=T,
                                               barwidth=lgd.width2,
                                               raster = raster2,
                                               ticks.linewidth
                                               = lgd.ticks.linewidth2,
                                               frame.colour = lgd.frame.col,
                                               frame.linewidth
                                               = lgd.frame.linewidth2),
                         na.value = diag.col,
                         limits=c(min_lower,max_lower),
                         name = set.lgd2.title)
  
  #Color Sidebar
  groups<-levels(as.factor(annot.grps))
  col_table<-data.frame("Grps" = groups,"Cols" = annot.pal)

  dframe.annot<-data.frame("Samples" = colnames(data),"Groups" = annot.grps)
  if(annot.pos == "left"){
    correlation.order<-rev(correlation.order)
  }
  dframe.annot$Samples<-factor(dframe.annot$Samples,
                               levels = levels(dframe.annot$Samples)
                               [correlation.order])
  dframe.annot$Groups<-factor(dframe.annot$Groups, levels = mixedsort(groups))
  
  #Plot Color Sidebar
  col_sidebar<-ggplot(dframe.annot, aes(Samples,"Groups")) +
    # theme_minimal() +
    geom_tile(aes(fill = Groups)) +
    theme(legend.justification = "center",
          legend.position = c(0.5,0.4),
          legend.text=element_text(size= 12),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          panel.grid = element_blank(),
          plot.margin=unit(c(0,0,0,0),"cm")) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values=as.character(col_table$Cols))

  if(annot.pos == "top"){
    col_sidebar<-col_sidebar +
      theme(axis.text.x = axis.text.x,
            axis.text.y = element_blank(),
            axis.ticks.x = axis.ticks,
            axis.ticks.y = element_blank()) +
      scale_x_discrete(position = "top",expand = c(0,0)) +
      xlab(set.x.title)
    if(dendro.pos !="top"){
      col_sidebar<-col_sidebar +
        theme(axis.title.x = axis.title.x, axis.title.y = element_blank())
    } else {
      col_sidebar<-col_sidebar +
        theme(axis.title = element_blank())
    }
  } else {
    col_sidebar<-col_sidebar +
      theme(axis.text.x = element_blank(),
            axis.text.y = axis.text.y,
            axis.ticks.x = element_blank(),
            axis.ticks.y = axis.ticks) +
      scale_x_discrete(expand = c(0,0)) +
      coord_flip() +
      xlab(set.y.title)
    if(dendro.pos !="left"){
      col_sidebar<-col_sidebar +
        theme(axis.title.x = element_blank(), axis.title.y = axis.title.y)
    } else {
      col_sidebar<-col_sidebar +
        theme(axis.title = element_blank())
    }
  }
  
  #Remove all legends
  upper.ggplot.nolgd<-upper.ggplot + theme(legend.position = "none")
  lower.ggplot.nolgd<-lower.ggplot + theme(legend.position = "none")
  sidebar.nolgd<-col_sidebar + theme(legend.position = "none")
  #Create grob for lower matrix
  lower.grob<-ggplotGrob(lower.ggplot.nolgd)
  main_grob<-upper.ggplot.nolgd+annotation_custom(lower.grob)

  main_grob<-ggplotGrob(main_grob)
  sidebar_grob<-ggplotGrob(sidebar.nolgd)
  if(dendro.pos!='none'){ dendro_grob<-ggplotGrob(ddgr_seg) }

  upper.legend <- get.lgd(upper.ggplot)
  lower.legend <- get.lgd(lower.ggplot)
  sidebar.legend <- get.lgd(col_sidebar)
  
  if(annot.pos == "top"){
    #Since all object have the same width you can align them on the same column
    if(dendro.pos == "top"){
      #Get common width of a list of objects widths
      maxWidth = unit.pmax(main_grob$widths[2:5], sidebar_grob$widths[2:5],
                           dendro_grob$widths[2:5])
      #Assign back the maximum width to all object of the list
      main_grob$widths[2:5] <- as.list(maxWidth)
      sidebar_grob$widths[2:5] <- as.list(maxWidth)
      dendro_grob$widths[2:5] <- as.list(maxWidth)
      
      main_grob<-arrangeGrob(dendro_grob, sidebar_grob, main_grob, ncol=1,
                              heights = c(10+dendro.size,9+annot.size,40))
    } else {
      #Get common width of a list of objects widths
      maxWidth = unit.pmax(main_grob$widths[2:5], sidebar_grob$widths[2:5])
      #Assign back the maximum width to all object of the list
      main_grob$widths[2:5] <- as.list(maxWidth)
      sidebar_grob$widths[2:5] <- as.list(maxWidth)
      
      main_grob<-arrangeGrob(sidebar_grob, main_grob, ncol=1,
                              heights = c(9+annot.size,40))
    }

    right.legends<-arrangeGrob(upper.legend,sidebar.legend,nrow = 1)
    
  } else{ #Annotation on the left    
    #Since all object have the same width you can align them on the same column
    if(dendro.pos == "left"){
      #Get common width of a list of objects widths
      maxHeight = unit.pmax(main_grob$heights[3:8], sidebar_grob$heights[3:8],
                            dendro_grob$heights[3:8])
      #Assign back the maximum width to all object of the list
      main_grob$heights[3:8] <- as.list(maxHeight)
      sidebar_grob$heights[3:8] <- as.list(maxHeight)
      dendro_grob$heights[3:8] <- as.list(maxHeight)
      
      main_grob<-arrangeGrob(dendro_grob, sidebar_grob, main_grob, nrow=1,
                              widths = c(10+dendro.size,9+annot.size,40))
    } else {
      #Get common width of a list of objects widths
      maxHeight = unit.pmax(main_grob$heights[3:8], sidebar_grob$heights[3:8])
      #Assign back the maximum width to all object of the list
      main_grob$heights[3:8] <- as.list(maxHeight)
      sidebar_grob$heights[3:8] <- as.list(maxHeight)
      
      main_grob<-arrangeGrob(sidebar_grob, main_grob, nrow=1,
                              widths = c(9+annot.size,40))
    }
    right.legends<-arrangeGrob(upper.legend,sidebar.legend,nrow = 1)
  }
  
  #Plot Final Figure
  grid.arrange(main_grob, right.legends, lower.legend, ncol = 2, nrow = 2,
               heights = c(40,3), widths = c(20,6))
  
  return(list("upper.res" = upper.res, "lower.res" = lower.res))
}
