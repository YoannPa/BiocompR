
##IMPORTS
Imports = c("ggplot2","data.table")
lapply(Imports, library, character.only = T)

#' Checks if a list's attributes has for class 'element_blank'.
#'
#' @param arg A \code{list}.
#' @return A \code{logical}.
#' @author Yoann Pageaud.
#' @keywords internal

is.elt_blank<-function(arg){
  bool<-attributes(arg)$class[1] == "element_blank"
  return(bool)
}

#' Extracts legend from a ggplot2 object.
#'
#' @param gg2.obj  A \code{gg} object with legends.
#' @return A \code{gg} object only containing the legends of the plot.
#' @author Yoann Pageaud.
#' @keywords internal

get.lgd<-function(gg2.obj){
  tmp <- ggplot_gtable(ggplot_build(gg2.obj))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' Creates a dendogram in ggplot2.
#'
#' @param df          A \code{data.frame} containing variables and values to be
#'                    used to create the dendrogram.
#' @param orientation A \code{character} specifying the orientation of the
#'                    dendrogram. Possible values are "top" and "left".
#' @return A \code{gg} plot of the dendrogram.
#' @author Yoann Pageaud.

ggdend <- function(df, orientation) {
  ddplot<- ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          axis.ticks.length = unit(0,"pt")) +
    expand_limits(x = c(0.5,max(df$x)+0.5), y = 0)
  if(orientation == "top"){
    ddplot <- ddplot +
      scale_x_continuous(expand = c(0,0)) +
      theme(plot.margin=unit(c(0.1,0,0,0),"cm")) +
      scale_y_continuous(expand = c(0, 0))
  } else if(orientation == "left"){
    ddplot <- ddplot +
      theme(plot.margin=unit(c(0,0,0,0.1),"cm")) +
      scale_y_reverse(expand = c(0,0)) +
      scale_x_reverse(expand = c(0,0)) +
      coord_flip()
  } else {stop("dendrogram's orientation value not supported by ggdend().")}
  return(ddplot)
}

#' Draws a ggplot for a basic upper or lower triangle.
#'
#' @param melt.tri       A \code{data.frame} melted triangle containing test
#'                       values.
#' @param grid.col       A \code{character} specifying the color of the grid.
#' @param grid.thickness A \code{double} value for the thickness of the grid.
#' @return A \code{gg} object of a basic triangle plot (a 'geom_tile()').
#' @author Yoann Pageaud.
#' @keywords internal

basic.ggplot.tri<-function(melt.tri, grid.col, grid.thickness, lgd.title,
                           lgd.text, lgd.pal, min_tri, max_tri, lgd.breaks,
                           lgd.round, lgd.ticks, lgd.nbin, lgd.height,lgd.width,
                           rasteri, lgd.ticks.linewidth, lgd.frame.col,
                           lgd.frame.linewidth, diag.col, set.lgd.title){
  ggplot() +
    geom_tile(data = melt.tri, aes(x=Var1, y=Var2, fill = value),
              color = grid.col, size=grid.thickness) +
    theme(legend.title = lgd.title,
          legend.text = lgd.text,
          legend.justification = c(1, 0),
          plot.margin = margin(0, 0, 0, 0),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent")) +
    scale_x_discrete(position = "top",expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradientn(colours = lgd.pal,
                         breaks = seq(min_tri,max_tri,length.out = lgd.breaks),
                         labels = round(seq(
                           min_tri,max_tri,length.out = lgd.breaks),lgd.round),
                         guide=guide_colourbar(
                           ticks=lgd.ticks, nbin=lgd.nbin, barheight=lgd.height,
                           label=T, barwidth=lgd.width, raster = rasteri,
                           ticks.linewidth = lgd.ticks.linewidth,
                           frame.colour = lgd.frame.col,
                           frame.linewidth = lgd.frame.linewidth),
                         na.value = diag.col, limits=c(min_tri,max_tri),
                         name = set.lgd.title)
}

#' Draws a ggplot2 of a basic sidebar.
#'
#' @param data    A \code{data.frame} with the column names 'Samples','.id' and
#'                'Groups' in this order.
#' @param palette A \code{character} vector containing R colors like a palette.
#' @return A \code{gg} object of the basic sidebar (a 'geom_tile()').
#' @author Yoann Pageaud.

basic.sidebar<-function(data, palette){
  ggplot(data = data) +
    geom_tile(aes(x = Samples, y = .id, fill = Groups)) +
    theme(legend.justification = c(0,1),
          legend.position = c(0.5,0.5),
          legend.text=element_text(size= 12),
          legend.title = element_text(size = 12, face = "bold"),
          legend.margin = margin(-35,0,0,-35), #Seems to fit in grid
          axis.text = element_text(size = 12, face = "bold"),
          panel.grid = element_blank(),
          plot.margin = margin(0,0,0,0),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    scale_fill_manual(values=as.character(palette))
}

#' Creates a colored side annotation bars in ggplot2.
#'
#' @param sample.names A \code{type} parameter description.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.

#TODO: Write documentation!
plot.col.sidebar<-function(
  sample.names, annot.grps, annot.pal, annot.pos, cor.order, split.annot = TRUE,
  merge.lgd = FALSE, lgd.title = "Legends", axis.text.x, axis.text.y,
  axis.ticks, axis.title.x, axis.title.y, set.x.title, set.y.title, dendro.pos){
  #Create list of groups
  groups<-lapply(lapply(annot.grps,as.factor),levels)
  #Create list of color tables
  if(is.list(annot.pal)) { #If a list of palettes is provided
    if(length(groups) == length(annot.pal)){ #if annotations match palettes
      col_table<-lapply(seq_along(groups), function(i){
        if(length(groups[[i]]) == length(annot.pal[[i]])){
          #if groups match colors
          data.frame("Grps"=groups[[i]],"Cols"=annot.pal[[i]])
        } else {
          stop(paste0("The length of annotation '",names(groups)[i],
                      "' levels do not match the length of the corresponding ",
                      "palette."))}
      })
    } else {
      stop("The number of annotations does not match the number of palettes provided.")
    }
  } else {
    if(is.vector(annot.pal)){ #if a single palette is provided
      col_table<-lapply(seq_along(groups), function(i){
        if(length(groups[[i]]) == length(annot.pal)){ #if groups match colors
          data.frame("Grps"=groups[[i]],"Cols"=annot.pal)
        } else {
          stop(paste0("The length of annotation '",names(groups)[i],
                      "' levels do not match the length of the corresponding ",
                      "palette."))}
      })
    } else { #If not a list or a vector
      stop("Unknown type for 'annot.pal'. 'annot.pal' should be either a list or a vector.")
    }
  }
  #Create list of annotation data.frames
  dframe.annot<-lapply(annot.grps, function(i){
    data.frame("Samples" = sample.names,"Groups" = i)
  })
  #Order samples following the correlation order provided
  # and categories by alphabetical order
  if(annot.pos == "left"){ cor.order<-rev(cor.order)}
  dframe.annot<-lapply(dframe.annot, function(i){
    i[["Samples"]]<-factor(
      i[["Samples"]],levels = levels(i[["Samples"]])[cor.order])
    i[["Groups"]]<-factor(
      i[["Groups"]], levels = sort(unique(i[["Groups"]])))
    i
  })
  #Rbind all annotations
  dframe.annot<-rbindlist(dframe.annot,idcol = TRUE)
  dframe.annot$.id<-as.factor(dframe.annot$.id)
  if(!split.annot){
    if(annot.pos == "top"){
      dframe.annot$.id<-factor(dframe.annot$.id,
                               levels = rev(levels(dframe.annot$.id)))
    }
  }
  col_table<-rbindlist(col_table,idcol = TRUE)
  if(any(duplicated(col_table$Grps))){
    stop("Duplicated group name provided. Make sure 'annot.grps' does not contain duplicated group names in the list.")
  }
  #Plot color sidebars
  col_sidebar<-basic.sidebar(data = dframe.annot, palette = col_table$Cols)
  #Modify base plot following its position
  if(annot.pos == "top"){
    col_sidebar<-col_sidebar +
      theme(axis.text.x.top = axis.text.x,axis.ticks.x = axis.ticks) +
      scale_y_discrete(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0),position = "top") + xlab(set.x.title)
    if(dendro.pos !="top"){
      col_sidebar<-col_sidebar +
        theme(axis.title.x = axis.title.x, axis.title.y = element_blank())
    } else { col_sidebar<-col_sidebar + theme(axis.title = element_blank()) }
    if(split.annot){
      col_sidebar<-col_sidebar +
        facet_grid(.id ~ ., scales = "free", space = "free_y")
    }
  } else if(annot.pos == "left"){
    col_sidebar<-col_sidebar +
      coord_flip() +
      theme(axis.text.y = axis.text.y, axis.ticks.y = axis.ticks,
            axis.text.x.top = element_text(angle = 90,hjust = 0, vjust = 0.5)) +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0), position = "right") +
      xlab(set.y.title)
    if(dendro.pos !="left"){
      col_sidebar<-col_sidebar +
        theme(axis.title.x = element_blank(), axis.title.y = axis.title.y)
    } else { col_sidebar<-col_sidebar + theme(axis.title = element_blank()) }
    if(split.annot){
      stop("A geom_tile vertically faceted in ggplot2_3.2.0 does not support heights redimensioning after being converted into a grob.")
      # col_sidebar<-col_sidebar +
      #   facet_grid(. ~ .id, scales = "free", space = "free_x")
    }
  }

  if(merge.lgd){# Do not split legends
    sidebar.lgd<-list(get.lgd(col_sidebar + labs(fill = lgd.title)))
  } else {# Split legends and return a list of legends
    if(!split.annot){
      if(annot.pos == "top"){
        dframe.annot$.id<-factor(
          dframe.annot$.id, levels = rev(levels(dframe.annot$.id)))
      }
    }
    sidebar.lgd<-lapply(seq_along(levels(dframe.annot$.id)), function(i){
      get.lgd(
        basic.sidebar(data = dframe.annot[.id == levels(dframe.annot$.id)[i]],
                      palette = col_table[.id == i]$Cols) +
          labs(fill = levels(dframe.annot$.id)[i])
      )
    })
  }
  return(list("sidebar" = col_sidebar + theme(legend.position = "none"),
              "legends" = sidebar.lgd))
}
