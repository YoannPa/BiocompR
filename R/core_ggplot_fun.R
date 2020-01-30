
##IMPORTS
Imports = c("ggplot2","data.table")
lapply(Imports, library, character.only = T)

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
#' @return An error message if something goes wrong during annotations checks.
#' @author Yoann Pageaud.
#' @export
#' @references
#' @keywords internal

check.annotations<-function(data, annot.grps, annot.pal){
  #Groups checking
  if(any(unlist(lapply(annot.grps,length))!= ncol(data))){
    stop("samples are not all assigned to a group.")
  } else{#Print groups values
    invisible(lapply(seq_along(annot.grps), function(i){
      cat(paste0(names(annot.grps)[i],": ",paste(unique(annot.grps[[i]]),
                                                 collapse=", "),".\n"))
    }))
  }
  #Color checking
  if(is.list(annot.pal)){
    if(length(annot.grps) == length(annot.pal)){
      invisible(lapply(seq_along(annot.pal), function(i){
        if(length(annot.pal[[i]])!=length(levels(as.factor(annot.grps[[i]])))){
          stop(paste0("The length of annotation '",names(annot.grps)[i],
                      "' levels do not match the length of the corresponding ",
                      "palette."))
        }
      }))
    } else {
      stop("The number of palettes does not match the number of annotations provided.")
    }
  } else if(is.vector(annot.pal)){ #if a single palette is provided
    invisible(lapply(seq_along(annot.grps), function(i){
      if(length(levels(as.factor(annot.grps[[i]]))) != length(annot.pal)){
        stop(paste0("The length of annotation '",names(annot.grps)[i],
                    "' levels do not match the length of the corresponding ",
                    "palette."))
      }
    }))
  } else { #If not a list or a vector
    stop("Unknown type for 'annot.pal'. 'annot.pal' should be either a list or a vector.")
  }
}

#' Checks if a list's attributes has for class 'element_blank'.
#'
#' @param arg A \code{list}.
#' @return A \code{logical}.
#' @author Yoann Pageaud.
#' @export
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
#' @export
#' @keywords internal
#' @references \href{https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs}{Share a legend between two ggplot2 graphs - Mara Averick}

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
#' @param plot.type   A \code{character} specifying whether the plot is a
#'                    correlation plot ('corrplot') or a heatmap ('heatmap').
#'                    This parameter is used to know if a dendrogram displayed
#'                    on the left of a plot should be vertically reversed (like
#'                    in a correlation plot) or not (like in a heatmap).
#' @return A \code{gg} plot of the dendrogram.
#' @author Yoann Pageaud.
#' @export

ggdend <- function(df, orientation, plot.type) {
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
      coord_flip()
    if(plot.type == "corrplot"){
      ddplot <- ddplot + scale_x_reverse(expand = c(0,0))
    } else if(plot.type == "heatmap"){
      ddplot <- ddplot + scale_x_continuous(expand = c(0,0))
    } else {
      stop("Unknown plot.type. Supported plot.type: c('corrplot', 'heatmap').")
    }
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
#' @export
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
#' @export

basic.sidebar<-function(data, palette){
  ggplot(data = data) +
    geom_tile(aes(x = Samples, y = .id, fill = Groups)) +
    theme(legend.justification = 'left',
          # legend.position = c(0.5,0.5),
          legend.text=element_text(size= 12),
          legend.title = element_text(size = 12),
          # legend.margin = margin(-35,0,0,-35), #Seems to fit in grid
          axis.text = element_text(size = 12),
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
#' @export

#TODO: Write documentation!
plot.col.sidebar<-function(
  sample.names, annot.grps, annot.pal, annot.pos, cor.order, split.annot = TRUE,
  merge.lgd = FALSE, right = FALSE, lgd.lab = "Legends", lgd.title = NULL,
  axis.text.x, axis.text.y, axis.ticks.x, axis.ticks.y, axis.title.x,
  axis.title.y, set.x.title, set.y.title, dendro.pos){
  #Create list of groups
  groups<-lapply(lapply(annot.grps,as.factor),levels)

  #Create list of color tables
  #TODO: try to merge with check.annotations()
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
  } else if(is.vector(annot.pal)){ #if a single palette is provided
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

  if(!is.null(lgd.title)){ #Add legend parameters if some
    col_sidebar<- col_sidebar + theme(legend.title = lgd.title)
  }
  #Modify base plot following its position
  if(annot.pos == "top"){
    col_sidebar<-col_sidebar +
      theme(axis.text.x.top = axis.text.x, axis.text.y = axis.text.y,
            axis.ticks.x = axis.ticks.x, axis.ticks.y = axis.ticks.y) +
      scale_x_discrete(expand = c(0,0),position = "top") + xlab(set.x.title)
    if(right){
      col_sidebar <- col_sidebar +
        scale_y_discrete(position = 'right', expand = c(0,0))
    } else { col_sidebar <- col_sidebar + scale_y_discrete(expand = c(0,0)) }
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
      theme(axis.text.y = axis.text.y, axis.ticks.y = axis.ticks.y,
            axis.text.x.top = axis.text.x) +
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
    sidebar.lgd<-list(get.lgd(col_sidebar + labs(fill = lgd.lab)))
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
#'         unit.pmax() function.
#' @author Yoann Pageaud.
#' @export
#' @references \href{https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page}{Align two plots on a page - Mara Averick}

resize.grobs<-function(ls.grobs, dimensions, start.unit, end.unit){
  #Get dimension units from the list of grobs to redimension
  ls.dim<-lapply(X = ls.grobs, FUN = function(i){
    i[[dimensions]][start.unit:end.unit]
  })
  #Calculate maxilum of all unit objects including the main grob.
  max.dim <-eval(parse(
    text = paste("unit.pmax(",paste(paste(
      rep("ls.dim[[",length(ls.dim)), seq(length(ls.dim)),"]]", sep = ""),
      collapse = ", "), ")", sep = "")))
  #Apply changes to grobs dimensions
  ls.grobs<-lapply(X = ls.grobs, FUN = function(i){
    i[[dimensions]][start.unit:end.unit]<-as.list(max.dim)
    i
  })
  return(ls.grobs)
}
