
#' Loads pre-defined palettes.
#'
#' @return Loads \code{character} vectors of colors to be used as palettes.
#' @author Yoann Pageaud.
#' @export
#' @examples load.palettes()

load.palettes <- function(){
  pal_sunset <<- c("red", "gold", "blue4")
  pal_westworld <<- c("sienna1", "lightgoldenrod", "skyblue3")
  pal_startrek <<- c("red", "goldenrod1", "dodgerblue")
  pal_margesimpson <<- c("lightgreen", "gold", "dodgerblue2")
  pal_eighties <<- c("darkviolet", "deeppink", "blue4")
  pal_insta <<- c("deeppink", "red", "goldenrod1")
}

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
#' @param verbose    A \code{logical} to specify wether the function should be
#'                   run on verbose mode (verbose = TRUE) or not
#'                   (Default: verbose = FALSE).
#' @return An error message if something goes wrong during annotations checks.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

check.annotations <- function(data, annot.grps, annot.pal, verbose = FALSE){
  #Groups checking
  if(is.list(annot.grps)){
    if(any(unlist(lapply(annot.grps, length))!= ncol(data))){
      stop("samples are not all assigned to a group.")
    } else{
      if(verbose){
        #Print groups values
        invisible(lapply(seq_along(annot.grps), function(i){
          cat(paste0(names(annot.grps)[i],": ", paste(unique(annot.grps[[i]]),
                                                      collapse = ", "), ".\n"))
        }))
      }
    }
  } else {
    stop(paste("'annot.grps' must be a named list.",
               "e.g. annot.grps = list('annotation_name' = annotation_vector)",
               sep = "\n"))
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
      if(length(annot.pal) == 0){ stop("'annot.pal' is empty.") }
      if(length(levels(as.factor(annot.grps[[i]]))) != length(annot.pal)){
        if(isTRUE(all.equal(target = annot.pal,
                            current = grDevices::rainbow(n = ncol(data))))){
          stop(paste("A specific palette must be defined in 'annot.pal' to",
                     "match the annotation provided."))
        } else {
          stop(paste0("The length of annotation '",names(annot.grps)[i],
                      "' levels do not match the length of the corresponding ",
                      "palette."))
        }
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

is.elt_blank <- function(arg){
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

get.lgd <- function(gg2.obj){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(gg2.obj))
  leg <- which(vapply(X = tmp$grobs, FUN = function(x) x$name,
                      FUN.VALUE = character(length = 1)) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' Creates a dendogram in ggplot2.
#'
#' @param df          A \code{data.frame} containing variables and values to be
#'                    used to create the dendrogram.
#' @param orientation A \code{character} specifying the orientation of the
#'                    dendrogram. Possible values are "top" and "left".
#' @param reverse.x   A \code{logical} specifying whether the X-axis should be
#'                    reversed (reverse.x = TRUE like in a correlation plot), or
#'                    kept unchanged (reverse.x = TRUE like in a heatmap), for a
#'                    dendrogram to be displayed on the left of a plot
#'                    (Default: reverse.x = FALSE).
#' @return A \code{gg} plot of the dendrogram.
#' @author Yoann Pageaud.
#' @export

ggdend <- function(df, orientation, reverse.x = FALSE) {
  ddplot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "transparent"),
      plot.background = ggplot2::element_rect(fill = "transparent"),
      axis.ticks.length = ggplot2::unit(0, "pt")) +
    ggplot2::expand_limits(x = c(0.5, max(df$x) + 0.5), y = 0)
  if(orientation == "top"){
    ddplot <- ddplot +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0, 0, 0), "cm")) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
  } else if(orientation == "left"){
    ddplot <- ddplot +
      ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0.1), "cm")) +
      ggplot2::scale_y_reverse(expand = c(0, 0)) +
      ggplot2::coord_flip()
    if(reverse.x){
      ddplot <- ddplot + ggplot2::scale_x_reverse(expand = c(0, 0))
    } else { ddplot <- ddplot + ggplot2::scale_x_continuous(expand = c(0, 0)) }
  } else { stop("dendrogram's orientation value not supported by ggdend().") }
  return(ddplot)
}

#' Draws a ggplot for a basic upper or lower triangle.
#'
#' @param melt.tri            A \code{data.frame} melted triangle containing
#'                            test values.
#' @param grid.col            A \code{character} specifying the color of the
#'                            grid.
#' @param grid.thickness      A \code{double} value for the thickness of the
#'                            grid.
#' @param lgd.title           An \code{element_text} object to setup the legend
#'                            title text.
#' @param lgd.text            An \code{element_text} object to setup the legend
#'                            labels text.
#' @param lgd.pal             A \code{character} vector of colors to use for the
#'                            palette of the triangle plot.
#' @param min_tri             A \code{numeric} defining the minimum limit for
#'                            the legend of the values in the triangle plot.
#' @param max_tri             A \code{numeric} defining the maximum limit for
#'                            the legend of the values in the triangle plot.
#' @param lgd.breaks          An \code{integer} defining the number of breaks
#'                            wanted in the legend.
#' @param lgd.round           An \code{integer} indicating the number of decimal
#'                            places to be used for the default legends labels.
#' @param lgd.ticks           A \code{logical} to specify wether ticks should be
#'                            diplayed on both axes.
#' @param lgd.nbin            An \code{integer} specifying the number of bins
#'                            for drawing both colorbars. A smoother colorbar
#'                            results from a larger value.
#' @param lgd.height          A \code{double} specifying the height of the
#'                            colorbar.
#' @param lgd.width           A \code{double} specifying the width of the
#'                            colorbar.
#' @param rasteri             A \code{logical} to specify whether or not colors
#'                            in the legend should be rasterized.
#' @param lgd.ticks.linewidth A \code{double} value to specify the thickness of
#'                            legends ticks.
#' @param lgd.frame.col       A \code{character} defining the color of the
#'                            frame of legends.
#' @param lgd.frame.linewidth A \code{double} defining the thickness of the
#'                            frame of both legends.
#' @param diag.col            A \code{character} defining the color of cells
#'                            with of the empty diagonal.
#' @param set.lgd.title       A \code{character} to used as the name of the
#'                            legend.
#' @return A \code{gg} object of a basic triangle plot (a 'geom_tile()').
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

basic.ggplot.tri<-function(melt.tri, grid.col, grid.thickness, lgd.title,
                           lgd.text, lgd.pal, min_tri, max_tri, lgd.breaks,
                           lgd.round, lgd.ticks, lgd.nbin, lgd.height,lgd.width,
                           rasteri, lgd.ticks.linewidth, lgd.frame.col,
                           lgd.frame.linewidth, diag.col, set.lgd.title){
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = melt.tri, ggplot2::aes(x = Var1, y = Var2, fill = value),
      color = grid.col, size=grid.thickness) +
    ggplot2::theme(
      legend.title = lgd.title, legend.text = lgd.text,
      legend.justification = c(1, 0), plot.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "transparent"),
      plot.background = ggplot2::element_rect(fill = "transparent")) +
    ggplot2::scale_x_discrete(position = "top",expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_fill_gradientn(
      colours = lgd.pal, breaks = seq(min_tri, max_tri,length.out = lgd.breaks),
      labels = round(seq(min_tri, max_tri, length.out = lgd.breaks), lgd.round),
      guide = ggplot2::guide_colourbar(
        ticks = lgd.ticks, nbin = lgd.nbin, barheight = lgd.height,
        label = TRUE, barwidth = lgd.width, raster = rasteri,
        ticks.linewidth = lgd.ticks.linewidth, frame.colour = lgd.frame.col,
        frame.linewidth = lgd.frame.linewidth), na.value = diag.col,
      limits = c(min_tri, max_tri), name = set.lgd.title)
}

#' Draws a ggplot2 of a basic sidebar.
#'
#' @param data      A \code{data.frame} with the column names 'Samples','.id'
#'                  and 'Groups' in this order.
#' @param palette   A \code{character} vector containing R colors like a
#'                  palette.
#' @param annot.sep A \code{double} specifying the width of the separation
#'                  spaces between annotations (Default: annot.sep = 0).
#' @param annot.cut A \code{double} specifying the width of cuts separating
#'                  annotation cells (Default: annot.cut = 0).
#' @param lgd.ncol  An \code{integer} specifying the number of columns to be
#'                  used to display a legend (Default: lgd.ncol = 1).
#' @return A \code{gg} object of the basic sidebar (a 'geom_tile()').
#' @author Yoann Pageaud.
#' @export

basic.sidebar <- function(
  data, palette, annot.sep = 0, annot.cut = 0, lgd.ncol = 1, facet = NULL){
  basic <- ggplot2::ggplot() +
    ggplot2::theme(legend.justification = c(0, 1),
                   legend.position = c(0, 1),
                   legend.spacing.y = ggplot2::unit(0.05, 'cm'),
                   # legend.text = element_text(size = 12),
                   # legend.title = element_text(size = 12),
                   # legend.margin = margin(-35,0,0,-35), #Seems to fit in grid
                   axis.text = ggplot2::element_text(size = 12),
                   panel.grid = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0, 0, 0, 0),
                   # strip.background = element_blank(),
                   # strip.text = element_blank()
    ) +
    ggplot2::scale_fill_manual(values = as.character(palette)) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = lgd.ncol, byrow = TRUE))

  if(!is.null(facet)){
    #Add faceting
    dt.facet <- data[.id == facet]
    dt.facet[, facet.annot := Groups]
    data <- merge(x = data, y = dt.facet[, c("Samples", "facet.annot"), ],
                  by = "Samples", all.x = TRUE)
    #Plot annotation bar with facet
    basic <- basic + ggplot2::geom_tile(data = data, mapping = ggplot2::aes(
      x = Samples, y = .id, fill = Groups, height = 1 - annot.sep,
      width = 1 - annot.cut)) +
      ggplot2::facet_grid(. ~ facet.annot, scales = "free", space = "free") +
      ggplot2::theme(
        panel.spacing = ggplot2::unit(0, "lines"),
        strip.background = ggplot2::element_rect(color = "black", size = 0.5))
  } else {
    basic <- basic + ggplot2::geom_tile(data = data, mapping = ggplot2::aes(
      x = Samples, y = .id, fill = Groups, height = 1 - annot.sep,
      width = 1 - annot.cut))
  }
  return(basic)
}

#' Creates a colored side annotation bars in ggplot2.
#'
#' @param sample.names A \code{character} vector of the labels to be used for
#'                     annotation.
#' @param annot.grps   A \code{list} of vectors of groups to which variables
#'                     belongs for the annotation sidebars. vectors lengths
#'                     have to match the number of variables.
#' @param annot.pal    A \code{vector} or a list of vectors containing colors
#'                     as characters for the annotation sidebars. The length of
#'                     vectors has to match the number of levels of vectors
#'                     listed in 'annot.grps'.
#'                     \itemize{
#'                      \item{If annot.pal is a list: its length must match the
#'                            length of the list provided to 'annot.grps'.}
#'                      \item{If annot.pal is a vector: make sure that the
#'                            levels content of annotations listed in
#'                            'annot.grps' is the same, and that no annotation
#'                            contains less or more levels than another one in
#'                            'annot.grps'.}
#'                     }
#' @param annot.pos    A \code{character} specifying the position of the
#'                     annotation sidebar.\cr Possible values are: 'top',
#'                     'left' or 'both'.
#' @param annot.sep    A \code{double} specifying the width of the separation
#'                     spaces between annotations (Default: annot.sep = 0).
#' @param annot.cut    A \code{double} specifying the width of cuts separating
#'                     annotation cells (Default: annot.cut = 0).
#' @param cor.order    An \code{integer} vector to be used to reorder the labels
#'                     in 'sample.names'.
#' @param merge.lgd    A \code{logical} to specify whether annotation legends
#'                     should be merged (annot.lgd.merge = TRUE) or remain
#'                     separated (annot.lgd.merge = FALSE)
#'                     (Default: annot.lgd.merge = FALSE).
#' @param right        A \code{logical} to specify that Y-axis should be
#'                     displayed on the right side of the plot
#'                     (Default: right = FALSE).
#' @param lgd.name      A \code{character} to specify a title to the legend of
#'                     the plot, only if 'merge.lgd' = TRUE
#'                     (Default: lgd.name = "Legends").
#' @param lgd.title    An \code{element_text} object to setup legend titles
#'                     (Default: lgd.title = ggplot2::element_blank()).
#' @param lgd.text     An \code{element_text} object to setup legend labels
#'                     (Default: lgd.text = ggplot2::element_blank()).
#' @param lgd.ncol     An \code{integer} specifying the number of columns to be
#'                     used to display a legend (Default: lgd.ncol = 1).
#' @param axis.text.x  An \code{element_text} object to setup X axis text
#'                     (Default: axis.text.x = ggplot2::element_text(size = 12)).
#' @param axis.text.y  An \code{element_text} object to setup Y axis text
#'                     (Default: axis.text.y = ggplot2::element_text(size = 12)).
#' @param axis.ticks.x An \code{element_line} object to setup X axis ticks.
#' @param axis.ticks.y An \code{element_line} object to setup Y axis ticks.
#' @param axis.title.x An \code{element_text} object to setup X axis title.
#'                     \itemize{Exceptions:
#'                      \item{If annot.pos == 'left', then axis.title.x =
#'                      ggplot2::element_blank().}
#'                      \item{If annot.pos == dendro.pos == 'top', then
#'                      axis.title.x = ggplot2::element_blank().}
#'                     }
#' @param axis.title.y An \code{element_text} object to setup Y axis title.
#'                     \itemize{Exceptions:
#'                      \item{If annot.pos == 'top', then axis.title.y =
#'                      ggplot2::element_blank().}
#'                      \item{If annot.pos == dendro.pos == 'left', then
#'                      axis.title.y = ggplot2::element_blank().}
#'                     }
#' @param set.x.title  A \code{character}to be used as the title for the X axis.
#' @param set.y.title  A \code{character}to be used as the title for the Y axis.
#' @param dendro.pos   A \code{character} specifying the position of the
#'                     dendrogram (Supported: dendro.pos = c('top', 'left')).
#' @return A \code{list} of length 2:
#'         \itemize{
#'          \item{'sidebar' contains the colored sidebar plot.}
#'          \item{'legends' lists legends of the matching sidebar as \code{gg}
#'                objects.}
#'         }
#' @author Yoann Pageaud.
#' @export

plot.col.sidebar <- function(
  #TODO: change name of parameter cor.order to see if it can be remove
  sample.names, annot.grps, annot.pal, annot.pos = 0, annot.sep = 0, annot.cut,
  cor.order, merge.lgd = FALSE, right = FALSE, lgd.name = "Legends",
  lgd.title = ggplot2::element_blank(), lgd.text = ggplot2::element_blank(),
  lgd.ncol = 1, axis.text.x = ggplot2::element_text(size = 12),
  axis.text.y = ggplot2::element_text(size = 12), axis.ticks.x, axis.ticks.y,
  axis.title.x, axis.title.y, set.x.title, set.y.title, dendro.pos,
  facet = NULL){

  #Create list of groups in their original order
  origin.grps <- lapply(X = annot.grps, FUN = function(i){
    if(is.factor(i)){ levels(i) } else { levels(as.factor(i)) }
  })
  #Update levels following the new order of the annotation
  groups <- lapply(X = lapply(X = annot.grps, FUN = function(i){
    factor(x = i, levels =  unique(i))}), FUN = levels)
  #Create list of color tables
  #TODO: try to merge with check.annotations()
  if(is.list(annot.pal)) { #If a list of palettes is provided
    if(length(groups) == length(annot.pal)){ #if annotations match palettes
      #Map groups to palettes
      ls.df.grp.pal <- Map(data.frame, "Grps" = origin.grps, "Cols" = annot.pal,
                           stringsAsFactors = FALSE)
      col_table <- lapply(seq_along(groups), function(i){
        if(length(groups[[i]]) == length(annot.pal[[i]])){
          #if groups match colors
          ls.df.grp.pal[[i]][match(groups[[i]], ls.df.grp.pal[[i]]$Grps), ]
          # data.frame("Grps"=groups[[i]],"Cols"=annot.pal[[i]])
        } else {
          stop(paste0(
            "The length of annotation '", names(groups)[i],
            "' levels does not match the length of the corresponding palette."))
        }
      })
    } else {
      stop("The number of annotations does not match the number of palettes provided.")
    }
  } else if(is.vector(annot.pal)){ #if a single palette is provided
    #Map groups to the same palette
    ls.df.grp.pal <- lapply(X = origin.grps, FUN = function(grp){
      data.frame("Grps" = grp, "Cols" = annot.pal, stringsAsFactors = FALSE)
    })
    col_table <- lapply(seq_along(groups), function(i){
      if(length(groups[[i]]) == length(annot.pal)){ #if groups match colors
        ls.df.grp.pal[[i]][match(groups[[i]], ls.df.grp.pal[[i]]$Grps), ]
        # data.frame("Grps" = groups[[i]], "Cols" = annot.pal)
      } else {
        stop(paste0("The length of annotation '", names(groups)[i],
                    "' levels do not match the length of the corresponding ",
                    "palette."))}
    })
  } else { #If not a list or a vector
    stop("Unknown type for 'annot.pal'. 'annot.pal' should be either a list or a vector.")
  }

  #Create list of annotation data.frames
  dframe.annot <- lapply(annot.grps, function(i){
    data.frame("Samples" = sample.names, "Groups" = i)
  })
  #Order samples following the correlation order provided
  # and categories by alphabetical order
  if(annot.pos == "left"){ cor.order <- rev(cor.order) }
  dframe.annot <- lapply(dframe.annot, function(i){
    i[["Samples"]] <- factor(i[["Samples"]], levels = i[["Samples"]][cor.order])
    i[["Groups"]] <- factor(i[["Groups"]], levels = unique(i[["Groups"]]))
    i
  })
  #Rbind all annotations
  dframe.annot <- data.table::rbindlist(dframe.annot, idcol = TRUE)
  #Convert .id as factors
  dframe.annot$.id <- factor(
    x = dframe.annot$.id, levels = unique(dframe.annot$.id))
  # if(!split.annot){
  if(annot.pos == "top"){ #Change order of levels
    dframe.annot$.id <- factor(
      x = dframe.annot$.id, levels = rev(levels(dframe.annot$.id)))
  }
  # }

  #Check color tables
  col_table <- lapply(X = col_table, FUN = function(tbl){
    if(any(duplicated(tbl$Cols))){ #Check palette consistency
      stop("1 color in a palette has been associated to more than 1 group.")
    }
    if(any(duplicated(tbl$Grps))){ #Check annotation consistency
      warning("Duplicated group name provided. Removing duplicated...")
      tbl <- tbl[!duplicated(Grps)]
    } else { tbl }
  })
  col_table <- data.table::rbindlist(col_table, idcol = TRUE)
  if(is.vector(annot.pal)){
    if(any(duplicated(col_table$Cols))){ #Check palette consistency
      # warning(paste(
      #   "Some colors have been assigned to more than 1 group in the palette.",
      #   "Removing duplicated occurences..."))
      col_table <- col_table[!duplicated(x = Cols)]
    }
  }
  # if(any(duplicated(col_table$Grps))){
  #   warning("Duplicated group name provided. Removing duplicated...")
  #   col_table <- col_table[!duplicated(Grps)]
  # }
  #Plot color sidebars
  col_sidebar <- BiocompR::basic.sidebar(
    data = dframe.annot, palette = col_table$Cols, annot.sep = annot.sep,
    annot.cut = annot.cut, lgd.ncol = lgd.ncol, facet = facet)
  #Add legend parameters if some
  col_sidebar <- col_sidebar + ggplot2::theme(legend.title = lgd.title)
  #Modify base plot following its position
  if(annot.pos == "top"){
    col_sidebar <- col_sidebar +
      ggplot2::theme(axis.text.x.top = axis.text.x, axis.text.y = axis.text.y,
                     axis.ticks.x = axis.ticks.x, axis.ticks.y = axis.ticks.y) +
      ggplot2::scale_x_discrete(expand = c(0, 0), position = "top") +
      ggplot2::xlab(set.x.title)
    if(right){
      col_sidebar <- col_sidebar +
        ggplot2::scale_y_discrete(position = 'right', expand = c(0, 0))
    } else {
      col_sidebar <- col_sidebar + ggplot2::scale_y_discrete(expand = c(0, 0))
    }
    if(dendro.pos != "top"){
      col_sidebar <- col_sidebar +
        ggplot2::theme(axis.title.x = axis.title.x,
                       axis.title.y = ggplot2::element_blank())
    } else {
      col_sidebar <- col_sidebar +
        ggplot2::theme(axis.title = ggplot2::element_blank())
    }
    # if(split.annot){
    #   col_sidebar <- col_sidebar +
    #     facet_grid(.id ~ ., scales = "free", space = "free_y")
    # }
  } else if(annot.pos == "left"){
    col_sidebar <- col_sidebar +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = axis.text.y, axis.ticks.y = axis.ticks.y,
                     axis.text.x.top = axis.text.x) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0), position = "right") +
      ggplot2::xlab(set.y.title)
    if(dendro.pos !="left"){
      col_sidebar <- col_sidebar +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = axis.title.y)
    } else {
      col_sidebar <- col_sidebar +
        ggplot2::theme(axis.title = ggplot2::element_blank())
    }
    # if(split.annot){
    #   stop("A geom_tile vertically faceted in ggplot2_3.2.0 does not support heights redimensioning after being converted into a grob.")
    # }
  }
  if(merge.lgd){ # Do not split legends
    sidebar.lgd <- list(
      BiocompR::get.lgd(col_sidebar + ggplot2::labs(fill = lgd.name)))
  } else { # Split legends and return a list of legends
    # if(!split.annot){
    if(annot.pos == "top"){
      dframe.annot$.id <- factor(
        dframe.annot$.id, levels = rev(levels(dframe.annot$.id)))
    }
    # }

    #Generate separate legends if more than 1 palette available
    # or if only 1 annotation is used
    if((is.list(annot.pal) & length(annot.pal) > 1) |
       length(levels(dframe.annot$.id)) == 1){
      #Get all legends separately
      sidebar.lgd <- lapply(seq_along(levels(dframe.annot$.id)), function(i){
        BiocompR::get.lgd(
          BiocompR::basic.sidebar(
            data = dframe.annot[.id == levels(dframe.annot$.id)[i]],
            palette = col_table[.id == i]$Cols, lgd.ncol = lgd.ncol) +
            ggplot2::theme(legend.title = lgd.title, legend.text = lgd.text) +
            ggplot2::labs(fill = levels(dframe.annot$.id)[i])
        )
      })
    } else {
      stop("Cannot generate separated legends if only one annotation palette is given.")
    }
  }
  return(list("sidebar" = col_sidebar +
                ggplot2::theme(legend.position = "none"),
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
#'         grid::unit.pmax() function.
#' @author Yoann Pageaud.
#' @export
#' @references \href{https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page}{Align two plots on a page - Mara Averick}

resize.grobs <- function(ls.grobs, dimensions, start.unit, end.unit){
  #Get dimension units from the list of grobs to redimension
  ls.dim <- lapply(X = ls.grobs, FUN = function(i){
    if(length(i[[dimensions]]) < end.unit){
      i[[dimensions]][start.unit:length(i[[dimensions]])]
    } else { i[[dimensions]][start.unit:end.unit] }
  })
  #Calculate maximum of all unit objects including the main grob.
  max.dim <- eval(parse(
    text = paste("grid::unit.pmax(",paste(paste(
      rep("ls.dim[[",length(ls.dim)), seq(length(ls.dim)),"]]", sep = ""),
      collapse = ", "), ")", sep = "")))
  #Apply changes to grobs dimensions
  ls.grobs <- lapply(X = ls.grobs, FUN = function(i){
    i[[dimensions]][start.unit:end.unit] <- as.list(max.dim)
    i
  })
  return(ls.grobs)
}

#' Resizes heights or widths of a grob based on the dimensions of another grob.
#'
#' @param grob1      A \code{grob} to be modified.
#' @param grob2      A \code{grob} to be used as reference for dimensions
#'                   modifications.
#' @param dimensions A \code{character} specifying the type of dimensions to
#'                   resize, either 'heights' or 'widths'.
#' @param positions  An \code{integer} vector specifying indexes of the
#'                   dimensions to change.
#' @return A \code{grob} for which the dimensions have been modified.
#' @author Yoann Pageaud.
#' @export

resize.grob.oneway <- function(grob1, grob2, dimensions, positions){
  #Get dimension units from the list of grobs to redimension
  if(max(positions) > length(grob2[[dimensions]]) | min(positions) < 1){
    stop("positions out of range in grob2.")
  }
  if(max(positions) > length(grob1[[dimensions]])){
    warning("some positions out of range in grob1 ignored.")
    positions <- positions[positions <= length(grob1[[dimensions]])]
  }
  grob1[[dimensions]][positions] <- grob2[[dimensions]][positions]
  return(grob1)
}

#' Stack grobs legends vertically in separate spaces of specific heights.
#'
#' @param grobs.list A \code{list} of legends as grid objects.
#' @param annot.grps A \code{factor} list mapped to the legends.
#' @param height.lgds.space An \code{integer} specifying the height of the total
#'                          space which is supposed to gather all legends in
#'                          grobs.list.
#' @return A \code{grob} containing all legends with their respective heights.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

stack.grobs.legends <- function(grobs.list, annot.grps, height.lgds.space){
  size_lgd <- unlist(lapply(X = annot.grps, FUN = function(i){
    length(unique(i)) + 1
  }))
  #Add a ending space in the legend to stack them to the top
  grobs.list <- c(grobs.list, list(grid::textGrob("")))
  #Calculate height of the ending space
  hghts <- c(size_lgd, height.lgds.space-sum(size_lgd))
  sidebar_legend <- gridExtra::arrangeGrob(grobs = grobs.list, heights = hghts)
  return(sidebar_legend)
}

#' Rasterize a gg plot into a raster grob.
#'
#' @param gg.plot A \code{gg} plot to be rasterized.
#' @param filter  A \code{character} to be used as a filter for ggplot
#'                rasterization. The list of the supported rasterization
#'                filters is available in magick::filter_types()
#'                (Default: raster = "Lanczos"). Warning: Be aware that
#'                rasterization may take several minutes to process the ggplot.
#' @return A \code{grob} of the rasterized ggplot.
#' @author Yoann Pageaud.
#' @export

raster.ggplot.to.grob <- function(gg.plot, filter = "Lanczos"){
  #Catch heatmap in magick::image_graph()
  fig <- magick::image_graph(width = 2160, height = 2160, res = 96)
  print(gg.plot)
  grDevices::dev.off()
  rastered <- magick::image_resize(
    image = fig, geometry = "1080x1080", filter = filter)
  #Create raster grob
  raster.grob <- grid::rasterGrob(
    rastered, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
    interpolate = TRUE)
  #Return grob annotation
  return(raster.grob)
}

#' This is a special geom intended for use as static annotations derived from
#' ggplot2::annotation_custom() matching a specific panel on a faceted ggplot.
#'
#' @param grob      A \code{grob} to display.
#' @param xmin,xmax x location (in data coordinates) giving horizontal location
#'                  of raster.
#' @param ymin,ymax y location (in data coordinates) giving vertical location
#'                  of raster.
#' @param data      A subset of a \code{data.table} matching the panel where the
#'                  grob annotation should be displayed.
#' @return A \code{type} object returned description.
#' @export
#' @keywords internal

annotation_custom2 <- function(
  grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
  ggplot2::layer(
    data = data, stat = StatIdentity, position = PositionIdentity,
    geom = ggplot2:::GeomCustomAnn, inherit.aes = TRUE, params = list(
      grob = grob, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
}
