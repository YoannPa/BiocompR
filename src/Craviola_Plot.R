##IMPORTS
library(ggplot2)
library(reshape2)
source("/media/yoann/3CAD87DD271F7BEC/PhD_2018-2021/Hematopoiesis_DNMT1_Hypomorph/analysis/Methylation_Violin_Plot/bin_polygons.R")

##DUMMY DATA
dummy_data = data.frame(
  Samples=c(rep('Sample1', 1000), rep('Sample2', 1000),
            rep('Sample3', 1000), rep('Sample4', 1000),
            rep('Sample5', 1000), rep('Sample6', 1000)),
  Cell.types=c(rep('a', 2000), rep('b', 2000), rep('c', 2000)),
  Genotypes=c(rep('i', 1000), rep('j', 1000), rep('i', 1000),rep('j', 1000),
          rep('i', 1000), rep('j', 1000)),
  Beta=c(rnorm(1000), rnorm(1000, 0.5), rnorm(1000, 1), rnorm(1000, 1.5),
         rnorm(1000,-1), rnorm(1000, -0.5)),
  Coverage = c(rep(x = c(60,50,40,30,20,10,30,40,50,60), each=100),
               rep(x = c(60,50,40,30,20,10,30,40,50,60), each=100),
               rep(x = c(60,50,40,30,20,10,30,40,50,60), each=100),
               rep(x = c(60,50,40,30,20,10,30,40,50,60), each=100),
               rep(x = c(60,50,40,30,20,10,30,40,50,60), each=100),
               rep(x = c(60,50,40,30,20,10,30,40,50,60), each=100))
)

dummy_data.bis = data.frame(
  Samples=c(rep('Sample1', 1000), rep('Sample2', 1000),
            rep('Sample3', 1000), rep('Sample4', 1000),
            rep('Sample5', 1000), rep('Sample6', 1000)),
  Cell.types=c(rep('a', 2000), rep('b', 2000), rep('c', 2000)),
  Genotypes=c(rep('i', 1000), rep('j', 1000), rep('i', 1000),rep('j', 1000),
              rep('i', 1000), rep('j', 1000)),
  Beta=c(rnorm(1000), rnorm(1000, 0.5), rnorm(1000, 1), rnorm(1000, 1.5),
         rnorm(1000,-1), rnorm(1000, -0.5)))

##FUNCTIONS

# ggcraviola ###################################################################

#' Read Chromosomes Methylation from a specific sample.
#'
#' @param data           A \code{data.frame} of 5 columns:
#'                       column 1 must be the samples column,
#'                       column 2 must be the "grouping" variable column,
#'                       column 3 must be the "filling color" variable column,
#'                       column 4 must be the "value" column,
#'                       column 5 can be the additionnal "opacity" variable
#'                       column.
#' @param fill.color     A \code{character} vector of length 2 containing colors
#'                       to use to fill the craviolas.
#' @param craviola.width A \code{double} value to specify the width of
#'                       craviolas.
#' @param extrema        A \code{logical} specifying if maxima and minima should
#'                       be displayed or not (Default: extrema = TRUE).
#' @param boxplots       A \code{logical} specifying if boxplots should be
#'                       displayed or not (Default: boxplots = TRUE).
#' @param boxplot.width  A \code{double} specifying the width of boxplots
#'                       (Default: boxplot.width = 0.04).
#' @param mean.value     A \code{logical} specifying if the red dot of the mean
#'                       value should be displayed or not
#'                       (Default: mean.value = TRUE).
#' @param bins           A \code{logical} specifying if the violin plos should
#'                       be binned following specific quantiles to be displayed
#'                       following different opacity thanks to the values in the
#'                       5th column of the data data.frame
#'                       (Default: mean.value = FALSE).
#' @param bins.quantiles A \code{double} vector to define the limits between the
#'                       bins used by the quantile function.
#'                       (Default: bins.quantiles = seq(0.1,0.9,0.1)).
#' @param lines.col      A \code{character} matching a color to use for the
#'                       border lines of the craviola's bins.
#' @return a \code{gg} craviola plot.
#' @author Yoann Pageaud.

#TODO : cut density at min and max !!

ggcraviola<-function(data, fill.color=c("blue","red"), craviola.width = 1,
                     extrema = TRUE, boxplots = TRUE, boxplot.width=0.04,
                     mean.value = TRUE, bins=FALSE,
                     bins.quantiles=seq(0.1,0.9,0.1), lines.col = NA){

  Annot.table<-data[!duplicated(data[[1]]),1:3]

  if(length(unique(Annot.table[[3]])) > 2){
    stop("More than 2 conditions inputed. Only 2 conditions tolerated.")
  } else {
    original.var.col<-levels(Annot.table[[3]])
    levels(Annot.table[[3]])[1]<-"1"
    levels(Annot.table[[3]])[2]<-"2"
  }
  amount.grp<-length(unique(Annot.table[[2]]))
  if(amount.grp > 1){
    lapply(seq_along(unique(Annot.table[[2]])), function(i){
      levels(Annot.table[[2]])[i]<<- i
    })
  }

  mylist_data<-split(data,f = data[[1]])
  list_val1<-lapply(mylist_data,subset, select = 4)
  list_vect.val1<-lapply(list_val1,unlist)

  #Create stats plots
  list.bp.stat<-lapply(seq_along(list_vect.val1),function(i){
    qiles<-quantile(list_vect.val1[[i]])
    means<-mean(list_vect.val1[[i]])
    if(Annot.table[Annot.table[[1]] == names(list_vect.val1)[i],3] == 1){
      x.pos <- -boxplot.width
    } else { x.pos <- boxplot.width }
    if(Annot.table[Annot.table[[1]] == names(list_vect.val1)[i],2] != 1){
      x.pos<-x.pos + (as.integer(Annot.table[Annot.table[[1]] ==
                                               names(list_vect.val1)[i],2])-1)
    }
    data.frame(Var.col = Annot.table[Annot.table[[1]] ==
                                       names(list_vect.val1)[i],3],
               x=x.pos, min = qiles[1], low=qiles[2], mid=qiles[3],
               top=qiles[4], max = qiles[5],mean = means,
               pos.crav = round(x.pos))
  })

  box.dframe<-do.call(rbind, list.bp.stat)

  if(extrema){
    min.max.dframe<-melt(data = box.dframe,
                         id.vars = c("Var.col","pos.crav","x"),
                         measure.vars = c("min","max"))
  }

  #Create Craviola plot
  list_dens.res<-lapply(list_vect.val1,density)

  list_dens.df<-lapply(list_dens.res, function(i){
    data.frame(y.pos=i$x, dens.curv = i$y*craviola.width)
  })

  if(amount.grp > 1){
    list_oriented_dens<-lapply(seq_along(list_dens.df),function(i){
      if(Annot.table[Annot.table[[1]] == names(list_dens.df)[i],3] == 1){
        list_dens.df[[i]]$dens.curv<<-list_dens.df[[i]]$dens.curv * -1
        list_dens.df[[i]]
      } else { list_dens.df[[i]] }

      if(Annot.table[Annot.table[[1]] == names(list_dens.df)[i],2] != 1){
        list_dens.df[[i]]$dens.curv<<-list_dens.df[[i]]$dens.curv +
          (as.integer(Annot.table[Annot.table[[1]] ==
                                    names(list_dens.df)[i],2]) - 1)
        list_dens.df[[i]]
      } else { list_dens.df[[i]] }
    })
  }

  #Remove density values outside the extrema
  xtrems<-ls.quantile(ls = list_vect.val1, qtiles = c(0,1))
  bined.xtrm.dens<-bin.polygons(list_oriented_dens = list_oriented_dens,
                 list.quant.lim = xtrems,Annot.table = Annot.table)
  #Keep only bin 1 for each sample
  list_oriented_dens<-lapply(seq_along(bined.xtrm.dens), function(i){
    df<-bined.xtrm.dens[[i]][bined.xtrm.dens[[i]]$bin == 1,c(2,3)]
    rownames(df)<-NULL
    df
  })
  ##Create Bins based on a third variable
  if(bins){ #Create bin polygons
    #Bin polygons
    list.quant.lim<-ls.quantile(ls = list_vect.val1, qtiles = bins.quantiles)
    list.dfs<-bin.polygons(list_oriented_dens = list_oriented_dens,
                           list.quant.lim = list.quant.lim,
                           Annot.table = Annot.table)
    #Calculate average value on 3rd variable for each bin
    list.av.val2<-lapply(seq_along(list.quant.lim), function(i){
      smpl.data<-data[data[[1]] == names(list.quant.lim)[i],]
      # Check if external quantiles are min and max
      if(min(list_vect.val1[[i]]) != list.quant.lim[[i]][1]){
        list.quant.lim[[i]]<<-c(min(list_vect.val1[[i]]),list.quant.lim[[i]])
      }
      if(max(list_vect.val1[[i]]) != rev(list.quant.lim[[i]])[1]){
        list.quant.lim[[i]]<<-c(list.quant.lim[[i]],max(list_vect.val1[[i]]))
      }
      smpl.data[["bin.groups"]]<-findInterval(smpl.data[[4]],
                                              list.quant.lim[[i]],
                                              all.inside = T)-1
      unlist(lapply(sort(unique(smpl.data$bin.groups)), function(j){
        mean(smpl.data[smpl.data$bin.groups == j,5])
      }))
    })
    vec_av.val2<-unlist(list.av.val2) #Make vector average val2
    #Map Bins average value on 3rd variable to the dataframe list
    list.dfs<-lapply(seq_along(list.dfs), function(i){
      list.bins<-split(x =list.dfs[[i]],f= list.dfs[[i]]$bin)
      list.bins<-Map(cbind,list.bins, bin.av.val2 = list.av.val2[[i]])
      do.call(rbind,list.bins)
    })
    #Add Sample IDs, Var.grp and Var.col
    list.dframes<-Map(cbind,Var1 = Annot.table[[1]], Var.grp = Annot.table[[2]],
                      Var.col = Annot.table[[3]], list.dfs)
  } else { #No bin polygons
    #Add Sample IDs, Var.grp and Var.col
    list.dframes<-Map(cbind,Var1 = Annot.table[[1]], Var.grp = Annot.table[[2]],
                      Var.col = Annot.table[[3]], list_oriented_dens)
  }
  #Make data.frame
  dframe<-do.call(rbind,list.dframes)
  
  #Plot
  craviola.plot<-ggplot() +
    scale_x_continuous(breaks = as.integer(levels(dframe$Var.grp))-1,
                       labels = levels(data[[2]])) +
    xlab(colnames(Annot.table)[2]) + ylab(colnames(data)[4]) +
    labs(fill = colnames(Annot.table)[3], color = "Extrema",
         alpha = colnames(data)[5]) +
    guides(fill = guide_legend(order = 1)) +
    scale_fill_manual(values = fill.color, labels = original.var.col)
  
  #Plot Options
  if(bins){ #bins TRUE
    craviola.plot<-craviola.plot +
      geom_polygon(data = dframe, mapping = aes(dens.curv, y.pos,fill = Var.col,
                                               group = interaction(Var1,bin),
                                               alpha = bin.av.val2),
                   colour = lines.col) +
      guides(alpha = guide_legend(order = 2)) +
      scale_alpha_continuous(limits=c(floor(min(vec_av.val2)),
                                      ceiling(max(vec_av.val2))),
                             breaks=round(seq(floor(min(vec_av.val2)),
                                              ceiling(max(vec_av.val2)),
                                              length.out = 5)))
  } else { #bins FALSE
    craviola.plot<-craviola.plot +
      geom_polygon(data = dframe,
                   mapping = aes(dens.curv, y.pos,fill = Var.col,
                                 group = interaction(Var.col,Var1)),
                   colour = lines.col)
  }
  if(boxplots){ #boxplots TRUE
    craviola.plot<-craviola.plot +
      geom_boxplot(data = box.dframe,
                   mapping = aes(x=x, ymin = low,lower = low, middle = mid,
                                 upper = top, ymax = top), stat = "identity")
  }
  if(extrema){ #extrema TRUE
    craviola.plot<-craviola.plot +
      geom_segment(data=min.max.dframe,
                   mapping = aes(x = pos.crav, xend = x, y = value,
                                 yend = value, color = fill.color), size = 1) +
      guides(color = guide_legend(order = 3))
  }
  if(mean.value){ #mean.value TRUE
    craviola.plot<-craviola.plot +
      geom_point(data = box.dframe,mapping = aes(x = x,y = mean),size=2,
                 color = "red")
  }
  craviola.plot
}
