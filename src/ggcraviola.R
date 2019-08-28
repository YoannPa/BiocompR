##IMPORTS
Imports = c("ggplot2","data.table")
lapply(Imports, library, character.only = T)
source("src/bin_polygons.R")

##DUMMY DATA
dummy.complete = data.frame(
  Samples=rep(paste0("Sample",c(1:6)),each = 1000),
  Groups=rep(c('A','B','C'),each = 2000),
  Conditions=rep(c('I','J'),each = 1000,3),
  Values=c(rnorm(1000,0), rnorm(1000, 0.5),
           rnorm(1000, 3), rnorm(1000, 3.5),
           rnorm(1000,-3), rnorm(1000, -3.5)),
  "Function(2ndVar)" = rep(rep(x =c(60,50,40,30,20,10,30,40,50,60),each=100),6))

dummy.complete.bis = data.frame(
  Samples=rep(paste0("Sample",c(1:6)),each = 1000),
  Groups=rep(c('A','B','C'),each = 2000),
  Conditions=rep(c('I','J'),each = 1000,3),
  Values=c(rnorm(1000,0), rnorm(1000, 0.5),
           rnorm(1000, 3), rnorm(1000, 3.5),
           rnorm(1000,-3), rnorm(1000, -3.5)))

dummy.minimal = data.frame(
  Samples = rep(paste0("Sample",c(1:2)),each = 1000),
  Values=c(rnorm(1000,0), rnorm(1000, 0.5)))

##FUNCTIONS

# ggcraviola ###################################################################

#' Read Chromosomes Methylation from a specific sample.
#'
#' @param data           A \code{data.frame}. 2 formats of data.frame are
#'                       supported.
#'                       A "complete" data.frame of 5 columns:
#'                       column 1 must be the samples column,
#'                       column 2 must be the "grouping" variable column,
#'                       column 3 must be the "filling color" variable column,
#'                       column 4 must be the "value" column,
#'                       column 5 can be the additionnal "opacity" variable
#'                       column to be used if bins = TRUE.
#'                       A "minimal" data.frame of 2 columns:
#'                       column 1 must be the samples column. This column can
#'                       take a maximum of 2 possible levels. Based on this
#'                       column ggcraviola will guess the "grouping" and the
#'                       "filling".
#'                       column 2 must be the "value" column.
#'                       The "minimal" data.frame is usefull to easily plot 1
#'                       craviola with 2 distributions.
#'                       You cannot plot more than 1 craviola with the "minimal"
#'                       data.frame format. the "minimal" format also do not
#'                       support the binning. 
#'                       
#' @param fill.color     A \code{character} vector of length 2 containing colors
#'                       to use to fill the craviolas
#'                       (Default: fill.color = c("blue","red")).
#' @param craviola.width A \code{double} value to specify the width of
#'                       craviolas (Default: craviola.width = 1).
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
#'                       (Default: bins = FALSE).
#' @param bins.quantiles A \code{double} vector to define the limits between the
#'                       bins as percentiles of distributions
#'                       (Default: bins.quantiles = seq(0.1,0.9,0.1)).
#' @param bin.fun        A \code{character} to specify the function to apply on
#'                       values of the additional variable for each bin of
#'                       distrubutions
#'                       (Supported: bin.fun = c("sd","mad","mean");
#'                       Default: bin.fun = "sd"). 
#' @param lines.col      A \code{character} matching a color to use for the
#'                       border lines of the craviola's bins
#'                       (Default: lines.col = NA).
#' @return a \code{gg} craviola plot.
#' @author Yoann Pageaud.

#TODO : cut density at min and max !!

ggcraviola<-function(data, fill.color=c("blue","red"), craviola.width = 1,
                     boxplots = TRUE, boxplot.width=0.04, mean.value = TRUE,
                     bins=FALSE, bins.quantiles=seq(0.1,0.9,0.1), bin.fun="sd",
                     lines.col = NA){
  #Check if data is a data.table and convert if not
  if (!is.data.table(data)){data<-as.data.table(data)}
  #Make annotation table
  if (ncol(data) < 4){
    if (nrow(data[!duplicated(data[[1]])]) == 2) {
      if(is.factor(data[!duplicated(data[[1]])][[1]])){
        Annot.table<-data.table("Samples"=data[!duplicated(data[[1]])][[1]],
                                "Groups"=as.factor(c(1,1)),
                                "Conditions"=data[!duplicated(data[[1]])][[1]],
                                data[!duplicated(data[[1]]),2])
        colnames(data)[1]<-"Samples"
        data<-merge(x = Annot.table[,-4],y = data,by="Samples",all.y=TRUE)
      } else {
        stop("Column 1 in data is not of type 'factor'.")
      }
    } else{
      stop("Missing columns in the data provided.")
    }
  } else {
    Annot.table<-data[!duplicated(data[[1]]),1:3]
  }
  if(length(unique(Annot.table[[3]])) > 2){
    stop("More than 2 conditions inputed. Only 2 conditions tolerated.")
  } else {
    original.var.col<-levels(Annot.table[[3]])
    if (length(levels(Annot.table[[3]])) > 2) {
      stop("More levels than possible values. Only 2 conditions tolerated. Remove the excess levels.")
    } else {
    levels(Annot.table[[3]])[1]<-"1"
    levels(Annot.table[[3]])[2]<-"2"
    }
  }
  if (length(Annot.table[[1]]) < length(levels(Annot.table[[1]]))){
    stop("More levels than matching values found in column 1 of data.")
  }
  amount.grp<-length(unique(Annot.table[[2]]))
  if(amount.grp > 1){
    invisible(lapply(seq_along(unique(Annot.table[[2]])), function(i){
      levels(Annot.table[[2]])[i]<<- i
    }))
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
  
  #Create Craviola plot
  list_dens.res<-lapply(list_vect.val1,density)
  list_dens.df<-lapply(list_dens.res, function(i){
    data.frame(y.pos=i$x, dens.curv = i$y*craviola.width)
  })
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
  names(list_oriented_dens)<-names(list_dens.df)
  #Create Bins based on a third variable
  if(bins){ #Create bin polygons
    #Bin polygons
    list.quant.lim<-ls.quantile(ls = list_vect.val1, qtiles = bins.quantiles)
    list.dfs<-bin.polygons(list_oriented_dens = list_oriented_dens,
                           list.quant.lim = list.quant.lim,
                           Annot.table = Annot.table)
    #Calculate average value on 3rd variable for each bin
    list.fun.val2<-lapply(seq_along(list.quant.lim), function(i){
      smpl.data<-data[data[[1]] == names(list.quant.lim)[i],]
      # Check if external quantiles are min and max
      if(min(list_vect.val1[[i]]) != list.quant.lim[[i]][1]){
        # Add minimum value at the beginning of the vector
        list.quant.lim[[i]]<<-c(min(list_vect.val1[[i]]),list.quant.lim[[i]])
      }
      if(max(list_vect.val1[[i]]) != rev(list.quant.lim[[i]])[1]){
        # Add maximum value at the end of the vector
        list.quant.lim[[i]]<<-c(list.quant.lim[[i]],max(list_vect.val1[[i]]))
      }
      smpl.data[["bin.groups"]]<-findInterval(smpl.data[[4]],
                                              list.quant.lim[[i]],
                                              all.inside = T)-1
      if(bin.fun == "mean"){
        unlist(lapply(sort(unique(smpl.data$bin.groups)), function(j){
          mean(smpl.data[smpl.data$bin.groups == j][[5]],na.rm = T)
        })) 
      } else if(bin.fun == "sd"){
        unlist(lapply(sort(unique(smpl.data$bin.groups)), function(j){
          sd(smpl.data[smpl.data$bin.groups == j][[5]],na.rm = T)
        })) 
      } else if(bin.fun == "mad"){
        unlist(lapply(sort(unique(smpl.data$bin.groups)), function(j){
          mad(smpl.data[smpl.data$bin.groups == j][[5]],na.rm = T)
        }))
      } else {
        stop("Unsupported function. Supported functions: bin.fun = c('mean','sd' and 'mad').")
      }
    })
    vec_av.val2<-unlist(list.fun.val2) #Make vector average val2
    #Map Bins average value on 3rd variable to the dataframe list
    list.dfs<-lapply(seq_along(list.dfs), function(i){
      list.bins<-split(x =list.dfs[[i]],f= list.dfs[[i]]$bin)
      list.bins<-Map(cbind,list.bins, bin.av.val2 = list.fun.val2[[i]])
      do.call(rbind,list.bins)
    })
    #Add Sample IDs, Var.grp and Var.col
    list.dframes<-Map(cbind,Var1 = Annot.table[[1]], Var.grp = Annot.table[[2]],
                      Var.col = Annot.table[[3]], list.dfs)
  } else { #No bin polygons
    #Add Sample IDs, Var.grp and Var.col
    list.dframes<-Map(cbind,Var1 = Annot.table[[1]], Var.grp = Annot.table[[2]],
                      Var.col = Annot.table[[3]],
                      list_oriented_dens[Annot.table[[1]]])
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
                                 upper = top, ymax = top, group = x),
                   stat = "identity")
  }
  if(mean.value){ #mean.value TRUE
    craviola.plot<-craviola.plot +
      geom_point(data = box.dframe,mapping = aes(x = x,y = mean),size=2,
                 color = "red")
  }
  craviola.plot
}
