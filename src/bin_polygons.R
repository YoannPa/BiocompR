
##FUNCTIONS

# ls.quantile ##################################################################

#' @description Applies the quantile function on a list of vectors.
#' 
#' @param ls     A \code{list} of vectors.
#' @param qtiles A \code{double} vector specifying the values of percentiles of
#'               interest. Values must be between 0 and 1.
#' @value a \code{list} of vectors containing percentiles values, one vector by
#' distribution. 
#' @author Yoann Pageaud.

ls.quantile<-function(ls,qtiles){
  lapply(ls,quantile, qtiles)
}

# bin.polygons #################################################################

#' @description Bins density object following specific percentiles.
#' 
#' @param list_oriented_dens A \code{list} of data.frames, each dataframes
#'                           describeing a density distribution.
#' @param list.quant.lim     A \code{list} of vectors, each vectors containing
#'                           values of specific percentiles.
#' @param Annot.table        A \code{data.frame} containing annotations about
#'                           the distributions.
#' @value a \code{list} of modified data.frames of the density distributions
#' containing breaks to be used for delimitating bins.
#' @author Yoann Pageaud.

bin.polygons<-function(list_oriented_dens,list.quant.lim,Annot.table){
  #Add quantile values to density positions vector and sort it.
  list_dens.pos<-lapply(list_oriented_dens, function(i){ i$y.pos })
  merged.pos<-Map(c,list.quant.lim,list_dens.pos)
  merged.pos<-lapply(merged.pos,sort)
  #Get density values from position-1 and position+1
  list_dens.val<-lapply(list_oriented_dens, function(i){ i$dens.curv })
  quant.pos<-lapply(seq_along(merged.pos), function(i){
    match(merged.pos[[i]][names(list.quant.lim[[i]])],merged.pos[[i]])
  })
  val.before<-lapply(seq_along(quant.pos), function(i){
    unlist(lapply(seq_along(quant.pos[[i]]),function(j){
      list_dens.val[[i]][quant.pos[[i]][j]-j]
    }))
  })
  val.after<-lapply(seq_along(quant.pos), function(i){
    unlist(lapply(seq_along(quant.pos[[i]]),function(j){
      list_dens.val[[i]][quant.pos[[i]][j]+1-j]
    }))
  })
  #Calculate mean density values
  mean.dens.val<-lapply(seq_along(val.before), function(i){
    rowMeans(data.frame(val.before = val.before[[i]],
                        val.after = val.after[[i]]))
  })
  names(mean.dens.val)<-names(merged.pos)
  #Create 1 dataframe per position with the vector c(mean,0,0,mean)
  # with 4 times the position val
  lim.vect<-lapply(seq_along(mean.dens.val), function(i){
    lapply(seq_along(mean.dens.val[[i]]), function(j){
      x.position<-as.integer(Annot.table[
        Annot.table[[1]] == names(mean.dens.val)[i],2]) - 1
      data.frame(y.pos = rep(merged.pos[[i]][names(list.quant.lim[[i]])[j]],
                             4),
                 dens.curv = c(mean.dens.val[[i]][j],x.position,x.position,
                               mean.dens.val[[i]][j]))
    })
  })
  #Split density dataframes by position+1
  split.positions<-lapply(seq_along(quant.pos),function(i){
    unlist(lapply(seq_along(quant.pos[[i]]),
                  function(j){quant.pos[[i]][j]+1-j}))
  })
  splitted.dens.df<-lapply(seq_along(list_oriented_dens), function(i){
    split(list_oriented_dens[[i]],
          findInterval(1:nrow(list_oriented_dens[[i]]),split.positions[[i]]))
  })
  #Define bin numbers and Map them to each list of dataframes
  splitted.dens.df<-lapply(seq_along(splitted.dens.df), function(i){
    bin.names<-names(splitted.dens.df[[i]])
    Map(cbind,bin = bin.names,splitted.dens.df[[i]])
  })
  lim.vect<-lapply(seq_along(lim.vect),function(i){
    bin.names<-lapply(seq(0,length(lim.vect[[i]])-1), function(bin){
      c(bin, bin, bin+1, bin+1)
    })
    Map(cbind,bin = bin.names,lim.vect[[i]])
  })
  #Map rbind() the list of daframes and the list of position dataframes
  suppressWarnings(list.mix.df<-Map(rbind,splitted.dens.df,lim.vect))
  list.dfs<-lapply(list.mix.df,function(i){
    mrg_df<-do.call(rbind,i)
    mrg_df[!duplicated(mrg_df),] # remove appended rogue dataframe
  })
  return(list.dfs)
}
