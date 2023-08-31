
#' Applies the quantile function on a list of vectors.
#'
#' @param ls     A \code{list} of vectors.
#' @param qtiles A \code{double} vector specifying the values of percentiles of
#'               interest. Values must be between 0 and 1.
#' @return A \code{list} of vectors containing percentiles values, one vector by
#' distribution.
#' @author Yoann Pageaud.
#' @keywords internal

ls.quantile <- function(ls, qtiles){
    lapply(X = ls, FUN = stats::quantile, qtiles)
}

#' Bins density object following specific percentiles.
#'
#' @param list_oriented_dens A \code{list} of data.frames, each data.frame
#'                           describing a density distribution.
#' @param list.quant.lim     A \code{list} of vectors, each vectors containing
#'                           values of specific percentiles.
#' @param annot_table        A \code{data.frame} containing annotations about
#'                           the distributions.
#' @return A \code{list} of modified data.frames of the density distributions
#' containing breaks to be used for delimiting bins.
#' @author Yoann Pageaud.
#' @keywords internal

bin.polygons <- function(list_oriented_dens, list.quant.lim, annot_table){
    # Remove quantiles with duplicated values
    list.quant.lim <- lapply(X = list.quant.lim, FUN = function(i){
        i[!duplicated(i) & !duplicated(i, fromLast = TRUE)]
    })
    #Add quantile values to density positions vector and sort it.
    list_dens.pos <- lapply(X = list_oriented_dens, FUN = function(i){
        i$y.pos })
    merged.pos <- Map(c, list.quant.lim, list_dens.pos)
    merged.pos <- lapply(X = merged.pos, FUN = sort)
    #Get density values from position-1 and position+1
    list_dens.val <- lapply(X = list_oriented_dens, FUN = function(i){
        i$dens.curv })
    quant.pos <- lapply(X = seq_along(merged.pos), FUN = function(i){
        match(merged.pos[[i]][names(list.quant.lim[[i]])], merged.pos[[i]])
    })
    val.before <- lapply(X = seq_along(quant.pos), FUN = function(i){
        unlist(lapply(X = seq_along(quant.pos[[i]]), FUN = function(j){
            list_dens.val[[i]][quant.pos[[i]][j]-j]
        }))
    })
    val.after <- lapply(X = seq_along(quant.pos), FUN = function(i){
        unlist(lapply(X = seq_along(quant.pos[[i]]), FUN = function(j){
            list_dens.val[[i]][quant.pos[[i]][j]+1-j]
        }))
    })
    #Calculate mean density values
    mean.dens.val <- lapply(X = seq_along(val.before), FUN = function(i){
        rowMeans(data.table::data.table(
            "val.before" = val.before[[i]], "val.after" = val.after[[i]]))
        # rowMeans(data.frame(val.before = val.before[[i]],
        #                     val.after = val.after[[i]]))
    })
    names(mean.dens.val) <- names(merged.pos)
    #Create 1 data.frame per position with the vector c(mean,0,0,mean)
    # with 4 times the position val
    lim.vect <- lapply(X = seq_along(mean.dens.val), FUN = function(i){
        lapply(X = seq_along(mean.dens.val[[i]]), FUN = function(j){
            x.position <- as.integer(
                annot_table[samples == names(mean.dens.val)[i]]$groups) - 1
            # x.position <- as.integer(Annot.table[
            #     Annot.table[[1]] == names(mean.dens.val)[i],2]) - 1
            data.table::data.table(
                "y.pos" = rep(
                    merged.pos[[i]][names(list.quant.lim[[i]])[j]], 4),
                "dens.curv" = c(
                    mean.dens.val[[i]][j],x.position,x.position,
                    mean.dens.val[[i]][j]))
            # data.frame(
            #     y.pos = rep(merged.pos[[i]][names(list.quant.lim[[i]])[j]], 4),
            #     dens.curv = c(
            #         mean.dens.val[[i]][j],x.position,x.position,
            #         mean.dens.val[[i]][j]))
        })
    })
    #Split density data.frames by position+1
    split.positions <- lapply(X = seq_along(quant.pos), FUN = function(i){
        unlist(lapply(X = seq_along(
            quant.pos[[i]]), FUN = function(j){ quant.pos[[i]][j] + 1 - j }))
    })
    splitted.dens.df <- lapply(
        X = seq_along(list_oriented_dens), FUN = function(i){
            split(list_oriented_dens[[i]], findInterval(x = seq(nrow(
                list_oriented_dens[[i]])), vec = split.positions[[i]]))
        })
    #Define bin numbers and Map them to each list of data.frames
    splitted.dens.df <- lapply(
        X = seq_along(splitted.dens.df), FUN = function(i){
        bin.names <- names(splitted.dens.df[[i]])
        Map(cbind, bin = bin.names, splitted.dens.df[[i]])
    })
    lim.vect <- lapply(X = seq_along(lim.vect), FUN = function(i){
        bin.names <- lapply(
            X = seq(0, length(lim.vect[[i]])-1), FUN = function(bin){
            c(bin, bin, bin+1, bin+1)
        })
        Map(cbind, bin = bin.names, lim.vect[[i]])
    })
    #Map rbind() the list of daframes and the list of position dataframes
    BiocompR::warn.handle(
        pattern = paste(
            "number of columns of result is not a multiple of vector length",
            "\\(arg 2\\)"),
        list.mix.df <- Map(rbind, splitted.dens.df, lim.vect))
    list.dfs <- lapply(X = list.mix.df, FUN = function(i){
        mrg_df <- do.call(rbind, i)
        mrg_df[!duplicated(mrg_df), ] # Remove appended rogue dataframe
    })
    return(list.dfs)
}
