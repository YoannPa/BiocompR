
#' Computes pairwise Kolmogorov-Smirnov tests on a matrix and display results in
#' a fused plot
#'
#' @param param1 A \code{type} parameter description.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export

#TODO Write help
#TODO Test function
ks.plot <- function(
  data, ncores = 1, order.select, order.method, axis.text.y, axis.text.x,
  annot.grps, annot.pal, annot.size, lgd.title, lgd.breaks1, lgd.breaks2,
  lgd.nbin1, lgd.nbin2, lgd.pal1, lgd.pal2){
  #Run KS test
  upper<-pairwise.ks(data=data, statistic = 'stat', ncores=ncores)$res.statistic
  lower<-pairwise.ks(data = data, statistic = 'p', ncores=ncores)$res.statistic
  #Kolmogorov Plot
  ks_plot<-fused.view(
    sample.names = colnames(data), upper.mat = upper, lower.mat = lower,
    order.select = order.select, order.method = order.method,
    axis.text.y = axis.text.y, axis.text.x = axis.text.x,
    annot.grps = annot.grps, annot.pal = annot.pal, annot.size = annot.size,
    lgd.title = lgd.title,
    set.lgd1.title = "K.S.\nD-Stat", set.lgd2.title = "K.S.\nP-value",
    lgd.limits = c(0,1), lgd.breaks1 = lgd.breaks1, lgd.breaks2 = lgd.breaks2,
    lgd.nbin1 = lgd.nbin1, lgd.nbin2 = lgd.nbin2,
    lgd.pal1 = lgd.pal1, lgd.pal2 = lgd.pal2, lgd.width2 = 35,
    raster1 = TRUE, raster2 = FALSE)
  return(ks_plot)
}
