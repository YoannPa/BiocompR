##IMPORTS
Imports = c("quantmod","ggplot2","ggrepel","parallel")
lapply(Imports, library, character.only = T)

#' Computes in parallel and plot an histogram using ggplot2 from a given vector
#' of values.
#'
#' @param x          A \code{numeric} vector to be used for plotting the
#'                   histogram.
#' @param xmax       A \code{numeric} specifying the maximum limit on the x axis
#'                   to display values. All values above this limit will be
#'                   agregated on the last histogram bin.
#' @param nbreaks    An \code{integer} specifying the number of delimitations to
#'                   be use for histogram bins\cr(Default: nbreaks = 10).
#' @param ngrad      An \code{integer} specifying the number of graduations to
#'                   display on the X-axis of the plot\cr(Default: ngrad = 10).
#' @param round.grad An \code{integer} specifying the number significant digits
#'                   to be considered when calculating graduations\cr
#'                   (Default: round.grad = 1).
#' @param ncores     An \code{integer} to specifying the number of cores to use
#'                   when parallel-running the computation of the histogram\cr
#'                   (Default: ncores = 1).
#' @param xlab       A \code{character} to be used as X-axis label\cr
#'                   (Default: xlab = 'values').
#' @param ylab       A \code{character} to be used as Y-axis label\cr
#'                   (Default: ylab = 'Frequency').
#' @param bin.col    A \code{character} matching a R color code to be use to
#'                   fill the histogram bins.
#' @return A \code{gg} plot of an histogram.
#' @author Yoann Pageaud.
#' @export

fancy.hist<-function(x, xmax, nbreaks = 10, ngrad = 10, round.grad = 1,
                     ncores = 1, xlab="values", ylab = "Frequency",
                     bin.col = "#0570b0"){
  #Check if is any NAs
  if(any(is.na(x))){ x<-x[!is.na(x)] } #rm NAs

  #Set breaks
  xbreaks <- seq(0, xmax, length.out = nbreaks)
  #Set graduations
  xgrads <- round(x = seq(0, xmax, length.out = ngrad), digits = round.grad)
  #Set X axis labels
  xlabs<-as.character(xgrads)
  if(max(x) > xmax){
    #Set Maximum
    x[x > xmax] <- xmax
    xlabs[length(xlabs)]<-paste(xlabs[length(xlabs)],"\nor more")
  }
  #Compute quantities for each bin
  histdata<-mclapply(seq(length(xbreaks)-1), mc.cores = ncores, function(i){
    if(i == length(xbreaks)-1){ #If last bin take values equal to maximum too
      qs<-length(x[x >= xbreaks[i] & x <= xbreaks[i+1]])
    } else { qs<-length(x[x >= xbreaks[i] & x < xbreaks[i+1]]) }
  })
  histdata<-unlist(histdata)
  #Get recommended cut-off value and median
  cat("Compute Median and Cutoff\n")
  cutoff.val <- quantmod::findValleys(histdata)[1]
  cutoff.pos <-cutoff.val*(length(histdata)/xmax) + 0.5
  median.val<-median(x,na.rm = TRUE)
  median.pos<-median.val*(length(histdata)/xmax) + 0.5
  histbreaks<-xgrads*(length(histdata)/xmax) + 0.5 #Scale graduations
  dframe<-data.frame(x= seq(histdata), y=histdata) #Create dataframe
  #Plot
  cat("Plotting\n")
  gghist<-ggplot(data = dframe, aes(x=x, y=y)) +
    geom_bar(stat = "identity", width=1, fill = bin.col, alpha = 0.7) +
    scale_x_continuous(breaks = histbreaks, labels = xlabs,
                       limits = c(histbreaks[1],histbreaks[length(histbreaks)]),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_vline(xintercept = median.pos, color = "#313695", size = 0.7) +
    ggrepel::geom_label_repel(
      data = data.frame(), aes(x = median.pos, y = Inf, fontface = 2,
                               label = paste0("median = ", round(x = median.val,
                                                                 digits = 2))),
      vjust = 1.1, color = "#313695") +
    labs(x = xlab, y = ylab) +
    theme(plot.margin = unit(c(0.1,1,0.1,0.1),"cm"),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          panel.grid.major = element_line(colour = "grey"),
          panel.grid.minor = element_line(colour = "grey"),
          panel.background = element_rect(fill = "white"))
  #If recommended cut-off inferior or equal to median, plot cut-off
  if(cutoff.val <= median.val){
    gghist<-gghist +
      geom_vline(xintercept=cutoff.pos, color="#d7191c", size=0.7) +
      ggrepel::geom_label_repel(
        data = data.frame(),
        aes(x = cutoff.pos, y = Inf, fontface = 2,
            label = paste0("cut-off = ", round(x = cutoff.val, digits = 2))),
        vjust = 2.4,color = "#d7191c")
  }
  gghist
}
