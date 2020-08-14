
#' Computes in parallel and plot an histogram using ggplot2 from a given vector
#' of values.
#'
#' @param x          A \code{numeric} vector to be used for plotting the
#'                   histogram.
#' @param xmax       A \code{numeric} specifying the maximum limit on the X-axis
#'                   to display values. All values above this limit will be
#'                   agregated on the last histogram bin
#'                   (Default: max(x, na.rm = TRUE)).
#' @param nbreaks    An \code{integer} specifying the number of delimitations to
#'                   be use for histogram bins (Default: nbreaks = 11).
#' @param ngrad      An \code{integer} specifying the number of graduations to
#'                   display on the X-axis of the plot (Default: ngrad = 11).
#' @param round.grad An \code{integer} specifying the number significant digits
#'                   to be considered when calculating graduations
#'                   (Default: round.grad = 1).
#' @param bin.col    A \code{character} matching a R color code to be use to
#'                   fill the histogram bins (Default: bin.col = "#0570b0").
#' @param show.annot A \code{logical} specifying whether the median and the
#'                   cut-off annotation should be displayed (show.annot = TRUE)
#'                   or not (Default: show.annot = FALSE).
#' @return A \code{gg} plot of an histogram.
#' @author Yoann Pageaud.
#' @examples
#' # Basic use of fancy.hist()
#' fancy.hist(x = rnorm(100, 5, 1))
#' # Set a specific maximum limit on the X-axis to agregate values above it in a
#' # single bin.
#' fancy.hist(x = rnorm(100, 5, 1), xmax = 7)
#' # Change the number of graduations on X-axis
#' fancy.hist(x = rnorm(100, 5, 1), xmax = 7, ngrad = 8)
#' # Change the number of histogram bins
#' fancy.hist(x = rnorm(100, 5, 1), xmax = 7, ngrad = 8, nbreaks = 8)
#' # Change filling color of the histogram
#' fancy.hist(
#'  x = rnorm(100, 5, 1), xmax = 7, ngrad = 8, nbreaks = 8, bin.col = "orange")
#' # Display annotations
#' fancy.hist(
#'  x = rnorm(100, 5, 1), xmax = 7, ngrad = 8, nbreaks = 8, bin.col = "orange",
#'  show.annot = TRUE)
#' # Rename X-axis and Y-axis titles
#' fancy.hist(
#'  x = rnorm(100, 5, 1), xmax = 7, ngrad = 8, nbreaks = 8, bin.col = "orange",
#'  show.annot = TRUE) +
#' labs(x = "Distribution of the values", y = "Number of values in each bin")
#' @export

fancy.hist<-function(
  x, xmax = max(x, na.rm = TRUE), nbreaks = 11, ngrad = 11, round.grad = 1,
  bin.col = "#0570b0", show.annot = FALSE){
  #Check if is any NAs
  if(any(is.na(x))){ x<-x[!is.na(x)] } #rm NAs

  #Set breaks
  xbreaks <- seq(0, xmax, length.out = nbreaks)
  #Set graduations
  xgrads <- round(x = seq(0, xmax, length.out = ngrad), digits = round.grad)
  #Set X axis labels
  xlabs <- as.character(xgrads)
  if(max(x) > xmax){
    #Set Maximum
    x[x > xmax] <- xmax
    xlabs[length(xlabs)] <- paste(xlabs[length(xlabs)], "\nor more")
  }
  #Compute quantities for each bin
  histdata <- lapply(X = seq(length(xbreaks)-1), FUN = function(i){
    if(i == length(xbreaks)-1){ #If last bin take values equal to maximum too
      qs<-length(x[x >= xbreaks[i] & x <= xbreaks[i+1]])
    } else { qs <- length(x[x >= xbreaks[i] & x < xbreaks[i+1]]) }
  })
  histdata<-unlist(histdata)
  #Get recommended cut-off value and median
  if(show.annot){
    cat("Compute Median and Cutoff\n")
    cutoff.val <- quantmod::findValleys(histdata)[1]
    cutoff.pos <- cutoff.val*(length(histdata)/xmax) + 0.5
    median.val <- median(x,na.rm = TRUE)
    median.pos <- median.val*(length(histdata)/xmax) + 0.5
  }
  histbreaks <- xgrads*(length(histdata)/xmax) + 0.5 #Scale graduations
  histlim <- seq(0, xmax, length.out = ngrad)*(length(histdata)/xmax) + 0.5
  histlim <- c(histlim[1], rev(histlim)[1])
  dframe <- data.frame(x = seq(histdata), y = histdata) #Create dataframe
  #Plot
  cat("Plotting\n")
  gghist <- ggplot(data = dframe, aes(x = x, y = y)) +
    geom_bar(stat = "identity", width = 1, fill = bin.col, alpha = 0.7) +
    scale_x_continuous(
      breaks = histbreaks, labels = xlabs, limits = histlim, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Values", y = "Frequency") +
    theme(plot.margin = unit(c(0.1,1,0.1,0.1),"cm"),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          panel.grid.major = element_line(colour = "grey"),
          panel.grid.minor = element_line(colour = "grey"),
          panel.background = element_rect(fill = "white"))
  if(show.annot){
    gghist <- gghist +
      geom_vline(xintercept = median.pos, color = "#313695", size = 0.7) +
      ggrepel::geom_label_repel(data = data.frame(), aes(
        x = median.pos, y = Inf, fontface = 2,
        label = paste0("median = ", round(x = median.val, digits = 2))),
        vjust = 1.1, color = "#313695")
    #If recommended cut-off inferior or equal to median, plot cut-off
    if(cutoff.val <= median.val){
      gghist <- gghist +
        geom_vline(xintercept = cutoff.pos, color = "#d7191c", size = 0.7) +
        ggrepel::geom_label_repel(
          data = data.frame(),
          aes(x = cutoff.pos, y = Inf, fontface = 2,
              label = paste0("cut-off = ", round(x = cutoff.val, digits = 2))),
          vjust = 2.4, color = "#d7191c")
    }
  }
  gghist
}
