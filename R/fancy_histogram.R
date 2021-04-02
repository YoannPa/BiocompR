
#' Computes in parallel and plot an histogram using ggplot2 from a given vector
#' of values.
#'
#' @param x          A \code{numeric} or \code{character} vector to be used for
#'                   plotting the histogram (if x is a character vector, 'xmax',
#'                   'ngrad' and 'round.grad' will not be available).
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
#' @param ncores     An \code{integer} to specifying the number of cores to use
#'                   when parallel-running the computation of the histogram
#'                   (Default: ncores = 1).
#' @param bin.col    A \code{character} matching a R color code to be use to
#'                   fill the histogram bins (Default: bin.col = "#0570b0").
#' @param show.annot A \code{logical} specifying whether the median and the
#'                   cut-off annotation should be displayed (show.annot = TRUE)
#'                   or not (Default: show.annot = FALSE).
#' @param facet      A \code{numeric} or \code{character} vector of n levels to
#'                   be used to generate n histograms in separate panels (Only
#'                   works if x is a character vector).
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

#TODO: Add facet option for x as a numeric vector
#TODO: Adapt show.annot to the faceting and when X is a character.
fancy.hist <- function(
  x, xmax = max(x, na.rm = TRUE), nbreaks = 11, ngrad = 11, round.grad = 1,
  ncores = 1, bin.col = "#0570b0", show.annot = FALSE, facet = NULL){
  #Check if is any NAs
  if(any(is.na(x))){ x <- x[!is.na(x)] } #rm NAs
  #Check x type
  if(is.numeric(x)){
    xbreaks <- seq(0, xmax, length.out = nbreaks) #Set breaks
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
    histdata <- parallel::mclapply(
      X = seq(length(xbreaks)-1), mc.cores = ncores, FUN = function(i){
        if(i == length(xbreaks)-1){
          #If last bin take values equal to maximum too
          length(x[x >= xbreaks[i] & x <= xbreaks[i+1]])
        } else { length(x[x >= xbreaks[i] & x < xbreaks[i+1]]) }
      })
    histdata <- unlist(histdata)

    histbreaks <- xgrads*(length(histdata)/xmax) + 0.5 #Scale graduations
    histlim <- seq(0, xmax, length.out = ngrad)*(length(histdata)/xmax) + 0.5
    histlim <- c(histlim[1], rev(histlim)[1])

    dframe <- data.table(x = seq(histdata), y = histdata) #Create data.table
  } else if(is.character(x)){
    uniq.x <- unique(x)
    if(!is.null(facet)){
      if(length(uniq.x) != length(facet)){
        stop("facet must be of the same length as unique(x) to be applied.")
      } else {
        #Make dt.facet
        dt.facet <- data.table(uniq.x = uniq.x, facet = facet)
        #Get facet sizes
        dt.facet[, size := .N, by = "facet"]
        #Calculate number of breaks in each facet
        dt.facet[, n.breaks := round((size*nbreaks)/length(uniq.x))]
        reduce.dt.facet <- unique(dt.facet, by = "facet")[, -c(1), ]
        #Breakdown proportionality of breaks
        ls.xbreaks <- lapply(X = reduce.dt.facet$facet, FUN = function(i){
          seq(0, reduce.dt.facet[facet == i]$size,
              length.out = reduce.dt.facet[facet == i]$n.breaks)
        })
        ls.histdata <- lapply(X = seq_along(ls.xbreaks), FUN = function(i){
          histdata <- parallel::mclapply(
            X = seq(length(ls.xbreaks[[i]])-1), mc.cores = ncores,
            FUN = function(j){
              if(j == length(ls.xbreaks[[i]])-1){
                #If last bin take values equal to maximum too
                sub.chars <- uniq.x[
                  round(ls.xbreaks[[i]][j]):round(ls.xbreaks[[i]][j+1])]
              } else {
                sub.chars <- uniq.x[
                  round(ls.xbreaks[[i]][j]):round(ls.xbreaks[[i]][j+1]-1)]
              }
              length(x[x %in% sub.chars])
            })
          data.table(x = seq(histdata), y = unlist(histdata))
        })
        names(ls.histdata) <- reduce.dt.facet$facet
        dframe <- rbindlist(l = ls.histdata, idcol = "facet")
        dframe[, facet := factor(
          x = facet, levels = levels(reduce.dt.facet$facet))]
        dframe[, x := .I] #Update the rank on all data
      }
    } else {
      xmax <- length(uniq.x)
      xbreaks <- seq(0, length(uniq.x), length.out = nbreaks) #Set breaks
      histdata <- parallel::mclapply(
        X = seq(length(xbreaks)-1), mc.cores = ncores, FUN = function(i){
          if(i == length(xbreaks)-1){
            #If last bin take values equal to maximum too
            sub.chars <- uniq.x[round(xbreaks[i]):round(xbreaks[i+1])]
          } else {
            sub.chars <- uniq.x[round(xbreaks[i]):round(xbreaks[i+1]-1)]
          }
          length(x[x %in% sub.chars])
        })
      histdata <- unlist(histdata)
      dframe <- data.table(x = seq(histdata), y = histdata) #Create data.table
    }
  } else { stop("format not supported for 'x'.") }

  #Get recommended cut-off value and median
  if(show.annot){
    cat("Compute Median and Cutoff\n")
    cutoff.val <- quantmod::findValleys(x = histdata)[1]
    cutoff.pos <- cutoff.val*(length(histdata)/xmax) + 0.5
    median.val <- median(x, na.rm = TRUE)
    median.pos <- median.val*(length(histdata)/xmax) + 0.5
  }
  #Default ggplot
  gghist <- ggplot(data = dframe, aes(x = x, y = y)) +
    geom_bar(stat = "identity", width = 1, fill = bin.col, alpha = 0.7) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Values", y = "Frequency") +
    theme(plot.margin = unit(c(0.1,1,0.1,0.1),"cm"),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          panel.background = element_rect(fill = "white"))

  if(is.character(x)){
    #Plot for x as a character vector
    if(is.null(facet)){
      #Without facet
      gghist <- gghist +
        scale_x_continuous(expand = c(0, 0)) +
        theme(panel.grid.major.y = element_line(colour = "grey"),
              panel.grid.minor.y = element_line(colour = "grey"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank())
    } else {
      #With facet
      gghist <- gghist +
        scale_x_continuous(expand = c(0, 0)) +
        facet_grid(. ~ facet, scales = "free", space = "free") +
        theme(panel.grid.major.y = element_line(colour = "grey"),
              panel.grid.minor.y = element_line(colour = "grey"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank())
    }
  } else if(is.numeric(x)){
    #Plot for x as a numeric vector
    gghist <- gghist +
      scale_x_continuous(breaks = histbreaks, labels = xlabs, limits = histlim,
                         expand = c(0, 0)) +
      theme(panel.grid.major = element_line(colour = "grey"),
            panel.grid.minor = element_line(colour = "grey"))
  }
  #Add annotations
  if(show.annot){
    gghist <- gghist +
      geom_vline(xintercept = median.pos, color = "#313695", size = 0.7) +
      ggrepel::geom_label_repel(data = data.frame(), aes(
        x = median.pos, y = Inf, fontface = 2,
        label = paste0("median = ", round(x = median.val, digits = 2))),
        vjust = 0.5, color = "#313695")
    #If recommended cut-off inferior or equal to median, plot cut-off
    if(cutoff.val <= median.val){
      gghist <- gghist +
        geom_vline(xintercept = cutoff.pos, color = "#d7191c", size = 0.7) +
        ggrepel::geom_label_repel(
          data = data.frame(),
          aes(x = cutoff.pos, y = Inf, fontface = 2,
              label = paste0("cut-off = ", round(x = cutoff.val, digits = 2))),
          vjust = 0.5, color = "#d7191c")
    }
  }
  gghist
}
