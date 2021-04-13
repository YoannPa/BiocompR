
##FUNCTIONS

#' Computes boxplots or violins from 1 variable values against ranges of a 2nd
#' one.
#'
#' @param data   A two-columns \code{data.frame}:
#'                 \itemize{
#'                  \item{column 1 contains values of the 1st variable. These
#'                  values will be used for computing boxplots.}
#'                  \item{column 2 contains values of the 2nd variable. These
#'                  values will be used to define ranges, one range per
#'                  boxplot.}
#'                 }
#' @param violin   A \code{logical} to specify whether violins of categories
#'                 should be displayed with the boxplots (Warning: setting
#'                 violin = TRUE significantly increase the computing time.
#'                 Default: violin = FALSE).
#' @param cat.step A \code{numeric} to specify the size of ranges.
#'                 (Default: cat.step = 10L).
#' @param cat.max  Maximum amount of ranges to create (Default: cat.max = 10L).
#'                 Values outside the last range are grouped into an extra bin
#'                 named "over <cat.step * cat.max>". In total the plot will
#'                 display 1 boxplot for each range, and one additional boxplot
#'                 for the values falling into the extra bin.
#' @param fill     A \code{character} to specify a color to be used for filling
#'                 boxplot. The color will be used to derive light shades from
#'                 it. The given input color will match the color of the boxplot
#'                 computed on the highest number of values (in other words: for
#'                 the largest category). Consequently, all the shades generated
#'                 will be lighter than the original one. Thus, it is advised to
#'                 use a dark color (Default: fill = "deepskyblue3").
#' @return A \code{gg} plot object displaying boxplots and/or violin plots.
#' @author Yoann Pageaud, Yassen Assenov.
#' @export

bivar.plot <- function(data, violin = FALSE, cat.step = 10L, cat.max = 10L,
                       fill = "deepskyblue3"){
  #Convert into a data.table
  if(!data.table::is.data.table(data)){
    data <- data.table::as.data.table(data)
  }
  #Create colnames
  oldcolnames<-colnames(data)
  colnames(data)<-c("Var1","Var2")
  #Create Categories from Var2
  data[, category :=.(as.integer(floor(Var2 / cat.step)))]
  #Assign values to last category if there category is above cat.max
  data[category > cat.max, category := cat.max]
  #Get Number of CpGs by categories
  data[, cat.sizes := .N, by = "category"]
  #Calculate intervals
  cat.breaks <- cumsum(rep(cat.step, cat.max))
  cat.first <- ifelse(is.integer(data[, 1]), "1", "0")
  #Create X Axis Labels
  cat.labels <- paste0("] ", c(cat.first, cat.breaks[-length(cat.breaks)]),
                       ", ", cat.breaks, " ]")
  #Close first interval
  cat.labels[1]<-gsub(pattern = '^\\]\\s', replacement = "[ ", x=cat.labels[1])
  #Add last interval
  cat.labels <- c(cat.labels, paste("over", utils::tail(cat.breaks, 1)))
  #Convert categories as factors
  data[, category := .(factor(category, levels = c(0:(length(cat.breaks)))))]
  #Update levels
  levels(data$category) <- cat.labels
  #Get boxplot statistics
  dfr <- data[, stats::quantile(x = Var1), by = c("category", "cat.sizes")]
  #Categories as Levels
  levels(dfr$category)<-dfr[
    match(x = levels(category), table = category)][, category:=.(paste0(
      category,"\n",formatC(x = cat.sizes, format = "e",digits = 2)))]$category
  #Change category also in data
  levels(data$category)<-levels(dfr$category)

  #Plot categories boxplot
  bivar<-ggplot2::ggplot()
  if(violin){
    bivar <- bivar +
      ggplot2::geom_violin(data = data, mapping = ggplot2::aes(
        x = category, y = Var1, alpha = cat.sizes), fill = fill) +
      ggplot2::geom_boxplot(data = dfr, mapping = ggplot2::aes(
        x = category, y = V1), fill = "white", outlier.size = 1.5, width = 0.1)
  } else {
    bivar <- bivar +
      ggplot2::geom_boxplot(data = dfr, mapping = ggplot2::aes(
        x = category, y = V1), outlier.size = 1.5, color = NA,
        outlier.color = "black") +
      ggplot2::geom_boxplot(data = dfr, mapping = ggplot2::aes(
        x = category, y = V1, alpha = cat.sizes), fill = fill,
        outlier.shape = NA)
  }
  bivar <- bivar +
    ggplot2::theme(
      legend.position="none",
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey"),
      axis.title = ggplot2::element_text(size = 13, hjust = 0.5),
      axis.text = ggplot2::element_text(size = 11, color = "black"),
      axis.ticks = ggplot2::element_blank()) +
    ggplot2::labs(x = oldcolnames[2], y = oldcolnames[1])
  return(bivar)
}
