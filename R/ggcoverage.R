
#' Plots an annotated stacked barplot.
#'
#' @param data             A \code{data.table} containing:
#'                         \itemize{
#'                          \item{labels in column 1}
#'                          \item{total amounts in column 2}
#'                          \item{subset amounts in column 3}
#'                         }
#' @param round.unit       An \code{integer} to specify how many decimals should
#'                         be kept when rounding percentages
#'                         (Default: round.unit = 2).
#' @param rev.stack        A \code{logical} to specify whether the stacking
#'                         order of bars should be reversed or not
#'                         (Warning: color display order will remain unchanged.
#'                         Default: rev.stack = FALSE).
#' @param invert.percent   A \code{logical} to specify whether the calculated
#'                         percentage should reflect the share of the subset
#'                         initially defined in the data.table
#'                         (Default: invert.percent = FALSE), or should reflect
#'                         the remaining share of the total amount minus the
#'                         subset (invert.percent = TRUE).
#' @param horizontal       A \code{logical} to specify whether the plot should
#'                         be plotted horizontally or vertically
#'                         (Default: horizontal = FALSE).
#' @param log.scaled       A \code{logical} to specify whether the plot should
#'                         be log-scaled or not (Default: log.scaled = FALSE).
#' @param decreasing.order A \code{logical} to specify how the bars should be
#'                         ordered. By default, bars are ordered by increasing
#'                         order of the sum of the values they display. This
#'                         order can be inverted with: decreasing.order = TRUE.
#' @return A \code{gg} stacked barplot with annotations.
#' @author Yoann Pageaud.
#' @export

#TODO: Add examples!
ggcoverage <- function(
  data, round.unit = 2, rev.stack = FALSE, invert.percent = FALSE,
  horizontal = FALSE, log.scaled = FALSE, decreasing.order = FALSE){
  colnames(data)[1:3] <- c("IDs", "Total", "Subset")
  #Replace NAs by zeros
  data[is.na(Total), c("Total", "Subset") := 0]
  data[is.na(Subset), Subset := 0]
  #Calculate remaining amounts
  data[, remainings := .(Total-Subset)]
  if(invert.percent){
    data[, percents := .(round((Subset/Total)*100, round.unit))]
  } else {
    data[, percents := .(round((remainings/Total)*100, round.unit))]
  }
  data[, percents := .(paste0(percents, "%"))]
  if(log.scaled){
    data[Subset == 0, Subset := 1]
    if(rev.stack){
      data[, c("logTotal", "logSubset") := .(log10(Total), log10(Subset))]
      data[, logremainings := .(logTotal-logSubset)]
    } else {
      data[, c("logTotal", "logremainings") := .(
        log10(Total), log10(remainings))]
      data[, logSubset := .(logTotal-logremainings)]
    }
    data <- melt(data, id.vars = c(
      "IDs", "Total", "Subset", "logTotal", "remainings", "percents"),
      measure.vars = c("logSubset", "logremainings"))
  } else {
    data <- melt(data, id.vars = c("IDs", "Total", "percents"),
                 measure.vars = c("Subset", "remainings"))
  }
  if(!rev.stack){ #Change order for value before cumsum
    data <- data[order(variable, decreasing = TRUE)]
  }
  data[, label_ypos := .(cumsum(value) - 0.5*value), by = IDs]
  if(log.scaled){
    data[variable == "logremainings", value.char := remainings]
    data[variable == "logSubset", value.char := Subset]
  } else { data[, value.char := .(as.character(value))] }
  #Set the orientation of interest
  if (horizontal) {
    display.count.cutoff <- 0.04
    coeff.max.margin <- 0.1
  } else {
    display.count.cutoff <- 0.02
    coeff.max.margin <- 0.05
  }
  if(log.scaled){
    data$IDs <-
      factor(data$IDs, levels = data[variable == "logremainings"][
        order(Total, -value, IDs, decreasing = decreasing.order)]$IDs)
  } else {
    data$IDs <-
      factor(data$IDs, levels = data[variable == "remainings"][
        order(Total, -value, IDs, decreasing = decreasing.order)]$IDs)
  }
  if(rev.stack){
    data$variable <- factor(data$variable, levels = rev(levels(data$variable)))
  }
  #Removing duplicated strings to not display it
  data[, filter.val := .(value - display.count.cutoff*max(data$value))]
  warn.handle(
    pattern = "Coercing 'character' RHS to 'double' to match the type of the target column",
    data[filter.val < 0, value.char := " "])
  data[variable == "Subset", percents := " "]
  if(log.scaled){ data[, "Total" := logTotal] }
  #Barplot
  ggcov <- ggplot(data = data, aes(x = IDs, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    geom_text(aes(y = label_ypos, label = value.char), vjust = 0.5,
              color = "white", size = 4, fontface = "bold")
  if (horizontal) {
    ggcov <- ggcov +
      geom_text(aes(y = Total, label = percents), hjust = -0.1) +
      coord_flip()
  } else {
    ggcov <- ggcov +
      geom_text(aes(y = Total, label = percents), vjust = -1, hjust = 0.5)
  }
  return(ggcov)
}
