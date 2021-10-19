
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
#' @importFrom data.table `:=`
#' @export
#' @examples
#' # Basic use of ggcoverage()
#' df <- data.frame(
#'  col1 = LETTERS[1:6], col2 = 100:105, col3 = 35:40,
#'  col4 = rep(x = paste0("Group", 1:3), each = 2))
#' ggcoverage(data = df)
#' # Set the number of decimal for the percentage rounding to 0
#' ggcoverage(data = df, round.unit = 0)
#' # Compute percentages on the Subset share instead of the Remaining share
#' ggcoverage(data = df, round.unit = 0, invert.percent = TRUE)
#' # Order bars in decreasing order
#' ggcoverage(
#'  data = df, round.unit = 0, invert.percent = TRUE, decreasing.order = TRUE)
#' # Reverse stackings (Warning: color mapping does not change)
#' ggcoverage(
#'  data = df, round.unit = 0, invert.percent = TRUE, decreasing.order = TRUE,
#'  rev.stack = TRUE)
#' # Apply log-transformation to the barplot
#' ggcoverage(
#'  data = df, round.unit = 0, invert.percent = TRUE, decreasing.order = TRUE,
#'  rev.stack = TRUE, horizontal = TRUE, log.scaled = TRUE)
#' # Add elements of customization
#' ggcoverage(
#'   data = df, round.unit = 0, invert.percent = TRUE, decreasing.order = TRUE,
#'   rev.stack = TRUE, horizontal = TRUE, log.scaled = TRUE) +
#'   ggtitle("This is a ggcoverage barplot!") + # Add title
#'   labs(x = "Samples", y = "Coverages") + #Set X and Y axis titles
#'  theme(
#'     plot.title = element_text(hjust = 0.5),
#'    axis.text = element_text(size = 14, color = "black"), # Custom axis text
#'     axis.title = element_text(size = 15),
#'     # Change legend appearance
#'     legend.title = element_text(size = 13),
#'     legend.text = element_text(size = 12),
#'     # Change panel appearance
#'     panel.background = element_rect(color = "black", fill = NA),
#'     panel.grid.major.y = element_blank(),
#'     panel.grid.minor.y = element_blank(),
#'     panel.grid.major.x = element_line(color = "grey"),
#'     panel.grid.minor.x = element_line(color = "grey"),
#'     # Change strips appearance
#'     strip.background = element_rect(color = "black"),
#'     strip.text.y = element_text(size = 14, angle = 0)) +
#'   scale_y_continuous(
#'    # Expand fully plot panel on X-axis (coordinates have been flipped)
#'     expand = c(0, 0),
#'     # Set new limits (useful if percentages do not appear on the default plot)
#'     limits = c(0, log10(max(df$col2, na.rm = TRUE)) +
#'                  0.1*log10(max(df$col2, na.rm = TRUE))),
#'     # Change number of breaks
#'    breaks = c(0, 1, 2),
#'     # Change X-axis labels
#'     labels = c("0" = 1, "1" = 10, "2.0" = 100)) +
#'   guides(fill = guide_legend(title = "Coverage legend")) + # Set legend title
#'   scale_fill_manual(labels = c("Remaining", "Subset"), # Rename conditions
#'                     values = c("#D6604D", "#4393C3")) + # Change colors
#'   facet_grid(df$col4 ~ ., scales = "free", space = "free_y") # Add grouping

ggcoverage <- function(
  data, round.unit = 2, rev.stack = FALSE, invert.percent = FALSE,
  horizontal = FALSE, log.scaled = FALSE, decreasing.order = FALSE){
  #Fix BiocCheck() complaining about these objects initialization
  Total <- NULL
  Subset <- NULL
  Remaining <- NULL
  . <- NULL
  percents <- NULL
  logRemaining <- NULL
  logTotal <- NULL
  logSubset <- NULL
  variable <- NULL
  label_ypos <- NULL
  value <- NULL
  IDs <- NULL
  value.char <- NULL
  filter.val <- NULL
  #Check if data is a data.table and convert if not
  if(!data.table::is.data.table(data)){
    data <- data.table::as.data.table(data) }
  colnames(data)[seq(3)] <- c("IDs", "Total", "Subset")
  #Replace NAs by zeros
  data[is.na(Total), c("Total", "Subset") := 0]
  data[is.na(Subset), Subset := 0]
  #Calculate remaining amounts
  data[, Remaining := .(Total-Subset)]
  if(invert.percent){
    data[, percents := .(round((Subset/Total)*100, round.unit))]
  } else {
    data[, percents := .(round((Remaining/Total)*100, round.unit))]
  }
  data[, percents := .(paste0(percents, "%"))]
  if(log.scaled){
    data[Subset == 0, Subset := 1]
    if(rev.stack){
      data[, c("logTotal", "logSubset") := .(log10(Total), log10(Subset))]
      data[, logRemaining := .(logTotal-logSubset)]
    } else {
      data[, c("logTotal", "logRemaining") := .(
        log10(Total), log10(Remaining))]
      data[, logSubset := .(logTotal-logRemaining)]
    }
    data <- data.table::melt(data, id.vars = c(
      "IDs", "Total", "Subset", "logTotal", "Remaining", "percents"),
      measure.vars = c("logSubset", "logRemaining"))
  } else {
    data <- data.table::melt(data, id.vars = c("IDs", "Total", "percents"),
                             measure.vars = c("Subset", "Remaining"))
  }
  if(!rev.stack){ #Change order for value before cumsum
    data <- data[order(variable, decreasing = TRUE)]
  }
  data[, label_ypos := .(cumsum(value) - 0.5*value), by = IDs]
  if(log.scaled){
    data[variable == "logRemaining", value.char := Remaining]
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
      factor(data$IDs, levels = data[variable == "logRemaining"][
        order(Total, -value, IDs, decreasing = decreasing.order)]$IDs)
  } else {
    data$IDs <-
      factor(data$IDs, levels = data[variable == "Remaining"][
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
  ggcov <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(
    x = IDs, y = value, fill = variable)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(mapping = ggplot2::aes(
      y = label_ypos, label = value.char), vjust = 0.5, color = "white",
      size = 4, fontface = "bold")
  if (horizontal) {
    ggcov <- ggcov +
      ggplot2::geom_text(mapping = ggplot2::aes(
        y = Total, label = percents), hjust = -0.1) +
      ggplot2::coord_flip()
  } else {
    ggcov <- ggcov +
      ggplot2::geom_text(mapping = ggplot2::aes(y = Total, label = percents),
                         vjust = -1, hjust = 0.5)
  }
  return(ggcov)
}
