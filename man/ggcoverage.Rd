% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggcoverage.R
\name{ggcoverage}
\alias{ggcoverage}
\title{Plots an annotated stacked barplot.}
\usage{
ggcoverage(
  data,
  round.unit = 2,
  rev.stack = FALSE,
  invert.percent = FALSE,
  horizontal = FALSE,
  log.scaled = FALSE,
  decreasing.order = FALSE
)
}
\arguments{
\item{data}{A \code{data.table} containing:
\itemize{
 \item{labels in column 1}
 \item{total amounts in column 2}
 \item{subset amounts in column 3}
}}

\item{round.unit}{An \code{integer} to specify how many decimals should
be kept when rounding percentages
(Default: round.unit = 2).}

\item{rev.stack}{A \code{logical} to specify whether the stacking
order of bars should be reversed or not
(Warning: color display order will remain unchanged.
Default: rev.stack = FALSE).}

\item{invert.percent}{A \code{logical} to specify whether the calculated
percentage should reflect the share of the subset
initially defined in the data.table
(Default: invert.percent = FALSE), or should reflect
the remaining share of the total amount minus the
subset (invert.percent = TRUE).}

\item{horizontal}{A \code{logical} to specify whether the plot should
be plotted horizontally or vertically
(Default: horizontal = FALSE).}

\item{log.scaled}{A \code{logical} to specify whether the plot should
be log-scaled or not (Default: log.scaled = FALSE).}

\item{decreasing.order}{A \code{logical} to specify how the bars should be
ordered. By default, bars are ordered by increasing
order of the sum of the values they display. This
order can be inverted with: decreasing.order = TRUE.}
}
\value{
A \code{gg} stacked barplot with annotations.
}
\description{
Plots an annotated stacked barplot.
}
\examples{
# Basic use of ggcoverage()
df <- data.frame(
 col1 = LETTERS[1:6], col2 = 100:105, col3 = 35:40,
 col4 = rep(x = paste0("Group", 1:3), each = 2))
ggcoverage(data = df)
# Set the number of decimal for the percentage rounding to 0
ggcoverage(data = df, round.unit = 0)
# Compute percentages on the Subset share instead of the Remaining share
ggcoverage(data = df, round.unit = 0, invert.percent = TRUE)
# Order bars in decreasing order
ggcoverage(
 data = df, round.unit = 0, invert.percent = TRUE, decreasing.order = TRUE)
# Reverse stackings (Warning: color mapping does not change)
ggcoverage(
 data = df, round.unit = 0, invert.percent = TRUE, decreasing.order = TRUE,
 rev.stack = TRUE)
# Apply log-transformation to the barplot
ggcoverage(
 data = df, round.unit = 0, invert.percent = TRUE, decreasing.order = TRUE,
 rev.stack = TRUE, horizontal = TRUE, log.scaled = TRUE)
# Add elements of customization
ggcoverage(
  data = df, round.unit = 0, invert.percent = TRUE, decreasing.order = TRUE,
  rev.stack = TRUE, horizontal = TRUE, log.scaled = TRUE) +
  ggtitle("This is a ggcoverage barplot!") + # Add title
  labs(x = "Samples", y = "Coverages") + # Set X and Y axis titles
 theme(
    plot.title = element_text(hjust = 0.5),
   axis.text = element_text(size = 14, color = "black"), # Custom axis text
    axis.title = element_text(size = 15),
    # Change legend appearance
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    # Change panel appearance
    panel.background = element_rect(color = "black", fill = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey"),
    panel.grid.minor.x = element_line(color = "grey"),
    # Change strips appearance
    strip.background = element_rect(color = "black"),
    strip.text.y = element_text(size = 14, angle = 0)) +
  scale_y_continuous(
   # Expand fully plot panel on X-axis (coordinates have been flipped)
    expand = c(0, 0),
    # Set new limits (useful if percentages do not appear on the default plot)
    limits = c(0, log10(max(df$col2, na.rm = TRUE)) +
                 0.1*log10(max(df$col2, na.rm = TRUE))),
    # Change number of breaks
   breaks = c(0, 1, 2),
    # Change X-axis labels
    labels = c("0" = 1, "1" = 10, "2.0" = 100)) +
  guides(fill = guide_legend(title = "Coverage legend")) + # Set legend title
  scale_fill_manual(
    labels = c("Remaining", "Subset"), # Rename conditions
    # Change colors
    values = rev(biopalette(name = "BiocompR_cond", mute = TRUE))) +
  facet_grid(df$col4 ~ ., scales = "free", space = "free_y") # Add grouping
}
\author{
Yoann Pageaud.
}
