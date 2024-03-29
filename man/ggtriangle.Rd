% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_ggplot_fun.R
\name{ggtriangle}
\alias{ggtriangle}
\title{Draws a triangle plot from a basic molten triangle matrix.}
\usage{
ggtriangle(
  melt_tri,
  grid_col = "white",
  grid_linewidth = 0.3,
  ggtri_theme = NULL,
  scale_fill_grad = ggplot2::scale_fill_gradientn(colors = BiocompR::biopalette(name =
    "viridis_C_plasma"), na.value = "grey"),
  guide_custom_bar = ggplot2::guide_colorbar(ticks.linewidth = 0.5, ticks.colour =
    "black", frame.colour = "black", frame.linewidth = 0.5),
  x_axis_pos = "top",
  y_axis_pos = "right"
)
}
\arguments{
\item{melt_tri}{A \code{data.frame} melted triangle containing a
statistical test values.}

\item{grid_col}{A \code{character} specifying the color of the grid.}

\item{grid_linewidth}{A \code{double} value for the thickness of the grid.}

\item{ggtri_theme}{A ggplot2 \code{theme} to specify any theme parameter
you wish to custom on the triangle plot
(Default: ggtri_theme = NULL). For more information
about how to define a theme, see
\link[ggplot2]{theme}.}

\item{scale_fill_grad}{A \code{ScaleContinous} object generated by ggplot2
functions such as \link[ggplot2]{scale_fill_gradient}
,\link[ggplot2]{scale_fill_gradient2} or
\link[ggplot2]{scale_fill_gradientn} to customize
heatmap colors and the associated color bar.}

\item{guide_custom_bar}{A \code{guide} object generated by the ggplot2
function \link[ggplot2]{guide_colorbar} to custom the
triangle plot color bar appearance
(see also 'scale_fill_grad' option).}

\item{x_axis_pos}{A \code{character} to specify the position of the X
axis on the plot (Default: x_axis_pos = "top";
Supported: c("top", "bottom")).}

\item{y_axis_pos}{A \code{character} to specify the position of the Y
axis on the plot (Default: y_axis_pos = "right";
Supported: c("right", "left")).}
}
\value{
A \code{gg} object of a basic triangle plot (a 'geom_tile()').
}
\description{
Draws a triangle plot from a basic molten triangle matrix.
}
\examples{
# Create matrix from mtcars dataset
mat <- as.matrix(t(scale(mtcars)))
# Compute Pearson's correlation between cars
corr.res <- psych::corr.test(
    mat, use = "pairwise", method = "pearson", adjust = "BH")$r
# Order samples following a complete hierarchical clustering
correlation.order <- corrplot::corrMatOrder(
    corr = corr.res, order = "hclust", hclust.method = "complete")
corr.res <- corr.res[correlation.order, correlation.order]
# Remove duplicated data
corr.res[upper.tri(corr.res)] <- NA
# Melt Correlation matrix
dt.corr <- data.table::as.data.table(x = corr.res, keep.rownames = "Var1")
molt.corr <- data.table::melt.data.table(
    data = dt.corr, id.vars = "Var1", variable.name = "Var2", na.rm = TRUE,
    measure.vars = colnames(dt.corr)[-c(1)])
molt.corr[, Var1 := as.factor(x = Var1)]
molt.corr[, Var1 := factor(
    x = Var1, levels = as.character(unique(molt.corr$Var1)))]
# Invert order of samples
molt.corr[, Var2 := factor(x = Var2, levels = rev(levels(Var2)))]
# Replace identical correlations by NAs
molt.corr[Var1 == Var2, value := NA]
# Draw the triangle plot
ggtriangle(melt_tri = molt.corr)
}
\author{
Yoann Pageaud.
}
