
#' Checks classes of a prcomp object to extract loadings data.
#'
#' @param prcomp.res A PCA result of classes \code{prcomp} or
#'                   \code{irlba_prcomp} resulting from stats::prcomp() or
#'                   irlba::prcomp_irlba().
#' @return A \code{data.table} containing data to plot loadings on a PCA biplot.
#' @author Yoann Pageaud.
#' @keywords internal

get_loadings <- function(prcomp.res, PCs){
  if(methods::is(prcomp.res, "prcomp")){
    if(methods::is(prcomp.res, "irlba_prcomp")){
      loadings.data <- cbind(
        "labels" = names(prcomp.res$center),
        data.table::as.data.table(prcomp.res$rotation))[
          , c("labels", PCs), with = FALSE]
    } else {
      loadings.data <- data.table::as.data.table(
        prcomp.res$rotation, keep.rownames = "labels")[
          , c("labels", PCs), with = FALSE]
    }
  } else { stop("'prcomp.res' is not an object of class 'prcomp'.") }
  return(loadings.data)
}

#' Change colnames in annotation data.table to color, fill, or shape data.
#'
#' @param data       A \code{data.table} containing annotations matching data in
#'                   'prcomp.res'.
#' @param color.data A \code{character} specifying the column name in 'data' to
#'                   be used to map colors to points.
#' @param fill.data  A \code{character} specifying the column name in 'data' to
#'                   be used to map filling colors to points (Important: this
#'                   parameter will only work if associated with 'shape.data',
#'                   and only if input values provided to scale_shape_manual()
#'                   are comprised between 21 and 25).
#' @param shape.data A \code{character} specifying the column name in 'data' to
#'                   be used to map shapes to points.
#' @return A \code{data.table} with all annotations, and the colnames necessary
#'         for color, fill, or shape, changed.
#' @author Yoann Pageaud.
#' @keywords internal

rename_dt_col <- function(data, color.data, fill.data, shape.data){
  #Duplicate data.table
  dt.annot <- data.table::as.data.table(as.data.frame(data))
  #Check if any of color.data, fill.data or shape.data colnames is a unique
  # match
  if(is.null(color.data)){ varval <- NA } else {varval <- color.data}
  if(is.null(fill.data)){ varval <- c(varval, NA)
  } else {varval <- c(varval, fill.data)}
  if(is.null(shape.data)){ varval <- c(varval, NA)
  } else {varval <- c(varval, shape.data)}
  dt.varia <- data.table::data.table(
    varnames = c("color.data", "fill.data", "shape.data"), varval)
  if(any(table(dt.varia$varval) < 2)){
    #Get variables that should be set as new colnames
    vars.to.set <- dt.varia[varval %in% names(table(dt.varia$varval)[
      table(dt.varia$varval) == 1])]$varnames
    #Set new colnames
    lapply(X = vars.to.set, FUN = function(v){
      if(v == "color.data"){
        if(color.data %in% colnames(dt.annot)){
          data.table::setnames(
            x = dt.annot, old = color.data, new = "color.data")
          if(!is.numeric(dt.annot$color.data)){
            if(is.null(levels(dt.annot$color.data))){
              dt.annot[, color.data := as.factor(color.data)]
            }
          }
        } else { stop(
          "'color.data' does not match any column name in 'data'.") }
      } else if(v == "fill.data"){
        if(fill.data %in% colnames(dt.annot)){
          data.table::setnames(
            x = dt.annot, old = fill.data, new = "fill.data")
          if(is.null(levels(dt.annot$fill.data))){
            dt.annot[, fill.data := as.factor(fill.data)]
          }
        } else { stop(
          "'fill.data' does not match any column name in 'data'.") }
      } else if(v == "shape.data"){
        if(shape.data %in% colnames(dt.annot)){
          data.table::setnames(
            x = dt.annot, old = shape.data, new = "shape.data")
          if(is.null(levels(dt.annot$shape.data))){
            dt.annot[, shape.data := as.factor(shape.data)]
          }
        } else { stop(
          "'fill.data' does not match any column name in 'data'.") }
      }
    })
  }
  #Set priorities over variables when there is no unique colname match for
  # each: color > fill > shape
  if(any(table(dt.varia$varval) >= 2)){
    vars.to.set <- dt.varia[varval %in% names(table(dt.varia$varval)[
      table(dt.varia$varval) >= 2])]$varnames
    if(isTRUE(all.equal(
      target = c("color.data", "shape.data"), current = vars.to.set))){
      if(color.data %in% colnames(dt.annot)){
        data.table::setnames(
          x = dt.annot, old = color.data, new = "color.data")
        if(is.null(levels(dt.annot$color.data))){
          dt.annot[, color.data := as.factor(color.data)]
        }
      } else {
        stop("'color.data' does not match any column name in 'data'.") }
    } else if(isTRUE(all.equal(
      target = c("color.data", "fill.data"), current = vars.to.set))){
      if(color.data %in% colnames(dt.annot)){
        data.table::setnames(
          x = dt.annot, old = color.data, new = "color.data")
        if(is.null(levels(dt.annot$color.data))){
          dt.annot[, color.data := as.factor(color.data)]
        }
      } else {
        stop("'color.data' does not match any column name in 'data'.") }
    } else if(isTRUE(all.equal(
      target = c("fill.data", "shape.data"), current = vars.to.set))){
      if(fill.data %in% colnames(dt.annot)){
        data.table::setnames(
          x = dt.annot, old = fill.data, new = "fill.data")
        if(is.null(levels(dt.annot$fill.data))){
          dt.annot[, fill.data := as.factor(fill.data)]
        }
      } else {
        stop("'fill.data' does not match any column name in 'data'.") }
    } else if(isTRUE(all.equal(
      target = c("color.data", "fill.data", "shape.data"),
      current = vars.to.set))){
      if(color.data %in% colnames(dt.annot)){
        data.table::setnames(
          x = dt.annot, old = color.data, new = "color.data")
        if(is.null(levels(dt.annot$color.data))){
          dt.annot[, color.data := as.factor(color.data)]
        }
      } else {
        stop("'color.data' does not match any column name in 'data'.") }
    } else { stop("Unsupported parameters combination.") }
  }
  return(dt.annot)
}

#' Collects and computes needed metrics for PCA biplot.
#'
#' @param prcomp.res A PCA result of classes \code{prcomp} or
#'                   \code{irlba_prcomp} resulting from stats::prcomp() or
#'                   irlba::prcomp_irlba().
#' @param dt.annot   A \code{data.table} containing annotations matching data in
#'                   'prcomp.res'.
#' @param PCs        An \code{integer} vector matching principal components to
#'                   be used to generate the PCA biplot.
#' @param scale      A \code{double} scaling parameter, disabled by 0.
#' @return A \code{list} containing PC IDs, scaled values of the PCA, and PCs
#'         labels to be displayed on the margins of each axis.
#' @author Yoann Pageaud.
#' @export

prepare_pca_data <- function(prcomp.res, dt.annot, PCs, scale){
  # Checking the number of PCs requested
  if(length(PCs) > length(prcomp.res$sdev)){
    warning(
      paste("There are more components specified than actually available.",
            "Reducing the number of PCs to all existing ones."))
    PCs <- seq_along(prcomp.res$sdev)
  }
  # Make PC names vector
  PC <- paste0("PC", PCs)
  # Get sdev from selected PCs
  lam <- prcomp.res$sdev[PCs]
  # Create scaling factor
  lam <- lam * sqrt(nrow(prcomp.res$x))
  lam <- lam^scale
  # Scale PCA data
  dt.scaled.pc <- data.table::as.data.table(t(t(prcomp.res$x[, PC])/lam))
  # dt.scaled.pc <- data.table::merge.data.table(
  #   x = dt.scaled.pc, y = dt.annot, by.x = "rn", by.y = colnames(dt.annot)[1])
  # # Reorder dt.scaled.pc
  # ordered.cols <- colnames(dt.scaled.pc)[
  #   c(2:(length(PCs) + 1), 1, (length(PCs) + 2):ncol(dt.scaled.pc))]
  # dt.scaled.pc <- dt.scaled.pc[, ..ordered.cols, ]
  # # Restore old key colname
  # setnames(x = dt.scaled.pc, old = "rn", new = colnames(dt.annot)[1])
  dt.scaled.pc <- cbind(dt.scaled.pc, dt.annot)
  
  # Calculate the percentage of variability explained by the principal component
  ve <- prcomp.res$sdev^2/sum(prcomp.res$sdev^2)
  lab.PC <- paste0(PC, " (", round(ve[PCs] * 100, 2), "%)")
  
  ls.res <- list(
    "PC" = PC, "scaled_PC" = dt.scaled.pc, "labels" = lab.PC,
    "var.explained" = ve[PCs])
  return(ls.res)
}

#' ggplot2 theme for PCA biplot functions.
#'
#' @return A \code{ggplot2} object themed for biplot functions.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

theme_biplot <- function(){
  ggplot2::ggplot() +  ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.grid = ggplot2::element_line(colour = "grey"),
    axis.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 13),
    legend.text = ggplot2::element_text(size = 12),
    legend.key = ggplot2::element_blank())
}

#' Draws loadings from a PCA.
#'
#' @param loadings.data A \code{data.table} containing data necessary to plot
#'                      loadings.
#' @param loadings.col  A \code{character} specifying a color to be used for
#'                      loadings.
#' @return A \code{ggplot2} plot of PCA loadings.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

ggloadings <- function(ggbiplot, loadings.data, loadings.col){
  #Fix BiocCheck() complaining about these objects initialization
  PCx <- NULL
  PCy <- NULL
  #Add loadings and labels to PCA scatter plot
  ggbiplot + ggplot2::geom_segment(
    data = loadings.data, mapping = ggplot2::aes(
      x = 0, y = 0, xend = PCx, yend = PCy), arrow = grid::arrow(
        length = grid::unit(8, "points")), colour = loadings.col) +
    ggrepel::geom_label_repel(
      data = loadings.data, mapping = ggplot2::aes(
        x = PCx, y = PCy, label = labels), size = 3)
}

#' Computes and draws a custom PCA biplot.
#'
#' @param prcomp.res         A PCA result of classes \code{prcomp} or
#'                           \code{irlba_prcomp} resulting from stats::prcomp()
#'                           or irlba::prcomp_irlba().
#' @param data               A \code{data.table} containing annotations matching
#'                           data in 'prcomp.res'.
#' @param PCx                An \code{integer} matching the principal component
#'                           values to display on X-axis.
#' @param PCy                An \code{integer} matching the principal component
#'                           values to display on Y-axis.
#' @param scale              A \code{double} scaling parameter, disabled by 0.
#' @param point.size         A \code{double} specifying the size of points.
#' @param color.data         A \code{character} specifying the column name in
#'                           'data' to be used to map colors to points.
#' @param fill.data          A \code{character} specifying the column name in
#'                           'data' to be used to map filling colors to points
#'                           (Important: this parameter will only work if
#'                           associated with 'shape.data', and only if input
#'                           values provided to scale_shape_manual() are
#'                           comprised between 21 and 25).
#' @param shape.data         A \code{character} specifying the column name in
#'                           'data' to be used to map shapes to points.
#' @param loadings           A \code{logical} specifying whether the loadings
#'                           should be displayed (TRUE) or not (FALSE).
#' @param loadings.col       A \code{character} specifying a color to be used
#'                           for loadings.
#' @param top.load.by.quad   An \code{integer} specifying the top n most
#'                           important loadings to be displayed in the four
#'                           quadrants of the biplot graph (by quadrants). This
#'                           parameters allows to display only the most
#'                           important loadings, and to hide the less important
#'                           ones, to improve visibility when there is too many
#'                           of them.
#' @param load.above.x       A \code{double} specifying the minimum value of the
#'                           X-coordinate for loadings to be displayed on the
#'                           plot.
#' @param load.above.y       A \code{double} specifying the minimum value of the
#'                           Y-coordinate for loadings to be displayed on the
#'                           plot.
#' @param load.below.x       A \code{double} specifying the maximum value of the
#'                           X-coordinate for loadings to be displayed on the
#'                           plot.
#' @param load.below.y       A \code{double} specifying the maximum value of the
#'                           Y-coordinate for loadings to be displayed on the
#'                           plot.
#' @return A \code{gg} plot of a PCA biplot.
#' @author Yoann Pageaud.
#' @importFrom data.table `:=` `.SD`
#' @export
#' @examples
#' #Get PCA results
#' pca.res <- prcomp(iris[,-5])
#' #Draw the simplest biplot of PC1 and PC2:
#' ggbipca(prcomp.res = pca.res, data = iris)
#' #Draw a biplot of PC2 and PC3:
#' ggbipca(prcomp.res = pca.res, data = iris, PCx = 2, PCy = 3)
#' #Change data scale:
#' ggbipca(prcomp.res = pca.res, data = iris, scale = 2)
#' #Change points size:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2)
#' #Color points following Species:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         color.data = "Species")
#' #Change points shape following Species:
#' # (could use any other categorical / ordinal value if others available)
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         color.data = "Species", shape.data = "Species")
#' #Map custom colors:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'        color.data = "Species", shape.data = "Species") +
#'   scale_color_manual(values = c("green", "red", "orange"))
#' #Map custom point shapes:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         color.data = "Species", shape.data = "Species") +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' #Map custom filling colors:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         fill.data = "Species", shape.data = "Species") +
#'   scale_fill_manual(values = c("green", "royalblue", "orange")) +
#'   scale_shape_manual(values = c(21, 22, 23)) # Needs values between 21 and 25
#' #Enable custom filling colors and fix point shapes:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         fill.data = "Species", shape.data = "Species") +
#'   scale_fill_manual(values = c("green", "royalblue", "orange")) +
#'   scale_shape_manual(values = rep(21, 3)) # Maps shape 21 to 3 categories
#' #Map filling colors, outline colors, and point shapes all together:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         fill.data = "Species", shape.data = "Species",
#'         color.data = "Species") +
#'   scale_fill_manual(values = c("green", "royalblue", "orange")) +
#'   scale_color_manual(values = c("black", "red", "blue")) +
#'   scale_shape_manual(values = c(21, 22, 23)) # Needs values between 21 and 25
#' #Show loadings:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         color.data = "Species", shape.data = "Species", loadings = TRUE) +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' #Change loadings color:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         color.data = "Species", shape.data = "Species", loadings = TRUE,
#'         loadings.col = "purple") +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' #Display the top 1 loading in each quadrant:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         color.data = "Species", shape.data = "Species", loadings = TRUE,
#'         loadings.col = "purple", top.load.by.quad = 1) +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' #Display loadings for which X-coordinates are above 0.05 and Y-coordinates
#' # are above 0:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         color.data = "Species", shape.data = "Species", loadings = TRUE,
#'         loadings.col = "purple", load.above.x = 0.05, load.above.y = 0) +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' #Display top 1 loading in each quadrant, for which X-coordinates are
#' # above 0.1:
#' ggbipca(prcomp.res = pca.res, data = iris, point.size = 2,
#'         color.data = "Species", shape.data = "Species", loadings = TRUE,
#'         loadings.col = "purple", top.load.by.quad = 1, load.above.x = 0.1) +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' @references
#' \itemize{
#'  \item{\href{https://cran.r-project.org/web/packages/ggfortify/index.html}{ggfortify: Data Visualization Tools for Statistical Analysis Results}}
#'  \item{\href{https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html}{ggfortify: Plotting PCA (Principal Component Analysis)}}
#'  \item{\href{https://stats.stackexchange.com/questions/143905/loadings-vs-eigenvectors-in-pca-when-to-use-one-or-another}{Loadings vs eigenvectors in PCA: when to use one or another?}}
#'  \item{\href{https://stats.stackexchange.com/questions/119746/what-is-the-proper-association-measure-of-a-variable-with-a-pca-component-on-a/}{What is the proper association measure of a variable with a PCA component?}}
#' }

ggbipca <- function(
  prcomp.res, data, PCx = 1, PCy = 2, scale = 1, point.size = 1,
  color.data = NULL, fill.data = NULL, shape.data = NULL, loadings = FALSE,
  loadings.col = "red", top.load.by.quad = NULL, load.above.x = NULL,
  load.above.y = NULL, load.below.x = NULL, load.below.y = NULL){
  #Fix BiocCheck() complaining about these objects initialization
  load.sqrd.length <- NULL
  quadrant <- NULL
  
  dt.annot <- BiocompR:::rename_dt_col(
    data = data, color.data = color.data, fill.data = fill.data,
    shape.data = shape.data)
  
  prepare.res <- BiocompR:::prepare_pca_data(
    prcomp.res = prcomp.res, dt.annot = dt.annot, PCs = c(PCx, PCy),
    scale = scale)
  PC <- prepare.res$PC
  dt.scaled.pc <- prepare.res$scaled_PC
  lab.PC <- prepare.res$labels
  
  data.table::setnames(x = dt.scaled.pc, old = PC[1], new = "PCx")
  data.table::setnames(x = dt.scaled.pc, old = PC[2], new = "PCy")
  #Get loadings data for PCx & PCy
  if(loadings){
    # Get loadings data
    loadings.data <- BiocompR:::get_loadings(prcomp.res = prcomp.res, PCs = PC)
    
    data.table::setnames(x = loadings.data, old = PC[1], new = "PCx")
    data.table::setnames(x = loadings.data, old = PC[2], new = "PCy")
    #Define scale for plot data and loadings based on selected PCs
    scaler <- min(
      max(abs(dt.scaled.pc[["PCx"]]))/max(abs(loadings.data[["PCx"]])),
      max(abs(dt.scaled.pc[["PCy"]]))/max(abs(loadings.data[["PCy"]])))
    loadings.data[, 2:3 := lapply(
      X = .SD, FUN = function(i){ i * scaler * 0.8 }), .SDcols = c(2, 3)]
    
    if(!is.null(top.load.by.quad)){
      #Compute loadings length
      loadings.data[, load.sqrd.length := PCx^2 + PCy^2]
      #Assign quadrants
      loadings.data[PCx >= 0 & PCy >= 0, quadrant := "top-right"]
      loadings.data[PCx >= 0 & PCy < 0, quadrant := "bottom-right"]
      loadings.data[PCx < 0 & PCy < 0, quadrant := "bottom-left"]
      loadings.data[PCx < 0 & PCy >= 0, quadrant := "top-left"]
      #Keep top N longest arrows by quadrant
      loadings.data <- loadings.data[order(quadrant, -load.sqrd.length)]
      loadings.data <- loadings.data[, utils::head(
        .SD, top.load.by.quad), by = quadrant]
    }
    #Subset the loadings displayed if some cut-off
    if(!is.null(load.above.x) | !is.null(load.below.x)){
      if(!is.null(load.above.x) & is.null(load.below.x)){
        loadings.data <- loadings.data[PCx >= load.above.x]
      } else if(!is.null(load.below.x) & is.null(load.above.x)){
        loadings.data <- loadings.data[PCx <= load.below.x]
      } else {
        stop("Cannot process contradictory conditions on loadings display.")
      }
    }
    if(!is.null(load.above.y) | !is.null(load.below.y)){
      if(!is.null(load.above.y) & is.null(load.below.y)){
        loadings.data <- loadings.data[PCy >= load.above.y]
      } else if(!is.null(load.below.y) & is.null(load.above.y)){
        loadings.data <- loadings.data[PCy <= load.below.y]
      } else {
        stop("Cannot process contradictory conditions on loadings display.")
      }
    }
  }
  #Make PCA ggplot
  biplt <- theme_biplot() + ggplot2::theme(
    axis.ticks = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 13))
  
  if(!is.null(color.data) & is.null(shape.data) & is.null(fill.data)){
    biplt <- biplt +
      #Draw sample distribution
      ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
        x = PCx, y = PCy, color = color.data), size = point.size) +
      ggplot2::labs(x = lab.PC[1], y = lab.PC[2], color = color.data)
  } else if(!is.null(color.data) & !is.null(shape.data) & is.null(fill.data)){
    if(color.data != shape.data){
      biplt <- biplt +
        #Draw sample distribution
        ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
          x = PCx, y = PCy, color = color.data, shape = shape.data),
          size = point.size) +
        ggplot2::labs(
          x = lab.PC[1], y = lab.PC[2], shape = shape.data,
          color = color.data)
    } else {
      biplt <- biplt +
        #Draw sample distribution
        ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
          x = PCx, y = PCy, color = color.data, shape = color.data),
          size = point.size) +
        ggplot2::labs(
          x = lab.PC[1], y = lab.PC[2], shape = color.data,
          color = color.data)
    }
  } else if(is.null(color.data) & is.null(shape.data) & is.null(fill.data)){
    biplt <- biplt +
      ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
        x = PCx, y = PCy), color = "black", size = point.size) +
      ggplot2::labs(x = lab.PC[1], y = lab.PC[2])
  } else if(is.null(color.data) & !is.null(shape.data) & is.null(fill.data)){
    biplt <- biplt +
      ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
        x = PCx, y = PCy, shape = shape.data),
        color = "black", size = point.size) +
      ggplot2::labs(x = lab.PC[1], y = lab.PC[2], shape = shape.data)
  } else if(!is.null(color.data) & is.null(shape.data) & !is.null(fill.data)){
    biplt <- biplt +
      #Draw sample distribution
      ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
        x = PCx, y = PCy, color = color.data), size = point.size) +
      ggplot2::labs(x = lab.PC[1], y = lab.PC[2], color = color.data)
  } else if(
    !is.null(color.data) & !is.null(shape.data) & !is.null(fill.data)){
    if(color.data != shape.data){
      if(color.data != fill.data){
        biplt <- biplt +
          #Draw sample distribution
          ggplot2::geom_point(
            data = dt.scaled.pc, mapping = ggplot2::aes(
              x = PCx, y = PCy, color = color.data,
              shape = shape.data, fill = fill.data),
            size = point.size) +
          ggplot2::labs(
            x = lab.PC[1], y = lab.PC[2], shape = shape.data,
            color = color.data, fill = fill.data)
      } else {
        biplt <- biplt +
          #Draw sample distribution
          ggplot2::geom_point(
            data = dt.scaled.pc, mapping = ggplot2::aes(
              x = PCx, y = PCy, color = color.data,
              shape = shape.data, fill = color.data),
            size = point.size) +
          ggplot2::labs(
            x = lab.PC[1], y = lab.PC[2], shape = shape.data,
            color = color.data, fill = color.data)
      }
    } else {
      if(color.data != fill.data){
        biplt <- biplt +
          #Draw sample distribution
          ggplot2::geom_point(
            data = dt.scaled.pc, mapping = ggplot2::aes(
              x = PCx, y = PCy, color = color.data,
              shape = color.data, fill = fill.data),
            size = point.size) +
          ggplot2::labs(
            x = lab.PC[1], y = lab.PC[2], shape = color.data,
            color = color.data, fill = fill.data)
      } else {
        biplt <- biplt +
          #Draw sample distribution
          ggplot2::geom_point(
            data = dt.scaled.pc, mapping = ggplot2::aes(
              x = PCx, y = PCy, color = color.data,
              shape = color.data, fill = color.data),
            size = point.size) +
          ggplot2::labs(
            x = lab.PC[1], y = lab.PC[2], shape = color.data,
            color = color.data, fill = color.data)
      }
    }
  } else if(is.null(color.data) & is.null(shape.data) & !is.null(fill.data)){
    stop("Cannot set filling color if point shapes are not set to values between 21 to 25.")
  } else if(is.null(color.data) & !is.null(shape.data) & !is.null(fill.data)){
    if(shape.data != fill.data){
      biplt <- biplt +
        ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
          x = PCx, y = PCy, shape = shape.data, fill = fill.data),
          color = "black", size = point.size) +
        ggplot2::labs(
          x = lab.PC[1], y = lab.PC[2], shape = shape.data,
          fill = fill.data)
    } else {
      biplt <- biplt +
        ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
          x = PCx, y = PCy, shape = fill.data, fill = fill.data),
          color = "black", size = point.size) +
        ggplot2::labs(
          x = lab.PC[1], y = lab.PC[2], shape = fill.data,
          fill = fill.data)
    }
  }
  #Draw loadings
  if(loadings){
    biplt <- ggloadings(
      ggbiplot = biplt, loadings.data = loadings.data,
      loadings.col = loadings.col)
  }
  #Return PCA biplot
  return(biplt)
}

#' Computes and draws biplots for multiple principal components at once.
#'
#' @param prcomp.res         A PCA result of classes \code{prcomp} or
#'                           \code{irlba_prcomp} resulting from stats::prcomp()
#'                           or irlba::prcomp_irlba().
#' @param data               A \code{data.table} containing annotations matching
#'                           data in 'prcomp.res'.
#' @param PCs                An \code{integer} vector matching principal
#'                           components to be used to generate the cross-biplot.
#' @param scale              A \code{double} scaling parameter, disabled by 0.
#' @param point.size         A \code{double} specifying the size of points.
#' @param color.data         A \code{character} specifying the column name in
#'                           'data' to be used to map colors to points.
#' @param shape.data         A \code{character} specifying the column name in
#'                           'data' to be used to map shapes to points.
#' @param loadings           A \code{logical} specifying whether the loadings
#'                           should be displayed (TRUE) or not (FALSE).
#' @param loadings.col       A \code{character} specifying a color to be used
#'                           for loadings.
#' @param top.load.by.quad   An \code{integer} specifying the top n most
#'                           important loadings to be displayed in the four
#'                           quadrants of the biplot graph (by quadrants). This
#'                           parameters allows to display only the most
#'                           important loadings, and to hide the less important
#'                           ones, to improve visibility when there is too many
#'                           of them.
#' @return A \code{gg} plot displaying all combinations of biplots for the
#'         selected PCs.
#' @author Yoann Pageaud.
#' @importFrom data.table `:=` `.SD`
#' @export
#' @examples
#' #Get PCA results
#' pca.res <- prcomp(iris[,-5])
#' #Draw the simplest cross-biplot:
#' cross.biplot(prcomp.res = pca.res, data = iris)
#' #Draw a biplot of PC1 to PC4:
#' cross.biplot(prcomp.res = pca.res, data = iris, PCs = c(1:4))
#' #Change data scale:
#' cross.biplot(prcomp.res = pca.res, data = iris, scale = 2)
#' #Change points size:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2)
#' #Color points following Species:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              color.data = "Species")
#' #Change points shape following Species:
#' # (could use any other categorical / ordinal value if others available)
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              color.data = "Species", shape.data = "Species")
#' #Map custom colors:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              color.data = "Species", shape.data = "Species") +
#'   scale_color_manual(values = c("green", "red", "orange"))
#' #Map custom point shapes:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              color.data = "Species", shape.data = "Species") +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' #Map custom filling colors:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              fill.data = "Species", shape.data = "Species") +
#'   scale_fill_manual(values = c("green", "royalblue", "orange")) +
#'   scale_shape_manual(values = c(21, 22, 23)) # Needs values between 21 and 25
#' #Enable custom filling colors and fix point shapes:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              fill.data = "Species", shape.data = "Species") +
#'   scale_fill_manual(values = c("green", "royalblue", "orange")) +
#'   scale_shape_manual(values = rep(21, 3)) # Maps shape N°21 to 3 categories
#' #Map filling colors, outline colors, and point shapes all together:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              fill.data = "Species", shape.data = "Species",
#'              color.data = "Species") +
#'   scale_fill_manual(values = c("green", "royalblue", "orange")) +
#'   scale_color_manual(values = c("black", "red", "blue")) +
#'   scale_shape_manual(values = c(21, 22, 23)) # Needs values between 21 and 25
#' #Show loadings:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              color.data = "Species", shape.data = "Species",
#'              loadings = TRUE) +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' #Change loadings color:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              color.data = "Species", shape.data = "Species", loadings = TRUE,
#'              loadings.col = "purple") +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))
#' #Display the top 1 loading in each quadrant:
#' cross.biplot(prcomp.res = pca.res, data = iris, point.size = 2,
#'              color.data = "Species", shape.data = "Species", loadings = TRUE,
#'              loadings.col = "purple", top.load.by.quad = 1) +
#'   scale_color_manual(values = c("green", "red", "orange")) +
#'   scale_shape_manual(values = c(83, 8, 25))

cross.biplot <- function(
  prcomp.res, data, PCs = c(1, 2, 3), scale = 1, point.size = 1,
  color.data = NULL, fill.data = NULL, shape.data = NULL, loadings = FALSE,
  loadings.col = "red", top.load.by.quad = NULL){
  #Fix BiocCheck() complaining about these objects initialization
  Var1 <- NULL
  Var2 <- NULL
  . <- NULL
  name.X <- NULL
  PCx <- NULL
  name.Y <- NULL
  PCy <- NULL
  load.sqrd.length <- NULL
  quadrant <- NULL
  
  dt.annot <- BiocompR:::rename_dt_col(
    data = data, color.data = color.data, fill.data = fill.data,
    shape.data = shape.data)
  
  prepare.res <- BiocompR:::prepare_pca_data(
    prcomp.res = prcomp.res, dt.annot = dt.annot, PCs = PCs, scale = scale)
  PC <- prepare.res$PC
  dt.scaled.pc <- prepare.res$scaled_PC
  lab.PC <- prepare.res$labels
  lab.PC <- data.table::data.table(PC, lab.PC)
  #Create possible PCs combinations
  comb.pcs <- data.table::as.data.table(
    expand.grid(PC, PC, stringsAsFactors = FALSE))[Var1 != Var2]
  #Get loadings data for PCx & PCy
  if(loadings){
    # Get loadings data
    loadings.data <- BiocompR:::get_loadings(prcomp.res = prcomp.res, PCs = PC)
    #Define scale for plot data and loadings based on selected PCs
    scaler <- min(unlist(lapply(X = PC, FUN = function(p){
      max(abs(dt.scaled.pc[[p]]))/max(abs(loadings.data[[p]]))
    })))
    loadings.data[, c(PC) := lapply(
      X = .SD, FUN = function(i){ i * scaler * 0.8 }), .SDcols = PC]
    #Recreate loadings.data
    ls.dt.loadings <- lapply(X = seq(nrow(comb.pcs)), FUN = function(i){
      pc.select <- unlist(c("labels", comb.pcs[i]))
      sub.dt <- loadings.data[, ..pc.select, ]
      sub.dt[, c("name.Y", "name.X") := .(pc.select[2], pc.select[3])]
    })
    loadings.data <- data.table::rbindlist(
      l = ls.dt.loadings, use.names = FALSE)
    data.table::setnames(x = loadings.data, old = 2, new = "PCy")
    data.table::setnames(x = loadings.data, old = 3, new = "PCx")
    loadings.data <- loadings.data[, .(labels, name.X, PCx, name.Y, PCy)]
    loadings.data[, c("name.X", "name.Y") := .(
      as.factor(name.X), as.factor(name.Y))]
    data.table::setattr(loadings.data$name.X, "levels", lab.PC[
      order(match(PC, levels(loadings.data$name.X)))]$lab.PC)
    data.table::setattr(loadings.data$name.Y, "levels", lab.PC[
      order(match(PC, levels(loadings.data$name.Y)))]$lab.PC)
    if(!is.null(top.load.by.quad)){
      #Compute loadings length
      loadings.data[, load.sqrd.length := PCx^2 + PCy^2]
      #Assign quadrants
      loadings.data[PCx >= 0 & PCy >= 0, quadrant := "top-right"]
      loadings.data[PCx >= 0 & PCy < 0, quadrant := "bottom-right"]
      loadings.data[PCx < 0 & PCy < 0, quadrant := "bottom-left"]
      loadings.data[PCx < 0 & PCy >= 0, quadrant := "top-left"]
      #Keep top N longest arrows by quadrant and by PC combinations
      loadings.data <- loadings.data[order(
        quadrant, -load.sqrd.length), .SD, by = .(name.X, name.Y)]
      loadings.data <- loadings.data[, utils::head(
        .SD, top.load.by.quad), by = .(name.X, name.Y, quadrant)]
    }
  }
  #Recreate dt.scaled.pc
  ls.dt.scaled <- lapply(X = seq(nrow(comb.pcs)), FUN = function(i){
    pc.select <- unlist(c(comb.pcs[i], colnames(dt.annot)))
    sub.dt <- dt.scaled.pc[, ..pc.select, ]
    sub.dt[, c("name.Y", "name.X") := .(pc.select[1], pc.select[2])]
  })
  dt.scaled.pc <- data.table::rbindlist(l = ls.dt.scaled, use.names = FALSE)
  data.table::setnames(x = dt.scaled.pc, old = 1, new = "PCy")
  data.table::setnames(x = dt.scaled.pc, old = 2, new = "PCx")
  cols <- c("name.X", "PCx", "name.Y", "PCy", colnames(dt.annot))
  dt.scaled.pc <- dt.scaled.pc[, ..cols, ]
  dt.scaled.pc[, c("name.X", "name.Y") := .(
    as.factor(name.X), as.factor(name.Y))]
  data.table::setattr(dt.scaled.pc$name.X, "levels", lab.PC[
    order(match(PC, levels(dt.scaled.pc$name.X)))]$lab.PC)
  data.table::setattr(dt.scaled.pc$name.Y, "levels", lab.PC[
    order(match(PC, levels(dt.scaled.pc$name.Y)))]$lab.PC)
  
  #Make PCA ggplot
  biplt <- theme_biplot() + ggplot2::theme(
    panel.spacing = ggplot2::unit(0, "lines"),
    panel.border = ggplot2::element_rect(
      color = "black", fill = NA, linewidth = 0.5),
    axis.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = -45, hjust = 0, vjust = 0.5),
    strip.text = ggplot2::element_text(size = 13),
    strip.background = ggplot2::element_rect(
      fill = "white", color = "black", linewidth = 0.5)) +
    ggplot2::facet_grid(name.Y ~ name.X, scales = "free", space = "fixed")
  
  if(!is.null(color.data) & is.null(shape.data) & is.null(fill.data)){
    biplt <- biplt +
      #Draw sample distribution
      ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
        x = PCx, y = PCy, color = color.data), size = point.size) +
      ggplot2::labs(color = color.data)
  } else if(!is.null(color.data) & !is.null(shape.data) & is.null(fill.data)){
    if(color.data != shape.data){
      biplt <- biplt +
        #Draw sample distribution
        ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
          x = PCx, y = PCy, color = color.data, shape = shape.data),
          size = point.size) +
        ggplot2::labs(shape = shape.data, color = color.data)
    } else {
      biplt <- biplt +
        #Draw sample distribution
        ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
          x = PCx, y = PCy, color = color.data, shape = color.data),
          size = point.size) +
        ggplot2::labs(shape = color.data, color = color.data)
    }
  } else if(is.null(color.data) & is.null(shape.data) & is.null(fill.data)){
    biplt <- biplt +
      ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
        x = PCx, y = PCy), color = "black", size = point.size)
  } else if(is.null(color.data) & !is.null(shape.data) & is.null(fill.data)){
    biplt <- biplt +
      ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
        x = PCx, y = PCy, shape = shape.data),
        color = "black", size = point.size) +
      ggplot2::labs(shape = shape.data)
  } else if(!is.null(color.data) & is.null(shape.data) & !is.null(fill.data)){
    biplt <- biplt +
      #Draw sample distribution
      ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
        x = PCx, y = PCy, color = color.data), size = point.size) +
      ggplot2::labs(color = color.data)
  } else if(
    !is.null(color.data) & !is.null(shape.data) & !is.null(fill.data)){
    if(color.data != shape.data){
      if(color.data != fill.data){
        biplt <- biplt +
          #Draw sample distribution
          ggplot2::geom_point(
            data = dt.scaled.pc, mapping = ggplot2::aes(
              x = PCx, y = PCy, color = color.data,
              shape = shape.data, fill = fill.data),
            size = point.size) +
          ggplot2::labs(
            shape = shape.data, color = color.data, fill = fill.data)
      } else {
        biplt <- biplt +
          #Draw sample distribution
          ggplot2::geom_point(
            data = dt.scaled.pc, mapping = ggplot2::aes(
              x = PCx, y = PCy, color = color.data,
              shape = shape.data, fill = color.data),
            size = point.size) +
          ggplot2::labs(
            shape = shape.data, color = color.data, fill = color.data)
      }
    } else {
      if(color.data != fill.data){
        biplt <- biplt +
          #Draw sample distribution
          ggplot2::geom_point(
            data = dt.scaled.pc, mapping = ggplot2::aes(
              x = PCx, y = PCy, color = color.data,
              shape = color.data, fill = fill.data),
            size = point.size) +
          ggplot2::labs(
            shape = color.data, color = color.data, fill = fill.data)
      } else {
        biplt <- biplt +
          #Draw sample distribution
          ggplot2::geom_point(
            data = dt.scaled.pc, mapping = ggplot2::aes(
              x = PCx, y = PCy, color = color.data,
              shape = color.data, fill = color.data),
            size = point.size) +
          ggplot2::labs(
            shape = color.data, color = color.data, fill = color.data)
      }
    }
  } else if(is.null(color.data) & is.null(shape.data) & !is.null(fill.data)){
    stop("Cannot set filling color if point shapes are not set to values between 21 to 25.")
  } else if(is.null(color.data) & !is.null(shape.data) & !is.null(fill.data)){
    if(shape.data != fill.data){
      biplt <- biplt +
        ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
          x = PCx, y = PCy, shape = shape.data, fill = fill.data),
          color = "black", size = point.size) +
        ggplot2::labs(shape = shape.data, fill = fill.data)
    } else {
      biplt <- biplt +
        ggplot2::geom_point(data = dt.scaled.pc, mapping = ggplot2::aes(
          x = PCx, y = PCy, shape = fill.data, fill = fill.data),
          color = "black", size = point.size) +
        ggplot2::labs(shape = fill.data, fill = fill.data)
    }
  }
  
  #Draw loadings
  if(loadings){
    biplt <- ggloadings(
      ggbiplot = biplt, loadings.data = loadings.data,
      loadings.col = loadings.col)
  }
  #Return PCA biplot
  return(biplt)
}
